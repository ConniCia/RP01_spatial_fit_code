#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "header.h"

int MCMC(int do_symmetric, int do_skipdiagonal, int do_rescalekappa,
         int nrPAR, int iterMCMC, int nr_travellers, int version,
         long double **oldPAR, long double *minPAR, long double *aPAR,
				 long double acceptance_rate, long double tweaking_factor,
				 char output_folder[], int outer_loop_nr,
				 int min_pop,
 				 int nr_admins, int nr_patches,
 				 int *cum_nr_patches, int *patch_population,
 				 int **patch_distance, int **admin_flow_data,
         int **patch_radius_population_L, int **patch_radius_population_S,
         int *admin_nr_travellers, int *admin_population,
 				 long double **admin_flow_model, char *model,
         int iter_flow_model
) {

  /* Define and initialise variables ---------------------------------------- */
  // two sets of parameters: verify whether log-likelihood is improved
  long double newLOGL, oldLOGL = computeLOGL(do_symmetric, do_skipdiagonal, do_rescalekappa,
                                             min_pop, *oldPAR,
                                             nr_admins, nr_patches, nrPAR, nr_travellers, version,
                                             cum_nr_patches, patch_population,
                                             patch_distance, admin_flow_data,
                                             patch_radius_population_L, patch_radius_population_S,
                                             admin_nr_travellers, admin_population,
                                             admin_flow_model, model,
                                             output_folder, 0, 0);
  long double *newPAR = (long double *) malloc(nrPAR * sizeof(long double));
  for (int p = 0; p < nrPAR; p++)
    newPAR[p] = (*oldPAR)[p];
  // maximum parameter-specific step size
  long double *deltaPAR = (long double *) malloc(nrPAR * sizeof(long double));
  for (int p = 0; p < nrPAR; p++) {
    deltaPAR[p] = fabsl((*oldPAR)[p]) / 10.;
    if (deltaPAR[p] < minPAR[p])
      deltaPAR[p] = minPAR[p];
  }
  // adaptive parameter-specific step size ensures "enough" changes get accepted
  int *nr_accepted = (int *) calloc(nrPAR, sizeof(int));
  int tweak_counter=0, max_tweaks=100, iter_tweak=1000;
  double iter_tweak_par = (double) iter_tweak / (double) nrPAR;
  if ( tweaking_factor < 1. ) {
    max_tweaks = 0;
    fprintf(stderr, "Attention: tweaking turned off, parameter step size will not be adjusted during MCMC!!\n");
    fflush(stderr);
    fflush(stderr);
  }


  /* Files and headers ------------------------------------------------------ */
  // MCMC chain
  FILE *file_MCMC;
  char str[500];
  sprintf(str, "%s/MCMC_%d.tsv", output_folder, outer_loop_nr);
  if ((file_MCMC = fopen(str, "w")) == NULL) {
		fprintf(stderr, "Error: file \"MCMC_%d\" not found.\n", outer_loop_nr);
		exit(1);
	}
  fprintf(file_MCMC, "step\tLOGL\t\t");
  for (int p = 0; p < nrPAR; p++)
    fprintf(file_MCMC, "PAR%d\t", p);
  fprintf(file_MCMC, "\n0\t%.2LF\t", oldLOGL);
  for (int p = 0; p < nrPAR; p++)
    fprintf(file_MCMC, "%.2LF\t", (*oldPAR)[p]);
  // adaptive step size: acceptance rate (ar) and step size (deltaPAR)
  FILE *file_tweak;
  sprintf(str, "%s/MCMC_tweak_%d.tsv", output_folder, outer_loop_nr);
  if ((file_tweak = fopen(str, "w")) == NULL) {
		fprintf(stderr, "Error: file \"MCMC_tweak_%d\" not found.\n", outer_loop_nr);
		exit(1);
	}
  fprintf(file_tweak, "step\t");
  for (int p = 0; p < nrPAR; p++)
    fprintf(file_tweak, "AR%d\tdelta%d\t", p, p);

  /* MCMC loop -------------------------------------------------------------- */
  // auxiliary variables
  int j;
  long double delta;
  double ar;
  // loop: change one parameter value at a time
  for (int i = 0; i < iterMCMC; i++) {

    // draw new value uniformly from interval [-delta, +delta] around current value
    delta = ( ((double)rand() / (double)RAND_MAX) - 0.5 ) * 2.;
    j = i % nrPAR;
    newPAR[j] = (*oldPAR)[j] + deltaPAR[j] * delta;

    // discard unadmissible values
    if (newPAR[j] < aPAR[j]) {
      newPAR[j] = (*oldPAR)[j];
    }
    // accept new values if the improve log-likelihood or worsen it only slightly
    else {
      newLOGL = computeLOGL(do_symmetric, do_skipdiagonal, do_rescalekappa,
                            min_pop, newPAR,
                            nr_admins, nr_patches, nrPAR, nr_travellers, version,
                            cum_nr_patches, patch_population,
                            patch_distance, admin_flow_data,
                            patch_radius_population_L, patch_radius_population_S,
                            admin_nr_travellers, admin_population,
                            admin_flow_model, model,
                            output_folder, 0, 0);
      delta = newLOGL - oldLOGL;
      if ( (delta > 0) || ( ((double)rand() / (double)RAND_MAX) < expl(delta) ) ) {
        oldLOGL = newLOGL;
        (*oldPAR)[j] = newPAR[j];
        nr_accepted[j]++;
      }
      else {
        newPAR[j] = (*oldPAR)[j];
      }
    }

    // print MCMC chain
    if ( (i+1) % 100 == 0 ) {
      fprintf(file_MCMC, "\n%d\t%.2LF\t", i+1, oldLOGL);
      for (int p = 0; p < nrPAR; p++)
        fprintf(file_MCMC, "%.2LF\t", (*oldPAR)[p]);
      fflush(file_MCMC); // force immediate printing
			fflush(file_MCMC);
    }

    // modify parameter-specific step size according to number of accepted changes
    // N.B. do this only at the beginning of each chain to ensure good mixing
    if ( (i+1) % iter_tweak == 0 && tweak_counter < max_tweaks ) {
      tweak_counter++;
      fprintf(file_tweak, "\n%d\t", i+1);
      fflush(file_tweak); // force immediate printing
			fflush(file_tweak);

      // 1. too many accepted changes: increase step size
      // 2. + 3. parameter value is stuck: decrease step size
      for (int p = 0; p < nrPAR; p++) {
        ar = (double)nr_accepted[p] / iter_tweak_par; // acceptance rate
        if ( ar > acceptance_rate )
          deltaPAR[p] *= tweaking_factor;
        else if ( deltaPAR[p] / tweaking_factor > minPAR[p] )
          deltaPAR[p] /= tweaking_factor;
        else
          deltaPAR[p] = minPAR[p];
        fprintf(file_tweak, "%.2F\t%.2LF\t", ar, deltaPAR[p]);
        // reset parameter-specific counter
        nr_accepted[p] = 0;
      }
    }

  } // end MCMC iteration

  fclose(file_MCMC);
  fclose(file_tweak);

  free(newPAR);
  free(deltaPAR);
  free(nr_accepted);

  return 0;
}

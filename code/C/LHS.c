#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "header.h"

int LHS(int do_symmetric, int do_skipdiagonal, int do_rescalekappa,
        int iterLHS, int nr_chains, int nrPAR, int nr_travellers, int version,
        long double *aPAR, long double *bPAR, long double *PAR, long double ***chainPAR,
				char output_folder[], int min_pop,
				int nr_admins, int nr_patches,
				int *cum_nr_patches, int *patch_population,
				int **patch_distance, int **admin_flow_data,
        int **patch_radius_population_L, int **patch_radius_population_S,
        int *admin_nr_travellers, int *admin_population,
				long double **admin_flow_model, char *model,
        int iter_flow_model
) {

  /* Define and initialise variables ---------------------------------------- */
  // determine length of sub-intervals
  long double *lenPAR = (long double *) malloc(nrPAR * sizeof(long double));
  for (int p = 0; p < nrPAR; p++)
    lenPAR[p] = ( bPAR[p] - aPAR[p] ) / (long double) iterLHS;
  // define random order of subintervals
  int **randPAR = (int **) malloc(nrPAR * sizeof(int *));
  for (int p = 0; p < nrPAR; p++) {
    randPAR[p] = (int *) malloc(iterLHS * sizeof(int));
    for (int j = 0; j < iterLHS; j++) randPAR[p][j] = j;
    randomise(&randPAR[p], iterLHS);
  }
  // initialise log-likelihood values for each chain to something I'm sure to improve
  long double LOGL = computeLOGL(do_symmetric, do_skipdiagonal, do_rescalekappa,
                     min_pop, PAR,
                     nr_admins, nr_patches, nrPAR, nr_travellers, version,
                     cum_nr_patches, patch_population,
                     patch_distance, admin_flow_data,
                     patch_radius_population_L, patch_radius_population_S,
                     admin_nr_travellers, admin_population,
                     admin_flow_model, model,
                     output_folder, 0, 0);
  long double *chainLOGL = (long double *) malloc(nr_chains * sizeof(long double));
  for (int i = 0; i < nr_chains; i++)
    chainLOGL[i] = LOGL * (1+i);

  /* File and header -------------------------------------------------------- */
  FILE *file;
  char str[500];
  sprintf(str, "%s/LHS.tsv", output_folder);
  if ((file = fopen(str, "w")) == NULL) {
		fprintf(stderr, "Error: file \"LHS\" not found.\n");
		exit(1);
	}
  fprintf(file, "step\tLOGL\t\t");
  for (int p = 0; p < nrPAR; p++)
    fprintf(file, "PAR%d\t", p);
  fprintf(file, "\n0\t%.2LF\t", LOGL);
  for (int p = 0; p < nrPAR; p++)
    fprintf(file, "%.2LF\t", PAR[p]);

  /* LHS loop --------------------------------------------------------------- */
  long double delta;
  for (int i = 0; i < iterLHS; i++) {

    // draw new values uniformly from chosen intervals
    for (int p = 0; p < nrPAR; p++) {
      delta = ((double) rand() / (double) RAND_MAX);
      PAR[p] = aPAR[p] + (randPAR[p][i] + delta) * lenPAR[p];
    }

    LOGL = computeLOGL(do_symmetric, do_skipdiagonal, do_rescalekappa,
                       min_pop, PAR,
                       nr_admins, nr_patches, nrPAR, nr_travellers, version,
                       cum_nr_patches, patch_population,
                       patch_distance, admin_flow_data,
                       patch_radius_population_L, patch_radius_population_S,
                       admin_nr_travellers, admin_population,
                       admin_flow_model, model,
                       output_folder, 0, 0);

    // check if the new log-likelihood is among best computed so far
    if( !isnan(LOGL) && LOGL > chainLOGL[nr_chains - 1] ) {
      chainLOGL[nr_chains - 1] = LOGL;
      for (int p = 0; p < nrPAR; p++)
        (*chainPAR)[nr_chains - 1][p] = PAR[p];
      resort(nr_chains, nrPAR, &chainLOGL, chainPAR);
    }

    // print
    fprintf(file, "\n%d\t%.2LF\t", i+1, LOGL);
    for (int p = 0; p < nrPAR; p++)
      fprintf(file, "%.2LF\t", PAR[p]);
    fflush(file); // force immediate printing
    fflush(file);

  } // end of LHS loop

  fclose(file);

  free(lenPAR);
  for (int p = 0; p < nrPAR; p++)
    free(randPAR[p]);
  free(randPAR);
  free(chainLOGL);

  return 0;
}

// // debug
// for (int j = 0; j < nr_chains; j++) {
//   printf("%LF\t", chainLOGL[j]);
//   for (int p = 0; p < nrPAR; p++)
//     printf("%LF\t", (*chainPAR)[j][p]);
//   printf("\n");
// }
// printf("\n");

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "header.h" // header is in same directory as this file
#ifdef _OPENMP
  #include <omp.h> // for openMP parallelisation
#endif

int *arr_sort;

int main(int argc, char *argv[]) {

	clock_t CPUstart = clock();

  /* ------------------------------------------------------------------------ */
  /* Read command line input                                                  */
  /* ------------------------------------------------------------------------ */
  char *input_folder, *country, *scale, *model, *output_path, *output_suffix;
  int nrPAR, nr_travellers, version; // nr of parameters to fit varies depending on the model, version and method to compute nr of travellers (RM only)
  long double *PAR, *minPAR, *aPAR, *bPAR;
  int do_symmetric, do_skipdiagonal, do_rescalekappa;
  int min_pop, iter_flow_model, iterLHS, iterMCMC, nr_chains, chain_nr, seed;
  long double acceptance_rate, tweaking_factor;

  parser(argc, argv,
         &input_folder, &country, &scale, &model, &output_path, &output_suffix,
         &nrPAR, &nr_travellers, &version, &PAR, &minPAR, &aPAR, &bPAR,
         &do_symmetric, &do_skipdiagonal, &do_rescalekappa,
         &min_pop, &iter_flow_model, &iterLHS, &iterMCMC, &nr_chains, &chain_nr, &seed,
         &acceptance_rate, &tweaking_factor
  );

  char output_folder[500];
  print_parser(input_folder, country, scale, model, output_path, output_suffix,
  			       nrPAR, nr_travellers, version, PAR, minPAR, aPAR, bPAR,
               do_symmetric, do_skipdiagonal, do_rescalekappa,
  			       min_pop, iterLHS, iterMCMC, nr_chains, chain_nr, seed,
  			       acceptance_rate, tweaking_factor,
               output_folder
  );


  /* ------------------------------------------------------------------------ */
  /* Read file input                                                          */
  /* ------------------------------------------------------------------------ */
  int nr_admins, nr_patches;
  int *cum_nr_patches, *patch_population, *admin_population;
  int **patch_distance, **admin_flow_data;
  long double **admin_flow_model;

  read_input(input_folder, scale,
             &nr_admins, &nr_patches,
             &cum_nr_patches, &patch_population, &admin_population,
             &patch_distance, &admin_flow_data,
             &admin_flow_model
  );


  /* ------------------------------------------------------------------------ */
  /* Symmetrise flow data and remove diagonal                                 */
  /* ------------------------------------------------------------------------ */
  if (do_symmetric > 0)
		symmetrise_matrix(admin_flow_data, nr_admins);

  if (do_skipdiagonal > 0)
    remove_diagonal(admin_flow_data, nr_admins);


  /* ------------------------------------------------------------------------ */
  /* Compute/read matrices for radiation model                                */
  /* ------------------------------------------------------------------------ */
  int **patch_radius_population_L=NULL, **patch_radius_population_S=NULL, *admin_nr_travellers=NULL;
  if (strncmp(model, "RM", 3) == 0) {
    prepare_RM(&patch_radius_population_L, &patch_radius_population_S, &admin_nr_travellers,
               nr_patches, nr_admins, patch_population, patch_distance, admin_flow_data,
               input_folder, country, scale);
  }


  /* ------------------------------------------------------------------------ */
  /* Initialise chains                                                        */
  /* ------------------------------------------------------------------------ */
  srand(seed);
  long double **chainPAR = (long double **) malloc(nr_chains * sizeof(long double *));
  for (int i = 0; i < nr_chains; i++) {
    chainPAR[i] = (long double *) malloc(nrPAR * sizeof(long double));
    for (int p = 0; p < nrPAR; p++)
      chainPAR[i][p] = PAR[p];
  }
  long double *PAR_flow = (long double *)malloc(nrPAR * sizeof(long double));


  /* ------------------------------------------------------------------------ */
  /* Print admin_flow_data matrix and exit                                    */
  /* ------------------------------------------------------------------------ */

  if (iter_flow_model > 0) {
    for (int i = 0; i < iter_flow_model; i++) {

			for (int p = 0; p < nrPAR; p++) {
				PAR_flow[p] = PAR[i * nrPAR + p];
      }

      computeLOGL(do_symmetric, do_skipdiagonal, do_rescalekappa,
                  min_pop, PAR_flow,
                  nr_admins, nr_patches, nrPAR, nr_travellers, version,
                  cum_nr_patches, patch_population,
                  patch_distance, admin_flow_data,
                  patch_radius_population_L, patch_radius_population_S,
                  admin_nr_travellers, admin_population,
                  admin_flow_model, model,
                  output_folder, iter_flow_model, i);

    }
  }

  /* ------------------------------------------------------------------------ */
  /* Latin Hypercube Sampler                                                  */
  /* ------------------------------------------------------------------------ */
  if( iterLHS > 0 )
    LHS(do_symmetric, do_skipdiagonal, do_rescalekappa,
        iterLHS, nr_chains, nrPAR, nr_travellers, version,
        aPAR, bPAR, PAR, &chainPAR,
        output_folder, min_pop, nr_admins, nr_patches,
        cum_nr_patches, patch_population, patch_distance, admin_flow_data,
        patch_radius_population_L, patch_radius_population_S,
        admin_nr_travellers, admin_population,
        admin_flow_model, model, 0);


  /* ------------------------------------------------------------------------ */
  /* MCMC                                                                     */
  /* ------------------------------------------------------------------------ */
  if( iterMCMC > 0 ) {
    for (int i = 0; i < nr_chains; i++)
      MCMC(do_symmetric, do_skipdiagonal, do_rescalekappa,
        nrPAR, iterMCMC, nr_travellers, version,
        &chainPAR[i], minPAR, aPAR,
        acceptance_rate, tweaking_factor, output_folder, i,
        min_pop, nr_admins, nr_patches, cum_nr_patches,
        patch_population, patch_distance, admin_flow_data,
        patch_radius_population_L, patch_radius_population_S,
        admin_nr_travellers, admin_population,
        admin_flow_model, model, 0);
  }


  /* ------------------------------------------------------------------------ */
  /* Tidy up                                                                  */
  /* ------------------------------------------------------------------------ */

  // allocated in parser()
  free(PAR);
  free(minPAR);
  free(aPAR);
  free(bPAR);
  free(PAR_flow);

  // allocated in read_input()
  free(patch_population);
  for (int i = 0; i < nr_patches; i++)
    free(patch_distance[i]);
  free(patch_distance);
  free(cum_nr_patches);
  free(admin_population);
  for (int i = 0; i < nr_admins; i++) {
    free(admin_flow_data[i]);
    free(admin_flow_model[i]);
  }
  free(admin_flow_data);
  free(admin_flow_model);

  // allocated in prepare_RM()
  if (strncmp(model, "RM", 3) == 0) {
    for (int i = 0; i < nr_patches; i++) {
      free(patch_radius_population_L[i]);
      free(patch_radius_population_S[i]);
    }
    free(patch_radius_population_L);
    free(patch_radius_population_S);
    free(admin_nr_travellers);
  }

  // allocated in main()
  for (int i = 0; i < nr_chains; i++)
    free(chainPAR[i]);
  free(chainPAR);


  /* ------------------------------------------------------------------------ */
  /* Print time                                                               */
  /* ------------------------------------------------------------------------ */

	double CPUdiff = (clock() - CPUstart) / (double)CLOCKS_PER_SEC;
  print_time(output_folder, CPUdiff);


  return 0;
}

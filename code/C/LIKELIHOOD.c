#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#ifdef _OPENMP
  #include <omp.h> // for openMP parallelisation
#endif
#include "header.h"

long double computeLOGL(int do_symmetric, int do_skipdiagonal, int do_rescalekappa,
                        int min_pop, long double *PAR,
                        int nr_admins, int nr_patches, int nrPAR, int nr_travellers, int version,
                        int *cum_nr_patches, int *patch_population,
                        int **patch_distance, int **admin_flow_data,
                        int **patch_radius_population_L, int **patch_radius_population_S,
                        int *admin_nr_travellers, int *admin_population,
                        long double **admin_flow_model, char *model,
                        char output_folder[], int iter_flow_model, int nr_flow_model
) {
  // Define and initialise variables
  int i=0, big_i, ii;
  int j=0, big_j=0;
  long double temp=0., temp2=0.;

  long double LOGL=0.;
  long double nbPAR = 1./PAR[nrPAR-1]; // negative binomial parameter


  /* GRAVITY MODEL ---------------------------------------------------------- */
  if (strncmp(model, "GM", 3) == 0) {

    // parallelisation
    #ifdef _OPENMP
    #pragma omp parallel for default(shared) private(big_i, big_j, i, j, ii) schedule(dynamic) reduction(+:LOGL)
    #endif

    for (big_i = 0; big_i < nr_admins; big_i++) {

      ii = 0;
      if (do_symmetric > 0)
        ii = big_i;

      for (big_j = ii; big_j < nr_admins; big_j++) {

        if (do_skipdiagonal > 0 && big_i == big_j)
          continue;

        // compute matrix elements of model flow
        admin_flow_model[big_i][big_j] = 0.;
        for (i = cum_nr_patches[big_i]; i < cum_nr_patches[big_i + 1]; i++) {

          if (patch_population[i] < min_pop)
            continue;

          for (j = cum_nr_patches[big_j]; j < cum_nr_patches[big_j + 1]; j++) {

            if (patch_population[j] < min_pop)
              continue;

            // compute element of patch matrix
            if (do_symmetric > 0) {
              admin_flow_model[big_i][big_j] +=
                  ( pow((double) patch_population[i], PAR[1]) * pow((double) patch_population[j], PAR[2])
                  + pow((double) patch_population[i], PAR[2]) * pow((double) patch_population[j], PAR[1]) )
                / ( 1 + pow( (double)patch_distance[i][j] / pow(10., PAR[3]) , PAR[4] ) );
            } else {
              admin_flow_model[big_i][big_j] +=
                  pow((double) patch_population[i], PAR[1]) * pow((double) patch_population[j], PAR[2])
                / ( 1 + pow( (double)patch_distance[i][j] / pow(10., PAR[3]) , PAR[4] ) );
            }

          }
        }

        if (do_rescalekappa > 0) {
          admin_flow_model[big_i][big_j] *= pow( 10., PAR[0] - PAR[3] * PAR[4] );
        } else {
          admin_flow_model[big_i][big_j] *= pow( 10., PAR[0] );
        }

        // compute log-likelihood
        LOGL +=
            lgammal((long double) admin_flow_data[big_i][big_j] + nbPAR)
          - lgammal((long double) admin_flow_data[big_i][big_j] + 1)
          - lgammal(nbPAR)
          + (long double) admin_flow_data[big_i][big_j]
            * logl( admin_flow_model[big_i][big_j] / (admin_flow_model[big_i][big_j] + nbPAR) )
          + nbPAR * logl( nbPAR / (admin_flow_model[big_i][big_j] + nbPAR) );

      }
    }

  } // END OF GRAVITY MODEL

  /* RADIATION MODEL -------------------------------------------------------- */
  else {

    // parallelisation
    #ifdef _OPENMP
    #pragma omp parallel for default(shared) private(big_i, big_j, i, j, ii, temp, temp2) schedule(dynamic) reduction(+:LOGL)
    #endif

    for (big_i = 0; big_i < nr_admins; big_i++) {

      ii = 0;
      if (do_symmetric > 0)
        ii = big_i;

      for (big_j = ii; big_j < nr_admins; big_j++) {

        // skip within-admin travel for radiation models
        if (big_i == big_j)
          continue;

        // compute matrix elements of model flow
        admin_flow_model[big_i][big_j] = 0.;
        for (i = cum_nr_patches[big_i]; i < cum_nr_patches[big_i + 1]; i++) {

          if (patch_population[i] < min_pop)
            continue;

          for (j = cum_nr_patches[big_j]; j < cum_nr_patches[big_j + 1]; j++) {

            if (patch_population[j] < min_pop)
              continue;

            // compute element of patch matrix
            temp2 = 0;
            if (version == 0) { // Wesolowski 2015 Plos Computational Biology
              temp =
                  ( (long double) patch_population[i] * (long double) patch_population[j] )
                / ( (long double) patch_radius_population_L[i][j] * (long double) patch_radius_population_S[i][j] );

              if (do_symmetric > 0)
                temp2 =
                    ( (long double) patch_population[j] * (long double) patch_population[i] )
                  / ( (long double) patch_radius_population_L[j][i] * (long double) patch_radius_population_S[j][i] );

            }
            else if (version == 1) { // Marshall 2018 Scientific Reports I
              temp =
                  ( PAR[1] * (long double) patch_population[i] * (long double) patch_population[j] )
                / ( ((long double) patch_radius_population_L[i][j] + (long double) patch_population[i] * (PAR[1] - 1.)) * ((long double) patch_radius_population_S[i][j] + (long double) patch_population[i] * (PAR[1] - 1.)) );

              if (do_symmetric > 0)
                temp2 =
                  ( PAR[1] * (long double) patch_population[j] * (long double) patch_population[i] )
                / ( ((long double) patch_radius_population_L[j][i] + (long double) patch_population[j] * (PAR[1] - 1.)) * ((long double) patch_radius_population_S[j][i] + (long double) patch_population[j] * (PAR[1] - 1.)) );

            }
            else { // Marshall 2018 Scientific Reports II
              temp =
                  ( ((long double) patch_population[i] + PAR[1]) * (long double) patch_population[j] )
                / ( ((long double) patch_radius_population_L[i][j] + PAR[1]) * ((long double) patch_radius_population_S[i][j] + PAR[1]) );

              if (do_symmetric > 0)
                temp2 =
                  ( ((long double) patch_population[j] + PAR[1]) *  (long double) patch_population[i] )
                / ( ((long double) patch_radius_population_L[j][i] + PAR[1]) * ((long double) patch_radius_population_S[j][i] + PAR[1]) );

            }


            // fit nr_travellers
            switch (nr_travellers) {
              case 0:
                admin_flow_model[big_i][big_j] += pow(10., PAR[0])                                                         * (temp + temp2);
                break;
              case 1:
                admin_flow_model[big_i][big_j] += pow(10., PAR[0]) * ( (long double) admin_nr_travellers[big_i]            * temp + (long double) admin_nr_travellers[big_j]            * temp2 );
                break;
              case 2:
                admin_flow_model[big_i][big_j] += pow(10., PAR[0]) * ( pow((double) admin_population[big_i], PAR[nrPAR-2]) * temp + pow((double) admin_population[big_j], PAR[nrPAR-2]) * temp2 );
                break;
              case 3:
                admin_flow_model[big_i][big_j] += pow(10., PAR[0]) * ( pow((double) patch_population[i],     PAR[nrPAR-2]) * temp + pow((double) patch_population[j],     PAR[nrPAR-2]) * temp2 );
                break;
              case 4:
                admin_flow_model[big_i][big_j] += pow(10., PAR[0]) * ( (long double) patch_population[i]                   * temp + (long double) patch_population[j]                   * temp2 );
                break;
              case 5:
                admin_flow_model[big_i][big_j] += PAR[0]           * ( (long double) patch_population[i]                   * temp + (long double) patch_population[j]                   * temp2 );
                break;
              default:
                fprintf(stderr, "Error in switch statement in log-likelihood computation!\n");
                exit(1);
            }

          } // end of compute element of patch matrix

        }

        // compute log-likelihood
        LOGL +=
            lgammal((long double) admin_flow_data[big_i][big_j] + nbPAR)
          - lgammal((long double) admin_flow_data[big_i][big_j] + 1)
          - lgammal(nbPAR)
          + (long double) admin_flow_data[big_i][big_j]
            * logl( admin_flow_model[big_i][big_j] / (admin_flow_model[big_i][big_j] + nbPAR) )
          + nbPAR * logl( nbPAR / (admin_flow_model[big_i][big_j] + nbPAR) );

      }
    }

  } // END OF RADIATION MODEL


  if (do_symmetric > 0)
    LOGL *= 2.;


  // print admin_flow_model to file
  if (iter_flow_model > 0)
    fprint_flow_model(LOGL, PAR, admin_flow_model,
                      nrPAR, nr_flow_model, iter_flow_model, nr_admins, do_symmetric, output_folder);

  // printf("%LF\n", admin_flow_model[nr_admins-2][nr_admins-1] );
  // printf("%LF\n", admin_flow_model[nr_admins-1][nr_admins-2] );
  // printf("%LF\n", LOGL );

  // exit(1);

  return LOGL;
}

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <sys/stat.h> // in order to use mkdir()
#if defined(_WIN32)
	#include <direct.h> // in order to use _mkdir on Windows
#endif
#include "header.h"

extern int *arr_sort;

static int compar(
  const void *a, const void *b
) {
  int aa = *((int *) a), bb = *((int *) b);

  if (arr_sort[aa] < arr_sort[bb])
    return -1;
  else if (arr_sort[aa] == arr_sort[bb])
    return 0;
  else if (arr_sort[aa] > arr_sort[bb])
    return 1;
  else {
    fprintf(stderr, "Error: comparison did not resolve order.\n");
    exit(1);
  }

  return 0;
}


int prepare_RM(
  int ***patch_radius_population_L, int ***patch_radius_population_S, int **admin_nr_travellers,
  int nr_patches, int nr_admins, int *patch_population, int **patch_distance, int **admin_flow_data,
  char *input_folder, char *country, char *scale
) {

  // for each origin i and destination j, compute the population contained in a radius r_ij around i
  *patch_radius_population_L = calloc(nr_patches, sizeof(int*));
  *patch_radius_population_S = calloc(nr_patches, sizeof(int*));
  for (int i = 0; i < nr_patches; i++) {
    (*patch_radius_population_L)[i] = calloc(nr_patches, sizeof(int));
    (*patch_radius_population_S)[i] = calloc(nr_patches, sizeof(int));
  }

  /* If existing, read in variables from files, else create and save to files */
  char str[500];
  sprintf(str, "%s/%s__patch_radius_population_L.tsv", input_folder, scale); // win: \\, mac /
  struct stat s;
  if (stat(str, &s) == 0)
    read_RM_variables(patch_radius_population_L, patch_radius_population_S, admin_nr_travellers,
                      nr_patches, nr_admins, patch_population, patch_distance, admin_flow_data,
                      input_folder, country, scale);
  else
    create_RM_variables(patch_radius_population_L, patch_radius_population_S, admin_nr_travellers,
                        nr_patches, nr_admins, patch_population, patch_distance, admin_flow_data,
                        input_folder, country, scale);


  return 0;
}


int read_RM_variables(
  int ***patch_radius_population_L, int ***patch_radius_population_S, int **admin_nr_travellers,
  int nr_patches, int nr_admins, int *patch_population, int **patch_distance, int **admin_flow_data,
  char *input_folder, char *country, char *scale
) {
  char str[300];
  FILE *file;


  /* Read radius population (L) --------------------------------------------- */
  sprintf(str, "%s/%s__patch_radius_population_L.tsv", input_folder, scale);
  if ((file = fopen(str, "r")) == NULL) {
    fprintf(stderr, "Error: file \"patch_radius_population_L.tsv\" not found.\n");
    exit(1);
  }
  for (int i = 0; i < nr_patches; i++) {
    for (int j = 0; j < nr_patches; j++) {
      if (!fscanf(file, "%d", &(*patch_radius_population_L)[i][j])) {
        fprintf(stderr, "Error reading in file \"patch_radius_population_L.tsv\"");
        exit(1);
      }
    }
  }
  fclose(file);


  /* Read radius population (S) --------------------------------------------- */
  sprintf(str, "%s/%s__patch_radius_population_S.tsv", input_folder, scale);
  if ((file = fopen(str, "r")) == NULL) {
    fprintf(stderr, "Error: file \"patch_radius_population_S.tsv\" not found.\n");
    exit(1);
  }
  for (int i = 0; i < nr_patches; i++) {
    for (int j = 0; j < nr_patches; j++) {
      if (!fscanf(file, "%d", &(*patch_radius_population_S)[i][j])) {
        fprintf(stderr, "Error reading in file \"patch_radius_population_S.tsv\"");
        exit(1);
      }
    }
  }
  fclose(file);


  /* Namibia: read estimated number of residents who travel ----------------- */
  if (strncmp(country, "Namibia", 8) == 0) {
    sprintf(str, "%s/%s__admin_nr_travellers.tsv", input_folder, scale);
    if ((file = fopen(str, "r")) == NULL) {
      fprintf(stderr, "Error: file \"admin_nr_travellers.tsv\" not found.\n");
      exit(1);
    }
		*admin_nr_travellers = calloc(nr_admins, sizeof(int));
    for (int i = 0; i < nr_admins; i++) {
      if (!fscanf(file, "%d", &(*admin_nr_travellers)[i])) {
        fprintf(stderr, "Error reading in file \"admin_nr_travellers.tsv\"");
        exit(1);
      }
    }
    fclose(file);
  }

  return 0;
}


int create_RM_variables(
  int ***patch_radius_population_L, int ***patch_radius_population_S, int **admin_nr_travellers,
  int nr_patches, int nr_admins, int *patch_population, int **patch_distance, int **admin_flow_data,
  char *input_folder, char *country, char *scale
) {

	/* Define and initialise variables ---------------------------------------- */
	// matrix of indices: for each origin, sort the indices by distance in ascending order
	int **patch_idx_mat = malloc(nr_patches * sizeof(int*));
	for (int i = 0; i < nr_patches; i++) {
		patch_idx_mat[i] = malloc(nr_patches * sizeof(int));
		for (int j = 0; j < nr_patches; j++)
			patch_idx_mat[i][j] = j;
	}

  /* Find sorting order of admin_flow_data matrix rows ---------------------- */
  for (int i = 0; i < nr_patches; i++) {
    arr_sort = patch_distance[i];
    qsort (patch_idx_mat[i], nr_patches, sizeof (int), compar);
  }

  /* Compute radius population ---------------------------------------------- */
  int temp;
  for (int i = 0; i < nr_patches; i++) {
    temp = 0;
    for (int j = 0; j < nr_patches; j++) {
      (*patch_radius_population_S)[i][patch_idx_mat[i][j]] = temp;
      temp += patch_population[ patch_idx_mat[i][j] ];
      (*patch_radius_population_L)[i][patch_idx_mat[i][j]] = temp;
    }
  }

  // verify (*patch_radius_population_L)[ : ][1] == patch_population[ : ]
  for (int i = 0; i < nr_patches; i++) {
    if ( (*patch_radius_population_L)[i][i] != patch_population[i] ) {
      fprintf(stderr, "Error: patch population in smallest radius does not correspond to origin population.\n");
      exit(1);
    }
  }


  /* Write matrices to file ------------------------------------------------- */
  // create files
  char str[500];
  FILE *file_radius_pop_L, *file_radius_pop_S, *file_trav;

  sprintf(str, "%s/%s__patch_radius_population_L.tsv", input_folder, scale);
  if ( (file_radius_pop_L = fopen(str, "w")) == 0 ) {
    fprintf(stderr, "Error: file \"patch_radius_population_L.tsv\" could not be opened.\n");
    exit(1);
  }
  sprintf(str, "%s/%s__patch_radius_population_S.tsv", input_folder, scale);
  if ( (file_radius_pop_S = fopen(str, "w")) == 0 ) {
    fprintf(stderr, "Error: file \"patch_radius_population_S.tsv\" could not be opened.\n");
    exit(1);
  }

  // write matrices to files
  for (int i = 0; i < nr_patches; i++) {
    for (int j = 0; j < nr_patches; j++) {
      fprintf(file_radius_pop_L, "%d\t", (*patch_radius_population_L)[i][j]);
      fprintf(file_radius_pop_S, "%d\t", (*patch_radius_population_S)[i][j]);
    }
    fprintf(file_radius_pop_L, "\n");
    fprintf(file_radius_pop_S, "\n");
  }

  fclose(file_radius_pop_L);
  fclose(file_radius_pop_S);
	for (int i = 0; i < nr_patches; i++)
		free(patch_idx_mat[i]);
  free(patch_idx_mat);

  /* Namibia: estimate number of residents who travel from nr of trips ------ */
  if (strncmp(country, "Namibia", 8) == 0) {
    // initialise array of nr of travellers per admin
    *admin_nr_travellers = calloc(nr_admins, sizeof(int));
    for (int i = 0; i < nr_admins; i++) {
      for (int j = 0; j < nr_admins; j++) {

        // skip within-admin travel
        if( j == i )
          continue;
        (*admin_nr_travellers)[i] += admin_flow_data[i][j];

      }
    }

    // create file
    sprintf(str, "%s/%s__admin_nr_travellers.tsv", input_folder, scale);
    if ( (file_trav = fopen(str, "w")) == 0 ) {
      fprintf(stderr, "Error: file \"admin_nr_travellers.tsv\" could not be opened.\n");
      exit(1);
    }

    // write array to file
    for (int i = 0; i < nr_admins; i++)
      fprintf(file_trav, "%d\n", (*admin_nr_travellers)[i]);

    fclose(file_trav);
  }

  return 0;
}

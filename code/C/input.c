#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "header.h"

int read_input(char *input_folder, char *scale,
  				     int *nr_admins, int *nr_patches,
  			 	     int **cum_nr_patches, int **patch_population, int **admin_population,
  			 	     int ***patch_distance, int ***admin_flow_data,
  			 	     long double ***admin_flow_model
) {
  char str[300], var_name[50], var_value[50];
  FILE *file;

  /* Read nr of patches and admins ------------------------------------------ */
  sprintf(str, "%s/%s__metadata.tsv", input_folder, scale); // win: \\, mac /
  if ((file = fopen(str, "r")) == NULL) {
    fprintf(stderr, "Error: file \"%s__metadata.tsv\" not found.\n", scale);
    exit(1);
  }
  while (!feof(file)) {
    if (fgets(str, 300, file)) {
      sscanf(str, "%s\t%s", var_name, var_value);

      if (strncmp(var_name, "nr_patches", 15) == 0) {
        *nr_patches = get_long(var_value, var_name);
      }
      else if (strncmp(var_name, "nr_admins", 15) == 0) {
        *nr_admins = get_long(var_value, var_name);
      }
      else {
        fprintf(stderr, "Unrecongized parameter name: \"%s\"\n", var_name);
        exit(1);
      }

    }
  }
  fclose(file);

  /* Allocate memory to dynamic arrays--------------------------------------- */
  *patch_population        = (int          *) malloc((*nr_patches) * sizeof(int));
  *patch_distance          = (int         **) malloc((*nr_patches) * sizeof(int *));
  for (int i = 0; i < (*nr_patches); i++) {
    (*patch_distance)[i]   = (int          *) malloc((*nr_patches) * sizeof(int));
  }

  *cum_nr_patches          = (int          *) malloc((*nr_admins + 1) * sizeof(int));
  *admin_population        = (int          *) malloc((*nr_admins) * sizeof(int));
  *admin_flow_data         = (int         **) malloc((*nr_admins) * sizeof(int *));
  *admin_flow_model        = (long double **) malloc((*nr_admins) * sizeof(long double *));
  for (int i = 0; i < (*nr_admins); i++) {
    (*admin_flow_data)[i]  = (int          *) malloc((*nr_admins) * sizeof(int));
    (*admin_flow_model)[i] = (long double  *) malloc((*nr_admins) * sizeof(long double));
  }

  /* Read patch population -------------------------------------------------- */
  sprintf(str, "%s/%s__patch_population.tsv", input_folder, scale);
  if ((file = fopen(str, "r")) == NULL) {
    fprintf(stderr, "Error: file \"%s__patch_population.tsv\" not found.\n", scale);
    exit(1);
  }
  for (int i = 0; i < (*nr_patches); i++) {
    if (!fscanf(file, "%d", &(*patch_population)[i])) {
      fprintf(stderr, "Error reading in file \"%s__patch_population.tsv\"", scale);
      exit(1);
    }
  }
  fclose(file);

  /* Read **admin** population ---------------------------------------------- */
  sprintf(str, "%s/%s__admin_population.tsv", input_folder, scale);
  if ((file = fopen(str, "r")) == NULL) {
    fprintf(stderr, "Error: file \"%s__admin_population.tsv\" not found.\n", scale);
    exit(1);
  }
  for (int i = 0; i < (*nr_admins); i++) {
    if (!fscanf(file, "%d", &(*admin_population)[i])) {
      fprintf(stderr, "Error reading in file \"%s__admin_population.tsv\"", scale);
      exit(1);
    }
  }
  fclose(file);

  /* Read patch distance (convert upper triangular to symmetric matrix) ----- */
  sprintf(str, "%s/%s__patch_distance.tsv", input_folder, scale);
  if ((file = fopen(str, "r")) == 0) {
    fprintf(stderr, "Error: file \"%s__patch_distance.tsv\" not found.\n", scale);
    exit(1);
  }
  for (int i = 0; i < (*nr_patches); i++) {
    (*patch_distance)[i][i] = 0;
    for (int j = i+1; j < (*nr_patches); j++) {
      if (!fscanf(file, "%d", &(*patch_distance)[i][j])) {
        fprintf(stderr, "Error reading in file \"%s__patch_distance.tsv\"", scale);
        exit(1);
      }
      (*patch_distance)[j][i] = (*patch_distance)[i][j];
    }
  }
  fclose(file);

  /* Read cumulative number of patches per admin ---------------------------- */
  sprintf(str, "%s/%s__patches_per_admin.tsv", input_folder, scale);
  if ((file = fopen(str, "r")) == NULL) {
    fprintf(stderr, "Error: file \"%s__patches_per_admin.tsv\" not found.\n", scale);
    exit(1);
  }
  int *temp = (int *) malloc((*nr_admins) * sizeof(int));
  for (int i = 0; i < (*nr_admins); i++) {
    if (!fscanf(file, "%d", &temp[i])) {
      fprintf(stderr, "Error reading in file \"%s__patches_per_admin.tsv\"", scale);
      exit(1);
    }
  }
  fclose(file);
  // convert to cumulative array - necessary for parallelisation
  (*cum_nr_patches)[0] = 0;
  long int cum = 0;
  for (int i = 0; i < (*nr_admins); i++) {
    cum += temp[i];
    (*cum_nr_patches)[i + 1] = cum;
  }
  if (cum != (*nr_patches)) {
    fprintf(stderr, "Error in read_input: sum does not correspond to number of patchs\n");
    exit(1);
  }
  free(temp);

  /* Read flow data between admins ------------------------------------------ */
  sprintf(str, "%s/%s__admin_flow_data.tsv", input_folder, scale);
  if ((file = fopen(str, "r")) == 0) {
    fprintf(stderr, "Error: file \"%s__flow_data.tsv\" not found.\n", scale);
    exit(1);
  }
  for (int i = 0; i < (*nr_admins); i++) {
    for (int j = 0; j < (*nr_admins); j++) {
      if (!fscanf(file, "%d", &(*admin_flow_data)[i][j])) {
        fprintf(stderr, "Error reading in file \"%s__flow_data.tsv\"", scale);
        exit(1);
      }
    }
  }
  fclose(file);

  return 0;
}

// // debug
// printf("%d %d\n", nr_patches, nr_admins);
// printf("patch_population [%d, %d]\n",   (*patch_population)[0],    (*patch_population)[*nr_patches-1]);
// printf("patch_distance   [%d, %d]\n",   (*patch_distance)  [0][1], (*patch_distance)  [*nr_patches-2][*nr_patches-1]);
// printf("cum_nr_patches   [%d, %d]\n",   (*cum_nr_patches)  [0],    (*cum_nr_patches)  [*nr_admins]);
// printf("admin_population [%d, %d]\n",   (*admin_population)[0],    (*admin_population)[*nr_admins-1]);
// printf("admin_flow_data  [%d, %d]\n",   (*admin_flow_data) [0][1], (*admin_flow_data) [*nr_admins-2][*nr_admins-1]);
// printf("admin_flow_model [%LF, %LF]\n", (*admin_flow_model)[0][0], (*admin_flow_model)[*nr_admins-1][*nr_admins-1]);

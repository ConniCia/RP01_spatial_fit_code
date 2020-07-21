#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "header.h"

/* -------------------------------------------------------------------------- */
/* Transform string to long int (and catch errors)                            */
/* -------------------------------------------------------------------------- */
long int get_long(char *str_value, char *str_flag) {
	char *end;
	long int res = strtol(str_value, &end, 10);

	if ((*end != '\0') && (end != NULL) && (strncmp(end, "", 2))) {
		fprintf(stderr, "error: string contains additional non-numeric characters\n");
		fprintf(stderr, "strtol() could not convert \"%s\" for parameter \"%s\"\n", end, str_flag);
		exit(1);
	}
	if (end == str_value) {
		fprintf(stderr, "error: no conversion has taken place due to wrong format\n");
		fprintf(stderr, "strtol() could not convert \"%s\" for parameter \"%s\"\n", end, str_flag);
		exit(1);
	}

	return res;
}


/* -------------------------------------------------------------------------- */
/* Transform string to long double (and catch errors)                         */
/* -------------------------------------------------------------------------- */
long double get_longdouble(char *str_value, char *str_flag) {
	char *end;
	long double res = strtold(str_value, &end);

	if ((*end != '\0') && (end != NULL) && (strncmp(end, "", 3))) {
		fprintf(stderr, "error: string contains additional non-numeric characters\n");
		fprintf(stderr, "strtol() could not convert \"%s\" for parameter \"%s\"\n", end, str_flag);
		exit(1);
	}
	if (end == str_value) {
		fprintf(stderr, "error: no conversion has taken place due to wrong format\n");
		fprintf(stderr, "strtol() could not convert \"%s\" for parameter \"%s\"\n", end, str_flag);
		exit(1);
	}

	return res;
}


/* -------------------------------------------------------------------------- */
/* Symmetrise matrix by summing its transpose                                 */
/* -------------------------------------------------------------------------- */
int symmetrise_matrix(int **M, long int n) {
	long int i, j;

	// make M symmetric (by summing M and M^t)
	for (i = 0; i < n; i++) {
		M[i][i] += M[i][i];
		for (j = i + 1; j < n; j++) { // using < we leave out the diagonal
			M[i][j] += M[j][i];         // take care of the upper diagonal
			M[j][i] = M[i][j];          // then of the lower diagonal
		}
	}

	return 0;
} // end of make_symmetric() // make input matrix symmetric by adding its transpose


/* -------------------------------------------------------------------------- */
/* Remove diagonal from flow matrix                                           */
/* -------------------------------------------------------------------------- */
int remove_diagonal(int **M, long int n) {
	long int i;

	// substitute zero on the diagonal
	for (i = 0; i < n; i++) {
		M[i][i] = 0;
	}

	return 0;
} // end of remove_diagonal() // Remove diagonal from flow matrix


/* -------------------------------------------------------------------------- */
/* Randomise order of array elements                                          */
/* -------------------------------------------------------------------------- */
int randomise(int **arr, int n) {

	int temp;
	// Start from the last element and swap one by one
	for (int i = n - 1; i > 0; i--)
	{
		// Pick a random index from 0 to i
		int j = rand() % (i + 1);
		// Swap arr[i] with arr[j]
		temp = (*arr)[i];
		(*arr)[i] = (*arr)[j];
		(*arr)[j] = temp;
	}

	return 0;
}

/* -------------------------------------------------------------------------- */
/* Resort chainPAR according to chainLOGL values                              */
/* -------------------------------------------------------------------------- */
int resort(int nr_chains, int nrPAR, long double **chainLOGL, long double ***chainPAR) {

	long double xLOGL = (*chainLOGL)[nr_chains - 1];
	long double *xPAR = (long double *) malloc(nrPAR * sizeof(long double));
	for (int i = 0; i < nrPAR; i++)
	  xPAR[i] = (*chainPAR)[nr_chains - 1][i];

	for (int i = nr_chains - 2; i >= 0; i--) {

	  if (xLOGL > (*chainLOGL)[i]) {
	    (*chainLOGL)[i+1] = (*chainLOGL)[i];
	    for (int j = 0; j < nrPAR; j++) (*chainPAR)[i+1][j] = (*chainPAR)[i][j];
	    if (i == 0) {
	      (*chainLOGL)[0] = xLOGL;
	      for (int j = 0; j < nrPAR; j++) (*chainPAR)[0][j] = xPAR[j];
	    }
	  }
	  else {
	    (*chainLOGL)[i+1] = xLOGL;
	    for (int j = 0; j < nrPAR; j++) (*chainPAR)[i+1][j] = xPAR[j];
	    break;
	  }
	}
	free(xPAR);

	return 0;
} // end of newbestorder() // sort bestLOGL and the relative sets of parameters

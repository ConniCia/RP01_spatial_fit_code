#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h> // in order to use mkdir()
#if defined(_WIN32)
	#include <direct.h> // in order to use _mkdir on Windows
#endif
#include "header.h"


int print_parser(char *input_folder, char *country, char *scale, char *model, char *output_path, char *output_suffix,
       					 int nrPAR, int nr_travellers, int version, long double *PAR, long double *minPAR, long double *aPAR, long double *bPAR,
								 int do_symmetric, int do_skipdiagonal, int do_rescalekappa,
			 				 	 int min_pop, int iterLHS, int iterMCMC, int nr_chains, int chain_nr, int seed,
			 				 	 long double acceptance_rate, long double tweaking_factor,
								 char output_folder[]
) {

  // create output folder if not already existing
  sprintf(output_folder, "%s/output_files_%s", output_path, output_suffix);

	struct stat s;
	if (stat(output_folder, &s) != 0) {
#if defined(_WIN32)
		_mkdir(output_folder);
#else
		mkdir(output_folder, S_IRWXU);
#endif
	}

  // create file
	char str[200];
	FILE *file_parser;
	sprintf(str, "%s/parser.tsv", output_folder);
	if ( (file_parser = fopen(str, "w")) == 0 ) {
		fprintf(stderr, "Error: file \"parser\" could not be opened.\n");
		exit(1);
	}

	// print to file
	fprintf(file_parser, "input_folder\t%s\n",     input_folder);
	fprintf(file_parser, "country\t\t%s\n",        country);
	fprintf(file_parser, "scale\t\t%s\n",          scale);
  fprintf(file_parser, "model\t\t%s\n",          model);
  fprintf(file_parser, "output_path\t%s\n",      output_path);
  fprintf(file_parser, "output_suffix\t%s\n",    output_suffix);
  fprintf(file_parser, "\nmin_pop\t\t%d\n",      min_pop);
	fprintf(file_parser, "iterLHS\t\t%d\n",        iterLHS);
	fprintf(file_parser, "iterMCMC\t%d\n",         iterMCMC);
  fprintf(file_parser, "nrPAR\t\t%d\n",          nrPAR);
	fprintf(file_parser, "nr_travellers\t%d\n",    nr_travellers);
	fprintf(file_parser, "version\t\t%d\n",        version);
	fprintf(file_parser, "PAR\t");
  for (int p = 0; p < nrPAR; p++)
    fprintf(file_parser, "\t%Lf",                PAR[p]);
  fprintf(file_parser, "\nminPAR\t");
  for (int p = 0; p < nrPAR; p++)
    fprintf(file_parser, "\t%Lf",                minPAR[p]);
  fprintf(file_parser, "\naPAR\t");
  for (int p = 0; p < nrPAR; p++)
    fprintf(file_parser, "\t%Lf",                aPAR[p]);
  fprintf(file_parser, "\nbPAR\t");
  for (int p = 0; p < nrPAR; p++)
    fprintf(file_parser, "\t%Lf",                bPAR[p]);
	fprintf(file_parser, "\ndo_symmetric\t%d\n",   do_symmetric);
	fprintf(file_parser, "do_skipdiagonal\t%d\n",  do_skipdiagonal);
	fprintf(file_parser, "do_rescalekappa\t%d\n",  do_rescalekappa);
	fprintf(file_parser, "nr_chains\t%d\n",        nr_chains);
	fprintf(file_parser, "chain_nr\t%d\n",         chain_nr);
	fprintf(file_parser, "seed\t\t%u\n",           seed);
	fprintf(file_parser, "acceptance_rate\t%Lf\n", acceptance_rate);
	fprintf(file_parser, "tweaking_factor\t%Lf\n", tweaking_factor);

  // close file
  fclose(file_parser);

  return 0;
}


int print_time(char output_folder[], double CPUdiff) {

	// create file
	char str[200];
	FILE *file_time;
	sprintf(str, "%s/time.tsv", output_folder);
	if ( (file_time = fopen(str, "w")) == 0 ) {
		fprintf(stderr, "Error: file \"time\" could not be opened.\n");
		exit(1);
	}

	if (CPUdiff < 60) {
    fprintf(file_time, "\nCPU time:\t%.2Fs\n", CPUdiff);
  }
  else if (CPUdiff < 3600) {
    fprintf(file_time, "\nCPU time:\t%.2Fm\n", CPUdiff / 60.);
  }
  else {
    fprintf(file_time, "\nCPU time:\t%.2Fh\n", CPUdiff / 3600.);
  }

	// close file
	fclose(file_time);

	return 0;
}


int fprint_flow_model (long double LOGL, long double *PAR, long double **admin_flow_model,
											 int nrPAR, int nr_flow_model, int iter_flow_model, int nr_admins, int do_symmetric,
											 char output_folder[]) {

  char str1[200], str2[200];
	FILE *file_LOGL, *file_matrix;
	sprintf(str1, "%s/flow_model_LOGL.tsv", output_folder);
	sprintf(str2, "%s/flow_model_matrix.tsv", output_folder);

	/* create files and headers ----------------------------------------------- */
	if (nr_flow_model < 0.5) {

		if ( (file_LOGL = fopen(str1, "w")) == 0 ) {
			fprintf(stderr, "Error: file \"flow_model_LOGL\" could not be opened.\n");
			exit(1);
		}
		fprintf(file_LOGL, "id_sample\tLOGL");
		for (int p = 0; p < nrPAR; p++)
	    fprintf(file_LOGL, "\tPAR%d", p);
		fprintf(file_LOGL, "\n");

		if ( (file_matrix = fopen(str2, "w")) == 0 ) {
			fprintf(stderr, "Error: file \"flow_model_matrix\" could not be opened.\n");
			exit(1);
		}
		fprintf(file_matrix, "id_sample\tid_patch_orig\tid_patch_dest\tflow_estimate\n");

	} else {
		file_LOGL   = fopen(str1, "a");
		file_matrix = fopen(str2, "a");
	}

	/* print to file ---------------------------------------------------------- */
	fprintf(file_LOGL, "%d\t%.2LF", nr_flow_model, LOGL);
	for (int p = 0; p < nrPAR; p++)
		fprintf(file_LOGL, "\t%.2LF", PAR[p]);
  fprintf(file_LOGL, "\n");
	fflush(file_LOGL); // force immediate printing
	fflush(file_LOGL);

	if (do_symmetric > 0) {
		for (int i = 0; i < nr_admins; i++) {
			for (int j = 0; j < nr_admins; j++) {
				if (j<i)
					fprintf(file_matrix, "%d\t%d\t%d\t%.2LF\n", nr_flow_model, i, j, admin_flow_model[j][i]);
				else
					fprintf(file_matrix, "%d\t%d\t%d\t%.2LF\n", nr_flow_model, i, j, admin_flow_model[i][j]);
			}
		}
	} else {
		for (int i = 0; i < nr_admins; i++) {
			for (int j = 0; j < nr_admins; j++)
				fprintf(file_matrix, "%d\t%d\t%d\t%.2LF\n", nr_flow_model, i, j, admin_flow_model[i][j]);
		}
	}
	fflush(file_matrix); // force immediate printing
	fflush(file_matrix);

	/* close file ------------------------------------------------------------- */
	if (nr_flow_model > iter_flow_model - 1.5) {
		fclose(file_LOGL);
		fclose(file_matrix);
	}

	return(0);
}

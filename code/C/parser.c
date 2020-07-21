#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "header.h"

int parser(int n, char *arr[],
	     		 char **input_folder, char **country, char **scale, char **model, char **output_path, char **output_suffix,
			 	 	 int *nrPAR, int *nr_travellers, int *version, long double **PAR, long double **minPAR, long double **aPAR, long double **bPAR,
					 int *do_symmetric, int *do_skipdiagonal, int *do_rescalekappa,
			 	 	 int *min_pop, int *iter_flow_model, int *iterLHS, int *iterMCMC, int *nr_chains, int *chain_nr, int *seed,
			 	 	 long double *acceptance_rate, long double *tweaking_factor
) {

  // help message
	char help[] = "Invoke programme \"spatialfit\" as (order not relevant as long as model is given before *PAR):\n\
	./*exec*\n\
	-input_folder      [full path to input folder]\n\
	-country           [country]\n\
	-scale             [scale]\n\
	-model             [GM or RM (if RM also: 0-2 version / 0-5 nr_travellers )]\n\
	-min_pop           [minimum patch population (else skip patch)]\n\
	-iter_flow_model   [number of flow matrices to print to file]\n\
	-iterLHS           [number of LHS iterations]\n\
	-iterMCMC          [number of MCMC iterations]\n\
	-PAR               [GM: kappa alpha beta gamma epsilon neg_bin, RM: par neg_bin]\n\
	-minPAR            [min step size for MCMC]\n\
	-aPAR              [left boundary of LHS interval]\n\
	-bPAR              [right boundary of LHS interval]\n\
	-do_symmetric      [symmetrise input flow matrix]\n\
	-do_skipdiagonal   [skip diagonal in input flow matrix]\n\
	-do_rescalekappa   [for GM: if TRUE use 10^{kappa - gamma * epsilon} instead of 10^kappa]\n\
	\n\
	-output_path       [path to output folder: .]\n\
	-output_suffix     [suffix for output folder: ""]\n\
	-nr_chains         [number of MCMC chains: 1]\n\
	-chain_nr          [chain_nr: 0]\n\
	-seed              [seed: 0]\n\
	-acceptance_rate   [acceptance rate: 0.25]\n\
	-tweaking_factor   [tweaking factor: 2]\n";

	// first group of flags is essential: no default values are given for these
	int ess_flags = 15, temp;

	// assign default values for other flags (or variables that are only used for RM)
	*nr_travellers     = 0;
  *version           = 0;
	*output_path       = ".";
	*output_suffix     = "";
	*nr_chains         = 1;
	*chain_nr          = 0;
	*seed              = 0;
	*acceptance_rate   = 0.25;
	*tweaking_factor   = 2.;


	// show help message if user asks for help
	if (strncmp(arr[1], "-h", 3) == 0 || strncmp(arr[1], "--help", 7) == 0 || strncmp(arr[1], "-help", 6) == 0) {
			fprintf(stderr, "%s\n", help);
			exit(1);
	}

	// read flag values and verify all essential flags are provided
	int j=1, flag_len=20, verify=0;
	for (int i = 1; i < n; i++) {
		j=i;

		// parameters whose values have to be passed from the command line
		if (strncmp(arr[i], "-input_folder", flag_len) == 0) {
			*input_folder = arr[++i];
			verify++;
		}
		else if (strncmp(arr[i], "-country", flag_len) == 0) {
			*country = arr[++i];
			verify++;
		}
		else if (strncmp(arr[i], "-scale", flag_len) == 0) {
			*scale = arr[++i];
			verify++;
		}
		else if (strncmp(arr[i], "-model", flag_len) == 0) {
			*model = arr[++i];
			verify++;

			// nr of parameters varies depending on model; define parameter arrays
			*nrPAR = 6;                          // model == "GM"
			if( strncmp(*model, "RM", 3) == 0 ) {// model == "RM"
				*version = get_long(arr[++i], arr[j]);
				*nr_travellers = get_long(arr[++i], arr[j]);

				if (*version == 0) { // Wesolowski 2015 Plos Computational Biology
					switch (*nr_travellers) {
						case 0:
						case 1:
						case 4:
						case 5:
							*nrPAR = 2;
							break;
						case 2:
						case 3:
							*nrPAR = 3;
							break;
						default:
							fprintf(stderr, "Unrecognized method to compute nr of travellers!\n");
							exit(1);
					}
				}
				else { // Marshall 2018 Scientific Reports
					switch (*nr_travellers) {
						case 0:
						case 1:
						case 4:
						case 5:
							*nrPAR = 3;
							break;
						case 2:
						case 3:
							*nrPAR = 4;
							break;
						default:
							fprintf(stderr, "Unrecognized method to compute nr of travellers!\n");
							exit(1);
					}
				}
			}

		}
		else if (strncmp(arr[i], "-PAR", flag_len) == 0) {
			temp = *nrPAR;
			if (*iter_flow_model > 0)
				temp = *nrPAR * *iter_flow_model;

			*PAR    = (long double *)malloc(temp * sizeof(long double));
			for (int p = 0; p < temp; p++) {
				(*PAR)[p] = get_longdouble(arr[++i], arr[j]);
			}
			verify++;
		}
		else if (strncmp(arr[i], "-minPAR", flag_len) == 0) {
			*minPAR = (long double *)malloc(*nrPAR * sizeof(long double));
			for (int p = 0; p < *nrPAR; p++) {
				(*minPAR)[p] = get_longdouble(arr[++i], arr[j]);
			}
			verify++;
		}
		else if (strncmp(arr[i], "-aPAR", flag_len) == 0) {
			*aPAR   = (long double *)malloc(*nrPAR * sizeof(long double));
			for (int p = 0; p < *nrPAR; p++) {
				(*aPAR)[p] = get_longdouble(arr[++i], arr[j]);
			}
			verify++;
		}
		else if (strncmp(arr[i], "-bPAR", flag_len) == 0) {
			*bPAR   = (long double *)malloc(*nrPAR * sizeof(long double));
			for (int p = 0; p < *nrPAR; p++) {
				(*bPAR)[p] = get_longdouble(arr[++i], arr[j]);
			}
			verify++;
		}
		else if (strncmp(arr[i], "-min_pop", flag_len) == 0) {
			*min_pop = get_long(arr[++i], arr[j]);
			verify++;
		}
		else if (strncmp(arr[i], "-iter_flow_model", flag_len) == 0) {
			*iter_flow_model = get_long(arr[++i], arr[j]);
			verify++;
		}
		else if (strncmp(arr[i], "-iterLHS", flag_len) == 0) {
			*iterLHS = get_long(arr[++i], arr[j]);
			verify++;
		}
		else if (strncmp(arr[i], "-iterMCMC", flag_len) == 0) {
			*iterMCMC = get_long(arr[++i], arr[j]);
			verify++;
		}
		else if (strncmp(arr[i], "-do_symmetric", flag_len) == 0) {
			*do_symmetric = get_long(arr[++i], arr[j]);
			verify++;
		}
		else if (strncmp(arr[i], "-do_skipdiagonal", flag_len) == 0) {
			*do_skipdiagonal = get_long(arr[++i], arr[j]);
			verify++;
		}
		else if (strncmp(arr[i], "-do_rescalekappa", flag_len) == 0) {
			*do_rescalekappa = get_long(arr[++i], arr[j]);
			verify++;
		}

		// parameters for which default values are given
		else if (strncmp(arr[i], "-output_path", flag_len) == 0) {
			*output_path = arr[++i];
		}
		else if (strncmp(arr[i], "-output_suffix", flag_len) == 0) {
			*output_suffix = arr[++i];
		}
		else if (strncmp(arr[i], "-nr_chains", flag_len) == 0) {
			*nr_chains = get_long(arr[++i], arr[j]);
		}
		else if (strncmp(arr[i], "-chain_nr", flag_len) == 0) {
			*chain_nr = get_long(arr[++i], arr[j]);
		}
		else if (strncmp(arr[i], "-seed", flag_len) == 0) {
			*seed = get_long(arr[++i], arr[j]);
		}
		else if (strncmp(arr[i], "-acceptance_rate", flag_len) == 0) {
			*acceptance_rate = get_longdouble(arr[++i], arr[j]);
		}
		else if (strncmp(arr[i], "-tweaking_factor", flag_len) == 0) {
			*tweaking_factor = get_longdouble(arr[++i], arr[j]);
		}
		else {
			fprintf(stderr, "Unrecongized flag: \"%s\"\n", arr[i]);
			exit(1);
		}
	} // end of for loop

	// verify that all the essential flags have been provided
	if (verify != ess_flags) {
		fprintf(stderr, "Not all of the essential flags have been provided! %d != %d\n", verify, ess_flags);
		exit(1);
	}

	// verify that flags provided won't break anything later on in the programme
	if ( strncmp(*country, "Kenya", flag_len) !=0 && strncmp(*country, "Namibia", flag_len) ) {
		fprintf(stderr, "Unrecognized country: \"%s\"\n", *country);
		exit(1);
	}
	if ( strncmp(*model, "GM", flag_len) !=0 && strncmp(*model, "RM", flag_len) != 0 ) {
		fprintf(stderr, "Unrecognized mobility model: \"%s\"\n", *model);
		exit(1);
	}
	if ( *version < 0 || *version > 2 ) {
		fprintf(stderr, "Unrecognized version of the radiation model: \"%d\"\n", *version);
		exit(1);
	}
	if ( *iterLHS < 0 || *iterMCMC < 0 || *iter_flow_model < 0) {
		fprintf(stderr, "Number of LHS, MCMC and flow_model iterations must be >= 0!\n");
		exit(1);
	}
	if ( *iterLHS == 0 && *iterMCMC == 0 && *iter_flow_model == 0) {
		fprintf(stderr, "LHS, MCMC and flow_model iterations can't all be null!\n");
		exit(1);
	}
	if ( *nr_chains < 1 ) {
		fprintf(stderr, "Minimum number of chains is 1!\n");
		exit(1);
	}


	return 0;
}

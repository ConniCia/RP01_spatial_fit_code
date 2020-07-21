/* -------------------------------------------------------------------------- */
/* Function prototypes                                                        */
/* -------------------------------------------------------------------------- */

int parser(int n, char *arr[],
	     		 char **input_folder, char **country, char **scale, char **model, char **output_path, char **output_suffix,
			 	 	 int *nrPAR, int *nr_travellers, int *version, long double **PAR, long double **minPAR, long double **aPAR, long double **bPAR,
					 int *do_symmetric, int *do_skipdiagonal, int *do_rescalekappa,
			 	 	 int *min_pop, int *iter_flow_model, int *iterLHS, int *iterMCMC, int *nr_chains, int *chain_nr, int *seed,
			 	 	 long double *acceptance_rate, long double *tweaking_factor
);

int read_input(char *input_folder, char *scale,
  				 		 int *nr_admins, int *nr_patches,
  			 	 		 int **cum_nr_patches, int **patch_population, int **admin_population,
  			 	 		 int ***patch_distance, int ***admin_flow_data,
  			 	 		 long double ***admin_flow_model
);

int prepare_RM(int ***patch_radius_population_L, int ***patch_radius_population_S, int **admin_nr_travellers,
  					 	 int nr_patches, int nr_admins, int *patch_population, int **patch_distance, int **admin_flow_data,
							 char *input_folder, char *country, char *scale
);

int read_RM_variables(int ***patch_radius_population_L, int ***patch_radius_population_S, int **admin_nr_travellers,
  										int nr_patches, int nr_admins, int *patch_population, int **patch_distance, int **admin_flow_data,
											char *input_folder, char *country, char *scale
);

int create_RM_variables(int ***patch_radius_population_L, int ***patch_radius_population_S, int **admin_nr_travellers,
  											int nr_patches, int nr_admins, int *patch_population, int **patch_distance, int **admin_flow_data,
												char *input_folder, char *country, char *scale
);

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
				int iter_flow_model // doesn't do anything in LHS
);

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
				 int iter_flow_model // doesn't do anything in MCMC
);

long double computeLOGL(int do_symmetric, int do_skipdiagonal, int do_rescalekappa, int min_pop,
												long double *PAR,
                        int nr_admins, int nr_patches, int nrPAR, int nr_travellers, int version,
                        int *cum_nr_patches, int *patch_population,
                        int **patch_distance, int **admin_flow_data,
												int **patch_radius_population_L, int **patch_radius_population_S,
												int *admin_nr_travellers, int *admin_population,
                        long double **admin_flow_model, char *model,
												char output_folder[], int iter_flow_model, int nr_flow_model
);

/* Print functions ---------------------------------------------------------- */
// parser
int print_parser(char *input_folder, char *country, char *scale, char *model, char *output_path, char *output_suffix,
			 					 int nrPAR, int nr_travellers, int version, long double *PAR, long double *minPAR, long double *aPAR, long double *bPAR,
								 int do_symmetric, int do_skipdiagonal, int do_rescalekappa,
			 				 	 int min_pop, int iterLHS, int iterMCMC, int nr_chains, int chain_nr, int seed,
			 				 	 long double acceptance_rate, long double tweaking_factor,
								 char output_folder[]
);
// time
int print_time(char output_folder[], double CPUdiff);
// admin flow model
int fprint_flow_model (long double LOGL, long double *PAR, long double **admin_flow_model,
											 int nrPAR, int nr_flow_model, int iter_flow_model, int nr_admins, int do_symmetric,
											 char output_folder[]);

/* Utility functions -------------------------------------------------------- */
// parser
long int get_long(char *str_value, char *str_flag);
long double get_longdouble(char *str_value, char *str_flag);
int symmetrise_matrix(int **M, long int n);
int remove_diagonal(int **M, long int n);
static int compar(const void *a, const void *b);
int randomise(int **arr, int n);
int resort(int nr_chains, int nrPAR, long double **chainLOGL, long double ***chainPAR);

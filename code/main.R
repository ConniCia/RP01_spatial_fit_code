#===============================================================================
#
# DISCLAIMERES
#
# This code assumes you have installed the R packages listed below as well as
# gcc and OpenMP (can be installed via homebrew on Mac).
# Creating the input files for 5km scales takes a long time, I recommend using
# a cluster.
# Running the C simulations takes an aweful long time, use a cluster! The admin
# and 20km scales are ok on a standard computer, but the 10 and 5km scales are
# demanding even on a cluster. 1 million iterations for the 5km scale on a 32
# core cluster take about 3 months.
# I know some of the code is not ideal. This paper was my first solo code
# development and I was still wrapping my head around the best coding practices.
# I've rewritten parts of the code in 2020: I collated all the different parts
# of the code in this single main script and have rewritten all of the figure
# and table functions. Code is not perfect, but reporducible. C code is as
# streamlined, and hence fast, as I could possibly make it.
#
#===============================================================================


# load packages ----------------------------------------------------------------

library(data.table)
library(dplyr)
library(tidyr)
library(geosphere)  # compute distances
library(tictoc)     # compute time in cleanRawData()
library(stringr)
library(xtable)     # print LaTeX tables
library(scales)     # format printed numbers
library(ggplot2)
library(gridExtra)  # arrange plot panels
library(grid)       # arrange plot panels
library(sf)         # handle and draw maps
library(stars)      # find boundaries from pixel data
library(latex2exp)  # Greek letters in figures


# source all R functions -------------------------------------------------------

list.files(path = "code/", pattern = "^R", full.names = TRUE)  %>%
  list.files(pattern = ".R$", recursive = TRUE, full.names = TRUE) %>%
  lapply(source) %>%
  invisible()


# Options ----------------------------------------------------------------------

countries      <- c("Kenya", "Namibia")
models         <- c("GM",
                    "RM-v0-t5", "RM-v0-t3",
                    "RM-v1-t0", "RM-v1-t1", "RM-v1-t2", "RM-v1-t3",
                    "RM-v2-t0", "RM-v2-t1", "RM-v2-t2", "RM-v2-t3", "RM-v2-t4")
grid_scales    <- c(20, 10, 5) # on top of admin scale
min_pop_limit  <- 10 # to apply at 5km scale

models_DEFAULT <- c("GM", "RM-v2-t3", "RM4") # main models (discussed in main text of paper)
scales_DEFAULT <- c("admin", "20") # we only apply more detailed grid to default tests

# country-model-scale combinations that give the best log-likelihood
dt_best <- data.frame(
  test_type = rep(c("DEFAULT_symmetric_skipdiagonal_rescalekappa", "DEFAULT"), 2),
  model   = c("GM",    "GM",      "RM-v2-t3", "RM-v2-t3"),
  country = c("Kenya", "Namibia", "Kenya",    "Namibia"),
  scale   = c("5",     "admin",   "20",       "admin"),
  notes   = ""
)


#===============================================================================

# clean raw data for simulation input --------------------------------------------
#
# wrapper_cleanRawData(countries            = countries,
#                      grid_scales          = grid_scales,
#                      scales_DEFAULT       = scales_DEFAULT,
#                      grid_changes         = grid_changes,
#                      grid_changes_Namibia = grid_changes_Namibia)
#
#
# # compile C code ---------------------------------------------------------------
#
# list.files(dir_C, ".c", full.names = TRUE) %>%
#   paste(collapse = " ") %>%
#   sprintf(fmt = "gcc-10 %s -g -fopenmp -lm -o %s", path_exec_mac) %>%
#   system()
#
#
# # create batch files to launch simulations -------------------------------------
#
# wrapper_createBatchFiles(grid_changes         = grid_changes,
#                          grid_changes_Namibia = grid_changes_Namibia,
#                          flag_changes         = flag_changes,
#                          models               = models,
#                          models_DEFAULT       = models_DEFAULT,
#                          countries            = countries,
#                          scales               = c("admin", grid_scales),
#                          scales_DEFAULT       = scales_DEFAULT,
#                          stages = c("LHS", "preMCMC", "MCMC", "flow_model",
#                                     "flow_model_5km_admin", "flow_model_default_data"))
#
#
# # clean simulation output ------------------------------------------------------
#
# wrapper_cleanSimulationOutput(
#   stages = c("input_data", "pixel_data",
#              "LHS", "preMCMC", "MCMC",
#              "flow_model", "_flow_model"), # first matrices, second LOGL
#   default_data = flag_changes$test_type)
#
#
# #===============================================================================
#
# # create tables (for reference and publication) --------------------------------
#
# tab_grid_cells()
# tab_models()
# tab_compare_flow_estimates()
# tab_parameter_interpretation()
#
#
# # create figures for validation and reference ----------------------------------
#
# fig_LHS()
# fig_MCMC(input_file = "preMCMC.rds", output_folder = "preMCMC_convergence")
# fig_MCMC(input_file = "MCMC.rds",    output_folder = "MCMC_mixing")
# fig_MCMC(input_file = "MCMC.rds",    output_folder = "MCMC_posterior_density")
#
#
# # create figures for publication -----------------------------------------------
#
# # main text
# do.paper_name <- TRUE
# fig_parameter_trends(do.paper_name)
# fig_heatmap(do.paper_name)
# fig_distance_profile(do.paper_name)
# map_population_density(do.paper_name)
#
# # SI figures
# do.paper_name <- FALSE
# fig_cell_population(do.paper_name)
# map_connectivity(do.paper_name)
# fig_heatmap(do.paper_name, filter = "5km_admin")
# fig_heatmap(do.paper_name, filter = "random50")
# fig_heatmap(do.paper_name, filter = "random15")
# fig_distance_profile(do.paper_name, filter = "5km_admin")
# fig_distance_profile(do.paper_name, filter = "random50")
# fig_distance_profile(do.paper_name, filter = "random15")
#
# # SI methods
# fig_symmetricity(do.paper_name)
# fig_parameter_correlation(do.paper_name)
# fig_parameter_interpretation(do.paper_name)

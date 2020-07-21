# Aim:
#   create batch files to launch simulations using the cluster web interface
# Steps:
#   - loop over all folders where to put batch files
#   - test whether to create or launch LHS, preMCMC, MCMC and flow_model

wrapper_createBatchFiles <- function (grid_changes, grid_changes_Namibia, flag_changes,
                              models, models_DEFAULT, countries,
                              scales, scales_DEFAULT, stages) {


  # combine grid and flag options ---------------------------------------------#

  test_data_Namibia <- expand.grid(
    country         = "Namibia",
    do_symmetric    = 0,
    do_skipdiagonal = 0,
    do_rescalekappa = 0,
    data_type       = c(grid_changes, grid_changes_Namibia)
  ) %>%
    mutate(test_type = data_type)
  test_data_Kenya <- expand.grid(
    country         = "Kenya",
    do_symmetric    = 1,
    do_skipdiagonal = 1,
    do_rescalekappa = 1,
    data_type       = c(grid_changes)
  ) %>%
    mutate(test_type = data_type)
  test_data_combos <- bind_rows(flag_changes, test_data_Namibia, test_data_Kenya)


  # loop ----------------------------------------------------------------------#

  for (tdc in seq_len(nrow(test_data_combos))) {
    for (model in models) {

      td <- test_data_combos[tdc, ]

      # DEFAULT TESTS ON DEFAULT MODELS!!
      if (grepl("DEFAULT", td$test_type) & model %in% models_DEFAULT) {

        # LHS
        if ("LHS" %in% stages & "admin" %in% scales)
          createBatch(td, model, "LHS", "admin")

        # preMCMC, MCMC, flow_model
        my.stages <- intersect(stages, c("preMCMC", "MCMC", "flow_model"))
        for (stage in my.stages) {
          for (scale in scales) {
            createBatch(td, model, stage, scale)
          }
        }

        # flow_model_5km_admin
        if ("flow_model_5km_admin" %in% stages & "5" %in% scales)
          createBatch(td, model, "flow_model_5km_admin", "5")


      # NON DEFAULT TESTS ON DEFAULT MODELS, DEFAULT TESTS ON NON DEFAULT MODELS
      } else if ((!grepl("DEFAULT", td$test_type) &  model %in% models_DEFAULT) |
                 ( grepl("DEFAULT", td$test_type) & !model %in% models_DEFAULT)) {

        if (td$country == "Kenya" & grepl("-t1", model))
          next

        # LHS
        if ("LHS" %in% stages & "admin" %in% scales)
          createBatch(td, model, "LHS", "admin")

        # preMCMC, MCMC, flow_model
        my.stages <- intersect(stages, c("preMCMC", "MCMC", "flow_model"))
        my.scales <- intersect(scales, c("admin", "20"))
        for (stage in my.stages) {
          for (scale in my.scales) {
            createBatch(td, model, stage, scale)
          }
        }

        # flow_model_default_data
        if (grepl("random", td$test_type)) {
          my.stages <- intersect(stages, "flow_model_default_data")
          my.scales <- intersect(scales, c("admin", "20"))
          for (stage in my.stages) {
            for (scale in my.scales) {
              createBatch(td, model, stage, scale)
            }
          }
        }

      }

    }
  }

  if (!do.cluster)
    system(paste0("chmod +x ", dir_proc, "*/*/*.sh"))

}


# function ---------------------------------------------------------------------

createBatch <- function (td, model, stage, scale) {

  # define paths --------------------------------------------------------------#

  data_type <- if (stage != "flow_model_default_data") td$data_type else "DEFAULT"

  input_folder_win    <- file.path(dir_cln_win,  paste_(td$country, data_type), fsep = "\\")
  input_folder_mac    <- file.path(dir_cln,      paste_(td$country, data_type))
  output_path_win     <- file.path(dir_proc_win, paste_("test", td$test_type), paste_(model, td$country, scale), fsep = "\\")
  output_path_mac     <- file.path(dir_proc,     paste_("test", td$test_type), paste_(model, td$country, scale))
  output_path_mac_adm <- file.path(dir_proc,     paste_("test", td$test_type), paste_(model, td$country, "admin"))
  extension <- if (do.cluster) "bat" else "sh"


  # initialise flags ----------------------------------------------------------#

  df_flags <- initialiseFlags(td, model, stage, scale,
                              input_folder_win, input_folder_mac,
                              output_path_win, output_path_mac)


  # modify initalisation of PAR -----------------------------------------------#

  ii <- if (grepl("MCMC", stage)) 1:4 else 1
  for (i in ii) {

    stage_print <- if (length(ii) > 1) paste_(stage, i) else stage
    df_flags$output_suffix <- stage_print

    path_mac <- file.path(output_path_mac, paste.(stage_print, extension))
    path_win <- file.path(output_path_win, paste.(stage_print, "bat"), fsep = "\\")
    suppressWarnings(dir.create(output_path_mac, recursive = TRUE))

    # important: if batch file already exists, don't waste tima and skip to next job
    if (file.exists(path_mac))
      next

    if (stage == "LHS") {
      # use default data
      newPAR <- df$PAR
    } else if (stage == "preMCMC" & scale == "admin") {

      # use i-th best admin LHS parameter combination
      input <- file.path(output_path_mac, "output_files_LHS/LHS.tsv")
      if (!file.exists(input))
        next

      newPAR <- input %>%
        read.table(header = TRUE) %>%
        arrange(-LOGL) %>%
        slice(i) %>%
        select(-(1:2)) %>%
        paste(collapse = " ")

    } else if (stage == "preMCMC" & scale != "admin") {

      # use end of i-th admin preMCMC chain
      input <- file.path(output_path_mac_adm, paste_("output_files", stage_print), "MCMC_0.tsv")
      if (!file.exists(input))
        next

      newPAR <- input %>%
        read.table(header = TRUE) %>%
        slice_tail() %>%
        select(-(1:2)) %>%
        paste(collapse = " ")

    } else if (stage == "MCMC") {

      # use end of i-th preMCMC chain
      input <- file.path(output_path_mac, paste0("output_files_pre", stage_print), "MCMC_0.tsv")
      if (!file.exists(input))
        next

      newPAR <- input %>%
        read.table(header = TRUE) %>%
        slice_tail() %>%
        select(-(1:2)) %>%
        paste(collapse = " ")

    } else if (stage %in% c("flow_model", "flow_model_5km_admin", "flow_model_default_data")) {

      # use 100 random samples of pooled MCMC posteriors
      input <- list.files(path = if (stage != "flow_model_5km_admin") output_path_mac else output_path_mac_adm,
                          pattern = "output_files_MCMC", full.names = TRUE)

      if (length(input) < 4 | !all(file.exists(input)))
        next

      set.seed(1) # for slice_sample() below
      newPAR <- input %>%
        lapply(function(x) {
          list.files(x, pattern = "MCMC_0", full.names = TRUE) %>%
            read.table(header = TRUE) %>%
            # remove burn-in
            removeBurnIn()
        }) %>%
        bind_rows() %>%
        slice_sample(n = 100) %>%
        select(-(1:2))

      meanPAR <- newPAR %>%
        summarise(across(everything(), ~round(mean(.x), 2)))

      newPAR <- bind_rows(newPAR, meanPAR) %>%
        t() %>%
        paste(collapse = " ")


    }

    df_flags$PAR <- newPAR

    # print to file -----------------------------------------------------------#
    line <- seq_along(df_flags) %>%
      sapply(function (x) paste0("-", names(df_flags[x]), " ", df_flags[x])) %>%
      paste(collapse = " ") %>%
      sub(".*? ", "", .) %>%
      writeLines(path_mac)

    if (do.cluster) {
      cat(path_win, "\n")
    } else {
      cat(path_mac, "\n")
    }

  }

}


# helper function --------------------------------------------------------------

initialiseFlags <- function (td, model, stage, scale,
                             input_folder_win, input_folder_mac,
                             output_path_win, output_path_mac) {

  # initialise most flags
  df_flags <- data.frame(
    path_exec       = if (do.cluster) path_exec_win else path_exec_mac,
    input_folder    = if (do.cluster) input_folder_win else input_folder_mac,
    output_path     = if (do.cluster) output_path_win else output_path_mac,
    output_suffix   = "", # to be defined below
    country         = td$country,
    scale           = scale,
    model           = gsub("-[a-z]", " ", model),
    min_pop         = if (scale == "5") min_pop_limit else 1,
    iter_flow_model = if (stage %in% c("flow_model", "flow_model_5km_admin", "flow_model_default_data")) 101 else 0,
    iterLHS         = if (stage == "LHS") 10000 else 0,
    iterMCMC        = if (stage == "preMCMC") {300000} else if (stage == "MCMC") {1000000} else {0},
    do_symmetric    = td$do_symmetric,
    do_skipdiagonal = td$do_skipdiagonal,
    do_rescalekappa = td$do_rescalekappa
  )

  # initialise vectors
  ## N.B. the value of PAR will get changed in case of preMCMC or MCMC
  if (model == "GM") {                                      # kappa, alpha, beta, gamma, epsilon, disp
    df <- data.frame(minPAR = "0.1  0.1  0.1  0.1  0.1  0.1",
                     aPAR   = " -5    0    0    0    0    0",
                     bPAR   = "  6    5    5    5    5    5",
                     PAR    = "  1    1    1    2    2    2")
  } else if (model == "RM-v0-t5") {                         # kappa, disp
    df <- data.frame(minPAR = "0.1  0.1",
                     aPAR   = "  0    0",
                     bPAR   = " 50    5",
                     PAR    = "  6    2")
  } else if ( model == "RM-v1-t0" | model == "RM-v1-t1" ) { # kappa, eta, disp
    df <- data.frame(minPAR = "0.1  0.1  0.1",
                     aPAR   = " -5    0    0",
                     bPAR   = "  6  500    5",
                     PAR    = "  3    5    3")
  } else if ( model == "RM-v1-t2" | model == "RM-v1-t3" ) { # kappa, eta, alpha, disp
    df <- data.frame(minPAR = "0.1  0.1  0.1  0.1",
                     aPAR   = " -5    0    0    0",
                     bPAR   = "  6  500    5    5",
                     PAR    = "  1    3    1    3")
  } else if ( model == "RM-v2-t0" | model == "RM-v2-t1" ) { # kappa, theta, disp
    df <- data.frame(minPAR = "0.1       0.1  0.1",
                     aPAR   = " -5         0    0",
                     bPAR   = "  6  10000000    5",
                     PAR    = "  3      1000    3")
  } else if ( model == "RM-v2-t2" | model == "RM-v2-t3" ) { # kappa, theta, alpha, disp
    df <- data.frame(minPAR = "0.1       0.1  0.1  0.1",
                     aPAR   = " -5         0    0    0",
                     bPAR   = "  6  10000000    5    5",
                     PAR    = " -2    100000    1    2")
  } else if (model == "RM-v0-t3") {                         # kappa, alpha, disp
    df <- data.frame(minPAR = "0.1  0.1  0.1",
                     aPAR   = " -5    0    0",
                     bPAR   = "  6    5    5",
                     PAR    = "  6    1    2")
  } else if (model == "RM-v2-t4") {                         # kappa, theta, disp
    df <- data.frame(minPAR = "0.1       0.1  0.1",
                     aPAR   = " -5         0    0",
                     bPAR   = "  6  10000000    5",
                     PAR    = " -2    100000    2")
  }

  # combine
  df_flags <- bind_cols(df_flags, df)

  return(df_flags)

}

tab_grid_cells <- function () {

  cat("Create table of nr of patches per country and scale.\n")

  readRDS(file.path(dir_res, "metadata.rds")) %>%
    filter(data_type == "DEFAULT") %>%
    select(-data_type) %>%
    pivot_wider(names_from = country, values_from = nr_patches) %>%
    rename(Scale = scale) %>%
    xtable(type = "latex") %>%
    print(file = file.path(dir_tbl, "patches_per_grid.tex"),
          include.rownames = FALSE,
          floating = FALSE)

}


tab_models <- function () {

  cat("Create tables of log-likelihood (and parameter estimate) means and 95% CrI.\n")


  # read in and transform -----------------------------------------------------#

  DT <- readRDS(file.path(dir_res, "MCMC.rds")) %>%
    # compute and format CrIs of log-likelihood and parameters
    cleanTable()


  # print tables by model -----------------------------------------------------#

  DT %>%
    # arrange
    arrange(Country) %>%
    # print by model
    group_split(Model) %>%
    lapply(function(dt) {

      model <- unique(dt$Model)
      print(model)

      dt <- dt %>%
        select(-Model) %>%
        # reshape table
        pivot_wider(
          names_from = variable,
          values_from = value
        )

      dt %>%
        select(-test_type) %>%
        xtable(type = "latex") %>%
        print(file = file.path(dir_tbl, paste.(model, "tex")),
              include.rownames = FALSE,
              # allows me to print in math mode
              sanitize.text.function = function(x){x},
              floating = FALSE)

      if (model %in% models_DEFAULT & any(grepl("DEFAULT", unique(dt$test_type)))) {
        dt %>%
          filter(grepl("DEFAULT", test_type)) %>%
          select(-test_type) %>%
          xtable(type = "latex") %>%
          print(file = file.path(dir_tbl, paste_(model, "DEFAULT.tex")),
                include.rownames = FALSE,
                # allows me to print in math mode
                sanitize.text.function = function(x){x},
                floating = FALSE)
      }

      return(NULL)

    })


  # print table by test_type --------------------------------------------------#

  DT %>%
    group_split(test_type) %>%
    lapply(function(dt) {

      test_type <- unique(dt$test_type)
      print(test_type)

      dt <- dt %>%
        # filter by log-likelihood
        # filter(variable == "log-likelihood") %>%
        # reshape table
        pivot_wider(
          names_from = variable,
          values_from = value
        )

      if (!grepl("random", test_type)) {

        dt %>%
          select(-test_type) %>%
          # arrange
          relocate(Country, Scale, Model, `log-likelihood`) %>%
          arrange(Country, Scale, Model, `log-likelihood`) %>%
          xtable(type = "latex") %>%
          print(file = file.path(dir_tbl, paste.(test_type, "tex")),
                include.rownames = FALSE,
                # allows me to print in math mode
                sanitize.text.function = function(x){x},
                floating = FALSE)

      } else {

        # add flow model matrix LOGL results on full input data
        dt2 <- readRDS(file.path(dir_res, "_flow_model.rds")) %>%
          # remove sampel produced by mean of posterior parameter distributions
          filter(id_sample < max(id_sample) - 0.5) %>%
          # compute and format CrIs of log-likelihood and parameters
          cleanTable(filter.notes = TRUE) %>%
          # filter by log-likelihood
          filter(variable == "log-likelihood") %>%
          # reshape table
          pivot_wider(
            names_from = variable,
            values_from = value
          )

        # join tables
        left_join(dt, dt2,
                  by = c("test_type", "Model", "Country", "Scale"),
                  suffix = c(" partial dataset", " full dataset")) %>%
          select(-test_type) %>%
          # arrange
          relocate(Country, Scale, Model, `log-likelihood full dataset`, `log-likelihood partial dataset`) %>%
          arrange(Country, Scale, Model, `log-likelihood full dataset`, `log-likelihood partial dataset`) %>%
          xtable(type = "latex") %>%
          print(file = file.path(dir_tbl, paste.(test_type, "tex")),
                include.rownames = FALSE,
                # allows me to print in math mode
                sanitize.text.function = function(x){x},
                floating = FALSE)

      }

      return(NULL)

    })

}


tab_compare_flow_estimates <- function (print_to_latex = TRUE) {

  cat("Create table to compare CrI of flow estimates with data.\n")

  # read in data ---------------------------------------------------------------

  # read in flow counts
  dt_data <- readRDS(file.path(dir_res, "admin_flow_data.rds")) %>%
    # filter default input data
    filter(grepl("DEFAULT", data_type)) %>%
    select(-scale)

  # read in simulation estimates
  dt_estim <- readRDS(file.path(dir_res, "flow_model.rds")) %>%
    # filter default input data
    filter(grepl("DEFAULT", test_type) & model %in% models_DEFAULT & notes == "") %>%
    select(-data_type)

  # combine tables
  dt <- compareFlowEstimates(dt_data, dt_estim) %>%
    # compute stats
    group_by(country, model, scale, notes) %>%
    mutate(nr_admin_patches = n()) %>%
    group_by(country, model, scale, notes, within95CrI) %>%
    summarise(n = round(100 * n() / unique(nr_admin_patches), 2)) %>%
    ungroup() %>%
    # reshape
    pivot_wider(names_from = within95CrI,
                values_from = n) %>%
    # formatting
    arrange(-within) %>%
    mutate(model = if_else(model == "RM-v2-t3", "RM4", model)) %>%
    rename(Country                   = country,
           Model                     = model,
           Scale                     = scale,
           Notes                     = notes,
           `Lower than 95% CrI (%)`  = lower,
           `Within 95% CrI (%)`      = within,
           `Higher than 95% CrI (%)` = higher)

  if (print_to_latex) {
    dt %>%
      xtable(type = "latex") %>%
      print(file = file.path(dir_tbl, "compare_flow_CrIs.tex"),
            include.rownames = FALSE,
            floating = FALSE)
  } else {
    return(dt)
  }

}


tab_parameter_interpretation <- function () {

  cat("Create table to explain interpretation of model parameters.\n")

  dt <- list()

  dt[[1]] <- data.frame(
    Parameter      = "$disp$",
    Formula        = "log-likelihood",
    Position       = "over-dispersion parameter of negative binomial distribution",
    Reference      = "-",
    Interpretation = ""
  )

  dt[[2]] <- data.frame(
    Parameter      = "$\\kappa$",
    Formula        = "GM and RM",
    Position       = "proportionality constant",
    Reference      = "-",
    Interpretation = ""
  )

  dt[[3]] <- data.frame(
    Parameter      = "$\\alpha$",
    Formula        = "GM, RM2 and RM4",
    Position       = "exponent of origin population",
    Reference      = "Fig S6A",
    Interpretation = ""
  )

  dt[[4]] <- data.frame(
    Parameter      = "$\\beta$",
    Formula        = "GM",
    Position       = "exponent of destination population",
    Reference      = "Fig S6A",
    Interpretation = ""
  )

  dt[[5]] <- data.frame(
    Parameter      = "$\\gamma$",
    Formula        = "spatial kernel of GM",
    Position       = "distance scale",
    Reference      = "Fig S6C",
    Interpretation = ""
  )

  dt[[6]] <- data.frame(
    Parameter      = "$\\varepsilon$",
    Formula        = "spatial kernel of GM",
    Position       = "exponent of scaled distance",
    Reference      = "Fig S6C",
    Interpretation = ""
  )

  dt[[7]] <- data.frame(
    Parameter      = "$\\theta$",
    Formula        = "RM3 and RM4",
    Position       = "increment of origin population",
    Reference      = "Fig S6B",
    Interpretation = ""
  )

  dt <- dt %>%
    bind_rows() %>%
    mutate(Parameter = as.character(Parameter))%>%
    xtable(type = "latex") %>%
    print(file = file.path(dir_tbl, "interpretation.tex"),
          include.rownames = FALSE,
          # allows me to print in math mode
          sanitize.text.function = function(x){x},
          floating = FALSE)


}


# Helper functions -------------------------------------------------------------

# compute CrI in print format
compute_CrI <- function(X, method) {

  if (is.na(sum(X))) {
    NA
  } else if (method == "signif3") {
    paste0(signif(    mean(X),        3), " (",
           signif(quantile(X, 0.025), 3), ", ",
           signif(quantile(X, 0.975), 3), ")"
    )
  } else if (method == "round") {
    paste0(round(    mean(X)),        " (",
           round(quantile(X, 0.025)), ", ",
           round(quantile(X, 0.975)), ")"
    )
  }

}


# compute and format CrIs of log-likelihood and parameters
cleanTable <- function (dt, filter.notes = FALSE) {

  if (filter.notes) {
    dt <- dt %>% filter(notes == "default_data")
  }

  dt <- dt %>%
    # sort variables
    mutate(variable = factor(variable, levels = getVariableNames("all"))) %>%
    # remove kappa
    # filter(variable != "kappa") %>%
    # model main text names
    mutate(model = case_when(
      model == "RM-v0-t5"  ~  "RM1",
      model == "RM-v0-t3"  ~  "RM2",
      model == "RM-v2-t4"  ~  "RM3",
      model == "RM-v2-t3"  ~  "RM4",
      TRUE                 ~  model
    ))  %>%
    # mutate(model = getModelName(model)) %>%
    # compute CrIs
    group_split(test_type, model, country, scale) %>% # pool chains
    lapply(function(dt) {

      dt1 <- dt %>%
        # compute CrI, rounded
        filter(variable == "log-likelihood") %>%
        group_by(test_type, model, country, scale, variable) %>%
        summarise(value = compute_CrI(value, method = "round")) %>%
        ungroup()

      dt2 <- dt %>%
        # compute CrI with 3 significant figures
        filter(variable != "log-likelihood") %>%
        group_by(test_type, model, country, scale, variable) %>%
        summarise(value = compute_CrI(value, method = "signif3")) %>%
        ungroup() %>%
        # modify variable names so they print well in LaTeX
        mutate(variable = if_else(variable == "disp",
                                  "$disp$", if_else(variable == "epsilon",
                                                    "$\\varepsilon$",
                                                    paste0("$\\", variable, "$"))))

      bind_rows(dt1, dt2)

    }) %>%
    bind_rows() %>%
    rename(Country = country,
           Scale   = scale,
           Model   = model)

  # rename admin
  levels(dt$Scale)[levels(dt$Scale) == "admin"] <- "Administrative unit"

  return(dt)

}


# Check whether input flows fall into 95% CrIs of estimated flows
compareFlowEstimates <- function (dt_data, dt_estim) {

  # merge
  ## Namibia GM
  ### input data came with within-unit travel (i.e. with diagonal) and was directional (i.e. no need to symmetrise)
  ### GM fits the diagonal
  ### -> input data and simulation estimates can be compared
  dt_Nam_GM <- inner_join(
    dt_data  %>% filter(country == "Namibia" & data_type == "DEFAULT"),
    dt_estim %>% filter(country == "Namibia" & model == "GM" & test_type == "DEFAULT")
  )

  ## Namibia RM
  ### input data came with within-unit travel (i.e. with diagonal) and was directional (i.e. no need to symmetrise)
  ### RM never fits the diagonal
  ### -> input data needs to have the diagonal removed to be compared with simulation estimates
  dt_Nam_RM <- inner_join(
    dt_data  %>% filter(country == "Namibia" & data_type == "DEFAULT_skipdiagonal"),
    dt_estim %>% filter(country == "Namibia" & model == "RM-v2-t3" & test_type == "DEFAULT")
  )

  ## Kenya GM
  ### input data came without diagonal and was not directional (i.e. had to be symmetrised)
  ### C simulation symmetrised the matrix and GM skipped to fit the diagonal
  ### -> input data has to be symmetrised before comparison with simulation estimates
  dt_Ken_GM <- inner_join(
    dt_data  %>% filter(country == "Kenya" & data_type == "DEFAULT_symmetric_skipdiagonal"),
    dt_estim %>% filter(country == "Kenya" & model == "GM" & test_type == "DEFAULT_symmetric_skipdiagonal_rescalekappa")
  )

  ## Kenya RM
  ### input data came without diagonal and was not directional (i.e. had to be symmetrised)
  ### C simulation symmetrised the matrix and RM never fits the diagonal
  ### -> input data has to be symmetrised before comparison with simulation estimates
  dt_Ken_RM <- inner_join(
    dt_data  %>% filter(country == "Kenya" & data_type == "DEFAULT_symmetric_skipdiagonal"),
    dt_estim %>% filter(country == "Kenya" & model == "RM-v2-t3" & test_type == "DEFAULT_symmetric_skipdiagonal_rescalekappa")
  )

  # combine tables
  dt <- bind_rows(dt_Nam_GM, dt_Nam_RM, dt_Ken_GM, dt_Ken_RM) %>%
    # compute whether data falls within 95% CrI of estimates
    mutate(within95CrI = if_else(flow_count < q1 - 0.001,
                                 "lower",
                                 if_else(flow_count < q2 + 0.001,
                                         "within",
                                         "higher"))) %>%
    mutate(within95CrI = factor(within95CrI,
                                levels = c("lower", "within", "higher")))

  return(dt)

}

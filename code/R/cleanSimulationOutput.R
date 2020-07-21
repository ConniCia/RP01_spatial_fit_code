# Aim:
#   clean simulation output data so that I can then create figures

wrapper_cleanSimulationOutput <- function (stages, default_data) {

  for (stage in stages) {

    switch(
      stage,
      "input_data"  = cleanInputData(),
      "pixel_data"  = cleanPixelData(),
      "LHS"         = cleanSimParameterEstimates(stage, "LHS.tsv", default_data),
      "preMCMC"     = ,
      "MCMC"        = cleanSimParameterEstimates(stage, "MCMC_0.tsv", default_data),
      "flow_model"  = cleanSimParameterEstimates(stage, "flow_model_matrix", default_data),
      "_flow_model" = cleanSimParameterEstimates(stage, "flow_model_LOGL", default_data)
    )

  }

}

# functions --------------------------------------------------------------------

cleanInputData <- function () {

  cat("# Clean input data\n")

  dirs <- list.files(path = dir_cln, full.names = TRUE,
                     include.dirs = TRUE, recursive = FALSE)


  # metadata ------------------------------------------------------------------#

  cat("## Metadata\n")

  dirs %>%
    list.files(pattern = "__metadata", full.names = TRUE) %>%
    lapply(function(x) {

      m <- read.table(x, header = FALSE)

      parse_input_path(x) %>%
        mutate(nr_patches = m[2, 2])

    }) %>%
    bind_rows() %>%
    as_tibble() %>%
    mutate(scale = factor(scale, levels = c("admin", "20", "10", "5"))) %>%
    arrange(country, scale) %>%
    saveRDS(file = file.path(dir_res, "metadata.rds"))


  # patch population ----------------------------------------------------------#

  cat("## Population data\n")

  dirs %>%
    list.files(pattern = "__patch_population", full.names = TRUE) %>%
    lapply(function(x) {

      read.table(x, header = FALSE, col.names = "population") %>%
        mutate(id_patch = row_number() - 1) %>%
        bind_cols(parse_input_path(x))

    }) %>%
    bind_rows() %>%
    mutate(scale = factor(scale, levels = c("admin", "20", "10", "5"))) %>%
    arrange(country, scale) %>%
    as_tibble() %>%
    saveRDS(file = file.path(dir_res, "patch_population.rds"))


  # admin distance ------------------------------------------------------------#

  cat("## Distance data\n")

  dirs %>%
    list.files(pattern = "admin__patch_distance", full.names = TRUE) %>%
    lapply(function(x) {

      # read in lower triangular (in vector form)
      v <- read.table(x, header = FALSE)[[1]]

      # create matrix of correct dimension
      m <- v %>%
        length() %>%
        matrix_dim() %>%
        matrix(0, ncol = ., nrow = .)

      # insert data into matrix
      m[lower.tri(m)] <- v
      m <- m + t(m)

      # transform into table
      matrix_to_table(m, x, "distance") %>%
        mutate(distance = round(distance / 1000))

    }) %>%
    bind_rows() %>%
    as_tibble() %>%
    saveRDS(file = file.path(dir_res, "admin_distance.rds"))


  # admin flow data -----------------------------------------------------------#

  cat("## Flow data\n")

  dirs %>%
    list.files(pattern = "admin__admin_flow_data", full.names = TRUE) %>%
    lapply(function(x) {

      # read in flow data in matrix form
      m <- read.table(x, header = FALSE) %>% as.matrix()
      if (sum(diag(m)) == 0)
        diag(m) <- 0

      # tranform into table
      dt <- matrix_to_table(m, x, "flow_count")

      # from DEFAULT: create symmetric, skipdiagonal and symmetric_skipdiagonal
      ## mimic changes applied in C simulation
      if (unique(dt$data_type) == "DEFAULT") {

        # symmetric
        dt1 <- matrix_to_table(m + t(m), x, "flow_count") %>%
          mutate(data_type = paste_(data_type, "symmetric"))

        # skipdiagonal
        diag(m) <- 0
        dt2 <- matrix_to_table(m, x, "flow_count") %>%
          mutate(data_type = paste_(data_type, "skipdiagonal"))

        # symmetric & skipdiagonal
        dt3 <- matrix_to_table(m + t(m), x, "flow_count") %>%
          mutate(data_type = paste_(data_type, "symmetric_skipdiagonal"))

        # bind tables
        dt <- bind_rows(dt, dt1, dt2, dt3)

      }

      return(dt)

    }) %>%
    bind_rows() %>%
    as_tibble() %>%
    saveRDS(file = file.path(dir_res, "admin_flow_data.rds"))

}


cleanPixelData <- function () {

  # pixel data of all pixels and admins ---------------------------------------#

  cat("## Admin boundaries with metadata\n")

  list.files(path = dir_cln, pattern = "my_total", full.names = TRUE,
             include.dirs = TRUE, recursive = FALSE) %>%
    lapply(function(x) {

      # read data
      id <- fread(x) %>%
        bind_cols(country = x %>% basename() %>% strsplit("_") %>% unlist() %>% nth(3))

      # create admin boundaries
      # Namibia:
      # correct for the 32 pixels that are wronly assigned to 330103 instead of 330102
      ## -> does not make a difference in the mobility model because total
      ##    population size there is only 3, but it's visible in the plot
      if (unique(id$country) == "Namibia")
        id[lat < -27.5 & lon > 19, code := 330102]


      # transform into sf file
      shp <-
        # transform into raster dataset
        raster::rasterFromXYZ(id[, .(lon, lat, code)]) %>%
        # transform into sf polygons and rename columns
        st_as_stars() %>%
        st_as_sf() %>%    # this is the raster to polygons part
        # dissolve borders by tower_original or tower_merged
        group_by(code) %>%
        summarise(temp = mean(code)) %>%
        st_cast() %>%
        select(-temp)

      # compute population-weighted centroids for each admin
      id <- id[, .(pwlon    = sum(pop * lon) / sum(pop),
                   pwlat    = sum(pop * lat) / sum(pop),
                   pop_dens = sum(pop) / .N),
               by = .(country, code)]

      # merge with shp
      shp <- merge(x = shp, y = id, by = "code")

      return(shp)

    }) %>%
    bind_rows() %>%
    arrange(country, code) %>%
    saveRDS(file = file.path(dir_res, "admin_total_coords_shp.rds"))


  # pixel data of admins and pixels used as DEFAULT ---------------------------#

  cat("## Admin boundaries with metadata\n")

  list.files(path = dir_cln, pattern = "my_C", full.names = TRUE,
             include.dirs = TRUE, recursive = FALSE) %>%
    lapply(function(x) {

      # read data
      id <- fread(x) %>%
        bind_cols(country = x %>% basename() %>% strsplit("_") %>% unlist() %>% nth(3))

      # create admin boundaries
      # Namibia:
      # correct for the 32 pixels that are wronly assigned to 330103 instead of 330102
      ## -> does not make a difference in the mobility model because total
      ##    population size there is only 3, but it's visible in the plot
      if (unique(id$country) == "Namibia")
        id[lat < -27.5 & lon > 19, code := 330102]

      # compute population-weighted centroids for each admin
      id <- id[, .(pwlon    = sum(pop * lon) / sum(pop),
                   pwlat    = sum(pop * lat) / sum(pop),
                   pop_dens = sum(pop) / .N),
               by = .(country, code)]

      return(id)

    }) %>%
    bind_rows() %>%
    arrange(country, code) %>%
    saveRDS(file = file.path(dir_res, "admin_clean_coords.rds"))

}


cleanSimParameterEstimates <- function (stage, input_file, default_data) {

  cat("# Clean", stage, "\n")

  list.files(pattern = paste_("output_files", stage),
             recursive = TRUE, full.names = TRUE, include.dirs = TRUE) %>%
    list.files(pattern = input_file, recursive = TRUE, full.names = TRUE) %>%
    lapply(function(x) {

      print(x)

      # read in and add metadata
      dt <- read.table(x, header = TRUE) %>%
        bind_cols(parse_sim_path(x, default_data))


      # additional steps for LHS ----------------------------------------------#
      if (grepl("LHS", stage)) {
        dt <- dt %>%
          # find best 4 LOGL values
          arrange(-LOGL) %>%
          mutate(is.MCMCinit = if_else(row_number() > 4.5, FALSE, TRUE))
      }


      # additional steps for preMCMC and MCMC chains --------------------------#
      if (grepl("MCMC", stage)) {

        ## add chain id
        dt <- x %>%
          dirname() %>%
          basename() %>%
          str_extract_all("\\d+") %>%
          unlist() %>%
          as.factor() %>%
          bind_cols(dt, id_chain = .)

        ## symmetric gravity model: if alpha < beta, swap alpha (PAR1) and beta (PAR2)!
        ## parameters are interchangeable, given symmetry, but one will be smaller than the other
        if (dt$model[1] == "GM" & grepl("symmetric", dt$test_type[1])) {
          dt2 <- dt %>% filter(step > max(step) - 50000)
          if (mean(dt2$PAR1) < mean(dt2$PAR2)) {
            temp <- dt$PAR1
            dt$PAR1 <- dt$PAR2
            dt$PAR2 <- temp
          }
        }

        # remove burn-in
        if (stage == "MCMC")
          dt <- dt %>% removeBurnIn()

      }


      # additional steps for flow_model ---------------------------------------#
      if (grepl("flow_model_LOGL", input_file)) {

        dt <- dt %>%
          mutate(id_sample = as.integer(id_sample))

      }


      # additional steps for flow_model ---------------------------------------#
      if (grepl("flow_model_matrix", input_file)) {

        dt <- dt %>%
          mutate(id_sample     = as.integer(id_sample),
                 id_patch_dest = as.integer(id_patch_dest),
                 id_patch_orig = as.integer(id_patch_orig)) %>%
          # compute summary statistics
          group_by(id_patch_orig, id_patch_dest,
                   data_type, test_type, model, country, scale, notes) %>% # pool samples
          summarise(flow = mean(flow_estimate),
                    q1   = quantile(flow_estimate, 0.025),
                    q2   = quantile(flow_estimate, 0.975)) %>%
          ungroup()

      }


      # additional steps ------------------------------------------------------#

      if (!grepl("flow_model_matrix", input_file)) {
        # melt table, rename variables with formula name and as factors
        dt <- dt %>%
          pivot_longer(
            cols = matches("LOGL|PAR", ignore.case = FALSE),
            names_to = "variable",
            values_to = "value",
            values_drop_na = TRUE
          ) %>%
          mutate(variable = factor(variable, labels = getVariableNames(unique(model))))
      }


      # return ----------------------------------------------------------------#

      return(dt)

    }) %>%
    bind_rows() %>%
    as_tibble() %>%
    mutate(scale = factor(scale, levels = c("admin", "20", "10", "5"))) %>%
    saveRDS(file = file.path(dir_res, paste.(stage, "rds")))

}


# Helper functions -------------------------------------------------------------

parse_input_path <- function (x) {

  subdir <- strsplit(basename(dirname(x)), "_")[[1]]

  data.frame(
    country   = subdir[1],
    data_type = paste(subdir[-1], collapse = "_"),
    scale     = strsplit(basename(x), "_")[[1]][1])

}


parse_sim_path <- function (path, default_data) {

  splitpath <- strsplit(path, "/")[[1]]

  test_type <- splitpath %>%
    grep(pattern = "test_", x = ., value = TRUE) %>%
    sub("test_", "", x = .)

  data_type <- if (test_type %in% default_data) "DEFAULT" else test_type

  subtest_type <- splitpath %>% grep(pattern = "GM|RM", x = ., value = TRUE) %>%
    strsplit("_") %>%
    unlist()

  dt <- data.frame(
    data_type = data_type,
    test_type = test_type,
    model     = subtest_type[1],
    country   = subtest_type[2],
    scale     = subtest_type[3]
  )

  if (grepl("flow_model_5km_admin", path)) {
    dt$notes <- "admin_pars"
  } else if (grepl("flow_model_default_data", path)) {
    dt$notes <- "default_data"
  } else {
    dt$notes <- ""
  }

  return(dt)

}


# Name parameters like in the formulas
getVariableNames <- function (model) {

  v <- switch(
    model,
    "GM"       = c("kappa", "alpha", "beta", "gamma", "epsilon"),
    "RM-v0-t5" = c("kappa"),
    "RM-v1-t0" = ,
    "RM-v1-t1" = c("kappa", "eta"),
    "RM-v1-t2" = ,
    "RM-v1-t3" = c("kappa", "eta", "alpha"),
    "RM-v2-t0" = ,
    "RM-v2-t1" = c("kappa", "theta"),
    "RM-v2-t2" = ,
    "RM-v2-t3" = c("kappa", "theta", "alpha"),
    "RM-v0-t3" = c("kappa", "alpha"),
    "RM-v2-t4" = c("kappa", "theta"),
    "all"      = c("kappa", "alpha", "beta", "gamma", "epsilon", "eta", "theta")
  )

  c("log-likelihood", v, "disp")

}


removeBurnIn <- function (dt) {

  dt %>%
    ## split into windows of 50k steps
    mutate(bin = findInterval(step, seq.int(from = 1, to = max(step), by = 50000))) %>%
    ## compute mean log-likelihood for each window
    group_by(bin) %>%
    mutate(mean_LOGL = mean(LOGL)) %>%
    ungroup() %>%
    ## compute the minimum log-likelihood across subsequent windows
    cbind(., future_min_mean_LOGL = sapply(1:nrow(.), function(x) min(.$mean_LOGL[x:nrow(.)]))) %>%
    ## keep all iterations after the first window such that there is a subsequent window with smaller mean log-likelihood
    mutate(convergence_condition = cumsum(mean_LOGL > future_min_mean_LOGL)) %>%
    ## throw away first good window as well
    ### this is to make sure I'm throwing the very first window away
    ### the mean of the very first window might be ok, but the first values might still be off!
    filter(convergence_condition > 1) %>%
    select(-bin, -mean_LOGL, -future_min_mean_LOGL, -convergence_condition)

}


thinChains <- function (dt, n = 1000) {

  unique_steps <- dt %>%
    select(step) %>%
    unlist() %>%
    unique()

  selected_idx <- seq.int(from = 1, to = length(unique_steps), length.out = n) %>%
    as.integer() %>%
    unique()

  dt %>% filter(step %in% unique_steps[selected_idx])

}


matrix_dim <- function (nl) {

  # nl = number of elements on lower triangular

  as.integer((1 + sqrt(1 + 8 * nl)) / 2)

}


matrix_to_table <- function (m, x, values_to) {
  m %>%
    as.data.frame() %>%
    `colnames<-`(1:ncol(.) - 1L) %>%
    mutate(id_patch_orig = row_number() - 1L) %>%
    pivot_longer(-id_patch_orig, names_to = "id_patch_dest", values_to = values_to) %>%
    mutate(id_patch_dest = as.integer(id_patch_dest),
           id_patch_orig = as.integer(id_patch_orig)) %>%
    bind_cols(parse_input_path(x))
}

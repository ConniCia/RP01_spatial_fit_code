# MCMC validation --------------------------------------------------------------

fig_LHS <- function () {

  cat("# Create LHS figures\n")

  dir_out <- file.path(dir_fig, "LHS")
  suppressWarnings(dir.create(dir_out))

  readRDS(file.path(dir_res, "LHS.rds")) %>%
    group_by(test_type, model, country, scale) %>%
    group_split() %>%
    lapply(function(dt) {

      file_name <- dt %>%
        select(test_type, model, country, scale) %>%
        mutate(scale = as.character(scale)) %>%
        summarise(across(everything(), unique)) %>%
        paste(collapse = "_")

      print(file_name)

      ggplot() +
        facet_wrap(vars(variable), scales = "free", labeller = label_parsed) +
        geom_point(data = dt, aes(x = step, y = value),
                   colour = my.blue, size = 0.3) +
        geom_point(data = dt %>% filter(is.MCMCinit), aes(x = step, y = value),
                   colour = my.red, size = 1) +
        scale_y_continuous(labels = number) +
        my.ggsave(path = dir_out, filename = file_name,
                  device = "tiff",
                  width = 183, height = 150)

    })

  return()

}


fig_MCMC <- function (input_file, output_folder) {

  cat("# Create", output_folder, "figures\n")

  dir_out <- file.path(dir_fig, output_folder)
  suppressWarnings(dir.create(dir_out))

  readRDS(file.path(dir_res, input_file)) %>%
    group_by(test_type, model, country, scale) %>%
    group_split() %>%
    lapply(function(dt) {

      file_name <- dt %>%
        select(test_type, model, country, scale) %>%
        mutate(scale = as.character(scale)) %>%
        summarise(across(everything(), unique)) %>%
        paste(collapse = "_")

      print(file_name)

      p <- dt %>%
        thinChains() %>%
        ggplot() +
        facet_wrap(vars(variable), scales = "free", labeller = label_parsed)

      if (output_folder == "MCMC_posterior_density") {
        p <- p +
          geom_density(aes(x = value, colour = id_chain, fill = id_chain),
                       alpha = 0.1)
      } else {
        p <- p +
          geom_line(aes(x = step, y = value, colour = id_chain))
      }

      p +
        scale_x_continuous(labels = number) +
        scale_y_continuous(labels = number) +
        my.ggsave(path = dir_out, filename = file_name,
                  device = "tiff",
                  width = 183, height = 150)

    })

  return()

}


# main text --------------------------------------------------------------------

fig_parameter_trends <- function (do.paper_name) {

  if (do.paper_name) {
    fig_name <- "Fig1"
    device <- "tiff"
  } else {
    fig_name <- "parameter_trends"
    device <- "png"
  }

  cat("# Create graph of parameter trends for main paper text\n")


  # read in and transform -----------------------------------------------------#

  dt <- readRDS(file.path(dir_res, "MCMC.rds")) %>%
    # sort variables
    mutate(variable = factor(variable, levels = getVariableNames("all"))) %>%
    # filter
    filter(grepl("DEFAULT", test_type) &
             model %in% models_DEFAULT &
             !variable == "kappa") %>%
    # compute stats
    group_by(model, country, scale, variable) %>% # pool chains
    summarise(
      m  = mean(value),
      q1 = quantile(value, 0.025),
      q2 = quantile(value, 0.975)
    ) %>%
    ungroup() %>%
    mutate(model = if_else(model == "RM-v2-t3", "RM4", model))


  # plot ----------------------------------------------------------------------#

  # LOGL
  p_A <- ggplot(data = dt %>% filter(variable == "log-likelihood"),
                aes(x = scale, y = m, ymin = q1, ymax = q2,
                    colour = model, fill = model, group = model)) +
    geom_line() +
    geom_ribbon(alpha = 0.4, colour = NA) +
    facet_wrap(vars(country), scales = "free", nrow = 1) +
    scale_y_continuous("", labels = number) +
    my.scale_model("fill") + my.scale_model("colour")

  # RM4
  p_B <- ggplot(data = dt %>% filter(model == "RM4" & variable != "log-likelihood"),
                aes(x = scale, y = m, ymin = q1, ymax = q2,
                    colour = country, fill = country, group = country)) +
    geom_line() +
    geom_ribbon(alpha = 0.4, colour = NA) +
    facet_wrap(vars(variable), scales = "free", nrow = 1, labeller = label_parsed) +
    scale_y_continuous("", labels = number) +
    my.scale_country("fill") + my.scale_country("colour")

  # GM
  p_C <- ggplot(data = dt %>% filter(model == "GM" & variable != "log-likelihood"),
                aes(x = scale, y = m, ymin = q1, ymax = q2,
                    colour = country, fill = country, group = country)) +
    geom_line() +
    geom_ribbon(alpha = 0.4, colour = NA) +
    facet_wrap(vars(variable), scales = "free", nrow = 1, labeller = label_parsed) +
    scale_y_continuous("", labels = number_format(accuracy = 0.1)) +
    my.scale_country("fill") + my.scale_country("colour")


  # arrange plots -------------------------------------------------------------#

  p_D <- my.ggplot.legend(p_A + theme(legend.position = "bottom"))
  p_E <- my.ggplot.legend(p_B + theme(legend.position = "bottom"))

  p_A <- p_A + labs(tag = "A") + theme(legend.position = "none", plot.margin = margin(0,0,0,-1))
  p_B <- p_B + labs(tag = "B") + theme(legend.position = "none", plot.margin = margin(0,0,0,-1))
  p_C <- p_C + labs(tag = "C") + theme(legend.position = "none", plot.margin = margin(0,0,0,-1))

  lay <- matrix(c(1,1,1,1,2,2,2,2,2,2,
                  3,3,3,3,3,3,3,3,3,3,
                  4,4,4,4,4,5,5,5,5,5), nrow = 3, byrow = TRUE)

  p <- arrangeGrob(p_A, p_B, p_C, p_D, p_E,
                   layout_matrix = lay,
                   heights = c(1,1,0.2))

  # print to file
  my.ggsave(plot = p,
            path = dir_fig, filename = fig_name,
            device = device,
            width = 7.5, height = 3.5, units = "in")

}


fig_heatmap <- function (do.paper_name, filter = "best") {

  device <- if (do.paper_name) "tiff" else "png"

  if (filter == "best") {
    fig_name <- if (do.paper_name) "Fig2" else "heatmaps"
    dt_filter <- dt_best
    cat("# Create graph of heatmaps for main paper text\n")
  } else if (filter == "5km_admin") {
    fig_name <- paste_("heatmaps", filter)
    dt_filter <- data.frame(scale = "5",
                            notes = "admin_pars")
    cat("# Create graph of heatmaps for paper SI\n")
  } else if (filter == "random50") {
    fig_name <- paste_("heatmaps", filter)
    dt_filter <- data.frame(test_type = "random50admins2",
                            scale = "admin",
                            notes = "default_data")
    cat("# Create graph of heatmaps for paper SI\n")
  } else if (filter == "random15") {
    fig_name <- paste_("heatmaps", filter)
    dt_filter <- data.frame(test_type = "random15admins2",
                            scale = "admin",
                            notes = "default_data")
    cat("# Create graph of heatmaps for paper SI\n")
  }

  # prepare data --------------------------------------------------------------#

  # read in flow counts: skipdiagonal as I can't map that anyway, symmetric as my map won't be directional
  dt_data <- readRDS(file.path(dir_res, "admin_flow_data.rds")) %>%
    # filter default input data
    filter((country == "Kenya" & data_type == "DEFAULT_symmetric_skipdiagonal") |
             (country == "Namibia" & data_type == "DEFAULT")) %>%
    select(-scale, -data_type) %>%
    mutate(model = "data") %>%
    rename(flow = flow_count)

  # read in simulation estimates
  dt_estim <- readRDS(file.path(dir_res, "flow_model.rds")) %>%
    # only get best scale per country and per model
    inner_join(dt_filter) %>%
    mutate(model = if_else(model == "RM-v2-t3", "RM4", model))

  # bind tables
  dt <- bind_rows(dt_data, dt_estim) %>%
    # remove diagonal
    mutate(flow = if_else(id_patch_orig == id_patch_dest, NA_real_, flow)) #%>%
  # # normalise
  # group_by(country, model) %>%
  # mutate(max_flow = max(flow, na.rm = TRUE)) %>%
  # ungroup() %>%
  # mutate(flow = 100 * flow / max_flow)


  # plot ----------------------------------------------------------------------#

  # find range for colour bar
  dt_A <- dt %>% filter(country == "Kenya")
  rng_A <- range(dt_A$flow, na.rm = TRUE)
  dt_B <- dt %>% filter(country == "Namibia")
  rng_B <- range(dt_B$flow, na.rm = TRUE)

  # plot
  p_A <- ggplot(dt_A,
                aes(x = id_patch_orig, y = id_patch_dest, fill = flow)) +
    geom_raster() +
    facet_grid(rows = vars(country), cols = vars(model), scales = "free") +
    labs(x = "origin", y = "destination") +
    theme(aspect.ratio = 1,
          panel.grid = element_blank(),
          axis.text = element_blank(), axis.ticks = element_blank()) +
    scale_fill_viridis_c(name = "",
                         limits = rng_A,
                         trans = "log1p",
                         breaks = my_log_breaks,
                         labels = number_format(accuracy = 1),
                         option = "D")

  p_B <- ggplot(dt_B,
                aes(x = id_patch_orig, y = id_patch_dest, fill = flow)) +
    geom_raster() +
    facet_grid(rows = vars(country), cols = vars(model), scales = "free") +
    labs(x = "origin", y = "destination") +
    theme(aspect.ratio = 1,
          panel.grid = element_blank(),
          axis.text = element_blank(), axis.ticks = element_blank(),
          strip.background.x = element_blank(), strip.text.x = element_blank()) +
    scale_fill_viridis_c(name = "",
                         limits = rng_B,
                         trans = "log1p",
                         breaks = my_log_breaks,
                         labels = number_format(accuracy = 1),
                         option = "D")


  # arrange plots -------------------------------------------------------------#

  lay <- matrix(c(1,1,
                  2,3), ncol = 2, byrow = TRUE)

  p <- arrangeGrob(p_A, p_B, ggplot(), layout_matrix = lay,
                   heights = c(1.1, 1), widths = c(12, 0.51))

  # print to file
  my.ggsave(plot = p,
            path = dir_fig, filename = fig_name,
            device = device,
            width = 7.5, height = 5, units = "in")

}


fig_distance_profile <- function (do.paper_name, filter = "best") {

  device <- if (do.paper_name) "tiff" else "png"

  if (filter == "best") {
    fig_name <- if (do.paper_name) "Fig3" else "distance_profiles"
    dt_filter <- dt_best
    cat("# Create graph of distance vs. trip frequency for main paper text\n")
  } else if (filter == "5km_admin") {
    fig_name <- paste_("distance_profiles", filter)
    dt_filter <- data.frame(scale = "5", notes = "admin_pars")
    cat("# Create graph of distance vs. trip frequency for paper SI\n")
  } else if (filter == "random50") {
    fig_name <- paste_("distance_profiles", filter)
    dt_filter <- data.frame(test_type = "random50admins2",
                            scale = "admin",
                            notes = "default_data")
    cat("# Create graph of distance vs. trip frequency for paper SI\n")
  } else if (filter == "random15") {
    fig_name <- paste_("distance_profiles", filter)
    dt_filter <- data.frame(test_type = "random15admins2",
                            scale = "admin",
                            notes = "default_data")
    cat("# Create graph of distance vs. trip frequency for paper SI\n")
  }

  # prepare data --------------------------------------------------------------#

  # read in distance table
  dt_dist <- readRDS(file.path(dir_res, "admin_distance.rds")) %>%
    filter(data_type == "DEFAULT") %>%
    select(-data_type, -scale)

  # read in flow counts
  dt_data <- readRDS(file.path(dir_res, "admin_flow_data.rds")) %>%
    # filter default input data: so that models are comparable
    filter((country == "Kenya" & data_type == "DEFAULT_symmetric_skipdiagonal") |
             (country == "Namibia" & data_type == "DEFAULT_skipdiagonal")) %>%
    mutate(model = "data") %>%
    rename(flow = flow_count)

  # read in simulation estimates
  dt_estim <- readRDS(file.path(dir_res, "flow_model.rds")) %>%
    # only get best scale per country and per model
    inner_join(dt_filter) %>%
    mutate(model = if_else(model == "RM-v2-t3", "RM4", model))

  # bind tables
  dt_AB <- get_distance_vs_frequency(dt_data, dt_estim, dt_dist, filter = "all")
  dt_C  <- get_distance_vs_frequency(dt_data, dt_estim, dt_dist, filter = "largest")
  dt_D  <- get_distance_vs_frequency(dt_data, dt_estim, dt_dist, filter = "smallest")


  # plot ----------------------------------------------------------------------#

  p_A <- ggplot(dt_AB, aes(x = distance, y = cumflow, colour = model)) +
    geom_path(lwd = 1) +
    facet_wrap(vars(country), scale = "free_x") +
    scale_x_log10("Distance (km)") +
    scale_y_continuous("Cumulative trip frequency (%)") +
    my.scale_model("colour", name = "")

  p_B <- ggplot(dt_AB, aes(x = distance, y = flow, colour = model, fill = model)) +
    geom_point(size = 0.5, alpha = 0.5) +
    facet_wrap(vars(country), scale = "free_x") +
    scale_x_log10("Distance (km)") +
    scale_y_log10("Trip frequency (%)  ") +
    my.scale_model("colour", name = "") +
    theme(strip.background.x = element_blank(), strip.text.x = element_blank())

  p_C <- ggplot(dt_C, aes(x = distance, y = flow, colour = model, fill = model)) +
    geom_point(size = 1) +
    facet_wrap(vars(country), scale = "free_x") +
    scale_x_log10("Distance (km)") +
    scale_y_log10("Trip frequency (%)  ") +
    my.scale_model("colour", name = "") +
    theme(strip.background.x = element_blank(), strip.text.x = element_blank())

  p_D <- ggplot(dt_D, aes(x = distance, y = flow, colour = model, fill = model)) +
    geom_point(size = 1) +
    facet_wrap(vars(country), scale = "free_x") +
    scale_x_log10("Distance (km)") +
    scale_y_log10("Trip frequency (%)  ") +
    my.scale_model("colour", name = "") +
    theme(strip.background.x = element_blank(), strip.text.x = element_blank())


  # arrange plots -------------------------------------------------------------#

  g <- my.arrange_panels(
    list(p_A + theme(plot.margin = margin(-1,4,-1,3)),
         p_B + theme(plot.margin = margin(-1,4,-1,3)),
         p_C + theme(plot.margin = margin(-1,4,-1,3)),
         p_D + theme(plot.margin = margin(-1,4,-1,3))),
    common_legend = "bottom")
  p <- arrangeGrob(grobs = g, ncol = 1, heights = c(1.6,1,1,1,0.3))

  # Print to file
  my.ggsave(plot = p,
            path = dir_fig, filename = fig_name,
            device = device,
            width = 6.5, height = 7.75, units = "in")

}


map_population_density <- function (do.paper_name) {

  if (do.paper_name) {
    fig_name <- "Fig4"
    device <- "tiff"
  } else {
    fig_name <- "map_admin_population"
    device <- "png"
  }

  cat("# Create map of population density and admin boundaries for main paper text\n")

  # read in
  shp <- readRDS(file.path(dir_res, "admin_total_coords_shp.rds"))

  # compute common colour range
  rng <- range(shp$pop_dens)


  # plots ---------------------------------------------------------------------#

  p <- shp %>%
    group_split(country) %>%
    lapply(function(x) {

      p_x <- ggplot(x) +
        geom_sf(aes(fill = pop_dens), colour = "grey", size = 0.1) +
        coord_sf() +
        geom_point(aes(x = pwlon, y = pwlat),
                   fill = "black", size = 0.7,
                   pch = 21, colour = "grey90", stroke = 0.2) +
        scale_fill_viridis_c(name = "",
                             trans = "log1p",
                             limits = rng,
                             breaks = my_log_breaks,
                             # labels = number,
                             option = "D") +
        theme_void()

      p_x

    })


  # arrange plots -------------------------------------------------------------#

  g <- my.arrange_panels(p, common_legend = "right")
  p <- arrangeGrob(grobs = g, nrow = 1, widths = c(1, 1.36, 0.3))

  # Print to file
  my.ggsave(plot = p,
            path = dir_fig, filename = fig_name,
            device = device,
            width = 7.5, height = 3.5, units = "in")

}


# SI figures -------------------------------------------------------------------

fig_cell_population <- function (do.paper_name) {

  fig_name <- "cell_population"
  device <- "png"

  cat("# Create graph of cell population frequency for paper SI\n")


  # read in and transform -----------------------------------------------------#

  dt <- readRDS(file.path(dir_res, "patch_population.rds")) %>%
    filter(data_type == "DEFAULT") %>%
    group_by(country, scale, population) %>%
    summarise(sum_cells      = n(),
              sum_population = sum(population)) %>%
    ungroup() %>%
    arrange(population) %>%
    group_by(country, scale) %>%
    mutate(prop_cells      = 100 * sum_cells      / sum(sum_cells),
           prop_population = 100 * sum_population / sum(sum_population),
           cumsum_prop_cells      = cumsum(prop_cells),
           cumsum_prop_population = cumsum(prop_population)) %>%
    ungroup()


  # plot ----------------------------------------------------------------------#

  p_A <- ggplot(dt, aes(x = population, y = cumsum_prop_cells,
                        colour = country)) +
    geom_point(alpha = 0.5) +
    facet_wrap(vars(scale), nrow = 1) +
    scale_x_log10("Cell population", labels = number) +
    ylab("Cumulative cell\nfrequency (%)") +
    my.scale_country("colour", name = "") +
    theme(plot.margin = margin(-2,2,-1,2))

  p_B <- ggplot(dt, aes(x = population, y = cumsum_prop_population,
                        colour = country)) +
    geom_point(alpha = 0.5) +
    facet_wrap(vars(scale), nrow = 1) +
    scale_x_log10("Cell population", labels = number) +
    ylab("Cumulative proportion\nof total population (%)") +
    my.scale_country("colour", name = "") +
    theme(strip.background.x = element_blank(), strip.text.x = element_blank(),
          plot.margin = margin(-1,2,-2,2))


  # arrange plots -------------------------------------------------------------#

  g <- my.arrange_panels(list(p_A, p_B), common_legend = "bottom")
  p <- arrangeGrob(grobs = g, ncol = 1, heights = c(1,1,0.3))

  # Print to file
  my.ggsave(plot = p,
            path = dir_fig, filename = fig_name,
            device = device,
            width = 7.5, height = 4, units = "in")

}


map_connectivity <- function (do.paper_name) {

  fig_name <- "map_connectivity"
  device <- "png"

  cat("# Create connectivity maps for paper SI\n")

  # read in and transform -----------------------------------------------------#

  # shp shapefile
  shp <- readRDS(file.path(dir_res, "admin_total_coords_shp.rds")) %>%
    # admin population-weighted centroid coordinates
    arrange(code) %>%
    group_by(country) %>%
    mutate(id_patch = row_number() - 1) %>%
    ungroup() %>%
    select(-code, -pop_dens)

  # read in flow counts: skipdiagonal as I can't map that anyway, symmetric as my map won't be directional
  dt_data <- readRDS(file.path(dir_res, "admin_flow_data.rds")) %>%
    # filter default input data
    filter(data_type == "DEFAULT_symmetric_skipdiagonal") %>%
    select(-scale, -data_type) %>%
    mutate(model = "data") %>%
    rename(flow = flow_count) %>%
    # filter only the upper triangular
    filter(id_patch_orig < id_patch_dest)

  # read in simulation estimates
  dt_estim <- readRDS(file.path(dir_res, "flow_model.rds")) %>%
    # filter default input data
    filter(grepl("DEFAULT", test_type)) %>%
    # only get best scale per country and per model
    inner_join(dt_best) %>%
    # make symmetric if it's not already
    group_split(model, country) %>%
    lapply(function(x) {

      # make symmetric
      if (!grepl("symmetric", unique(x$test_type))) {
        x <- inner_join(x, x,
                        by = c("id_patch_orig" = "id_patch_dest",
                               "id_patch_dest" = "id_patch_orig",
                               "country", "model", "test_type")) %>%
          mutate(flow = flow.x + flow.y)
      }

      # filter only the upper triangular
      x %>%
        filter(id_patch_orig < id_patch_dest) %>%
        select(country, model, id_patch_orig, id_patch_dest, flow)

    }) %>%
    bind_rows() %>%
    mutate(model = if_else(model == "RM-v2-t3", "RM4", model))


  # combine tables ------------------------------------------------------------#

  # merge with shp
  dt <- bind_rows(dt_data, dt_estim) %>%
    group_by(country) %>%
    mutate(max_flow = max(flow)) %>%
    ungroup() %>%
    mutate(flow = 100 * flow / max_flow) %>%
    # add admin coordinates
    inner_join(y = st_drop_geometry(shp),
               by = c("country", "id_patch_orig" = "id_patch")) %>%
    inner_join(y = st_drop_geometry(shp),
               by = c("country", "id_patch_dest" = "id_patch"),
               suffix = c("_orig", "_dest")) %>%
    # make sure stronger links are plotted on top
    arrange(flow)


  # plot ----------------------------------------------------------------------#

  min_flow <- 0.01
  p_A <- ggplot() +
    geom_sf(data = shp %>% filter(country == "Kenya"), size = 0.4) +
    geom_segment(data = dt %>% filter(country == "Kenya" & flow > min_flow),
                 aes(x    = pwlon_orig, y    = pwlat_orig,
                     xend = pwlon_dest, yend = pwlat_dest,
                     color = flow),
                 alpha = 1) +
    facet_grid(rows = vars(country), cols = vars(model)) +
    scale_colour_viridis_c(name = "",
                           trans = "log1p",
                           breaks = my_log_breaks,
                           limits = c(0, 100),
                           # labels = number,
                           direction = 1,
                           option = "D") +
    scale_size_identity() +
    theme(panel.grid = element_blank(), axis.title = element_blank(),
          axis.text = element_blank(), axis.ticks = element_blank())

  min_flow <- 0.1
  p_B <- ggplot() +
    geom_sf(data = shp %>% filter(country == "Namibia"), size = 0.4) +
    geom_segment(data = dt %>% filter(country == "Namibia" & flow > min_flow),
                 aes(x    = pwlon_orig, y    = pwlat_orig,
                     xend = pwlon_dest, yend = pwlat_dest,
                     color = flow),
                 alpha = 1) +
    facet_grid(rows = vars(country), cols = vars(model)) +
    scale_colour_viridis_c(name = "",
                           trans = "log1p",
                           breaks = my_log_breaks,
                           limits = c(0, 100),
                           # labels = number,
                           direction = 1,
                           option = "D") +
    scale_size_identity() +
    theme(panel.grid = element_blank(), axis.title = element_blank(),
          axis.text = element_blank(), axis.ticks = element_blank(),
          strip.background.x = element_blank(), strip.text.x = element_blank())


  # arrange plots -------------------------------------------------------------#

  g <- my.arrange_panels(list(p_A, p_B), common_legend = "right", add_tag = FALSE)

  lay <- matrix(c(1, 3,
                  2, 3), nrow = 2, byrow = TRUE)
  p <- arrangeGrob(grobs = g, layout_matrix = lay,
                   widths = c(3, 0.3), heights = c(1.425, 1))

  # Print to file
  my.ggsave(plot = p,
            path = dir_fig, filename = fig_name,
            device = device,
            width = 7.5, height = 5, units = "in")

}


# SI methods -------------------------------------------------------------------

fig_symmetricity <- function (do.paper_name) {

  fig_name <- "symmetricity"
  device <- "png"

  cat("# Create graph of data symetricity for paper SI\n")

  # read in and transform -----------------------------------------------------#
  dt <- readRDS(file.path(dir_res, "admin_flow_data.rds")) %>%
    # filter default input data
    filter(data_type == "DEFAULT") %>%
    select(-data_type, -scale) %>%
    # join with itself but symmetric
    inner_join(x = ., y = .,
               by = c("id_patch_orig" = "id_patch_dest",
                      "id_patch_dest" = "id_patch_orig",
                      "country")) %>%
    group_by(country, id_patch_orig, id_patch_dest, flow_count.x, flow_count.y) %>%
    summarise(max_count = max(flow_count.x, flow_count.y),
              residual  = flow_count.x - flow_count.y,
              prop_res  = 100 * residual / max_count) %>%
    ungroup() %>%
    mutate(prop_res = if_else(is.na(prop_res), 0, prop_res))


  # plot ----------------------------------------------------------------------#

  p_A <- ggplot(dt %>% filter(country == "Kenya"),
                aes(x = id_patch_orig, y = id_patch_dest, fill = prop_res)) +
    geom_raster() +
    labs(x = "origin", y = "destination") +
    scale_colour_gradient2(name = "",
                           low = "#313695",
                           mid = "#f5f095",
                           high = "#a62f26",
                           midpoint = 0,
                           limits = c(-100, 100),
                           space = "Lab",
                           na.value = "grey50",
                           guide = "colourbar",
                           aesthetics = "fill") +
    theme(aspect.ratio = 1,
          panel.grid = element_blank(),
          axis.text = element_blank(), axis.ticks = element_blank())

  p_B <- ggplot(dt %>% filter(country == "Namibia"),
                aes(x = id_patch_orig, y = id_patch_dest, fill = prop_res)) +
    geom_raster() +
    labs(x = "origin", y = "destination") +
    scale_colour_gradient2(name = "",
                           low = "#313695",
                           mid = "#f5f095",
                           high = "#a62f26",
                           midpoint = 0,
                           limits = c(-100, 100),
                           space = "Lab",
                           na.value = "grey50",
                           guide = "colourbar",
                           aesthetics = "fill") +
    theme(aspect.ratio = 1,
          panel.grid = element_blank(),
          axis.text = element_blank(), axis.ticks = element_blank())


  # arrange panels ------------------------------------------------------------#

  g <- my.arrange_panels(list(p_A, p_B), common_legend = "right")
  p <- arrangeGrob(grobs = g, nrow = 1, widths = c(1, 1, 0.3))

  # print to file
  my.ggsave(plot = p,
            path = dir_fig, filename = fig_name,
            device = device,
            width = 7.5, height = 3.5, units = "in")

  dt %>%
    group_by(country) %>%
    summarise(mean_abs_prop = mean(abs(prop_res)),
              q1_abs_prop   = quantile(abs(prop_res), 0.025),
              q2_abs_prop   = quantile(abs(prop_res), 0.975)) %>%
    ungroup() %>%
    xtable(type = "latex") %>%
    print(file = file.path(dir_tbl, "symmetricity.tex"))

}


fig_parameter_correlation <- function (do.paper_name) {

  fig_name <- "parameter_correlation"
  device <- "png"

  cat("# Create graph of correlation between GM parameters for Kenya paper SI\n")

  # read in
  dt <- readRDS(file.path(dir_res, "MCMC.rds")) %>%
    filter(country == "Kenya" & model == "GM" &
             test_type %in% c("DEFAULT_symmetric_skipdiagonal_rescalekappa",
                              "symmetric_skipdiagonal")) %>%
    thinChains(n = 100) %>% # variable number of points due to burnin removal
    pivot_wider(names_from = variable, values_from = value)

  DT <- dt %>%
    mutate(test_type = factor(test_type,
                              levels = c("DEFAULT_symmetric_skipdiagonal_rescalekappa",
                                         "symmetric_skipdiagonal"),
                              labels = c(TeX("GM factor $10^{\\kappa - \\gamma \\epsilon}$ (default)"),
                                         TeX("GM factor $10^{\\kappa}$"))))

  # plot
  p <- ggplot(DT, aes(x = kappa, y = gamma)) +
    geom_point(size = 0.8) +
    facet_grid(rows = vars(scale), cols = vars(test_type),
               scales = "free", labeller = label_parsed) +
    labs(x = TeX("$\\kappa$"), y = TeX("$\\gamma$"))

  # print to file
  my.ggsave(plot = p,
            path = dir_fig, filename = fig_name,
            device = device,
            width = 4, height = 3.5, units = "in")

}


fig_parameter_interpretation <- function (do.paper_name) {

  fig_name <- "parameter_interpretation"
  device <- "png"

  cat("# Create graph of parameter interpretation for paper SI\n")


  # define functions ----------------------------------------------------------#

  # exponentiation
  fun_A <- function (alpha, population) {
    data.frame(value = population^alpha)
  }

  # augmented population in RM
  fun_B <- function (theta, radius_pop, dest_pop) {
    data.frame(value = theta / ((theta + radius_pop + dest_pop) * (theta + radius_pop)))
  }

  # distance kernel
  fun_C <- function (gamma, epsilon, distance) {
    epsilon <- as.numeric(epsilon)
    data.frame(value = 1 / (1 + (distance / (10^gamma))^epsilon))
  }

  # # negative binomial distribution
  # fun_D <- function (r, k, p) {
  #   data.frame(value = gamma(k + r) / (gamma(k + 1) * gamma(r)) * p^k * (1 - p)^r)
  # }
  #
  # # poisson distribution
  # fun_E <- function (lambda, k) {
  #   data.frame(value = lambda^k * exp(-lambda) / gamma(k + 1))
  # }


  my_alpha <- function(x) {
    TeX(paste("$\\alpha$ =", x))
  }

  my_theta <- function(x) {
    TeX(paste("$\\theta$ =", x))
  }

  my_gamma <- function(x) {
    TeX(paste("$\\gamma$ =", x))
  }


  # plot ----------------------------------------------------------------------#

  p_A <- tidyr::crossing(alpha = c(0.5, 0.75, 1, 1.5, 2),
                         population = seq(0, 20, length.out = 100)) %>%
    bind_cols(purrr::pmap_df(., fun_A)) %>%
    ggplot(aes(population, value, colour = factor(alpha))) +
    geom_line() +
    ylim(0, 20) +
    scale_x_continuous("Population") +
    scale_color_discrete("", labels = my_alpha) +
    theme(legend.text.align = 0)

  p_B <- tidyr::crossing(theta = c(25, 50, 100),
                         radius_pop = seq(0, 300, length.out = 100),
                         dest_pop   = 50) %>%
    bind_cols(purrr::pmap_df(., fun_B)) %>%
    ggplot(aes(radius_pop, value, colour = factor(theta))) +
    geom_line() +
    scale_x_continuous("Radius population") +
    scale_color_discrete("", labels = my_theta) +
    theme(legend.text.align = 0)

  p_C <- tidyr::crossing(gamma = c(0.5, 1, 2, 3),
                         epsilon = factor(c("0.5", "1", "2"),
                                          c(TeX("$\\epsilon$ = 0.5"),
                                            TeX("$\\epsilon$ = 1"),
                                            TeX("$\\epsilon$ = 2"))),
                         distance = seq(0, 500, length.out = 100)) %>%
    bind_cols(purrr::pmap_df(., fun_C)) %>%
    ggplot(aes(distance, value, colour = factor(gamma))) +
    facet_wrap(vars(epsilon), scales = "free", ncol = 1, labeller = label_parsed) +
    geom_line() +
    scale_x_continuous("Distance") +
    scale_color_discrete("", labels = my_gamma) +
    theme(legend.text.align = 0)

  # p_D <- tidyr::crossing(r = c(0.5, 2, 5, 10),
  #                        k = seq(0, 50, length.out = 100),
  #                        p = c(0, 0.25, 0.5, 0.75)) %>%
  #   bind_cols(purrr::pmap_df(., fun_D)) %>%
  #   ggplot(aes(k, value, colour = factor(p))) +
  #   facet_wrap(vars(r), scales = "free") +
  #   geom_line()
  #
  # p_E <- tidyr::crossing(lambda = c(0.5, 2, 5, 10),
  #                        k = seq(0, 50, length.out = 100)) %>%
  #   bind_cols(purrr::pmap_df(., fun_E)) %>%
  #   ggplot(aes(k, value, colour = factor(lambda))) +
  #   geom_line()


  # arrange panels ------------------------------------------------------------#

  p_A <- p_A + labs(tag = "A")
  p_B <- p_B + labs(tag = "B")
  p_C <- p_C + labs(tag = "C")

  lay <- matrix(c(1,3,
                  2,3), nrow = 2, byrow = TRUE)
  p <- arrangeGrob(p_A, p_B, p_C, layout_matrix = lay)

  # print to file
  my.ggsave(plot = p,
            path = dir_fig, filename = fig_name,
            device = device,
            width = 7.5, height = 4, units = "in")

}

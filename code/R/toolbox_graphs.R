
# so-called magic number, multiply with font size to get real size (like in Word)
pt_size_factor <- 0.352777778


# Global graph options ---------------------------------------------------------

theme_set(
  theme_light() +
    theme(
      strip.background = element_rect(fill = "lightgrey", colour = NULL),
      panel.border = element_rect(fill = NA, colour = NA),
      axis.ticks = element_line(colour = "black"),
      strip.text = element_text(colour = "black")
    )
)


# Colours ----------------------------------------------------------------------

my.black  <- "#000000"
my.red    <- "#F8766D"
my.blue   <- "#00BFC4"
my.orange <- "#E7861B"
my.purple <- "#C77CFF"


# black, blue, red
data_colour <- c("data" = my.black,     "GM" = my.blue, "RM4" = my.red)
data_labels <- c("data" = "Input data", "GM" = "GM",    "RM4" = "RM4")

# purple, orange
country_colour <- c("Kenya" = my.orange, "Namibia" = my.purple)
country_labels <- c("Kenya" = "Kenya",   "Namibia" = "Namibia")


# Discrete colour and/or fill scales -------------------------------------------

my.scale_country <- function (aesthetics, name = "Country", ...) {
  scale_colour_manual(
    aesthetics = aesthetics,
    name = name,
    values = country_colour, labels = country_labels,
    ...
  )
}

my.scale_model <- function (aesthetics, name = "Model", ...) {
  scale_colour_manual(
    aesthetics = aesthetics,
    name = name,
    values = data_colour, labels = data_labels,
    ...
  )
}


# graph helper functions -------------------------------------------------------

my.heatmap <- function (dt, rng, scale, my.fill) {

  p <- ggplot(
    dt, aes_string(x = "id_patch_orig", y = "id_patch_dest", fill = my.fill)) +
    geom_raster() +
    labs(x = "origin", y = "destination") +
    theme(aspect.ratio = 1,
          panel.grid = element_blank(),
          axis.text=element_blank(), axis.ticks=element_blank())

  if (scale == "discrete") {
    p <- p +
      scale_fill_viridis_d(name = "",
                           direction = 1,
                           option = "A")
  } else if (scale == "log1p") {
    p <- p +
      scale_fill_viridis_c(name = "",
                           limits = rng,
                           trans = "log1p",
                           breaks = my_log_breaks,
                           labels = number,
                           direction = 1,
                           option = "A")
  } else {
    p <- p +
      scale_colour_gradient2(
        name = "",
        low = "#313695",
        mid = "#f5f095",
        high = "#a62f26",
        midpoint = 0,
        limits = rng,
        space = "Lab",
        na.value = "grey50",
        guide = "colourbar",
        aesthetics = "fill"
      )
  }

  return(p)

}


my.violin <- function (dt, my.y) {

  ggplot(
    dt,
    aes_string(x = "within95CrI", y = my.y,
               fill = "within95CrI", group = "within95CrI")
  ) +
    geom_violin() +
    geom_boxplot(outlier.shape=NA, width = 0.2, show.legend = FALSE) +
    scale_y_continuous(labels = number_format)

}


get_distance_vs_frequency <- function (dt_data, dt_estim, dt_dist, filter = NULL) {

  # find origin largest or smallest flow (same result if we look at destination instead!)
  if (filter == "largest") {
    dt_filter <- dt_data %>%
      group_by(country, id_patch_orig) %>%
      summarise(sum_flow = sum(flow)) %>%
      slice_max(sum_flow) %>%
      ungroup() %>%
      select(-sum_flow)
  } else if (filter == "smallest") {
    dt_filter <- dt_data %>%
      group_by(country, id_patch_orig) %>%
      summarise(sum_flow = sum(flow)) %>%
      slice_min(sum_flow) %>%
      ungroup() %>%
      select(-sum_flow)
  } else {
    dt_filter <- dt_data %>% select(country, id_patch_orig) %>% unique()
  }


  dt <- bind_rows(dt_data, dt_estim) %>%
    inner_join(dt_dist) %>%
    # filter by origin
    inner_join(dt_filter) %>%
    # skip when distance is 0
    filter(distance > 0) %>%
    # if symmetric: take only upper triangular
    group_split(data_type, test_type) %>%
    lapply(function(x) {
      type <- x %>% summarise(paste0(unique(data_type), unique(test_type))) %>% unlist()
      if (grepl("symmetric", type)) {
        x <- x %>% filter(id_patch_orig < id_patch_dest)
      }
      return(x)
    }) %>%
    bind_rows() %>%
    # mean by distance
    group_by(country, model, distance) %>%
    summarise(flow = mean(flow)) %>%
    ungroup()

  if (!filter %in% c("largest", "smallest")) {
    dt <- dt %>%
      # compute cumsum and normalise
      arrange(distance) %>%
      group_by(country, model) %>%
      mutate(cumflow = 100 * cumsum(flow) / sum(flow),
             flow    = 100 * flow / sum(flow)) %>%
      ungroup()
  } else {
    dt <- dt %>%
      # normalise
      arrange(distance) %>%
      group_by(country, model) %>%
      mutate(flow = 100 * flow / sum(flow)) %>%
      ungroup()
  }

  dt <- dt %>%
    arrange(model)

  return(dt)

}

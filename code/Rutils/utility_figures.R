# Extract ggplot legend out of plot --------------------------------------------

my.ggplot.legend <- function(a.gplot) {

  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend

}


# Arrange plots in regular grid ------------------------------------------------

my.arrange_panels <- function (plotlist, add_tag = TRUE,
                               remove_panel_legend = TRUE,
                               common_legend = "none") {

  # get legend
  if (common_legend != "none")
    p.L <- my.ggplot.legend(plotlist[[1]] + theme(legend.position = common_legend))

  # add tag and remove legends
  plotlist <- lapply(seq_along(plotlist), function(i) {
    p <- plotlist[[i]]

    if (add_tag)
      p <- p + labs(tag = LETTERS[i])

    if (remove_panel_legend)
      p <- p + theme(legend.position = "none")

    return(p)
  })

  # convert to gtables
  grobs <- lapply(plotlist, ggplotGrob)

  # get width
  maxWidth <- lapply(grobs, function(g) g$widths[2:5])
  maxWidth <- do.call(unit.pmax, maxWidth)

  # apply width
  for(i in 1:length(grobs))
    grobs[[i]]$widths[2:5] <- maxWidth

  # add legend
  if (common_legend != "none")
    grobs[[length(grobs) + 1]] <- p.L

  return(grobs)

}


# Simplified ggsave function ---------------------------------------------------

my.ggsave <- function (plot = last_plot(), filename, path = NULL, width, height,
                       units = "mm", device = "tiff", dpi = 600,
                       type = "cairo", compression = "lzw") {

  if (device == "tiff") {
    ggsave(plot,
           path = path, filename = paste(filename, device, sep = "."),
           width = width, height = height, units = units,
           device = device, dpi = dpi,
           type = "cairo", compression = "lzw")
  } else {
    ggsave(plot,
           path = path, filename = paste(filename, device, sep = "."),
           width = width, height = height, units = units,
           device = device, dpi = dpi)
  }

}


# Convenient breaks for colour scales ------------------------------------------

my_log_breaks <- function(limits) {
  as.integer(10^(0:9))
}

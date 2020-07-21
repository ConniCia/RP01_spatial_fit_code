# useful to drop columns that only contain NAs
# not_all_na <- function(x) {!all(is.na(x))}
#
#
# fig_CrI <- function () {
#
#   # read in flow data and estimates and filter by above models and scales
#   DT <- compareFlowEstimates() %>%
#     inner_join(dt_best) %>%
#     inner_join(dt_dist) %>%
#     inner_join(dt_pop, by = c("country", "id_patch_orig" = "id_patch")) %>%
#     rename(population_orig = population) %>%
#     inner_join(dt_pop, by = c("country", "id_patch_dest" = "id_patch")) %>%
#     rename(population_dest = population)
#
#   # heatmap IQR -------------------------------------------------------------#
#
#   p_CrI_GM  <- my.heatmap(dt %>% filter(model == "GM"),
#                           rng = NULL, scale = "discrete",
#                           my.fill = "within95CrI")
#   p_CrI_RM4 <- my.heatmap(dt %>% filter(model == "RM-v2-t3"),
#                           rng = NULL, scale = "discrete",
#                           my.fill = "within95CrI")
#
#
#   # boxplot with scatter diagram --------------------------------------------#
#
#   p_dist_GM  <- my.violin(dt %>% filter(model == "GM" & !is.na(within95CrI)),
#                           my.y = "distance")
#   p_dist_RM4 <- my.violin(dt %>% filter(model == "RM-v2-t3" & !is.na(within95CrI)),
#                           my.y = "distance")
#
#   p_pop_orig_GM  <- my.violin(dt %>% filter(model == "GM" & !is.na(within95CrI)),
#                               my.y = "population_orig")
#   p_pop_orig_RM4 <- my.violin(dt %>% filter(model == "RM-v2-t3" & !is.na(within95CrI)),
#                               my.y = "population_orig")
#
#   p_pop_dest_GM  <- my.violin(dt %>% filter(model == "GM" & !is.na(within95CrI)),
#                               my.y = "population_dest")
#   p_pop_dest_RM4 <- my.violin(dt %>% filter(model == "RM-v2-t3" & !is.na(within95CrI)),
#                               my.y = "population_dest")
#
# }

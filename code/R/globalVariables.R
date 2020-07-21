# Global options and variables =================================================

# global options
options(stringsAsFactors = FALSE, scipen = 5)


# standard paths ---------------------------------------------------------------

## work on Windows cluster or on personal Mac?
do.cluster <- if (file.exists("/Volumes/Mobility/")) TRUE else FALSE

root <- if (do.cluster) "/Volumes/Mobility" else "."
## code/
dir_code <- "code/"
dir_C    <- "code/C/"
dir_R    <- "code/R/"
dir_Ruti <- "code/Rutils/"
## data/
dir_raw  <- file.path(root, "data/raw/")
dir_cln  <- file.path(root, "data/clean/")
dir_proc <- file.path(root, "data/processing/")
dir_res  <- file.path(root, "data/results/")
dir_tbl  <- file.path(root, "data/tables/")
dir_fig  <- file.path(root, "data/figures/")


# Other paths ------------------------------------------------------------------

root_win <- "\\\\fi--dideclusthn.dide.ic.ac.uk\\Mobility\\"
dir_cln_win  <- paste0(root_win, "data\\clean\\")
dir_proc_win <- paste0(root_win, "data\\processing\\")
path_exec_win <- paste0(dir_proc_win, "execRP01.exe")
path_exec_mac <- paste0("./", dir_proc, "/execRP01")


# Analyses ---------------------------------------------------------------------

grid_changes <- c(
  "shiftlon", "shiftlat", "shiftlonlat",
  "random50admins1", "random50admins2", # take two separate samples
  "random30admins1", "random30admins2",
  "random15admins1", "random15admins2"
)
grid_changes_Namibia <- c(
  "hole", "flow", "noNE", "onlydenseN",
  "hole_flow", "hole_ne", "hole_onlydenseN", "flow_ne", "flow_onlydenseN",
  "hole_flow_ne", "hole_flow_onlydenseN"
)
flag_changes <- read.table(header = TRUE, text = "
  country do_symmetric do_skipdiagonal do_rescalekappa data_type test_type
  Kenya              1               1               1 DEFAULT   DEFAULT_symmetric_skipdiagonal_rescalekappa
  Kenya              1               1               0 DEFAULT   symmetric_skipdiagonal
  Namibia            0               0               0 DEFAULT   DEFAULT
  Namibia            0               0               1 DEFAULT   rescalekappa
  Namibia            0               1               0 DEFAULT   skipdiagonal
  Namibia            0               1               1 DEFAULT   skipdiagnal_rescalekappa
  Namibia            1               0               0 DEFAULT   symmetric
  Namibia            1               0               1 DEFAULT   symmetric_rescalekappa
  Namibia            1               1               0 DEFAULT   symmetric_skipdiagonal
  Namibia            1               1               1 DEFAULT   symmetric_skipdiagonal_rescalekappa")

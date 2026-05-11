# R/globals.R
# Suppress R CMD check notes for non-standard evaluation (dplyr/ggplot2)
# and global environment variables used by the LinkDeplete workflow

utils::globalVariables(c(
  # ggplot2 / dplyr NSE columns
  "season", "month", "year", "catch", "season_step", "series",
  "obs", "pred", "se", "scale_mu", "in_nll",
  "est", "lo", "hi",
  "xmin", "xmax", "ymin", "ymax",
  "D", "D_lo", "D_hi",
  "B_ratio", "F_ratio",
  "B0avg", "Bavg", "Bend",
  "R", "R_next",

  # Global environment objects (<<- workflow)
  "mod", "mout", "Datain",
  "dat_monthly", "seasons", "S", "TS",
  "N_cpue", "catch_mat",
  "cpue_mats", "cpue_se_mats", "cpue_specs",

  # RTMB getAll() unpacked parameters and data
  "log_B1", "log_R", "log_sigma_step",
  "M", "rec_ts", "ref_ts", "min_pts"
))

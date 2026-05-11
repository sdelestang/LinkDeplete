#' Set season structure from monthly data
#'
#' Assigns each row of a monthly data frame to a fishing season based on
#' a user-defined season start month, computes the within-season time step
#' index (\code{season_step}), and converts catch from kg to tonnes. Sets
#' \code{seasons}, \code{S}, and \code{TS} in the global environment.
#'
#' @param Dat Data frame with columns \code{year}, \code{month}, and
#'   \code{catch} (in kg).
#' @param SeasonStart Integer month (1--12) at which the fishing season
#'   begins. Default 1 (January).
#' @param FirstSeason Integer. Earliest season to retain after assigning
#'   season labels. Default 1970.
#'
#' @return A data frame (\code{dat_monthly}) with added columns
#'   \code{season} and \code{season_step}, and \code{catch} in tonnes.
#'   Also assigns \code{seasons} (sorted unique season vector),
#'   \code{S} (number of seasons), and \code{TS} (12) to the global
#'   environment via \code{<<-}.
#'
#' @import dplyr
#' @importFrom magrittr %>% %<>%
#'
#' @export
SetSeasons <- function(Dat = NULL, SeasonStart = 1, FirstSeason = 1970) {
  dat_monthly <- Dat %>%
    mutate(
      season      = if_else(month >= SeasonStart, year, year - (12-SeasonStart)),
      season_step = if_else(month >= SeasonStart, month - (SeasonStart-1), month + (12-((SeasonStart-1)))),
      catch       = catch / 1000
    ) %>% filter(season >= FirstSeason)
  seasons <<- sort(unique(dat_monthly$season))
  S  <<- length(seasons)
  TS <<- 12
  return(dat_monthly)
  }


#' Build the catch matrix
#'
#' Converts the \code{catch} column of \code{dat_monthly} into an
#' \code{S x TS} matrix using \code{\link{make_mat}} and replaces
#' \code{NA} cells with zero. Assigns \code{catch_mat} to the global
#' environment.
#'
#' @return Called for side effect. Assigns \code{catch_mat} (\code{S x TS}
#'   numeric matrix) to the global environment via \code{<<-}.
#'
#' @details
#' Requires \code{dat_monthly}, \code{seasons}, \code{S}, and \code{TS}
#' in the calling environment (typically set by \code{\link{SetSeasons}}).
#'
#' @seealso \code{\link{make_mat}}, \code{\link{SetSeasons}}
#'
#' @export
SetCatchMatrix <- function() {
  catch_mat <<- make_mat("catch")
  catch_mat[is.na(catch_mat)] <<- 0
}


#' Build CPUE and standard error matrices
#'
#' Detects the number of CPUE series from columns in \code{dat_monthly},
#' then builds \code{S x TS} matrices of observed CPUE and associated
#' standard errors for each series. Assigns \code{N_cpue},
#' \code{cpue_mats}, and \code{cpue_se_mats} to the global environment.
#'
#' @return Called for side effect. Assigns to the global environment:
#'   \describe{
#'     \item{\code{N_cpue}}{Integer number of CPUE series detected.}
#'     \item{\code{cpue_mats}}{List of \code{S x TS} CPUE matrices.}
#'     \item{\code{cpue_se_mats}}{List of \code{S x TS} standard error
#'       matrices.}
#'   }
#'
#' @details
#' Requires \code{dat_monthly} and \code{cpue_specs} in the calling
#' environment. The number of series is inferred by counting columns
#' matching \code{"cpue"} in \code{dat_monthly} and dividing by two
#' (one value column and one SE column per series).
#'
#' @seealso \code{\link{make_mat}}, \code{\link{SetSeasons}}
#'
#' @export
SetCPUEMatrix <- function() {
  N_cpue <<- sum(grepl("cpue", names(dat_monthly)))/2
  cpue_specs <- cpue_specs[1:N_cpue]
  cpue_mats  <<- lapply(cpue_specs, function(s) make_mat(s$col))
  cpue_se_mats <<- lapply(cpue_specs, function(s) make_mat(s$se_col))
}


#' Assemble the TMB/RTMB data input list
#'
#' Combines catch, CPUE, and natural mortality information into the
#' named list required by \code{\link{DepleteModel}}. CPUE matrices,
#' standard errors, depletion step vectors, and depletion flags are
#' appended dynamically for each series.
#'
#' @param M Numeric. Annual instantaneous natural mortality. Divided
#'   internally by \code{TS} to give per-time-step mortality. Default 0.5.
#' @param rec_ts Integer. Within-season time step at which recruitment
#'   enters. Default 3.
#' @param ref_ts Integer. Reference time step for mid-season biomass
#'   reporting (e.g. 6 for June in a Nov-start season). Default 6.
#' @param ref_levels Numeric vector of length 3 giving target, threshold,
#'   and limit depletion reference points as proportions of B0.
#'   Default \code{c(0.4, 0.3, 0.2)}.
#'
#' @return A named list suitable for passing to \code{\link{DepleteModel}}
#'   via \code{getAll()}.
#'
#' @details
#' Requires \code{catch_mat}, \code{N_cpue}, \code{S}, \code{TS},
#' \code{cpue_mats}, \code{cpue_se_mats}, and \code{cpue_specs} in the
#' calling environment.
#'
#' @seealso \code{\link{DepleteModel}}, \code{\link{SetCatchMatrix}},
#'   \code{\link{SetCPUEMatrix}}
#'
#' @export
MakeDataIn <- function(M = 0.5, rec_ts = 3, ref_ts = 6, ref_levels = c(0.4, 0.3, 0.2)) {
  ## ── Data input list ──────────────────────────────────────────────────────────
  dat <- list(
    catch   = catch_mat,           # catch matrix
    N_cpue  = N_cpue,              # number of cpue matrix
    M       = M / TS,            # Natural mort / number of timesteps - assumes even temporal space in time between steps
    min_pts = 3,                   # Dont do deplete if less than this many observations
    S       = S,                   # Number of seasons
    TS      = TS,                  # Number of Timesteps
    rec_ts = rec_ts,                    # Recruitment timestep
    ref_ts = ref_ts,                    # Timestep for reference comparison e.g. mid season 8 = June
    ref_levels = ref_levels  # Reference levels for Target, Threshold and Limit
  )

  ## Auto load the variuse number of cpues into dat ##
  for (k in seq_len(N_cpue)) {
    dat[[paste0("c", k)]]             <- cpue_mats[[k]]
    dat[[paste0("c", k, "se")]]       <- cpue_se_mats[[k]]
    dat[[paste0("cpue", k, "_depl")]] <- cpue_specs[[k]]$depl
    dat[[paste0("cpue_is_depl_", k)]] <- as.integer(cpue_specs[[k]]$is_depl)
  }

  return(dat)
}

#' Build starting parameter list from Leslie regressions
#'
#' Runs \code{\link{RunLeslie}} to obtain per-season B0 estimates,
#' interpolates gaps, and constructs the named parameter list required
#' by \code{\link{DepleteModel}}. Initial recruitment values are set
#' from Leslie B0 estimates, and catchability deviations are initialised
#' from mean Leslie slopes for depletion series.
#'
#' @return A named list of starting parameter values:
#'   \describe{
#'     \item{\code{log_B1}}{Log initial biomass (first Leslie B0).}
#'     \item{\code{log_R}}{Log recruitment vector (length \code{S - 1}).}
#'     \item{\code{log_sigma_step}}{Log recruitment random walk SD
#'       (starting value 0.5).}
#'     \item{\code{log_q_cpueK_dev}}{Log catchability deviation vector
#'       for each depletion series \emph{K}.}
#'     \item{\code{log_sigma_qK}}{Log catchability random walk SD for
#'       each depletion series \emph{K} (starting value 0.2).}
#'     \item{\code{log_q_cpueK}}{Log fixed catchability for each
#'       non-depletion series \emph{K} (starting value 0.003).}
#'   }
#'
#' @details
#' Requires \code{dat}, \code{seasons}, \code{S}, \code{N_cpue}, and
#' \code{cpue_specs} in the calling environment. Missing Leslie B0
#' values are filled by linear interpolation (\code{zoo::na.approx})
#' with remaining gaps set to the median.
#'
#' @seealso \code{\link{RunLeslie}}, \code{\link{DepleteModel}}
#'
#' @importFrom zoo na.approx
#'
#' @export
BuildPars <- function() {
  ##  Build starting log_R from Leslie B0 - helps to set parameters based ff individual estimates - use median for unknown years ##
  depl_series <- which(sapply(cpue_specs, `[[`, "is_depl"))

  ## RunLeslie function
  leslie_diag <- RunLeslie()
  leslie_B0 <- rep(NA, S)
  for (i in seq_along(seasons)) {
    for (k in depl_series) {
      row <- leslie_diag[leslie_diag$season == seasons[i] & leslie_diag$series == k, ]
      if (nrow(row) > 0 && row$flag == "ok") {
        leslie_B0[i] <- row$B0
        break
      }
    }
  }

  leslie_B0 <- zoo::na.approx(leslie_B0, na.rm = FALSE)
  leslie_B0[is.na(leslie_B0)] <- median(leslie_B0, na.rm = TRUE)
  ## Parameters ##
  par <- list(
    log_B1         = log(leslie_B0[1]),                      # Initial (model start) biomass
    log_R          = log(pmax(leslie_B0[-1], 10)),           # Recruitment based off individual depletions
    log_sigma_step = log(0.5)                                # Sigma for re applied to recruitment random walk
  )

  ## Add starting parameter values for catchability deviations from mean based on mean of good depletion for each index (log_q_cpueX_dev) and a sigma for each catchability ("log_q_cpueX") ##
  for (k in seq_len(N_cpue)) {
    if (cpue_specs[[k]]$is_depl) {
      ok_rows <- leslie_diag[leslie_diag$series == k & leslie_diag$flag == "ok", ]
      q_start <- if (nrow(ok_rows) > 0) mean(ok_rows$b) else 0.001
      par[[paste0("log_q_cpue", k, "_dev")]] <- rep(log(q_start), S)
      par[[paste0("log_sigma_q", k)]]        <- log(0.2)
    } else {
      par[[paste0("log_q_cpue", k)]] <- log(0.003)
    }
  }

  return(par)
}

#' Define random effects parameter names
#'
#' Returns the character vector of parameter names to be treated as
#' random effects by \code{RTMB::MakeADFun}. Includes recruitment
#' and, for each depletion CPUE series, the catchability deviation
#' vector.
#'
#' @param names Character vector of baseline random effect names.
#'   Default \code{"log_R"}.
#'
#' @return Character vector of random effect parameter names.
#'
#' @details
#' Requires \code{N_cpue} and \code{cpue_specs} in the calling
#' environment.
#'
#' @seealso \code{\link{BuildPars}}, \code{\link{DepleteModel}}
#'
#' @export
RandomPars <- function(names = "log_R") {
  ## Both recruitment and catchability (if depletion) are random effects ##
  random_pars <- names
  for (k in seq_len(N_cpue)) {
    if (cpue_specs[[k]]$is_depl) {
      random_pars <- c(random_pars, paste0("log_q_cpue", k, "_dev"))
    }
  }
  return(random_pars)
}

#' Save example input file to working directory
#'
#' Copies the template \code{DepletionDataIn.csv} from the package
#' installation to the user's current working directory (or a
#' specified path).
#'
#' @param dest Character string giving the destination directory.
#'   Defaults to the current working directory.
#' @param overwrite Logical. Overwrite an existing file? Default
#'   \code{FALSE}.
#'
#' @return The destination file path (invisibly).
#'
#' @export
save_example_input <- function(dest = getwd(), overwrite = FALSE) {
  src  <- system.file("extdata", "DepletionDataIn.csv", package = "LinkDeplete")
  out  <- file.path(dest, "DepletionDataIn.csv")

  if (file.exists(out) && !overwrite) {
    stop("File already exists. Use overwrite = TRUE to replace.", call. = FALSE)
  }

  file.copy(src, out, overwrite = overwrite)
  message("Saved to: ", out)
  invisible(out)
}

#' Build TMB model from DepleteModel
#'
#' Wrapper that sets the correct environment for RTMB AD taping
#' before calling \code{MakeADFun}.
#'
#' @param Pin Starting parameter list from \code{\link{BuildPars}}.
#' @param random Character vector of random effect names from
#'   \code{\link{RandomPars}}.
#' @param silent Logical. Suppress output? Default \code{FALSE}.
#'
#' @return An RTMB model object as returned by \code{MakeADFun}.
#'
#' @export
BuildModel <- function(Pin, random, silent = FALSE) {
  fn <- DepleteModel
  environment(fn) <- globalenv()
  MakeADFun(fn, Pin, random = random, silent = silent)
}

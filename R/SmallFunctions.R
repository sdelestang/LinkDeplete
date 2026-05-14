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
#' @param Model Model version to use. Only currently 1 version so defaults to DepleteModel.
#' @param Pin Starting parameter list from \code{\link{BuildPars}}.
#' @param random Character vector of random effect names from
#'   \code{\link{RandomPars}}.
#' @param silent Logical. Suppress output? Default \code{FALSE}.
#'
#' @return An RTMB model object as returned by \code{MakeADFun}.
#'
#' @export
BuildModel <- function(Model=DepleteModel, Pin, random, silent = FALSE) {
  fn <- Model
  environment(fn) <- globalenv()
  MakeADFun(fn, Pin, random = random, silent = silent)
}


#' Jitter test for convergence diagnostics
#'
#' Refits the model from multiple randomly perturbed starting values to
#' check whether the optimiser consistently finds the same minimum.
#' Starting parameters are jittered on the log scale using additive
#' normal perturbations with standard deviation \code{v}.
#'
#' @param model Model function to fit (e.g. \code{DepleteModel} or
#'   \code{DepleteModelv2}). Default \code{DepleteModel}.
#' @param n Integer. Number of jitter replicates. Default 10.
#' @param v Numeric. Standard deviation of additive jitter on log-scale
#'   parameters (equivalent to CV on natural scale). Default 0.1.
#'
#' @return A data frame with columns:
#'   \describe{
#'     \item{\code{run}}{Jitter replicate number (0 = base fit).}
#'     \item{\code{NLL}}{Negative log-likelihood at convergence.}
#'     \item{\code{conv}}{Convergence code from \code{nlminb} (0 = success).}
#'     \item{\code{grad}}{Maximum absolute gradient component.}
#'     \item{Parameter columns}{One column per fixed-effect parameter.}
#'   }
#'   Also produces a three-panel diagnostic plot and prints a summary
#'   to the console.
#'
#' @details
#' Requires \code{mod}, \code{mout}, \code{Pin}, \code{random_pars},
#' and the model function in the calling environment. The base fit
#' (run 0) is taken from \code{mout}. Each jitter run rebuilds the
#' model from scratch with perturbed starting values.
#'
#' Parameters on the log scale are jittered additively:
#' \code{p_jittered = p_original + rnorm(0, v)}, which is equivalent
#' to multiplying the natural-scale parameter by \code{exp(rnorm(0, v))}.
#'
#' Random effect starting values are also jittered with the same CV.
#'
#' @examples
#' \dontrun{
#' jit <- RunJitter(n = 10, v = 0.1)
#' jit <- RunJitter(model = DepleteModelv2, n = 20, v = 0.2)
#' }
#'
#' @import ggplot2
#' @import patchwork
#' @importFrom tidyr pivot_longer
#'
#' @export
RunJitter <- function(model = DepleteModel, n = 10, v = 0.1) {

  # ── Base fit (run 0) ─────────────────────────────────────────
  base_nll  <- mout$objective
  base_par  <- mout$par
  base_grad <- tryCatch(max(abs(mod$gr())), error = function(e) NA)

  par_names <- names(base_par)
  n_par     <- length(base_par)

  # Storage
  results <- data.frame(
    run  = 0L,
    NLL  = base_nll,
    conv = mout$convergence,
    grad = base_grad
  )
  results[par_names] <- as.list(base_par)

  cat("Run  0 (base): NLL =", round(base_nll, 4), "\n")

  # ── Jitter runs ──────────────────────────────────────────────
  for (i in seq_len(n)) {
    tryCatch({
      # Jitter fixed-effect starting values
      Pin_jit <- Pin
      for (nm in names(Pin_jit)) {
        val <- Pin_jit[[nm]]
        Pin_jit[[nm]] <- val + rnorm(length(val), 0, v)
      }

      mod_jit  <- BuildModel(model, Pin_jit, random = random_pars, silent = TRUE)
      mout_jit <- nlminb(mod_jit$par, mod_jit$fn, mod_jit$gr,
                         control = list(iter.max = 2000, eval.max = 4000))
      grad_jit <- tryCatch(max(abs(mod_jit$gr())), error = function(e) NA)

      row <- data.frame(
        run  = i,
        NLL  = mout_jit$objective,
        conv = mout_jit$convergence,
        grad = grad_jit
      )
      row[par_names] <- as.list(mout_jit$par)
      results <- rbind(results, row)

      flag <- if (mout_jit$convergence != 0) " *NO CONV*" else ""
      cat("Run ", sprintf("%2d", i), ": NLL =", round(mout_jit$objective, 4),
          " conv =", mout_jit$convergence, flag, "\n")

    }, error = function(e) {
      cat("Run ", sprintf("%2d", i), ": FAILED -", e$message, "\n")
      row <- data.frame(run = i, NLL = NA, conv = NA, grad = NA)
      row[par_names] <- NA
      results <<- rbind(results, row)
    })
  }

  # ── Console summary ─────────────────────────────────────────
  converged <- results[!is.na(results$conv) & results$conv == 0, ]
  cat("\n── Summary ──────────────────────────────────────\n")
  cat("Converged:", nrow(converged), "/", n + 1, "(including base)\n")
  cat("Base NLL: ", round(base_nll, 4), "\n")

  if (nrow(converged) > 1) {
    cat("NLL range: ", round(min(converged$NLL), 4), " – ",
        round(max(converged$NLL), 4), "\n")
    cat("NLL spread:", round(max(converged$NLL) - min(converged$NLL), 4), "\n")

    best <- converged[which.min(converged$NLL), ]
    if (best$NLL < base_nll - 0.01) {
      warning("Jitter run ", best$run, " found a better NLL (",
              round(best$NLL, 4), " vs base ", round(base_nll, 4),
              "). Consider refitting from those starting values.",
              call. = FALSE)
    } else {
      cat("Base fit is at or near the global minimum.\n")
    }
  }

  # ── Plots ───────────────────────────────────────────────────
  th <- theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(size = 10, face = "bold"))

  # 1. NLL dot plot
  plot_df <- results[!is.na(results$NLL), ]
  plot_df$is_base <- plot_df$run == 0
  plot_df$converged <- !is.na(plot_df$conv) & plot_df$conv == 0

  p_nll <- ggplot(plot_df, aes(x = factor(run), y = NLL)) +
    geom_hline(yintercept = base_nll, linetype = "dashed",
               colour = "grey50", linewidth = 0.5) +
    geom_point(aes(colour = converged, shape = is_base), size = 3) +
    scale_colour_manual(values = c("TRUE" = "#2166ac", "FALSE" = "#d6604d"),
                        labels = c("Failed", "Converged"), name = NULL) +
    scale_shape_manual(values = c("TRUE" = 17, "FALSE" = 16),
                       labels = c("Jitter", "Base"), name = NULL) +
    labs(x = "Run", y = "NLL", title = "NLL by jitter run") + th

  # 2. Parameter boxplots (natural scale)
  par_df <- converged[, c("run", par_names), drop = FALSE]

  # Transform to natural scale for interpretability
  par_natural <- par_df
  for (nm in par_names) {
    if (grepl("^log_", nm)) {
      par_natural[[nm]] <- exp(par_df[[nm]])
    }
  }
  # Clean names for display
  display_names <- gsub("^log_", "", par_names)
  names(par_natural)[match(par_names, names(par_natural))] <- display_names

  par_long <- tidyr::pivot_longer(par_natural, cols = all_of(display_names),
                                  names_to = "parameter", values_to = "value")

  # Base values on natural scale
  base_natural <- data.frame(parameter = display_names)
  base_natural$value <- sapply(par_names, function(nm) {
    val <- base_par[nm]
    if (grepl("^log_", nm)) exp(val) else val
  })

  p_par <- ggplot(par_long, aes(x = parameter, y = value)) +
    geom_boxplot(fill = "#2166ac", alpha = 0.3, outlier.shape = 1) +
    geom_jitter(width = 0.15, size = 1.5, alpha = 0.6, colour = "#2166ac") +
    geom_point(data = base_natural, colour = "red", size = 3, shape = 17) +
    facet_wrap(~parameter, scales = "free", nrow = 1) +
    labs(x = NULL, y = "Estimate (natural scale)",
         title = "Parameter estimates across jitters",
         subtitle = "Red triangle = base fit") +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + th

  # 3. NLL vs each parameter (converged runs only)
  nll_par_long <- tidyr::pivot_longer(converged[, c("NLL", par_names)],
                                      cols = all_of(par_names),
                                      names_to = "parameter",
                                      values_to = "value")

  p_nll_par <- ggplot(nll_par_long, aes(x = value, y = NLL)) +
    geom_point(colour = "#2166ac", size = 2, alpha = 0.7) +
    geom_point(data = data.frame(
      parameter = par_names,
      value     = as.numeric(base_par),
      NLL       = base_nll
    ), colour = "red", size = 3, shape = 17) +
    facet_wrap(~parameter, scales = "free_x", nrow = 1) +
    labs(x = "Parameter value (log scale)", y = "NLL",
         title = "NLL vs parameter value",
         subtitle = "Red triangle = base fit") + th

  print(p_nll / p_par / p_nll_par +
          plot_layout(heights = c(1, 1, 1)) +
          plot_annotation(
            title   = paste0("Jitter diagnostics (n = ", n, ", CV = ", v, ")"),
            theme   = theme(plot.title = element_text(face = "bold", size = 12))
          ))

  invisible(results)
}

#' Retrospective analysis
#'
#' Sequentially peels years from the end of the time series and refits
#' the model to diagnose retrospective bias. Produces a three-panel
#' plot of biomass, recruitment, and catchability trajectories for each
#' peel, and computes Mohn's rho for key quantities.
#'
#' @param model Model function to fit (e.g. \code{DepleteModel} or
#'   \code{DepleteModelv2}). Default \code{DepleteModel}.
#' @param n_peel Integer. Number of years to peel back. Default 5.
#'
#' @return A list (invisibly) with elements:
#'   \describe{
#'     \item{\code{mohns_rho}}{Named vector of Mohn's rho for biomass,
#'       recruitment, and depletion.}
#'     \item{\code{fits}}{List of per-peel data frames containing
#'       estimated trajectories.}
#'   }
#'   Also prints Mohn's rho values and produces a diagnostic plot.
#'
#' @details
#' Requires the original (unpeeled) data objects in the calling
#' environment: \code{dat_monthly}, \code{cpue_specs}, \code{Pin},
#' \code{random_pars}, and model configuration parameters passed to
#' \code{\link{MakeDataIn}} (accessed via \code{Datain}).
#'
#' Each peel rebuilds season structure, catch/CPUE matrices, and
#' starting parameters from scratch using \code{\link{SetSeasons}},
#' \code{\link{SetCatchMatrix}}, \code{\link{SetCPUEMatrix}},
#' \code{\link{BuildPars}}, and \code{\link{RandomPars}}.
#'
#' Mohn's rho is calculated as:
#' \deqn{\rho = \frac{1}{n} \sum_{p=1}^{n}
#'   \frac{\hat{X}_{peel,T_p} - \hat{X}_{full,T_p}}{\hat{X}_{full,T_p}}}
#' where \eqn{T_p} is the terminal year of peel \eqn{p}.
#'
#' @examples
#' \dontrun{
#' retro <- RunRetro(model = DepleteModelv2, n_peel = 5)
#' }
#'
#' @import ggplot2
#' @import patchwork
#'
#' @export
RunRetro <- function(model = DepleteModel, n_peel = 5) {

  # ── Save current global state ───────────────────────────────
  orig_dat_monthly  <- dat_monthly
  orig_seasons      <- seasons
  orig_S            <- S
  orig_TS           <- TS
  orig_catch_mat    <- catch_mat
  orig_cpue_mats    <- cpue_mats
  orig_cpue_se_mats <- cpue_se_mats
  orig_N_cpue       <- N_cpue
  orig_Datain       <- Datain
  orig_mod          <- mod
  orig_mout         <- mout

  # Extract settings from current Datain
  M_annual   <- Datain$M * Datain$TS
  rec_ts     <- Datain$rec_ts
  ref_ts     <- Datain$ref_ts
  ref_levels <- Datain$ref_levels
  min_pts    <- Datain$min_pts

  # Detect season start from dat_monthly
  first_month <- orig_dat_monthly$month[orig_dat_monthly$season_step == 1][1]
  first_season <- min(orig_seasons)

  on.exit({
    dat_monthly  <<- orig_dat_monthly
    seasons      <<- orig_seasons
    S            <<- orig_S
    TS           <<- orig_TS
    catch_mat    <<- orig_catch_mat
    cpue_mats    <<- orig_cpue_mats
    cpue_se_mats <<- orig_cpue_se_mats
    N_cpue       <<- orig_N_cpue
    Datain       <<- orig_Datain
    mod          <<- orig_mod
    mout         <<- orig_mout
  })

  # ── Full model sdreport ─────────────────────────────────────
  sdr_full  <- sdreport(orig_mod)
  est_full  <- summary(sdr_full, "report")

  extract_vec <- function(est, name) {
    rows <- rownames(est) == name
    exp(est[rows, "Estimate"])
  }

  B_full   <- extract_vec(est_full, "log_B_vec")
  R_full   <- extract_vec(est_full, "log_Brec_vec")
  D_full   <- extract_vec(est_full, "log_status_avg_vec")

  # Get q from full model
  bp_full    <- orig_mod$env$last.par.best
  q_full_all <- list()
  for (k in seq_len(orig_N_cpue)) {
    if (cpue_specs[[k]]$is_depl) {
      q_full_all[[k]] <- exp(bp_full[names(bp_full) == paste0("log_q_cpue", k, "_dev")])
    } else {
      q_val <- exp(orig_mout$par[paste0("log_q_cpue", k)])
      q_full_all[[k]] <- rep(q_val, orig_S)
    }
  }

  # Store full fit
  fits <- list()
  fits[[1]] <- data.frame(
    season = orig_seasons,
    peel   = 0,
    B      = B_full,
    R      = R_full,
    D      = D_full,
    q1     = q_full_all[[1]]
  )

  cat("Peel  0 (full): NLL =", round(orig_mout$objective, 4),
      " seasons =", orig_S, "\n")

  # ── Peel loop ───────────────────────────────────────────────
  for (p in seq_len(n_peel)) {
    tryCatch({
      last_season <- orig_seasons[orig_S - p]

      # Subset data
      dat_monthly <<- orig_dat_monthly[orig_dat_monthly$season <= last_season, ]
      seasons     <<- sort(unique(dat_monthly$season))
      S           <<- length(seasons)
      TS          <<- 12

      SetCatchMatrix()
      SetCPUEMatrix()

      Datain <<- MakeDataIn(M = M_annual, rec_ts = rec_ts,
                            ref_ts = ref_ts, ref_levels = ref_levels)

      Pin_p       <- BuildPars()
      random_p    <- RandomPars()
      mod_p       <- BuildModel(model, Pin_p, random = random_p, silent = TRUE)
      mout_p      <- nlminb(mod_p$par, mod_p$fn, mod_p$gr,
                            control = list(iter.max = 2000, eval.max = 4000))

      # Temporarily assign for sdreport
      mod  <<- mod_p
      mout <<- mout_p

      sdr_p <- sdreport(mod_p)
      est_p <- summary(sdr_p, "report")

      B_p <- extract_vec(est_p, "log_B_vec")
      R_p <- extract_vec(est_p, "log_Brec_vec")
      D_p <- extract_vec(est_p, "log_status_avg_vec")

      bp_p <- mod_p$env$last.par.best
      if (cpue_specs[[1]]$is_depl) {
        q1_p <- exp(bp_p[names(bp_p) == "log_q_cpue1_dev"])
      } else {
        q1_p <- rep(exp(mout_p$par["log_q_cpue1"]), S)
      }

      fits[[p + 1]] <- data.frame(
        season = seasons,
        peel   = p,
        B      = B_p,
        R      = R_p,
        D      = D_p,
        q1     = q1_p
      )

      cat("Peel ", sprintf("%2d", p), ": NLL =", round(mout_p$objective, 4),
          " conv =", mout_p$convergence,
          " seasons =", S, " (to ", last_season, ")\n")

    }, error = function(e) {
      cat("Peel ", sprintf("%2d", p), ": FAILED -", e$message, "\n")
    })
  }

  # ── Mohn's rho ──────────────────────────────────────────────
  rho_B <- numeric(0)
  rho_R <- numeric(0)
  rho_D <- numeric(0)

  for (p in seq_len(n_peel)) {
    if (length(fits) < p + 1) next
    peel_df <- fits[[p + 1]]
    term_season <- max(peel_df$season)
    idx_full <- which(orig_seasons == term_season)

    if (length(idx_full) == 1) {
      term_idx <- nrow(peel_df)
      rho_B <- c(rho_B, (peel_df$B[term_idx] - B_full[idx_full]) / B_full[idx_full])
      rho_R <- c(rho_R, (peel_df$R[term_idx] - R_full[idx_full]) / R_full[idx_full])
      rho_D <- c(rho_D, (peel_df$D[term_idx] - D_full[idx_full]) / D_full[idx_full])
    }
  }

  mohns_rho <- c(
    Biomass   = round(mean(rho_B), 4),
    Recruit   = round(mean(rho_R), 4),
    Depletion = round(mean(rho_D), 4)
  )

  cat("\n── Mohn's rho ──────────────────────────────────\n")
  cat("  Biomass:  ", mohns_rho["Biomass"], "\n")
  cat("  Recruit:  ", mohns_rho["Recruit"], "\n")
  cat("  Depletion:", mohns_rho["Depletion"], "\n")
  cat("  (|rho| < 0.15 generally acceptable for short-lived species)\n")

  # ── Plotting ────────────────────────────────────────────────
  all_df <- do.call(rbind, fits)
  all_df$peel_f <- factor(all_df$peel)

  th <- theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(size = 10, face = "bold"),
          legend.position = "bottom")

  peel_colours <- c("0" = "black",
                    setNames(scales::hue_pal()(n_peel), as.character(1:n_peel)))
  peel_sizes <- c("0" = 1.2, setNames(rep(0.7, n_peel), as.character(1:n_peel)))

  p_B <- ggplot(all_df, aes(season, B, colour = peel_f, linewidth = peel_f)) +
    geom_line() +
    scale_colour_manual(values = peel_colours, name = "Peel") +
    scale_linewidth_manual(values = peel_sizes, guide = "none") +
    labs(x = NULL, y = "Biomass (t)", title = "Jan biomass (post-recruit)") + th

  p_R <- ggplot(all_df, aes(season, R, colour = peel_f, linewidth = peel_f)) +
    geom_line() +
    scale_colour_manual(values = peel_colours, name = "Peel") +
    scale_linewidth_manual(values = peel_sizes, guide = "none") +
    labs(x = NULL, y = "Recruitment (t)", title = "Annual recruitment") + th

  p_q <- ggplot(all_df, aes(season, q1, colour = peel_f, linewidth = peel_f)) +
    geom_line() +
    scale_colour_manual(values = peel_colours, name = "Peel") +
    scale_linewidth_manual(values = peel_sizes, guide = "none") +
    labs(x = NULL, y = "Catchability (q)", title = "Time-varying q (series 1)") + th

  print(p_B / p_R / p_q +
          plot_layout(guides = "collect") +
          plot_annotation(
            title   = paste0("Retrospective analysis (", n_peel, " peels)"),
            caption = paste0("Mohn's rho:  Biomass = ", mohns_rho["Biomass"],
                             "   Recruitment = ", mohns_rho["Recruit"],
                             "   Depletion = ", mohns_rho["Depletion"]),
            theme   = theme(plot.title   = element_text(face = "bold", size = 12),
                            plot.caption = element_text(size = 9))
          ))

  invisible(list(mohns_rho = mohns_rho, fits = fits))
}



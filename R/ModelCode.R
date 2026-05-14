#' Depletion model negative log-likelihood
#'
#' RTMB objective function implementing a within-season Leslie depletion
#' model linked across seasons by a random walk on log-recruitment.
#' Biomass is tracked at monthly time steps with natural mortality applied
#' after catch removal. CPUE series are treated as either depletion
#' indices (fitted via closed-form weighted least squares within each
#' season) or biomass indices (fitted via lognormal likelihood against
#' mid-time-step biomass).
#'
#' @param Pin Named list of parameters unpacked by \code{getAll()}.
#'   Expected elements:
#'   \describe{
#'     \item{\code{log_B1}}{Log initial biomass.}
#'     \item{\code{log_R}}{Log recruitment vector (length \code{S - 1}).}
#'     \item{\code{log_sigma_step}}{Log standard deviation of the
#'       recruitment random walk.}
#'     \item{\code{log_q_cpueK}}{Log catchability for non-depletion
#'       series \emph{K} (one per series where \code{cpue_is_depl_K == 0}).}
#'     \item{\code{log_q_cpueK_dev}}{Log catchability deviations vector
#'       (length \code{S}) for depletion series \emph{K}.}
#'     \item{\code{log_sigma_qK}}{Log standard deviation of the
#'       catchability random walk for depletion series \emph{K}.}
#'   }
#'
#' @return Scalar negative log-likelihood. Side effects via
#'   \code{REPORT()} and \code{ADREPORT()} place the following on the
#'   AD tape:
#'   \describe{
#'     \item{\code{cpue_hat}}{List of predicted CPUE matrices
#'       (\code{S x TS}), one per series.}
#'     \item{\code{log_B_vec}}{Log post-recruitment biomass by season.}
#'     \item{\code{log_Brec_vec}}{Log recruitment biomass by season.}
#'     \item{\code{log_Bmid_vec}, \code{log_B0mid_vec}}{Log mid-season
#'       fished and unfished biomass at reference time step.}
#'     \item{\code{log_Bend_vec}, \code{log_B0end_vec}}{Log end-of-season
#'       fished and unfished biomass.}
#'     \item{\code{log_Bavg_vec}, \code{log_B0avg_vec}}{Log mean annual
#'       fished and unfished biomass.}
#'     \item{\code{log_status_vec}}{Log mid-season depletion (B/B0).}
#'     \item{\code{log_status_avg_vec}}{Log mean-annual depletion.}
#'     \item{\code{log_F_vec}}{Log annual instantaneous fishing mortality.}
#'     \item{\code{log_Bspm_vec}}{Log biomass comparable to surplus
#'       production model estimates.}
#'   }
#'
#' @details
#' The data list \code{dat} (accessed via \code{getAll()}) must contain:
#' \describe{
#'   \item{\code{S}}{Number of seasons.}
#'   \item{\code{TS}}{Number of monthly time steps per season.}
#'   \item{\code{N_cpue}}{Number of CPUE series.}
#'   \item{\code{M}}{Instantaneous natural mortality per time step.}
#'   \item{\code{rec_ts}}{Time step within season at which recruitment
#'     enters.}
#'   \item{\code{ref_ts}}{Reference time step for mid-season biomass
#'     reporting.}
#'   \item{\code{min_pts}}{Minimum data points required for a Leslie
#'     regression within a season.}
#'   \item{\code{catch}}{Catch matrix (\code{S x TS}).}
#'   \item{\code{cpue_is_depl_K}}{Integer flag (0/1) indicating whether
#'     series \emph{K} is a depletion index.}
#'   \item{\code{cpueK_depl}}{Integer vector of time steps used in the
#'     depletion regression for series \emph{K}.}
#'   \item{\code{cK}, \code{cKse}}{CPUE and standard error matrices
#'       (\code{S x TS}) for series \emph{K}.}
#' }
#'
#' @seealso \code{\link{RunLeslie}} for standalone Leslie regression
#'   diagnostics, \code{\link{Plot_fit}} and \code{\link{Plot_summary}}
#'   for visualising model output.
#'
#' @export
DepleteModel <- function(Pin) {
  getAll(Pin, Datain)

  sigma_step <- exp(log_sigma_step)
  nll <- 0
  nll <- nll - dlnorm(sigma_step, log(0.3), 0.5, log = TRUE)

  # ── Unpack parameters ─────────────────────────────────────────
  q_fixed <- list()
  sigma_q <- list()

  ## Put q into vectors - annual or fixed
  for (k in seq_len(N_cpue)) {
    is_depl_k <- Datain[[paste0("cpue_is_depl_", k)]]
    if (is_depl_k == 0) {
      q_fixed[[k]] <- exp(get(paste0("log_q_cpue", k)))
    } else {
      sigma_q[[k]] <- exp(get(paste0("log_sigma_q", k)))
    }
  }

  #  Random walk prior on log_Recruitment #
  for (s in 2:(S - 1)) {
    nll <- nll - dnorm(log_R[s] - log_R[s - 1], 0, sigma_step, log = TRUE)
  }

  #  Random walk priors on q deviations #
  for (k in seq_len(N_cpue)) {
    is_depl_k <- Datain[[paste0("cpue_is_depl_", k)]]
    if (is_depl_k == 1) {
      qdev <- get(paste0("log_q_cpue", k, "_dev"))
      for (s in 2:S) {
        nll <- nll - dnorm(qdev[s] - qdev[s - 1], 0, sigma_q[[k]], log = TRUE)
      }
    }
  }

  # objects for storage #
  cpue_hat <- list()
  for (k in seq_len(N_cpue)) cpue_hat[[k]] <- matrix(0, nrow=S, ncol=TS)

  log_B_vec      <- numeric(S)    # Biomass vector
  log_Brec_vec   <- numeric(S)    # Biomass after recruitment
  log_Bmid_vec   <- numeric(S)    # Mid years biomass - fished
  log_B0mid_vec  <- numeric(S)    # B0 is unfished version
  log_Bend_vec   <- numeric(S)    # Biomass at end of season
  log_B0end_vec  <- numeric(S)
  log_Bavg_vec   <- numeric(S)    # Mean annual fished biomass
  log_B0avg_vec  <- numeric(S)    # Mean annual unfished biomass
  log_status_vec <- numeric(S)     # B/B0 version
  log_status_avg_vec <- numeric(S) # B/B0 based on mean annual biomass
  log_Bspm_vec <- numeric(S)       # replicate SPM estimates for comparison
  log_F_vec      <- numeric(S)     # make annual F by compunding

  B_now       <- exp(log_B1)      # Start biomass
  B0_now      <- exp(log_B1)      # Second biomass which will not be fished
  B_end_prev  <- exp(log_B1)      # End of prevuious timesteps (TS) biomass
  B0_end_prev <- exp(log_B1)

  # Start Season loop #
  for (s in seq_len(S)) {

    q_year <- list()   # Store qs
    for (k in seq_len(N_cpue)) {
      is_depl_k <- Datain[[paste0("cpue_is_depl_", k)]]
      if (is_depl_k == 1) {
        q_year[[k]] <- exp(get(paste0("log_q_cpue", k, "_dev"))[s])
      } else {
        q_year[[k]] <- q_fixed[[k]]
      }
    }

    if (s > 1) {
      B_now  <- B_end_prev
      B0_now <- B0_end_prev
    }

    log_Brec_vec[s] <- if (s == 1) log_B1 else log_R[s - 1] # No recruitment in years 1 TS 1

    # Conduct Leslie depletion analysis using closed form linear model - faster
    # Weighted least squares (WLS) solved in closed form so everything stays on the RTMB AD tape.
    Sw   <- numeric(N_cpue)
    SwK  <- numeric(N_cpue)
    SwK2 <- numeric(N_cpue)
    SwU  <- numeric(N_cpue)
    SwUK <- numeric(N_cpue)
    n_k  <- integer(N_cpue)

    K_cum      <- 0
    B_mean_acc <- 0     # Mean biomass
    B0_mean_acc <- 0    # Unfished mean biomass accumulator
    C_total    <- 0

    if (s == 1) log_B_vec[1] <- log_B1

    ## Timestep loop ##
    for (m in seq_len(TS)) {

      if (m == rec_ts && s > 1) {              # Add recruitment
        B_now  <- B_now  + exp(log_R[s - 1])
        B0_now <- B0_now + exp(log_R[s - 1])
        log_B_vec[s] <- log(B_now)
      }

      C_m    <- catch[s, m]                     # Catch
      B_mid_raw <- B_now - C_m / 2
      B_mid  <- (B_mid_raw + sqrt(B_mid_raw^2 + 1e-8)) / 2   # AD-safe soft floor
      B0_mid <- B0_now

      for (k in seq_len(N_cpue)) {
        cpue_hat[[k]][s, m] <- q_year[[k]] * B_mid   # Model Est cpue in middle TS.
      }

      # Adding v small offset whenver thnings are logged to stop blow up
      if (m == ref_ts) {
        log_Bmid_vec[s]  <- log(B_mid  + 1e-6)
        log_B0mid_vec[s] <- log(B0_mid + 1e-6)
      }

      # Accumulate for annual mean biomass
      B_mean_acc  <- B_mean_acc  + B_mid
      B0_mean_acc <- B0_mean_acc + B0_mid
      C_total     <- C_total + C_m

      # Depletion accumulators if doing depletion
      for (k in seq_len(N_cpue)) {
        is_depl_k <- Datain[[paste0("cpue_is_depl_", k)]]
        if (is_depl_k == 0) next       # Move on

        depl_steps <- Datain[[paste0("cpue", k, "_depl")]]
        if (!(m %in% depl_steps)) next       # Move on

        ck   <- Datain[[paste0("c", k)]]        # Grab cpues + cpue se
        ckse <- Datain[[paste0("c", k, "se")]]
        # If cpue and se > Zero is not NA - add to accumulators to do analytic lm
        if (!is.na(ck[s, m]) && ck[s, m] > 0 &&
            !is.na(ckse[s, m]) && ckse[s, m] > 0) {
          K_mid <- K_cum + C_m / 2
          w <- 1 / ckse[s, m]^2
          Sw[k]   <- Sw[k]   + w
          SwK[k]  <- SwK[k]  + w * K_mid
          SwK2[k] <- SwK2[k] + w * K_mid^2
          SwU[k]  <- SwU[k]  + w * ck[s, m]
          SwUK[k] <- SwUK[k] + w * ck[s, m] * K_mid
          n_k[k]  <- n_k[k]  + 1
        }
      }

      # Non-depletion series: direct biomass index likelihood
      for (k in seq_len(N_cpue)) {
        is_depl_k <- Datain[[paste0("cpue_is_depl_", k)]]
        if (is_depl_k == 1) next

        ck   <- Datain[[paste0("c", k)]]
        ckse <- Datain[[paste0("c", k, "se")]]

        if (!is.na(ck[s, m]) && ck[s, m] > 0 &&
            !is.na(ckse[s, m]) && ckse[s, m] > 0) {
          cv  <- ckse[s, m] / ck[s, m]
          ## NLL for fit between observed and estimated catch rate
          nll <- nll - dnorm(log(ck[s, m]), log(q_fixed[[k]] * B_mid), cv, log = TRUE)
        }
      }

      K_cum  <- K_cum + C_m
      B_rem  <- B_now - C_m
      B_now  <- (B_rem + sqrt(B_rem^2 + 1e-8)) / 2 * exp(-M)
      B0_now <-  B0_now        * exp(-M)

    } # end timestep loop

    B_end_prev  <- B_now
    B0_end_prev <- B0_now

    ## Make ehaps of different vectors of indices for plotting asnd reporting
    ## Make biomass vector to compare with SPM
    log_Bspm_vec[s] <- log(exp(log_B_vec[s]) - sum(catch[s, ]) / 2 + 1e-6)

    # Mean annual biomass (fished and unfished)
    log_Bavg_vec[s]  <- log(B_mean_acc  / TS + 1e-6)
    log_B0avg_vec[s] <- log(B0_mean_acc / TS + 1e-6)

    # Depletion based on mean annual biomass
    log_status_avg_vec[s] <- log_Bavg_vec[s] - log_B0avg_vec[s]

    # Annual F from mean biomass
    B_mean  <- B_mean_acc / TS
    exploit <- min(C_total / (B_mean + 1e-6), 0.99)

    # Annual F from mid-year biomass
    C_total    <- sum(catch[s, ])
    B_post_rec <- exp(log_B_vec[s])
    B_mid_yr   <- B_post_rec * exp(-M * (TS / 2 - rec_ts))
    exploit    <- min(C_total / (B_mid_yr + 1e-6), 0.99)
    log_F_vec[s] <- log(-log(1 - exploit) + 1e-6)
    log_Bend_vec[s]   <- log(B_now + 1e-6)
    log_B0end_vec[s]  <- log(B0_now + 1e-6)
    log_status_vec[s] <- log_Bmid_vec[s] - log_B0mid_vec[s]

    # Use accumulators to compute Leslie likelihood for each depletion series
    for (k in seq_len(N_cpue)) {
      is_depl_k <- Datain[[paste0("cpue_is_depl_", k)]]
      if (is_depl_k == 0 || n_k[k] < min_pts) next

      denom_k <- Sw[k] * SwK2[k] - SwK[k]^2
      b_k     <- -(Sw[k] * SwUK[k] - SwU[k] * SwK[k]) / denom_k
      a_k     <- (SwU[k] - (-b_k) * SwK[k]) / Sw[k]
      B0_k    <- a_k / b_k

      if (b_k > 0 && B0_k > 0) {
        sig_k <- sigma_q[[k]]
        qdev  <- get(paste0("log_q_cpue", k, "_dev"))
        nll   <- nll - dnorm(log(b_k), qdev[s], sig_k, log = TRUE)

        var_b  <- Sw[k]   / denom_k
        var_a  <- SwK2[k] / denom_k
        cov_ab <- -SwK[k] / denom_k
        var_B0 <- (1/b_k)^2 * var_a + (a_k/b_k^2)^2 * var_b -
          2*(a_k/b_k^3)*cov_ab
        if (var_B0 > 0) {
          cv_B0 <- min(sqrt(var_B0) / B0_k, 5.0)
          ## NLL for fit between Leslise  and Model initial biomass after recruitment
          nll   <- nll - dnorm(log(B0_k), log_B_vec[s], cv_B0, log = TRUE)
        }
      }
    }

  } # end season loop

  # ── REPORT / ADREPORT ──────────────────────────────────────────────────
  REPORT(cpue_hat)
  ADREPORT(log_Bspm_vec)
  ADREPORT(log_B_vec)
  ADREPORT(log_Brec_vec)
  ADREPORT(log_Bmid_vec)
  ADREPORT(log_B0mid_vec)
  ADREPORT(log_Bend_vec)
  ADREPORT(log_B0end_vec)
  ADREPORT(log_Bavg_vec)
  ADREPORT(log_B0avg_vec)
  ADREPORT(log_status_vec)
  ADREPORT(log_status_avg_vec)
  ADREPORT(log_F_vec)

  nll
}

#' Depletion model v2 — direct CPUE likelihood
#'
#' Variant of \code{\link{DepleteModel}} that replaces the closed-form
#' Leslie regression likelihood with direct lognormal fits of observed
#' CPUE against \code{q * B_mid} for all series (depletion and
#' non-depletion alike). This avoids collapsing within-season
#' observations into summary statistics and lets each monthly data
#' point independently inform biomass and catchability.
#'
#' @inheritParams DepleteModel
#'
#' @return Scalar negative log-likelihood with the same REPORT/ADREPORT
#'   side effects as \code{\link{DepleteModel}}.
#'
#' @details
#' The data and parameter lists are identical to \code{\link{DepleteModel}}.
#' The only structural change is in how depletion-series CPUE enters the
#' likelihood: individual observations within the depletion window are
#' fitted directly rather than via Leslie slope and intercept summaries.
#'
#' @seealso \code{\link{DepleteModel}}, \code{\link{Plot_fit}}
#'
#' @export
DepleteModelv2 <- function(Pin) {
  getAll(Pin, Datain)

  sigma_step <- exp(log_sigma_step)
  nll <- 0

  # ── Prior on sigma_step ────────────────────────────────────────
  nll <- nll - dlnorm(sigma_step, log(0.3), 0.5, log = TRUE)

  # ── Unpack parameters ─────────────────────────────────────────
  q_fixed <- list()
  sigma_q <- list()

  for (k in seq_len(N_cpue)) {
    is_depl_k <- Datain[[paste0("cpue_is_depl_", k)]]
    if (is_depl_k == 0) {
      q_fixed[[k]] <- exp(get(paste0("log_q_cpue", k)))
    } else {
      sigma_q[[k]] <- exp(get(paste0("log_sigma_q", k)))
    }
  }

  # ── Random walk prior on log_Recruitment ───────────────────────
  for (s in 2:(S - 1)) {
    nll <- nll - dnorm(log_R[s] - log_R[s - 1], 0, sigma_step, log = TRUE)
  }

  # ── Random walk priors on q deviations ─────────────────────────
  for (k in seq_len(N_cpue)) {
    is_depl_k <- Datain[[paste0("cpue_is_depl_", k)]]
    if (is_depl_k == 1) {
      qdev <- get(paste0("log_q_cpue", k, "_dev"))
      for (s in 2:S) {
        nll <- nll - dnorm(qdev[s] - qdev[s - 1], 0, sigma_q[[k]], log = TRUE)
      }
    }
  }

  # ── Storage objects ────────────────────────────────────────────
  cpue_hat <- list()
  for (k in seq_len(N_cpue)) cpue_hat[[k]] <- matrix(0, nrow = S, ncol = TS)

  log_B_vec          <- numeric(S)
  log_Brec_vec       <- numeric(S)
  log_Bmid_vec       <- numeric(S)
  log_B0mid_vec      <- numeric(S)
  log_Bend_vec       <- numeric(S)
  log_B0end_vec      <- numeric(S)
  log_Bavg_vec       <- numeric(S)
  log_B0avg_vec      <- numeric(S)
  log_status_vec     <- numeric(S)
  log_status_avg_vec <- numeric(S)
  log_Bspm_vec       <- numeric(S)
  log_F_vec          <- numeric(S)

  B_now       <- exp(log_B1)
  B0_now      <- exp(log_B1)
  B_end_prev  <- exp(log_B1)
  B0_end_prev <- exp(log_B1)

  # ── Season loop ────────────────────────────────────────────────
  for (s in seq_len(S)) {

    q_year <- list()
    for (k in seq_len(N_cpue)) {
      is_depl_k <- Datain[[paste0("cpue_is_depl_", k)]]
      if (is_depl_k == 1) {
        q_year[[k]] <- exp(get(paste0("log_q_cpue", k, "_dev"))[s])
      } else {
        q_year[[k]] <- q_fixed[[k]]
      }
    }

    if (s > 1) {
      B_now  <- B_end_prev
      B0_now <- B0_end_prev
    }

    log_Brec_vec[s] <- if (s == 1) log_B1 else log_R[s - 1]

    K_cum       <- 0
    B_mean_acc  <- 0
    B0_mean_acc <- 0
    C_total     <- 0

    if (s == 1) log_B_vec[1] <- log_B1

    # ── Timestep loop ──────────────────────────────────────────
    for (m in seq_len(TS)) {

      if (m == rec_ts && s > 1) {
        B_now  <- B_now  + exp(log_R[s - 1])
        B0_now <- B0_now + exp(log_R[s - 1])
        log_B_vec[s] <- log(B_now)
      }

      C_m       <- catch[s, m]
      B_mid_raw <- B_now - C_m / 2
      B_mid     <- (B_mid_raw + sqrt(B_mid_raw^2 + 1e-8)) / 2   # AD-safe soft floor
      B0_mid    <- B0_now

      for (k in seq_len(N_cpue)) {
        cpue_hat[[k]][s, m] <- q_year[[k]] * B_mid
      }

      if (m == ref_ts) {
        log_Bmid_vec[s]  <- log(B_mid  + 1e-6)
        log_B0mid_vec[s] <- log(B0_mid + 1e-6)
      }

      B_mean_acc  <- B_mean_acc  + B_mid
      B0_mean_acc <- B0_mean_acc + B0_mid
      C_total     <- C_total + C_m

      # ── Direct CPUE likelihood for ALL series ──────────────
      for (k in seq_len(N_cpue)) {
        is_depl_k <- Datain[[paste0("cpue_is_depl_", k)]]

        # For depletion series, only fit within the depletion window
        if (is_depl_k == 1) {
          depl_steps <- Datain[[paste0("cpue", k, "_depl")]]
          if (!(m %in% depl_steps)) next
        }

        ck   <- Datain[[paste0("c", k)]]
        ckse <- Datain[[paste0("c", k, "se")]]

        if (!is.na(ck[s, m]) && ck[s, m] > 0 &&
            !is.na(ckse[s, m]) && ckse[s, m] > 0) {
          cv  <- ckse[s, m] / ck[s, m]
          nll <- nll - dnorm(log(ck[s, m]), log(q_year[[k]] * B_mid + 1e-8),
                             cv, log = TRUE)
        }
      }

      K_cum  <- K_cum + C_m
      B_rem  <- B_now - C_m
      B_now  <- (B_rem + sqrt(B_rem^2 + 1e-8)) / 2 * exp(-M)
      B0_now <-  B0_now * exp(-M)

    } # end timestep loop

    B_end_prev  <- B_now
    B0_end_prev <- B0_now

    # ── Reporting vectors ──────────────────────────────────────
    log_Bspm_vec[s] <- log(exp(log_B_vec[s]) - sum(catch[s, ]) / 2 + 1e-6)

    log_Bavg_vec[s]  <- log(B_mean_acc  / TS + 1e-6)
    log_B0avg_vec[s] <- log(B0_mean_acc / TS + 1e-6)

    log_status_avg_vec[s] <- log_Bavg_vec[s] - log_B0avg_vec[s]

    # Annual F from mid-year biomass
    B_mean     <- B_mean_acc / TS
    C_total    <- sum(catch[s, ])
    B_post_rec <- exp(log_B_vec[s])
    B_mid_yr   <- B_post_rec * exp(-M * (TS / 2 - rec_ts))
    exploit    <- min(C_total / (B_mid_yr + 1e-6), 0.99)
    log_F_vec[s] <- log(-log(1 - exploit) + 1e-6)

    log_Bend_vec[s]   <- log(B_now + 1e-6)
    log_B0end_vec[s]  <- log(B0_now + 1e-6)
    log_status_vec[s] <- log_Bmid_vec[s] - log_B0mid_vec[s]

  } # end season loop

  # ── REPORT / ADREPORT ──────────────────────────────────────────
  REPORT(cpue_hat)
  ADREPORT(log_Bspm_vec)
  ADREPORT(log_B_vec)
  ADREPORT(log_Brec_vec)
  ADREPORT(log_Bmid_vec)
  ADREPORT(log_B0mid_vec)
  ADREPORT(log_Bend_vec)
  ADREPORT(log_B0end_vec)
  ADREPORT(log_Bavg_vec)
  ADREPORT(log_B0avg_vec)
  ADREPORT(log_status_vec)
  ADREPORT(log_status_avg_vec)
  ADREPORT(log_F_vec)

  nll
}

#' Leslie regression diagnostic
#'
#' Runs standalone weighted least squares Leslie depletion regressions
#' for each season and each depletion CPUE series, independent of the
#' TMB model. Useful for diagnosing data support, identifying seasons
#' with negative slopes or insufficient data, and cross-checking
#' model-estimated B0 and catchability.
#'
#' @param x Unused placeholder argument (function relies on objects in
#'   the calling environment).
#'
#' @return A data frame with one row per season--series combination and
#'   columns:
#'   \describe{
#'     \item{\code{season}}{Season identifier.}
#'     \item{\code{series}}{Integer CPUE series index.}
#'     \item{\code{n}}{Number of data points in the regression.}
#'     \item{\code{b}}{Estimated catchability (Leslie slope, sign-corrected).}
#'     \item{\code{B0}}{Estimated initial biomass from the Leslie intercept
#'       and slope (\code{a / b}).}
#'     \item{\code{cv_B0}}{Coefficient of variation of B0 derived via the
#'       delta method.}
#'     \item{\code{flag}}{Diagnostic flag: \code{"ok"}, \code{"NEG_SLOPE"},
#'       \code{"NEG_B0"}, \code{"NEG_VAR"}, \code{"TOO_FEW_PTS"}, or
#'       \code{"NOT_DEPL"}.}
#'   }
#'
#' @details
#' Requires the following objects in the calling environment:
#' \describe{
#'   \item{\code{seasons}}{Vector of season identifiers.}
#'   \item{\code{N_cpue}}{Number of CPUE series.}
#'   \item{\code{TS}}{Number of time steps per season.}
#'   \item{\code{cpue_specs}}{List of per-series specification lists with
#'     elements \code{is_depl} and \code{depl}.}
#'   \item{\code{cpue_mats}, \code{cpue_se_mats}}{Lists of CPUE and
#'     standard error matrices (\code{S x TS}).}
#'   \item{\code{catch_mat}}{Catch matrix (\code{S x TS}).}
#'   \item{\code{dat}}{Data list containing \code{min_pts}.}
#' }
#'
#' @importFrom dplyr case_when
#'
#' @export
RunLeslie <- function(x) {
  leslie_diag <- do.call(rbind, lapply(seq_along(seasons), function(i) {
  do.call(rbind, lapply(seq_len(N_cpue), function(k) {

    if (!cpue_specs[[k]]$is_depl) {
      return(data.frame(season = seasons[i], series = k, n = 0,
                        b = NA, B0 = NA, cv_B0 = NA, flag = "NOT_DEPL"))
    }

    cpue_depl <- cpue_specs[[k]]$depl
    cpue_mat  <- cpue_mats[[k]]
    cpue_se   <- cpue_se_mats[[k]]

    K_cum <- 0
    Sw <- 0; SwK <- 0; SwK2 <- 0; SwU <- 0; SwUK <- 0; n <- 0

    for (m in seq_len(TS)) {
      C_m <- catch_mat[i, m]
      if (m %in% cpue_depl) {
        u  <- cpue_mat[i, m]
        se <- cpue_se[i, m]
        if (!is.na(u) && u > 0 && !is.na(se) && se > 0) {
          w    <- 1 / se^2
          Sw   <- Sw   + w
          SwK  <- SwK  + w * K_cum
          SwK2 <- SwK2 + w * K_cum^2
          SwU  <- SwU  + w * u
          SwUK <- SwUK + w * u * K_cum
          n    <- n + 1L
        }
      }
      K_cum <- K_cum + C_m
    }

    if (n >= Datain$min_pts) {
      denom <- Sw * SwK2 - SwK^2
      b_raw <- (Sw * SwUK - SwU * SwK) / denom
      b     <- -b_raw
      a     <- (SwU - b_raw * SwK) / Sw
      B0    <- a / b
      var_b  <- Sw   / denom
      var_a  <- SwK2 / denom
      cov_ab <- -SwK / denom
      var_B0 <- (1/b)^2 * var_a + (a/b^2)^2 * var_b - 2*(a/b^3)*cov_ab

      flag <- dplyr::case_when(
        b   <= 0    ~ "NEG_SLOPE",
        B0  <= 0    ~ "NEG_B0",
        var_B0 <= 0 ~ "NEG_VAR",
        TRUE        ~ "ok"
      )

      data.frame(season = seasons[i], series = k, n = n,
                 b = round(b, 6), B0 = round(B0, 1),
                 cv_B0 = if (var_B0 > 0 && B0 > 0) round(sqrt(var_B0)/B0, 3) else NA,
                 flag = flag)
    } else {
      data.frame(season = seasons[i], series = k, n = n,
                 b = NA, B0 = NA, cv_B0 = NA, flag = "TOO_FEW_PTS")
    }
  }))
}))}

#' Build a season-by-time-step matrix from monthly data
#'
#' Reshapes a column from the \code{dat_monthly} data frame into an
#' \code{S x TS} matrix indexed by season (rows) and within-season
#' time step (columns).
#'
#' @param col Character string naming the column in \code{dat_monthly}
#'   to extract.
#'
#' @return A numeric matrix of dimension \code{S x TS}. Cells without
#'   corresponding observations are \code{NA}.
#'
#' @details
#' Requires \code{dat_monthly} (with columns \code{season} and
#' \code{season_step}), \code{seasons}, \code{S}, and \code{TS} in the
#' calling environment.
#'
#' @export
make_mat <- function(col) {
  mat <- matrix(NA, S, TS)
  for (i in seq_along(seasons)) {
    sub <- dat_monthly[dat_monthly$season == seasons[i], ]
    mat[i, sub$season_step] <- sub[[col]]
  }
  mat
}


#' Project one season forward under a specified catch
#'
#' Deterministic forward projection of biomass through a single season
#' at monthly resolution, applying catch according to a user-specified
#' allocation across fishing time steps. Returns the end-of-season
#' depletion level (B/B0).
#'
#' @param total_catch Scalar total catch to distribute across the season.
#' @param B_end_prev End-of-season fished biomass from the previous season.
#' @param B0_end_prev End-of-season unfished biomass from the previous season.
#' @param R Recruitment biomass to be added at time step \code{rec_ts}.
#' @param M Instantaneous natural mortality per time step.
#' @param catch_split Numeric vector of proportional catch allocation
#'   across fishing time steps (must sum to 1 and have the same length
#'   as \code{fish_steps}).
#' @param fish_steps Integer vector of time steps in which fishing occurs.
#' @param rec_ts Integer time step at which recruitment enters.
#' @param TS Total number of time steps per season. Default 12.
#'
#' @return Scalar end-of-season depletion ratio (B/B0).
#'
#' @seealso \code{\link{solve_target_catch}}
#'
#' @export
project_depletion_end <- function(total_catch, B_end_prev, B0_end_prev, R, M,
                                  catch_split, fish_steps, rec_ts, TS = 12) {

  B_now  <- B_end_prev
  B0_now <- B0_end_prev

  for (m in seq_len(TS)) {
    if (m == rec_ts) {
      B_now  <- B_now  + R
      B0_now <- B0_now + R
    }

    C_m <- if (m %in% fish_steps) {
      total_catch * catch_split[which(fish_steps == m)]
    } else {
      0
    }

    B_now  <- (B_now  - C_m) * exp(-M)
    B0_now <-  B0_now         * exp(-M)
  }

  B_now / B0_now
}


#' Solve for the catch achieving a target depletion level
#'
#' Wrapper around \code{\link{project_depletion_end}} that uses
#' \code{\link[stats]{uniroot}} to find the total seasonal catch
#' producing a specified end-of-season depletion ratio (B/B0).
#'
#' @param target_depletion Target depletion ratio (e.g. 0.40 for 40\%
#'   of unfished biomass).
#' @param B_end_prev End-of-season fished biomass from the previous season.
#' @param B0_end_prev End-of-season unfished biomass from the previous season.
#' @param R Recruitment biomass to be added at time step \code{rec_ts}.
#' @param M Instantaneous natural mortality per time step.
#' @param catch_split Numeric vector of proportional catch allocation
#'   across fishing time steps.
#' @param fish_steps Integer vector of time steps in which fishing occurs.
#' @param rec_ts Integer time step at which recruitment enters.
#' @param TS Total number of time steps per season. Default 12.
#' @param tol Tolerance passed to \code{\link[stats]{uniroot}}.
#'   Default 0.01.
#'
#' @return A list as returned by \code{\link[stats]{uniroot}}, where
#'   \code{$root} is the total catch that achieves the target depletion.
#'
#' @details
#' The search interval runs from zero catch to \code{B_end_prev + R},
#' which represents the theoretical maximum removable biomass. If the
#' target depletion is unachievable within this range (e.g. the stock
#' is already below the target with zero catch), \code{uniroot} will
#' error.
#'
#' @seealso \code{\link{project_depletion_end}}, \code{\link{DepleteModel}}
#'
#'
#' @export
solve_target_catch <- function(target_depletion, B_end_prev, B0_end_prev, R, M,
                               catch_split, fish_steps, rec_ts, TS = 12,
                               tol = 0.01) {

  stats::uniroot(
    function(C) {
      project_depletion_end(C, B_end_prev, B0_end_prev, R, M,
                            catch_split, fish_steps, rec_ts, TS) - target_depletion
    },
    interval = c(0, B_end_prev + R),
    tol      = tol
  )
}

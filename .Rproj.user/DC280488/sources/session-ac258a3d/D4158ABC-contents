
privatise_y <- function(gg, is_private) {

  if (!is_private) return(gg)
  gg +
    scale_y_continuous(labels = NULL) +
    labs(y = paste0(gg$labels$y, " (relative)"))
}

rescale_private <- function(x) {
  # Index to series mean so shape is preserved

  mu <- mean(x, na.rm = TRUE)
  if (is.na(mu) || mu == 0) return(x)
  x / mu
}

#' Plot within-season CPUE fits by year
#'
#' Plots observed vs predicted CPUE across monthly time steps within each
#' fishing season, faceted by year.
#'
#' @param IsPrivate Logical. If \code{TRUE}, CPUE values are rescaled to the
#'   series mean and y-axis tick labels are suppressed. Default \code{FALSE}.
#'
#' @return A \code{ggplot} object.
#'
#' @details
#' Requires \code{mod}, \code{dat_monthly}, \code{Datain}, \code{seasons},
#' \code{N_cpue}, and \code{cpue_specs} in the calling environment.
#'
#' @seealso \code{\link{Plot_summary}}, \code{\link{Plot_stock_status}}
#'
#' @import ggplot2
#' @import dplyr
#' @import patchwork
#' @importFrom scales percent_format
#' @importFrom stats lm quantile median
#'
#' @export
Plot_fit <- function(IsPrivate = FALSE) {
  rep_det    <- mod$report()
  season_idx <- match(dat_monthly$season, seasons)

  # ── Build prediction vectors for each series ─────────────────────────────
  cpue_pred <- list()
  for (k in seq_len(N_cpue)) {
    hat_mat <- rep_det$cpue_hat[[k]]
    cpue_pred[[k]] <- mapply(function(i, s) hat_mat[i, s],
                             season_idx, dat_monthly$season_step)
  }

  # ── Build obs_long generically ───────────────────────────────────────────
  series_labels  <- paste0("CPUE ", seq_len(N_cpue))
  series_colours <- c("#7b2d8b", "steelblue", "goldenrod", "#d6604d")[seq_len(N_cpue)]
  names(series_colours) <- series_labels

  obs_long <- bind_rows(lapply(seq_len(N_cpue), function(k) {
    spec       <- cpue_specs[[k]]
    depl_steps <- Datain[[paste0("cpue", k, "_depl")]]

    pred_df      <- dat_monthly
    pred_df$pred <- cpue_pred[[k]]

    pred_df |>
      transmute(
        season, month, season_step,
        series = series_labels[k],
        obs    = .data[[spec$col]],
        se     = .data[[spec$se_col]],
        pred   = pred,
        in_nll = season_step %in% depl_steps
      )
  })) |>
    filter(!is.na(obs) | !is.na(pred))

  # ── Privacy rescaling ────────────────────────────────────────────────────
  if (IsPrivate) {
    obs_long <- obs_long |>
      group_by(series) |>
      mutate(
        scale_mu = mean(obs, na.rm = TRUE),
        scale_mu = ifelse(is.na(scale_mu) | scale_mu == 0, 1, scale_mu),
        obs      = obs  / scale_mu,
        pred     = pred / scale_mu,
        se       = se   / scale_mu
      ) |>
      select(-scale_mu) |>
      ungroup()
  }

  # ── Only show pred lines where obs exist in that season ──────────────────
  series_by_season <- obs_long |>
    filter(!is.na(obs)) |>
    distinct(season, series)

  pred_long <- obs_long |>
    inner_join(series_by_season, by = c("season", "series"))

  # ── Shading: use first depletion series for window ───────────────────────
  first_depl <- which(sapply(cpue_specs, `[[`, "is_depl"))[1]
  shade_steps <- Datain[[paste0("cpue", first_depl, "_depl")]]
  shade_df <- tibble(
    xmin = min(shade_steps) - 0.5,
    xmax = max(shade_steps) + 0.5,
    ymin = -Inf, ymax = Inf
  )

  step_labels <- c("N","D","J","F","M","A","M","J","J","A","S","O")

  p <- ggplot() +
    geom_rect(data = shade_df,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              inherit.aes = FALSE, fill = "grey80", alpha = 0.3) +
    geom_line(data = pred_long |> filter(!is.na(pred)),
              aes(x = season_step, y = pred, colour = series),
              linewidth = 0.7) +
    geom_point(data = obs_long |> filter(!is.na(obs) & in_nll),
               aes(x = season_step, y = obs, colour = series),
               size = 1.8) +
    geom_point(data = obs_long |> filter(!is.na(obs) & !in_nll),
               aes(x = season_step, y = obs, colour = series),
               size = 1.8, shape = 1) +
    geom_errorbar(data = obs_long |> filter(!is.na(obs) & !is.na(se)),
                  aes(x = season_step, y = obs,
                      ymin = obs - se, ymax = obs + se,
                      colour = series),
                  width = 0.3, linewidth = 0.4) +
    scale_colour_manual(values = series_colours) +
    scale_x_continuous(breaks = 1:12, labels = step_labels) +
    facet_wrap(~season, ncol = 6, scales = "free_y") +
    labs(x        = "Month (season start = Nov)",
         y        = if (IsPrivate) "CPUE (relative)" else "CPUE",
         colour   = NULL,
         title    = "Within-season CPUE fit by year",
         subtitle = "Filled = in NLL  |  Open = not in NLL  |  Shaded = fitting window") +
    theme_bw(base_size = 9) +
    theme(legend.position  = "bottom",
          strip.background = element_rect(fill = "grey90"),
          strip.text       = element_text(size = 7),
          axis.text.x      = element_text(size = 6))

  if (IsPrivate) {
    p <- p + scale_y_continuous(labels = NULL)
  }

  p
}

#' Plot depletion model summary
#'
#' Produces a six-panel summary of model outputs: January biomass,
#' annual recruitment, annual catch, time-varying catchability by CPUE
#' series, and two scatter plots examining the recruitment–biomass
#' relationship (recruitment vs mean annual biomass, and end-season
#' biomass vs following-year recruitment).
#'
#' @param IsPrivate Logical. If \code{TRUE}, annual catch is rescaled to its
#'   mean and y-axis tick labels on the catch panel are suppressed. All other
#'   panels (biomass, recruitment, catchability, scatter plots) are left
#'   unchanged as they represent model outputs rather than raw data.
#'   Default \code{FALSE}.
#'
#' @return A \code{patchwork} composite \code{ggplot} object.
#'
#' @details
#' Requires the following objects in the calling environment:
#' \describe{
#'   \item{\code{mod}}{Fitted TMB/RTMB model object with \code{par} and
#'     \code{env$last.par.best}.}
#'   \item{\code{mout}}{Optimiser output list with element \code{par}.}
#'   \item{\code{dat_monthly}}{Data frame of monthly observations.}
#'   \item{\code{dat}}{Named list containing depletion step vectors.}
#'   \item{\code{catch_mat}}{Matrix of catch by season (rows) and month (columns).}
#'   \item{\code{seasons}}{Vector of season labels.}
#'   \item{\code{S}}{Integer number of seasons.}
#'   \item{\code{N_cpue}}{Integer number of CPUE series.}
#'   \item{\code{cpue_specs}}{List of per-series specification lists.}
#' }
#'
#' The caption reports estimated process-error and catchability parameters
#' alongside R-squared values for the recruitment–biomass relationships.
#'
#' @seealso \code{\link{Plot_fit}}, \code{\link{Plot_stock_status}}
#' @export
Plot_summary <- function(IsPrivate=FALSE) {

  sdr     <- sdreport(mod)
  sdr_est <- summary(sdr, "report")

  get_adreport <- function(name) {
    rows <- rownames(sdr_est) == name
    data.frame(
      season = seasons,
      est    = exp(sdr_est[rows, "Estimate"]),
      lo     = exp(sdr_est[rows, "Estimate"] - 1.96 * sdr_est[rows, "Std. Error"]),
      hi     = exp(sdr_est[rows, "Estimate"] + 1.96 * sdr_est[rows, "Std. Error"])
    )
  }

  B_df    <- get_adreport("log_B_vec")
  Rec_df  <- get_adreport("log_Brec_vec")
  Bavg_df <- get_adreport("log_Bavg_vec")
  Bend_df <- get_adreport("log_Bend_vec")

  catch_annual <- data.frame(season = seasons, catch = rowSums(catch_mat))

  if (IsPrivate) {
    catch_mu <- mean(catch_annual$catch, na.rm = TRUE)
    if (is.na(catch_mu) || catch_mu == 0) catch_mu <- 1
    catch_annual$catch <- catch_annual$catch / catch_mu
  }


  # ── Build q data frame generically ───────────────────────────────────────
  bp <- mod$env$last.par.best
  series_labels  <- paste0("CPUE ", seq_len(N_cpue))
  series_colours <- c("#7b2d8b", "#1a9988", "#c8860a", "#d6604d")[seq_len(N_cpue)]
  names(series_colours) <- series_labels

  q_list <- lapply(seq_len(N_cpue), function(k) {
    spec <- cpue_specs[[k]]
    depl_steps <- Datain[[paste0("cpue", k, "_depl")]]

    # Get q values
    if (spec$is_depl) {
      q_vals <- exp(bp[names(bp) == paste0("log_q_cpue", k, "_dev")])
    } else {
      q_vals <- rep(exp(mout$par[paste0("log_q_cpue", k)]), S)
    }

    # Identify seasons with data
    seasons_with <- dat_monthly |>
      filter(season_step %in% depl_steps | !spec$is_depl,
             !is.na(.data[[spec$col]])) |>
      pull(season) |> unique()

    data.frame(
      season = seasons,
      q      = ifelse(seasons %in% seasons_with, q_vals, NA),
      series = series_labels[k]
    )
  })

  q_df <- bind_rows(q_list) |>
    mutate(series = factor(series, levels = series_labels))

  # ── Rec vs mean annual biomass ───────────────────────────────────────────
  RB_df <- data.frame(season = seasons, R = Rec_df$est, Bavg = Bavg_df$est)
  BR_df <- data.frame(
    season = seasons[-S],
    Bend   = Bend_df$est[-S],
    R_next = Rec_df$est[-1]
  )

  th <- theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          strip.text       = element_text(face = "bold"),
          plot.title       = element_text(size = 9))

  col_B <- "#2166ac"; col_R <- "#4dac26"; col_C <- "#d6604d"

  p_bio <- ggplot(B_df, aes(season, est)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), fill = col_B, alpha = 0.2) +
    geom_line(colour = col_B, linewidth = 0.9) +
    labs(x = NULL, y = "Biomass (t)", title = "Jan biomass (post-recruit)") + th

  p_rec <- ggplot(Rec_df, aes(season, est)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), fill = col_R, alpha = 0.2) +
    geom_line(colour = col_R, linewidth = 0.9) +
    labs(x = NULL, y = "Recruitment (t)", title = "Annual recruitment") + th

  p_cat <- ggplot(catch_annual, aes(season, catch)) +
    geom_col(fill = col_C, alpha = 0.8, width = 0.7) +
    labs(x = NULL,
         y = if (IsPrivate) "Catch (relative)" else "Catch (t)",
         title = "Annual catch") + th

  if (IsPrivate) {
    p_cat <- p_cat + scale_y_continuous(labels = NULL)
  }

  p_q <- ggplot(q_df |> filter(!is.na(q)), aes(season, q, colour = series)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 1.5) +
    scale_colour_manual(values = series_colours) +
    facet_wrap(~series, ncol = 1, scales = "free_y") +
    labs(x = NULL, y = "Catchability (q)",
         title = "Time-varying q\n(where data exist)") +
    theme(legend.position = "none") + th

  p_RB <- ggplot(RB_df, aes(R, Bavg)) +
    geom_point(aes(colour = season), size = 2.5) +
    geom_smooth(method = "lm", se = TRUE,
                colour = "grey30", linewidth = 0.8, linetype = "dashed") +
    geom_text(aes(label = season, colour = season), size = 2.3, vjust = -0.7) +
    scale_colour_viridis_c(option = "plasma", name = "Season") +
    labs(x = "Recruitment (t, Jan)", y = "Mean annual biomass (t)",
         title = "Recruitment vs mean annual biomass") + th

  p_BR <- ggplot(BR_df, aes(Bend, R_next)) +
    geom_point(aes(colour = season), size = 2.5) +
    geom_smooth(method = "lm", se = TRUE,
                colour = "grey30", linewidth = 0.8, linetype = "dashed") +
    geom_text(aes(label = season, colour = season), size = 2.3, vjust = -0.7) +
    scale_colour_viridis_c(option = "plasma", name = "Season") +
    labs(x = "End-season biomass (t)", y = "Following year recruitment (t)",
         title = "End-season B vs next recruitment") + th

  r2_RB <- summary(lm(Bavg ~ R,      data = RB_df))$r.squared
  r2_BR <- summary(lm(R_next ~ Bend, data = BR_df))$r.squared

  # ── Caption: build sigma/q strings generically ──────────────────────────
  caption_parts <- paste0("sigma_step = ", round(exp(mod$par["log_sigma_step"]), 3))
  for (k in seq_len(N_cpue)) {
    if (cpue_specs[[k]]$is_depl) {
      caption_parts <- c(caption_parts,
                         paste0("sigma_q", k, " = ",
                                round(exp(mod$par[paste0("log_sigma_q", k)]), 3)))
    } else {
      caption_parts <- c(caption_parts,
                         paste0("q_cpue", k, " = ",
                                round(exp(mout$par[paste0("log_q_cpue", k)]), 5)))
    }
  }
  caption_parts <- c(caption_parts,
                     paste0("R2 Rec-Bio = ", round(r2_RB, 3)),
                     paste0("R2 Bio-Rec = ", round(r2_BR, 3)))

  ((p_bio + p_rec + p_cat) / (p_q + p_RB + p_BR)) +
    plot_layout(heights = c(2, 3)) +
    plot_annotation(
      title   = "Depletion model Summary",
      caption = paste(caption_parts, collapse = "   "),
      theme   = theme(plot.title = element_text(face = "bold"))
    )
}

#' Plot stock status and Kobe diagram
#'
#' Produces a four-panel stock status summary: depletion trajectory
#' (mean annual B/B0) with target/threshold/limit reference bands,
#' fishing mortality with an empirical F reference, mean annual fished
#' vs unfished biomass, and a Kobe phase plot of B/B_target against
#' F/F_ref.
#'
#' @param IsPrivate Logical. Accepted for interface consistency with
#'   \code{\link{Plot_fit}} and \code{\link{Plot_summary}} but currently
#'   unused, as all panels display dimensionless ratios or model-derived
#'   quantities that do not expose confidential catch or CPUE data.
#'   Default \code{FALSE}.
#'
#' @return A \code{patchwork} composite \code{ggplot} object.
#'
#' @details
#' Requires the following objects in the calling environment:
#' \describe{
#'   \item{\code{mod}}{Fitted TMB/RTMB model object.}
#'   \item{\code{Datain}}{Named list containing \code{ref_levels}, a length-3
#'     numeric vector of target, threshold, and limit depletion reference
#'     points (as proportions of B0).}
#'   \item{\code{seasons}}{Vector of season labels.}
#' }
#'
#' The F reference level is computed as the mean F across seasons where
#' depletion was between the threshold and target (plus 5\%), falling back
#' to the median F if fewer than three qualifying seasons exist.
#'
#' @seealso \code{\link{Plot_fit}}, \code{\link{Plot_summary}}
#' @export
Plot_stock_status <- function(IsPrivate=FALSE) {

  sdr     <- sdreport(mod)
  sdr_est <- summary(sdr, "report")

  get_adreport <- function(name) {
    rows <- rownames(sdr_est) == name
    data.frame(
      season = seasons,
      est    = exp(sdr_est[rows, "Estimate"]),
      lo     = exp(sdr_est[rows, "Estimate"] - 1.96 * sdr_est[rows, "Std. Error"]),
      hi     = exp(sdr_est[rows, "Estimate"] + 1.96 * sdr_est[rows, "Std. Error"])
    )
  }

  get_adreport_log <- function(name) {
    rows <- rownames(sdr_est) == name
    data.frame(
      season = seasons,
      est    = sdr_est[rows, "Estimate"],
      se     = sdr_est[rows, "Std. Error"]
    )
  }

  # Unpack reference levels
  tgt   <- Datain$ref_levels[1]   # Target
  thr   <- Datain$ref_levels[2]   # Threshold
  lim   <- Datain$ref_levels[3]   # Limit

  Bavg_df    <- get_adreport("log_Bavg_vec")
  B0avg_df   <- get_adreport("log_B0avg_vec")
  F_df       <- get_adreport("log_F_vec")
  status_log <- get_adreport_log("log_status_avg_vec")

  status_df <- data.frame(
    season = seasons,
    D      = exp(status_log$est),
    D_lo   = exp(status_log$est - 1.96 * status_log$se),
    D_hi   = exp(status_log$est + 1.96 * status_log$se)
  )

  F_ref_years <- status_df$season[status_df$D >= thr & status_df$D <= (tgt + 0.05)]
  F_ref <- if (length(F_ref_years) >= 3) {
    mean(F_df$est[seasons %in% F_ref_years])
  } else {
    median(F_df$est)
  }

  kobe_df <- data.frame(
    season  = seasons,
    B_ratio = status_df$D / tgt,
    F_ratio = F_df$est / F_ref
  )

  th <- theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          plot.title       = element_text(size = 9))

  p_dep <- ggplot(status_df, aes(season, D)) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = tgt, ymax = Inf,
             fill = "green3", alpha = 0.10) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = thr, ymax = tgt,
             fill = "yellow2", alpha = 0.20) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = lim, ymax = thr,
             fill = "orange", alpha = 0.20) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = lim,
             fill = "red", alpha = 0.15) +
    geom_hline(yintercept = c(lim, thr, tgt),
               linetype = "dashed", colour = c("red","orange","darkgreen"),
               linewidth = 0.6) +
    geom_ribbon(aes(ymin = D_lo, ymax = D_hi), fill = "#2166ac", alpha = 0.25) +
    geom_line(colour = "#2166ac", linewidth = 0.9) +
    annotate("text", x = max(seasons), y = tgt + 0.01,
             label = paste0("Target (", tgt*100, "%)"),
             hjust = 1, size = 3, colour = "darkgreen") +
    annotate("text", x = max(seasons), y = thr + 0.01,
             label = paste0("Threshold (", thr*100, "%)"),
             hjust = 1, size = 3, colour = "orange") +
    annotate("text", x = max(seasons), y = lim + 0.01,
             label = paste0("Limit (", lim*100, "%)"),
             hjust = 1, size = 3, colour = "red") +
    scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
    labs(x = NULL, y = expression(bar(B) / bar(B)[0]),
         title = "Depletion status (mean annual biomass)") + th

  p_F <- ggplot(F_df, aes(season, est)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), fill = "#d6604d", alpha = 0.25) +
    geom_line(colour = "#d6604d", linewidth = 0.9) +
    geom_hline(yintercept = F_ref, linetype = "dashed",
               colour = "grey30", linewidth = 0.6) +
    annotate("text", x = min(seasons) + 1, y = F_ref * 1.05,
             label = paste0("F_ref (F", tgt*100, "%)"), hjust = 0, size = 3, colour = "grey30") +
    labs(x = NULL, y = "Annual instantaneous F",
         title = "Fishing mortality") + th

  F_ylim <- quantile(F_df$est, 0.95, na.rm = TRUE) * 1.5

  p_F <- p_F +
    coord_cartesian(ylim = c(0, F_ylim))

  bio_df <- data.frame(
    season = seasons,
    B0avg  = B0avg_df$est,
    Bavg   = Bavg_df$est
  )

  p_bio <- ggplot(bio_df, aes(season)) +
    geom_line(aes(y = B0avg), colour = "grey50", linewidth = 0.8,
              linetype = "dashed") +
    geom_line(aes(y = Bavg), colour = "#2166ac", linewidth = 0.9) +
    geom_ribbon(data = Bavg_df,
                aes(x = season, ymin = lo, ymax = hi), fill = "#2166ac", alpha = 0.20) +
    annotate("text", x = max(seasons), y = max(bio_df$B0avg),
             label = "Mean annual B0 (unfished)", hjust = 1, size = 3, colour = "grey40") +
    annotate("text", x = max(seasons), y = max(bio_df$Bavg) * 0.85,
             label = "Mean annual B (fished)", hjust = 1, size = 3, colour = "#2166ac") +
    labs(x = NULL, y = "Biomass (t)",
         title = "Mean annual fished vs unfished biomass") + th

  p_kobe <- ggplot() +
    annotate("rect", xmin = 1, xmax = Inf, ymin = -Inf, ymax = 1,
             fill = "green3", alpha = 0.12) +
    annotate("rect", xmin = -Inf, xmax = 1, ymin = -Inf, ymax = 1,
             fill = "yellow2", alpha = 0.15) +
    annotate("rect", xmin = 1, xmax = Inf, ymin = 1, ymax = Inf,
             fill = "orange", alpha = 0.15) +
    annotate("rect", xmin = -Inf, xmax = 1, ymin = 1, ymax = Inf,
             fill = "red", alpha = 0.12) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey30") +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey30") +
    geom_path(data = kobe_df, aes(B_ratio, F_ratio),
              colour = "grey40", linewidth = 0.5, alpha = 0.6) +
    geom_point(data = kobe_df, aes(B_ratio, F_ratio, colour = season),
               size = 2.5) +
    geom_text(data = kobe_df |>
                filter(season %in% c(min(seasons), max(seasons),
                                     seq(1985, 2020, by = 5))),
              aes(B_ratio, F_ratio, label = season, colour = season),
              size = 2.5, vjust = -0.7) +
    scale_colour_viridis_c(option = "plasma", name = "Season") +
    labs(x = bquote(bar(B) / B[target] ~ "(" * bar(B) / .(tgt*100) * "%" ~ bar(B)[0] * ")"),
         y = bquote(F / F[ref] ~ "(" * F / F[.(tgt*100) * "%"] * ")"),
         title = "Kobe plot") + th

  (p_dep + p_F ) / (p_bio + p_kobe) +
    plot_layout(heights = c(1, 1.5)) +
    plot_annotation(
      title   = "Assessment stock status",
      caption = paste0(
        "Depletion = mean(B) / mean(B0)   |   ",
        "Target = ", tgt*100, "%   Threshold = ", thr*100, "%   Limit = ", lim*100, "%   |   F = -log(1 - C/Bmean)"
      ),
      theme = theme(plot.title = element_text(face = "bold"))
    )
}

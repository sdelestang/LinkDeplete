# LinkDeplete

Within-season Leslie depletion model linked across seasons via a random walk on recruitment, fitted using RTMB. Designed for crustacean fisheries where catch and CPUE data are available at monthly resolution.

## Installation

```r
# install.packages("devtools")
devtools::install_github("your-username/LinkDeplete")
```

## Overview

LinkDeplete fits a biomass dynamics model that tracks depletion within each fishing season using the Leslie method, while linking seasons through a random walk on log-recruitment. The framework supports:

- Multiple CPUE series (depletion indices and/or biomass indices)
- Time-varying catchability via random walk
- Dynamic depletion-based reference points (target, threshold, limit)
- Forward projection to solve for a target catch
- Privacy mode for presenting results without exposing absolute catch or CPUE values

## Quick start

```r
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(LinkDeplete)

# ── Save and load example data ───────────────────────────────────────────
save_example_input()
dat <- read.csv("DepletionDataIn.csv")

# ── Set season structure ─────────────────────────────────────────────────
# SeasonStart = 11 means the fishing season runs November to October
dat_monthly <- SetSeasons(dat, SeasonStart = 11, FirstSeason = 1978)

# ── Define CPUE series ───────────────────────────────────────────────────
# depl:    time steps used in the within-season depletion regression
# is_depl: TRUE = depletion index (Leslie fit), FALSE = biomass index only
cpue_specs <- list(
  list(col = "cpue1", se_col = "cpue1_se", depl = 3:7,  is_depl = TRUE),
  list(col = "cpue2", se_col = "cpue2_se", depl = 3:7,  is_depl = TRUE),
  list(col = "cpue3", se_col = "cpue3_se", depl = 1:12, is_depl = FALSE),
  list(col = "cpue4", se_col = "cpue4_se", depl = 3:7,  is_depl = FALSE)
)

# ── Build catch and CPUE matrices ────────────────────────────────────────
SetCatchMatrix()
SetCPUEMatrix()

# ── Assemble model inputs ────────────────────────────────────────────────
Datain      <- MakeDataIn(M = 0.5, rec_ts = 3, ref_ts = 6,
                          ref_levels = c(0.4, 0.3, 0.2))
Pin         <- BuildPars()
random_pars <- RandomPars()

# ── Fit the model ────────────────────────────────────────────────────────
mod  <- MakeADFun(DepleteModel, Pin, random = random_pars, silent = FALSE)
mout <- nlminb(mod$par, mod$fn, mod$gr,
               control = list(iter.max = 2000, eval.max = 4000))

# ── Check convergence ────────────────────────────────────────────────────
mout                    # nlminb result (convergence == 0 is good)
max(abs(mod$gr()))      # final gradient — want < 0.001

# ── Plot results ─────────────────────────────────────────────────────────
Plot_fit(IsPrivate = FALSE)           # within-season CPUE fits
Plot_summary(IsPrivate = TRUE)        # model summary (catch axis hidden)
Plot_stock_status()                   # depletion, F, Kobe
```

## Privacy mode

All three plotting functions accept an `IsPrivate` argument. When `TRUE`:

- **Plot_fit**: CPUE values are rescaled to the series mean and y-axis tick labels are removed.
- **Plot_summary**: Annual catch is rescaled to its mean and y-axis tick labels on the catch panel are removed. Biomass, recruitment, and catchability panels are unchanged (model outputs).
- **Plot_stock_status**: No change needed — all panels already display dimensionless ratios or model-derived quantities.

This allows figures to be shown in presentations or reviews without exposing confidential catch or CPUE data.

## Input data format

The input CSV (`DepletionDataIn.csv`) should contain monthly records with at minimum:

| Column | Description |
|--------|-------------|
| `year` | Calendar year |
| `month` | Calendar month (1–12) |
| `catch` | Monthly catch in kg (converted to tonnes internally) |
| `cpue1` | CPUE series 1 (NA where no data) |
| `cpue1_se` | Standard error for CPUE series 1 |
| `cpue2` | CPUE series 2 (optional) |
| `cpue2_se` | Standard error for CPUE series 2 |
| ... | Additional CPUE series as needed |

Use `save_example_input()` to obtain a template.

## Key functions

| Function | Purpose |
|----------|---------|
| `SetSeasons()` | Assign monthly data to fishing seasons |
| `SetCatchMatrix()` | Build S × TS catch matrix |
| `SetCPUEMatrix()` | Build CPUE and SE matrices for all series |
| `MakeDataIn()` | Assemble the TMB data list |
| `BuildPars()` | Generate starting parameters from Leslie regressions |
| `RandomPars()` | Define random effects for RTMB |
| `DepleteModel()` | RTMB objective function |
| `RunLeslie()` | Standalone Leslie regression diagnostics |
| `Plot_fit()` | Within-season CPUE fit plots |
| `Plot_summary()` | Six-panel model summary |
| `Plot_stock_status()` | Stock status and Kobe plot |
| `project_depletion_end()` | Forward-project one season |
| `solve_target_catch()` | Find catch for a target depletion |
| `save_example_input()` | Copy example CSV to working directory |

## License

MIT



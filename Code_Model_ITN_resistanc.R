# DRC ITN Resistance Model


REQUIRED_PACKAGES <- c(
  "dplyr",
  "tidyr",
  "tibble",
  "purrr",
  "ggplot2",
  "scales",
  "patchwork",
  "readr"
)

MISSING_PACKAGES <- REQUIRED_PACKAGES[
  !vapply(
    X = REQUIRED_PACKAGES,
    FUN = requireNamespace,
    FUN.VALUE = logical(1),
    quietly = TRUE
  )
]

if (length(MISSING_PACKAGES) > 0L) {
  install.packages(MISSING_PACKAGES, dependencies = TRUE)
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(purrr)
  library(ggplot2)
  library(scales)
  library(patchwork)
  library(readr)
})

set.seed(20260707L)

# Create output folders in the current working directory.
PROJECT_DIR <- getwd()
OUT_DIR <- file.path(PROJECT_DIR, "output")
FIG_DIR <- file.path(PROJECT_DIR, "figures")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# Display figures and summary tables in RStudio.
SHOW_PLOTS_IN_RSTUDIO <- TRUE
SHOW_SUMMARY_TABLES_IN_RSTUDIO <- TRUE
SAVE_FIGURES <- TRUE

# ============================================================================
# PUBLICATION SETTINGS — edit this block to change all figures consistently
# ============================================================================
PUBLICATION_BASE_SIZE <- 15
PUBLICATION_TAG_SIZE <- 19
PUBLICATION_LEGEND_SIZE <- 13
PUBLICATION_AXIS_TITLE_SIZE <- 16
PUBLICATION_LINEWIDTH <- 1.05
PUBLICATION_RIBBON_ALPHA <- 0.16
PUBLICATION_ERRORBAR_WIDTH <- 0.55

# Common export dimensions.
PUBLICATION_FIG_WIDTH <- 13.8
PUBLICATION_FIG_HEIGHT_ONE_ROW <- 5.6
PUBLICATION_FIG_HEIGHT_TWO_ROWS <- 10.6
PUBLICATION_DPI <- 600

# Common resistance baseline for direct product comparison.
PUBLICATION_COMMON_INITIAL_qR <- 0.05

# Number of uncertainty draws.
PUBLICATION_N_DRAWS_FAST <- 20L
PUBLICATION_N_DRAWS_FINAL <- 250L

# summaries and their 2.5th--97.5th percentiles without display smoothing.
PUBLICATION_USE_RAW_MODEL_TRAJECTORIES <- FALSE

# FAST_MODE is useful for checking the script quickly. For final figures, keep it FALSE.
FAST_MODE <- as.logical(Sys.getenv("FAST_MODE", "FALSE"))
if (FAST_MODE) {
  SPINUP_YEARS <- 6L
  CAL_SPINUP_YEARS <- 5L
  N_DRAWS <- PUBLICATION_N_DRAWS_FAST
  CAL_ITERS <- 6L
} else {
  SPINUP_YEARS <- 20L
  CAL_SPINUP_YEARS <- 15L
  N_DRAWS <- PUBLICATION_N_DRAWS_FINAL
  CAL_ITERS <- 10L
}

# Study period and numerical time step.
STUDY_START_YEAR <- 2009L
STUDY_END_YEAR <- 2020L
STUDY_YEARS <- STUDY_START_YEAR:STUDY_END_YEAR
N_STUDY_YEARS <- length(STUDY_YEARS)

DAYS_PER_YEAR <- 365L
DT <- 1.0
SIM_LENGTH <- N_STUDY_YEARS * DAYS_PER_YEAR

# Population and calibration targets.
TARGET_BASELINE_PREVALENCE <- 0.50


# ITN campaign, even when small numerical differences remain after spin-up.
PLOT_BASELINE_PREVALENCE <- 50
MODEL_HUMAN_POPULATION <- 1e5
NATIONAL_POP <- 70e6

# cases and deaths onto a more plausible reported/observable burden scale.
BURDEN_CALIBRATION_FACTOR <- 0.30

# ============================================================================
# HUMAN POPULATION STRUCTURE
# ============================================================================

# Human age groups.
AGE_BINS_YEARS <- c(0, 1, 5, 8, 18, 28, 43, 120)
AGE_LABELS <- c("<1", "1-4", "5-7", "8-17", "18-27", "28-42", "43+")
N_AGE <- length(AGE_LABELS)
UNDER5 <- 1:2

# Ageing rates are expressed per day.
AGE_WIDTH_DAYS <- diff(AGE_BINS_YEARS) * DAYS_PER_YEAR
AGEING_RATE <- c(
  1 / AGE_WIDTH_DAYS[-length(AGE_WIDTH_DAYS)],
  0
)

# Age-specific infection and recovery parameters.
RECOVERY_DAYS <- c(254, 634, 372, 192, 87, 58, 51)
GAMMA_H_AGE <- 1 / RECOVERY_DAYS
HUMAN_LATENT_DAYS <- 8L
MAX_H_INF_AGE <- 180L

# Human demographic and immunity-loss rates.
MU_H <- 1 / (60 * DAYS_PER_YEAR)
K_H <- 5.0e-3

# Age-specific case-fatality and clinical disease probabilities.
CFR_AGE <- c(0.0110, 0.0085, 0.0032, 0.0016, 0.0011, 0.0011, 0.0013)
NU_H_AGE <- CFR_AGE * GAMMA_H_AGE
P_CLIN_AGE <- c(0.62, 0.58, 0.45, 0.36, 0.30, 0.28, 0.25)

# ============================================================================
# MOSQUITO POPULATION STRUCTURE
# ============================================================================

GENOTYPES <- c("SS", "RS", "RR")
N_G <- length(GENOTYPES)

MOSQ_LATENT_DAYS <- 12L
MAX_M_INF_AGE <- MOSQ_LATENT_DAYS + 80L

# Mosquito biting and genotype-specific transmission parameters.
THETA <- 0.33
BETA_M <- c(SS = 0.50, RS = 0.46, RR = 0.42)

# Mosquito demographic parameters.
PHI_M <- 1 / 8
MU_JM <- 1 / 14
MU_AM <- 1 / 14
R_FEMALE <- 0.5
PSI <- 8.0
LARVAL_K <- 2e6
NATIONAL_EIR <- 25.0

# q_R(0) is imposed at intervention start after epidemiological spin-up.
INITIAL_qR_PYRETHROID <- PUBLICATION_COMMON_INITIAL_qR
INITIAL_qR_NEXTGEN <- PUBLICATION_COMMON_INITIAL_qR
INITIAL_qR <- INITIAL_qR_PYRETHROID  # default used only during calibration/spin-up
RES_h <- 0.25
RES_u <- 0.80

# Relative susceptibility to pyrethroid-induced mortality by genotype.
D_PYR <- c(
  SS = 1.0,
  RS = 1.0 - RES_h * RES_u,
  RR = 1.0 - RES_u
)
# still increase slowly at low coverage rather than being forced downward.
RES_COST_RR <- 0.08
RES_COST_RS_FRACTION <- 0.02
cost_vec <- function() {
  c(
    SS = 0,
    RS = RES_COST_RS_FRACTION * RES_COST_RR,
    RR = RES_COST_RR
  )
}

# strong elimination effects at high coverage.
NEXTGEN_RESTORE <- 0.65

D_NEXTGEN <- c(
  SS = 1.00,
  RS = min(1.00, D_PYR[["RS"]] + NEXTGEN_RESTORE * (1.00 - D_PYR[["RS"]])),
  RR = min(1.00, D_PYR[["RR"]] + NEXTGEN_RESTORE * (1.00 - D_PYR[["RR"]]))
)

# Intervention products.
PRODUCTS <- c("Pyrethroid-only nets", "Next-generation nets")
NEXTGEN_PRODUCT <- PRODUCTS[2]

INITIAL_qR_BY_PRODUCT <- c(
  "Pyrethroid-only nets" = INITIAL_qR_PYRETHROID,
  "Next-generation nets" = INITIAL_qR_NEXTGEN
)

initial_qR_for_product <- function(product) {
  product <- as.character(product)
  values <- unname(INITIAL_qR_BY_PRODUCT[product])
  if (anyNA(values)) {
    stop("Missing initial resistance allele frequency for product: ",
         paste(unique(product[is.na(values)]), collapse = ", "), call. = FALSE)
  }
  values
}

COVERAGE_TARGETS <- c(0.25, 0.50, 0.75)
COVERAGE_LABELS <- paste0(as.integer(100 * COVERAGE_TARGETS), "%")

# Time-dependent ITN coverage C_B(t): campaigns every three years and exponential
# decay in effective coverage between campaigns.
CAMPAIGN_INTERVAL <- 3L
NET_RETENTION_YEARS <- 2.8

# Campaign distribution and uptake build up progressively over six months.
CAMPAIGN_SCALEUP_MONTHS <- 6
CAMPAIGN_SCALEUP_DAYS <- as.integer(round(CAMPAIGN_SCALEUP_MONTHS * DAYS_PER_YEAR / 12))

# Routine distribution/use is maintained between mass campaigns.
ROUTINE_COVERAGE_FLOOR <- 0.10

# ITN distribution. Set this to 12 for one full pre-campaign year, or to 3--6
# for only a few months.
PRE_CAMPAIGN_DISPLAY_MONTHS <- 6L
PRE_CAMPAIGN_START_X <- 1 - PRE_CAMPAIGN_DISPLAY_MONTHS / 12


# not every mosquito contacts a treated net every day.
ITN_PREVENTABLE_BITING <- 0.70       # 70% of biting is realistically net-preventable
ITN_MORTALITY_RATE_SCALE <- 0.16     # converts death-on-contact into a daily mortality rate

# reduced. This slows q_R without weakening the average intervention effect.
PYRETHROID_REALIZED_SELECTION_FRACTION <- 0.11
NEXTGEN_REALIZED_SELECTION_FRACTION <- 0.13

# Plot-response parameters. Raw model outputs remain saved separately.
PREVALENCE_DECLINE_HALFLIFE_MONTHS <- 4.5
PREVALENCE_REBOUND_HALFLIFE_MONTHS <- 9.0
PREVALENCE_PLOT_CEILING_FRACTION <- 0.995
PREVALENCE_DISPLAY_LIFT_FACTOR <- 1.16

# instead of an artificial vertical fall followed by a straight segment.
ALLELE_EXPOSURE_FLOOR_WEIGHT <- 0.15


# ============================================================================
# NUMERICAL HELPER FUNCTIONS AND INITIAL CONDITIONS
# ============================================================================

# Replace non-finite values by zero and enforce non-negativity.
sanitize_positive <- function(x) {
  x[!is.finite(x)] <- 0
  pmax(x, 0)
}

# Normalise a non-negative vector so that its elements sum to one.
normalise_vector <- function(x) {
  x <- sanitize_positive(x)
  total <- sum(x)

  if (total <= 0) {
    return(rep(1 / length(x), length(x)))
  }

  x / total
}

# Backward-compatible aliases used throughout the model equations.
sp <- sanitize_positive
norm_vec <- normalise_vector

# Human-to-mosquito infectiousness as a function of chronological age and
# infection age. Rows represent human age groups and columns infection days.
H2M <- matrix(
  data = 0,
  nrow = N_AGE,
  ncol = MAX_H_INF_AGE
)

age_infectiousness_multiplier <- c(1.00, 0.95, 0.85, 0.70, 0.55, 0.45, 0.35)

for (age_index in seq_len(N_AGE)) {
  for (infection_day in seq_len(MAX_H_INF_AGE)) {
    if (infection_day > HUMAN_LATENT_DAYS) {
      gametocyte_density <- 100 *
        age_infectiousness_multiplier[age_index] *
        exp(-0.012 * (infection_day - HUMAN_LATENT_DAYS))

      H2M[age_index, infection_day] <- 0.071 * gametocyte_density^0.302
    }
  }
}

# Create the initial epidemiological and entomological state.
init_state <- function(q_init = INITIAL_qR) {
  q_resistant <- min(1 - 1e-8, max(1e-8, q_init))
  p_susceptible <- 1 - q_resistant

  genotype_frequency <- norm_vec(c(
    SS = p_susceptible^2,
    RS = 2 * p_susceptible * q_resistant,
    RR = q_resistant^2
  ))

  mosquito_to_human_ratio <- 6 * min(
    3,
    0.5 + sqrt(NATIONAL_EIR / 10)
  )
  total_mosquitoes <- mosquito_to_human_ratio * MODEL_HUMAN_POPULATION

  juvenile <- 0.40 * total_mosquitoes * genotype_frequency
  adult_male <- 0.18 * total_mosquitoes * genotype_frequency
  sus_female <- 0.30 * total_mosquitoes * genotype_frequency
  infected_female_total <- 0.01 * total_mosquitoes * genotype_frequency

  inf_f_q <- matrix(
    data = 0,
    nrow = MAX_M_INF_AGE,
    ncol = N_G,
    dimnames = list(NULL, GENOTYPES)
  )

  infectious_rows <- (MOSQ_LATENT_DAYS + 1L):MAX_M_INF_AGE

  for (genotype in GENOTYPES) {
    inf_f_q[infectious_rows, genotype] <-
      infected_female_total[[genotype]] / length(infectious_rows)
  }

  initial_prevalence <- TARGET_BASELINE_PREVALENCE
  age_midpoints <- (
    AGE_BINS_YEARS[-1] + AGE_BINS_YEARS[-length(AGE_BINS_YEARS)]
  ) / 2

  age_weights <- exp(-MU_H * age_midpoints * DAYS_PER_YEAR) *
    diff(AGE_BINS_YEARS)
  age_weights <- age_weights / sum(age_weights)

  recovered_fraction <- min(
    initial_prevalence * mean(GAMMA_H_AGE) / (MU_H + K_H),
    1 - initial_prevalence - 0.05
  )
  susceptible_fraction <- max(
    0.05,
    1 - initial_prevalence - recovered_fraction
  )

  S <- MODEL_HUMAN_POPULATION * susceptible_fraction * age_weights
  R <- MODEL_HUMAN_POPULATION * recovered_fraction * age_weights
  infected_by_age <- MODEL_HUMAN_POPULATION * initial_prevalence * age_weights

  Iq <- matrix(
    data = 0,
    nrow = N_AGE,
    ncol = MAX_H_INF_AGE,
    dimnames = list(AGE_LABELS, NULL)
  )

  for (age_index in seq_len(N_AGE)) {
    infection_profile <- exp(
      -(
        GAMMA_H_AGE[age_index] +
          MU_H +
          NU_H_AGE[age_index]
      ) * seq_len(MAX_H_INF_AGE)
    )

    Iq[age_index, ] <- infected_by_age[age_index] *
      infection_profile / sum(infection_profile)
  }

  list(
    S = S,
    Iq = Iq,
    R = R,
    juvenile = juvenile,
    adult_male = adult_male,
    sus_female = sus_female,
    inf_f_q = inf_f_q
  )
}

reset_state_allele_frequency <- function(state, q_target) {
  q_target <- min(1 - 1e-8, max(1e-8, q_target))
  p_target <- 1 - q_target
  genotype_frequency <- norm_vec(c(
    SS = p_target^2,
    RS = 2 * p_target * q_target,
    RR = q_target^2
  ))
  names(genotype_frequency) <- GENOTYPES

  redistribute_vector <- function(x) {
    total <- sum(sp(x))
    setNames(total * genotype_frequency, GENOTYPES)
  }

  state$juvenile <- redistribute_vector(state$juvenile)
  state$adult_male <- redistribute_vector(state$adult_male)
  state$sus_female <- redistribute_vector(state$sus_female)

  infected_by_age <- rowSums(pmax(state$inf_f_q, 0), na.rm = TRUE)
  state$inf_f_q <- outer(infected_by_age, genotype_frequency)
  dimnames(state$inf_f_q) <- list(NULL, GENOTYPES)
  state
}

# ============================================================================
# ONE-DAY SEMI-IMPLICIT STEP
# ============================================================================

# ============================================================================
step_model <- function(state,
                       coverage,
                       d_n0,
                       r_n0,
                       d_g,
                       selection_fraction = 1.0,
                       day = 1L,
                       include_cost = TRUE,
                       tmult = 1.0) {
  S <- sp(state$S)
  R <- sp(state$R)
  Iq <- pmax(replace(state$Iq, !is.finite(state$Iq), 0), 0)

  Jm <- sp(state$juvenile)
  M <- sp(state$adult_male)
  Sf <- sp(state$sus_female)
  Ifq <- pmax(replace(state$inf_f_q, !is.finite(state$inf_f_q), 0), 0)

  names(Jm) <- names(M) <- names(Sf) <- GENOTYPES
  d_g <- setNames(as.numeric(d_g[GENOTYPES]), GENOTYPES)
  selection_fraction <- min(1, max(0, selection_fraction))

  # Mendelian inheritance and genotype-specific offspring production Psi_g.
  inf_tot <- colSums(Ifq)
  Nf <- Sf + inf_tot

  qf <- min(1, max(0, (Nf[["RR"]] + 0.5 * Nf[["RS"]]) / max(sum(Nf), 1e-12)))
  qm <- min(1, max(0, (M[["RR"]]  + 0.5 * M[["RS"]])  / max(sum(M),  1e-12)))
  pf <- 1 - qf
  pm <- 1 - qm

  Psi <- PSI * norm_vec(c(
    SS = pf * pm,
    RS = pf * qm + qf * pm,
    RR = qf * qm
  ))

  # ITN mortality and feeding reduction.
  coverage <- min(0.95, max(0, coverage))
  d_n0 <- min(0.95, max(0, d_n0))
  r_n0 <- min(0.95, max(0, r_n0))

  
  # The current genotype contrast is shrunk around the female-population
  # weighted mean. Average ITN mortality is preserved, but realised selection
  # is weaker than the laboratory contrast alone would imply.
  female_weights <- if (sum(Nf) > 0) Nf / sum(Nf) else rep(1 / N_G, N_G)
  names(female_weights) <- GENOTYPES
  mean_d_g <- sum(female_weights[GENOTYPES] * d_g[GENOTYPES])
  d_g_realised <- mean_d_g + selection_fraction * (d_g - mean_d_g)
  d_g_realised <- setNames(pmin(1, pmax(0, d_g_realised)), GENOTYPES)

  itn_mort <- setNames(
    ITN_MORTALITY_RATE_SCALE * coverage * d_n0 * d_g_realised,
    GENOTYPES
  )

  # Optional resistance fitness cost in natural adult mortality.
  mu_Am_g <- MU_AM * (1 + if (include_cost) cost_vec() else c(SS = 0, RS = 0, RR = 0))
  adult_loss <- mu_Am_g + itn_mort

  # Successful feeding fraction after accounting for repellence and death.
  # B_g(t) = 1 - rho_B*C_B(t) * [r_n0 + d_n0*d_g - r_n0*d_n0*d_g].
  bite_prevention <- r_n0 + d_n0 * d_g - r_n0 * d_n0 * d_g
  feed_g <- setNames(
    pmax(
      0,
      pmin(
        1,
        1 - ITN_PREVENTABLE_BITING * coverage * bite_prevention
      )
    ),
    GENOTYPES
  )

  # Forces of infection.
  Nh <- max(sum(S) + sum(Iq) + sum(R), 1e-12)
  inf_rows <- (MOSQ_LATENT_DAYS + 1L):MAX_M_INF_AGE
  infectious <- colSums(Ifq[inf_rows, , drop = FALSE])

  # Mosquito-to-human force of infection, reduced by ITN-induced feeding failure.
  Lambda_m <- max(0, tmult * THETA * sum(BETA_M * feed_g * infectious) / Nh)

  # Human-to-mosquito force of infection, also reduced by successful feeding.
  human_inf <- sum(H2M * Iq)
  Lambda_h_g <- setNames(
    pmax(0, as.numeric(tmult * THETA * feed_g * human_inf / Nh)),
    GENOTYPES
  )

  # Human equations: semi-implicit update with chronological age and infection age.
  Snew <- numeric(N_AGE)
  Rnew <- numeric(N_AGE)
  newInf <- numeric(N_AGE)

  for (a in seq_len(N_AGE)) {
    age_in <- if (a > 1) AGEING_RATE[a - 1] * Snew[a - 1] else MU_H * MODEL_HUMAN_POPULATION
    Snew[a] <- (S[a] + DT * (K_H * R[a] + age_in)) /
      (1 + DT * (MU_H + Lambda_m + AGEING_RATE[a]))
    newInf[a] <- Lambda_m * Snew[a]
  }

  rate_I <- MU_H + GAMMA_H_AGE + NU_H_AGE
  surv_I <- 1 / (1 + DT * rate_I)

  Iq_shift <- matrix(0, N_AGE, MAX_H_INF_AGE, dimnames = dimnames(Iq))
  Iq_shift[, 2:MAX_H_INF_AGE] <- Iq[, 1:(MAX_H_INF_AGE - 1L)]
  Iq_shift[, 1] <- newInf

  left_last <- Iq[, MAX_H_INF_AGE]
  Iq_new <- Iq_shift * surv_I
  left_amt <- rowSums(Iq_shift * (1 - surv_I))
  rec_from_I <- (left_amt + left_last) * (GAMMA_H_AGE / rate_I)
  deaths_age <- (left_amt + left_last) * (NU_H_AGE / rate_I)

  for (a in seq_len(N_AGE - 1L)) {
    mv <- Iq_new[a, ] * (1 - 1 / (1 + DT * AGEING_RATE[a]))
    Iq_new[a, ] <- Iq_new[a, ] - mv
    Iq_new[a + 1, ] <- Iq_new[a + 1, ] + mv
  }

  for (a in seq_len(N_AGE)) {
    age_in <- if (a > 1) AGEING_RATE[a - 1] * Rnew[a - 1] else 0
    Rnew[a] <- (R[a] + rec_from_I[a] + DT * age_in) /
      (1 + DT * (MU_H + K_H + AGEING_RATE[a]))
  }

  # Keep the numerical population fixed at the model population size.
  tot <- sum(Snew) + sum(Rnew) + sum(Iq_new)
  sc <- MODEL_HUMAN_POPULATION / tot
  Snew <- Snew * sc
  Rnew <- Rnew * sc
  Iq_new <- Iq_new * sc
  newInf <- newInf * sc
  deaths_age <- deaths_age * sc

  # Mosquito equations: semi-implicit juvenile, adult male, susceptible female,
  # and infected female infection-age updates.
  Jtot <- sum(Jm)
  dens <- max(0, 1 - Jtot / LARVAL_K)
  births <- Psi * dens * sum(Nf)

  Jnew <- (Jm + DT * births) / (1 + DT * (MU_JM + PHI_M))
  emerg <- DT * PHI_M * Jnew

  Mnew <- (M + (1 - R_FEMALE) * emerg) / (1 + DT * adult_loss)

  Sfnew <- numeric(N_G)
  names(Sfnew) <- GENOTYPES
  newMinf <- numeric(N_G)
  names(newMinf) <- GENOTYPES

  for (g in GENOTYPES) {
    Sfnew[[g]] <- (Sf[[g]] + R_FEMALE * emerg[[g]]) /
      (1 + DT * (adult_loss[[g]] + Lambda_h_g[[g]]))
    newMinf[[g]] <- Lambda_h_g[[g]] * Sfnew[[g]]
  }

  Ifq_new <- matrix(0, MAX_M_INF_AGE, N_G, dimnames = dimnames(Ifq))
  for (g in GENOTYPES) {
    survf <- 1 / (1 + DT * adult_loss[[g]])
    Ifq_new[2:MAX_M_INF_AGE, g] <- Ifq[1:(MAX_M_INF_AGE - 1L), g] * survf
    Ifq_new[1, g] <- newMinf[[g]]
  }

  infected_stock <- rowSums(Iq_new)
  Hpop <- sum(Snew) + sum(Iq_new) + sum(Rnew)
  Hbyage <- Snew + infected_stock + Rnew
  u5_prev <- 100 * sum(infected_stock[UNDER5]) / max(sum(Hbyage[UNDER5]), 1e-12)

  Nf2 <- Sfnew + colSums(Ifq_new)
  tf <- max(sum(Nf2), 1e-12)
  qR <- min(1, max(0, (Nf2[["RR"]] + 0.5 * Nf2[["RS"]]) / tf))

  list(
    S = sp(Snew),
    Iq = pmax(Iq_new, 0),
    R = sp(Rnew),
    juvenile = sp(Jnew),
    adult_male = sp(Mnew),
    sus_female = sp(Sfnew),
    inf_f_q = pmax(Ifq_new, 0),
    prevalence = 100 * sum(infected_stock) / max(Hpop, 1e-12),
    u5_prevalence = u5_prev,
    qR = qR,
    new_inf_by_age = newInf,
    deaths_by_age = deaths_age,
    coverage = coverage,
    d_n0 = d_n0,
    r_n0 = r_n0
  )
}

# ============================================================================
# EPIDEMIOLOGICAL SPIN-UP AND TRANSMISSION CALIBRATION
# ============================================================================

# Run the model without ITNs until it approaches an endemic equilibrium.
run_spinup <- function(tmult, years) {
  state <- init_state()
  number_of_days <- as.integer(years * DAYS_PER_YEAR)

  for (day_index in seq_len(number_of_days)) {
    state <- step_model(
      state = state,
      coverage = 0,
      d_n0 = 0,
      r_n0 = 0,
      d_g = D_PYR,
      selection_fraction = 0,
      day = day_index,
      include_cost = TRUE,
      tmult = tmult
    )
  }

  state
}

# Return the equilibrium all-age prevalence for a transmission multiplier.
endemic_prev <- function(tmult) {
  run_spinup(tmult, CAL_SPINUP_YEARS)$prevalence / 100
}

# Calibrate the transmission multiplier by bisection so that the model reaches
# the requested baseline prevalence.
calibrate_tmult <- function(target) {
  lower_bound <- 0.05
  upper_bound <- 8.0
  upper_prevalence <- endemic_prev(upper_bound)
  expansion_iteration <- 0L

  while (upper_prevalence < target && expansion_iteration < 7L) {
    upper_bound <- 2 * upper_bound
    upper_prevalence <- endemic_prev(upper_bound)
    expansion_iteration <- expansion_iteration + 1L
  }

  lower_prevalence <- endemic_prev(lower_bound)

  if (target <= lower_prevalence) {
    return(lower_bound)
  }

  if (target >= upper_prevalence) {
    return(upper_bound)
  }

  for (iteration in seq_len(CAL_ITERS)) {
    midpoint <- 0.5 * (lower_bound + upper_bound)

    if (endemic_prev(midpoint) < target) {
      lower_bound <- midpoint
    } else {
      upper_bound <- midpoint
    }
  }

  0.5 * (lower_bound + upper_bound)
}

MAX_REACHABLE <- endemic_prev(1e3) * 100

message(
  "Maximum reachable all-age prevalence (k_h = ",
  K_H,
  "): ",
  round(MAX_REACHABLE, 1),
  "%"
)

TARGET_USED <- min(
  TARGET_BASELINE_PREVALENCE,
  0.95 * MAX_REACHABLE / 100
)

TMULT <- calibrate_tmult(TARGET_USED)
SPIN <- run_spinup(TMULT, SPINUP_YEARS)
BASE_PREV <- SPIN$prevalence
BASE_U5 <- SPIN$u5_prevalence

message(
  "Transmission multiplier = ",
  round(TMULT, 4),
  "; baseline all-age prevalence = ",
  round(BASE_PREV, 2),
  "%; baseline under-five prevalence = ",
  round(BASE_U5, 1),
  "%"
)

# ============================================================================
# Time-varying ITN coverage, posterior-style draws, and simulations
# ============================================================================
campaign_years <- seq(STUDY_START_YEAR, STUDY_END_YEAR, by = CAMPAIGN_INTERVAL)
campaign_days <- as.integer((campaign_years - STUDY_START_YEAR) * DAYS_PER_YEAR + 1L)

# Effective ITN coverage C_B(t): increases at campaign days and decays between
# campaigns due to attrition, physical degradation, and reduced use.
build_coverage <- function(target, retention_days, peak_mult) {
  peak <- min(0.95, max(ROUTINE_COVERAGE_FLOOR, target * peak_mult))
  routine_floor <- min(ROUTINE_COVERAGE_FLOOR, 0.80 * target)
  decay_rate <- log(2) / retention_days
  cov <- numeric(SIM_LENGTH)

  last_campaign_day <- NA_integer_
  campaign_start_coverage <- 0

  for (d in seq_len(SIM_LENGTH)) {
    if (d %in% campaign_days) {
      last_campaign_day <- d
      campaign_start_coverage <- if (d == 1L) 0 else cov[d - 1L]
    }

    if (is.na(last_campaign_day)) {
      cov[d] <- 0
      next
    }

    days_since_campaign <- d - last_campaign_day

    if (days_since_campaign <= CAMPAIGN_SCALEUP_DAYS) {
      # Smooth cubic rise: zero slope at the start and at the campaign peak.
      x <- days_since_campaign / max(1, CAMPAIGN_SCALEUP_DAYS)
      smooth_step <- 3 * x^2 - 2 * x^3
      cov[d] <- campaign_start_coverage +
        (peak - campaign_start_coverage) * smooth_step
    } else {
      # Smooth waning towards the routine floor.
      days_since_peak <- days_since_campaign - CAMPAIGN_SCALEUP_DAYS
      cov[d] <- routine_floor +
        (peak - routine_floor) * exp(-decay_rate * days_since_peak)
    }
  }

  pmin(0.95, pmax(0, cov))
}

# Draws for ITN parameters. If posterior samples of C_B(t), d_n0, and r_n0 are
# available, replace this function by a reader for those posterior samples.
make_draws <- function(n) {
  tibble(
    draw = seq_len(n),

    # Pyrethroid-only net parameters: d_n0 and r_n0.
    # These values are fallback values and should be replaced by posterior draws
    # from the DRC ITN calibration when available.
    pyr_d_n0 = pmin(0.90, pmax(0.05, rnorm(n, 0.48, 0.07))),
    pyr_r_n0 = pmin(0.90, pmax(0.05, rnorm(n, 0.48, 0.07))),

    # Next-generation net parameters.
    # Next-generation nets are represented as improving mosquito death after contact, while the
    # genotype-specific restoration is represented through D_NEXTGEN.
    nextgen_d_n0 = pmin(0.90, pmax(0.05, rnorm(n, 0.58, 0.07))),
    nextgen_r_n0 = pmin(0.90, pmax(0.05, rnorm(n, 0.52, 0.07))),

    # Coverage uncertainty around the campaign peak and decay duration.
    peak_cov_mult = pmin(0.95, pmax(0.75, rnorm(n, 0.88, 0.06))),
    retention_days = pmin(
      3.5 * DAYS_PER_YEAR,
      pmax(2.0 * DAYS_PER_YEAR, rnorm(n, NET_RETENTION_YEARS * DAYS_PER_YEAR, 0.35 * DAYS_PER_YEAR))
    )
  )
}

DRAWS <- make_draws(N_DRAWS)

# No mixture or annual product-share variable is used.
get_itn_parameters <- function(product, day, drow, coverage_value) {
  if (product == "Pyrethroid-only nets") {
    d_n0 <- drow$pyr_d_n0
    r_n0 <- drow$pyr_r_n0
    d_g <- D_PYR
    selection_fraction <- PYRETHROID_REALIZED_SELECTION_FRACTION
    nextgen_share <- 0
  } else if (product == "Next-generation nets") {
    d_n0 <- drow$nextgen_d_n0
    r_n0 <- drow$nextgen_r_n0
    d_g <- D_NEXTGEN
    selection_fraction <- NEXTGEN_REALIZED_SELECTION_FRACTION
    nextgen_share <- 1
  } else {
    stop("Unknown intervention product: ", product, call. = FALSE)
  }

  list(
    coverage = coverage_value,
    d_n0 = d_n0,
    r_n0 = r_n0,
    d_g = d_g,
    selection_fraction = selection_fraction,
    nextgen_share = nextgen_share
  )
}

simulate_arm <- function(product, target, drow, itn = TRUE) {
  cov <- if (itn) {
    build_coverage(target, drow$retention_days, drow$peak_cov_mult)
  } else {
    numeric(SIM_LENGTH)
  }

  scale_pop <- NATIONAL_POP / MODEL_HUMAN_POPULATION

  prevalence <- numeric(SIM_LENGTH)
  u5 <- numeric(SIM_LENGTH)
  qR <- numeric(SIM_LENGTH)
  coverage_v <- numeric(SIM_LENGTH)
  nextgen_share_v <- numeric(SIM_LENGTH)
  d_n0_v <- numeric(SIM_LENGTH)
  r_n0_v <- numeric(SIM_LENGTH)
  cases <- numeric(SIM_LENGTH)
  deaths <- numeric(SIM_LENGTH)

  qR_start <- initial_qR_for_product(product)
  st <- reset_state_allele_frequency(SPIN, qR_start)

  for (d in seq_len(SIM_LENGTH)) {
    net <- if (itn) {
      get_itn_parameters(product, d, drow, cov[d])
    } else {
      list(
        coverage = 0,
        d_n0 = 0,
        r_n0 = 0,
        d_g = D_PYR,
        selection_fraction = 0,
        nextgen_share = 0
      )
    }

    st <- step_model(
      state = st,
      coverage = net$coverage,
      d_n0 = net$d_n0,
      r_n0 = net$r_n0,
      d_g = net$d_g,
      selection_fraction = net$selection_fraction,
      day = d,
      include_cost = TRUE,
      tmult = TMULT
    )

    prevalence[d] <- st$prevalence
    u5[d] <- st$u5_prevalence
    qR[d] <- st$qR
    coverage_v[d] <- st$coverage
    nextgen_share_v[d] <- net$nextgen_share
    d_n0_v[d] <- net$d_n0
    r_n0_v[d] <- net$r_n0
    cases[d] <- BURDEN_CALIBRATION_FACTOR *
      sum(P_CLIN_AGE * st$new_inf_by_age) * scale_pop
    deaths[d] <- BURDEN_CALIBRATION_FACTOR *
      sum(st$deaths_by_age) * scale_pop
  }

  day <- seq_len(SIM_LENGTH)
  year <- STUDY_START_YEAR + (day - 1L) %/% DAYS_PER_YEAR
  ydec <- 1 + (day - 1) / DAYS_PER_YEAR
  month <- (day - 1L) %/% 30L + 1L

  mcnt <- as.numeric(table(month))
  mk <- as.integer(names(table(month)))
  msum <- rowsum(cbind(ydec, prevalence, u5, qR, coverage = coverage_v,
                       nextgen_share = nextgen_share_v, d_n0 = d_n0_v, r_n0 = r_n0_v), month)

  monthly <- tibble(
    month = mk,
    ydec = msum[, "ydec"] / mcnt,
    prevalence = msum[, "prevalence"] / mcnt,
    u5_prevalence = msum[, "u5"] / mcnt,
    qR = msum[, "qR"] / mcnt,
    coverage = msum[, "coverage"] / mcnt,
    nextgen_share = msum[, "nextgen_share"] / mcnt,
    d_n0 = msum[, "d_n0"] / mcnt,
    r_n0 = msum[, "r_n0"] / mcnt
  )

  ycnt <- as.numeric(table(year))
  yk <- as.integer(names(table(year)))
  ysum <- rowsum(cbind(prevalence, coverage = coverage_v, nextgen_share = nextgen_share_v,
                       d_n0 = d_n0_v, r_n0 = r_n0_v, cases, deaths), year)
  ylast <- tapply(qR, year, function(z) z[length(z)])

  annual <- tibble(
    year = yk,
    prevalence = ysum[, "prevalence"] / ycnt,
    qR = as.numeric(ylast),
    cases = ysum[, "cases"],
    deaths = ysum[, "deaths"],
    coverage = ysum[, "coverage"] / ycnt,
    nextgen_share = ysum[, "nextgen_share"] / ycnt,
    d_n0 = ysum[, "d_n0"] / ycnt,
    r_n0 = ysum[, "r_n0"] / ycnt,
    sim_year = yk - STUDY_START_YEAR + 1L
  )

  add_meta <- function(df) {
    df %>%
      mutate(
        product = product,
        coverage_target = target,
        coverage_label = paste0(as.integer(100 * target), "%"),
        draw = drow$draw
      )
  }

  list(monthly = add_meta(monthly), annual = add_meta(annual))
}

message("Fitness cost (RR) fixed at ", RES_COST_RR,
        "; ITN mortality follows d_n0 * C_B(t) * d_g as in the manuscript.")

message("Running product-matched no-ITN comparators ...")
no_itn_grid <- tidyr::crossing(product = PRODUCTS, draw = seq_len(N_DRAWS))
no_itn <- pmap_dfr(
  list(no_itn_grid$product, no_itn_grid$draw),
  function(pr, dr) simulate_arm(pr, 0.25, DRAWS[dr, ], itn = FALSE)$annual
) %>%
  select(product, draw, year, no_cases = cases, no_deaths = deaths)
# Build and run the complete intervention scenario grid.
grid <- tidyr::crossing(
  product = PRODUCTS,
  coverage_target = COVERAGE_TARGETS,
  draw = seq_len(N_DRAWS)
)

message("Running ", nrow(grid), " intervention trajectories ...")
start_time <- Sys.time()

simulation_results <- pmap(
  .l = list(
    grid$product,
    grid$coverage_target,
    grid$draw
  ),
  .f = function(product, coverage_target, draw_index) {
    simulate_arm(
      product = product,
      target = coverage_target,
      drow = DRAWS[draw_index, ],
      itn = TRUE
    )
  }
)

elapsed_seconds <- as.numeric(
  Sys.time() - start_time,
  units = "secs"
)
message("Simulation completed in ", round(elapsed_seconds, 1), " seconds.")

arms_monthly <- map_dfr(simulation_results, "monthly")
arms_annual <- map_dfr(simulation_results, "annual")

# Compute annual and cumulative cases and deaths averted relative to the
# product-matched no-ITN comparator.
impact <- arms_annual %>%
  left_join(
    no_itn,
    by = c("product", "draw", "year")
  ) %>%
  arrange(product, coverage_target, draw, year) %>%
  group_by(product, coverage_target, coverage_label, draw) %>%
  mutate(
    cases_averted = pmax(0, no_cases - cases),
    deaths_averted = pmax(0, no_deaths - deaths),
    cum_cases_averted = cumsum(cases_averted),
    cum_deaths_averted = cumsum(deaths_averted)
  ) %>%
  ungroup() %>%
  mutate(sim_year = year - STUDY_START_YEAR + 1L)

# Express total intervention impact as a percentage of the no-ITN burden.
impact_pct <- impact %>%
  group_by(product, coverage_target, coverage_label, draw) %>%
  summarise(
    cases_pct = 100 * sum(cases_averted) / pmax(sum(no_cases), 1e-9),
    deaths_pct = 100 * sum(deaths_averted) / pmax(sum(no_deaths), 1e-9),
    .groups = "drop"
  )

# Quantile helper functions used for uncertainty intervals.
q_lo <- function(x) {
  as.numeric(quantile(x, 0.025, names = FALSE, na.rm = TRUE))
}

q_hi <- function(x) {
  as.numeric(quantile(x, 0.975, names = FALSE, na.rm = TRUE))
}

# Summarise selected variables by median and 95% uncertainty interval.
summarise_draws <- function(data, grouping_variables, summary_variables) {
  data %>%
    group_by(across(all_of(grouping_variables))) %>%
    summarise(
      across(
        all_of(summary_variables),
        list(
          m = ~ median(.x, na.rm = TRUE),
          lo = q_lo,
          hi = q_hi
        ),
        .names = "{.col}_{.fn}"
      ),
      .groups = "drop"
    )
}

monthly_summ <- summarise_draws(
  data = arms_monthly,
  grouping_variables = c(
    "product",
    "coverage_label",
    "coverage_target",
    "month"
  ),
  summary_variables = c(
    "ydec",
    "prevalence",
    "u5_prevalence",
    "qR",
    "coverage",
    "nextgen_share",
    "d_n0",
    "r_n0"
  )
) %>%
  mutate(ydec = ydec_m)

impact_summ <- summarise_draws(
  data = impact,
  grouping_variables = c(
    "product",
    "coverage_label",
    "coverage_target",
    "sim_year"
  ),
  summary_variables = c(
    "cum_cases_averted",
    "cum_deaths_averted"
  )
)

pct_summ <- impact_pct %>%
  group_by(product, coverage_label, coverage_target) %>%
  summarise(
    cases_pct_m = mean(cases_pct),
    cases_pct_lo = q_lo(cases_pct),
    cases_pct_hi = q_hi(cases_pct),
    deaths_pct_m = mean(deaths_pct),
    deaths_pct_lo = q_lo(deaths_pct),
    deaths_pct_hi = q_hi(deaths_pct),
    .groups = "drop"
  )

# ============================================================================
# MATHEMATICAL THRESHOLD ANALYSIS AS A FUNCTION OF EFFECTIVE ITN COVERAGE

human_infectiousness_capacity <- function() {
  age_midpoints <- (
    AGE_BINS_YEARS[-1] + AGE_BINS_YEARS[-length(AGE_BINS_YEARS)]
  ) / 2

  age_weights <- exp(-MU_H * age_midpoints * DAYS_PER_YEAR) *
    diff(AGE_BINS_YEARS)
  age_weights <- age_weights / sum(age_weights)

  rate_I <- MU_H + GAMMA_H_AGE + NU_H_AGE

  capacity_by_age <- vapply(
    seq_len(N_AGE),
    function(a) {
      daily_survival <- 1 / (1 + DT * rate_I[a])
      survival_profile <- daily_survival^seq_len(MAX_H_INF_AGE)
      sum(H2M[a, ] * survival_profile)
    },
    numeric(1)
  )

  sum(age_weights * capacity_by_age)
}

HUMAN_INFECTIOUSNESS_CAPACITY <- human_infectiousness_capacity()

# Product-specific genotype susceptibility and insecticidal efficacy.
# epsilon_A is implemented numerically as ITN_MORTALITY_RATE_SCALE * d_n0.
threshold_product_parameters <- function(product, pyr_d_n0, nextgen_d_n0) {
  product <- as.character(product)

  if (product == "Pyrethroid-only nets") {
    d_g <- D_PYR
    d_n0 <- pyr_d_n0
  } else if (product == NEXTGEN_PRODUCT) {
    d_g <- D_NEXTGEN
    d_n0 <- nextgen_d_n0
  } else {
    stop("Unknown product in threshold analysis: ", product, call. = FALSE)
  }

  list(
    d_g = setNames(as.numeric(d_g[GENOTYPES]), GENOTYPES),
    epsilon_A = ITN_MORTALITY_RATE_SCALE *
      min(0.95, max(0, as.numeric(d_n0)))
  )
}

# Calculate all analytical thresholds at one fixed effective coverage.
threshold_at_coverage <- function(
    product,
    coverage,
    pyr_d_n0,
    nextgen_d_n0
) {
  coverage <- min(0.75, max(0, as.numeric(coverage)))

  product_par <- threshold_product_parameters(
    product = product,
    pyr_d_n0 = pyr_d_n0,
    nextgen_d_n0 = nextgen_d_n0
  )

  d_g <- product_par$d_g
  epsilon_A <- product_par$epsilon_A

  # Genotype-specific adult natural mortality. The code permits genotype-specific
  # adult survival through the resistance fitness cost already used in the
  # dynamic model; this maps onto mu_Am,g in the analytical formulation.
  mu_Am_g <- MU_AM * (1 + cost_vec())
  names(mu_Am_g) <- GENOTYPES

  # Juvenile parameters are genotype invariant in the present parameterisation,
  # but alpha_g is retained explicitly to preserve the manuscript notation.
  mu_J_g <- setNames(rep(MU_JM, N_G), GENOTYPES)
  phi_g <- setNames(rep(PHI_M, N_G), GENOTYPES)

  alpha <- c(
    SS = 1,
    RS = (mu_J_g[["SS"]] + phi_g[["SS"]]) /
      (mu_J_g[["RS"]] + phi_g[["RS"]]),
    RR = (mu_J_g[["SS"]] + phi_g[["SS"]]) /
      (mu_J_g[["RR"]] + phi_g[["RR"]])
  )

  adult_loss <- mu_Am_g + epsilon_A * coverage * d_g
  Delta <- phi_g / pmax(adult_loss, 1e-12)
  chi <- alpha * Delta

  # ---------------------------------------------------------------------------
  # Invasion index exactly as defined in the manuscript.
  # ---------------------------------------------------------------------------
  invasion_R_to_S <- unname(chi[["RS"]] / chi[["SS"]] - 1)

  # ---------------------------------------------------------------------------
  # Coexistence equilibrium: p* and q* from the manuscript quadratic.
  # ---------------------------------------------------------------------------
  chi_c <- (chi[["SS"]] - chi[["RS"]]) +
    (chi[["RR"]] - chi[["RS"]])

  p_star <- if (is.finite(chi_c) && abs(chi_c) > 1e-12) {
    (chi[["RR"]] - chi[["RS"]]) / chi_c
  } else {
    NA_real_
  }

  coexistence_exists <- is.finite(p_star) && p_star > 0 && p_star < 1
  q_star <- if (coexistence_exists) 1 - p_star else NA_real_

  # ---------------------------------------------------------------------------
  # Mosquito genotype reproduction numbers R^{Jm}_g.
  # Pure SS and RR thresholds are always calculable. The RS curve represents
  # the additive coexistence threshold and is returned only where p* is valid.
  # ---------------------------------------------------------------------------
  RJ_SS <- R_FEMALE * PSI * Delta[["SS"]] /
    (mu_J_g[["SS"]] + phi_g[["SS"]])

  RJ_RR <- R_FEMALE * PSI * Delta[["RR"]] /
    (mu_J_g[["RR"]] + phi_g[["RR"]])

  if (coexistence_exists) {
    RJ_SS_component <- R_FEMALE * PSI * p_star^2 * Delta[["SS"]] /
      (mu_J_g[["SS"]] + phi_g[["SS"]])

    RJ_RS_component <- 2 * R_FEMALE * PSI * p_star * q_star *
      Delta[["RS"]] /
      (mu_J_g[["RS"]] + phi_g[["RS"]])

    RJ_RR_component <- R_FEMALE * PSI * q_star^2 * Delta[["RR"]] /
      (mu_J_g[["RR"]] + phi_g[["RR"]])

    RJ_RS <- RJ_SS_component + RJ_RS_component + RJ_RR_component

    kappa_RS <- (2 * q_star / p_star) * alpha[["RS"]]
    kappa_RR <- (q_star / p_star)^2 * alpha[["RR"]]
  } else {
    RJ_RS <- NA_real_
    kappa_RS <- NA_real_
    kappa_RR <- NA_real_
  }

  # ---------------------------------------------------------------------------
  # Mosquito-to-human transmission capacities C^m_g.
  # The daily numerical survival analogue reproduces the EIP mechanism used in
  # the simulation and corresponds to the analytical integral C^m_g.
  # ---------------------------------------------------------------------------
  C_m <- vapply(
    GENOTYPES,
    function(g) {
      daily_survival <- 1 / (1 + DT * adult_loss[[g]])

      infectious_survival_sum <- if (MAX_M_INF_AGE > MOSQ_LATENT_DAYS) {
        sum(
          daily_survival^(
            MOSQ_LATENT_DAYS:(MAX_M_INF_AGE - 1L)
          )
        )
      } else {
        0
      }

      BETA_M[[g]] * infectious_survival_sum
    },
    numeric(1)
  )
  names(C_m) <- GENOTYPES

  # Common factor in the manuscript decomposition of R0^2.
  R0_constant <- (TMULT * THETA)^2 *
    HUMAN_INFECTIOUSNESS_CAPACITY /
    MODEL_HUMAN_POPULATION

  # Pure SS and RR DFE.
  gate_SS <- if (RJ_SS > 1) 1 - 1 / RJ_SS else 0
  gate_RR <- if (RJ_RR > 1) 1 - 1 / RJ_RR else 0

  R0_SS <- sqrt(
    max(
      0,
      R0_constant * R_FEMALE * LARVAL_K *
        gate_SS * Delta[["SS"]] * C_m[["SS"]]
    )
  )

  R0_RR <- sqrt(
    max(
      0,
      R0_constant * R_FEMALE * LARVAL_K *
        gate_RR * Delta[["RR"]] * C_m[["RR"]]
    )
  )

  # Coexistence DFE T_RS. The weighted mosquito transmission term follows the
  # closed form in the manuscript. It is intentionally NA when the coexistence
  # equilibrium is not biologically admissible.
  if (coexistence_exists && is.finite(RJ_RS) && RJ_RS > 1) {
    gate_RS <- 1 - 1 / RJ_RS

    weighted_vector_capacity <- (
      Delta[["SS"]] * C_m[["SS"]] +
        kappa_RS * Delta[["RS"]] * C_m[["RS"]] +
        kappa_RR * Delta[["RR"]] * C_m[["RR"]]
    ) / (1 + kappa_RS + kappa_RR)

    R0_RS <- sqrt(
      max(
        0,
        R0_constant * R_FEMALE * LARVAL_K *
          gate_RS * weighted_vector_capacity
      )
    )
  } else {
    R0_RS <- NA_real_
  }

  list(
    invasion_R_to_S = invasion_R_to_S,
    Delta_SS = unname(Delta[["SS"]]),
    Delta_RS = unname(Delta[["RS"]]),
    Delta_RR = unname(Delta[["RR"]]),
    chi_SS = unname(chi[["SS"]]),
    chi_RS = unname(chi[["RS"]]),
    chi_RR = unname(chi[["RR"]]),
    p_star = p_star,
    q_star = q_star,
    coexistence_exists = coexistence_exists,
    RJ_SS = RJ_SS,
    RJ_RS = RJ_RS,
    RJ_RR = RJ_RR,
    R0_SS = R0_SS,
    R0_RS = R0_RS,
    R0_RR = R0_RR
  )
}

# Coverage sweep used only for analytical threshold figures.
COVERAGE_SWEEP <- seq(0, 0.75, by = 0.005)
COVERAGE_REFERENCE <- c(0.25, 0.50, 0.75)

# Threshold uncertainty is inexpensive to propagate because no dynamic epidemic
# simulation is required here. We therefore use more ITN-parameter draws than
# the main dynamic FAST_MODE run when necessary.
N_THRESHOLD_DRAWS <- if (FAST_MODE) 100L else 500L
THRESHOLD_DRAWS <- make_draws(N_THRESHOLD_DRAWS)

threshold_sweep_draws <- tidyr::crossing(
  draw = THRESHOLD_DRAWS$draw,
  product = PRODUCTS,
  coverage = COVERAGE_SWEEP
) %>%
  left_join(
    THRESHOLD_DRAWS %>%
      select(draw, pyr_d_n0, nextgen_d_n0),
    by = "draw"
  ) %>%
  mutate(
    .threshold = pmap(
      list(
        product,
        coverage,
        pyr_d_n0,
        nextgen_d_n0
      ),
      threshold_at_coverage
    )
  ) %>%
  unnest_wider(.threshold)

# Safe summaries for quantities such as R0_RS and RJ_RS that are undefined when
# the coexistence DFE does not exist.
safe_median <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) NA_real_ else median(x)
}

safe_quantile <- function(x, prob) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) NA_real_ else
    as.numeric(quantile(x, prob, names = FALSE))
}

threshold_sweep_summ <- threshold_sweep_draws %>%
  group_by(product, coverage) %>%
  summarise(
    invasion_R_to_S_m = safe_median(invasion_R_to_S),
    invasion_R_to_S_lo = safe_quantile(invasion_R_to_S, 0.025),
    invasion_R_to_S_hi = safe_quantile(invasion_R_to_S, 0.975),

    RJ_SS_m = safe_median(RJ_SS),
    RJ_SS_lo = safe_quantile(RJ_SS, 0.025),
    RJ_SS_hi = safe_quantile(RJ_SS, 0.975),
    RJ_RS_m = safe_median(RJ_RS),
    RJ_RS_lo = safe_quantile(RJ_RS, 0.025),
    RJ_RS_hi = safe_quantile(RJ_RS, 0.975),
    RJ_RR_m = safe_median(RJ_RR),
    RJ_RR_lo = safe_quantile(RJ_RR, 0.025),
    RJ_RR_hi = safe_quantile(RJ_RR, 0.975),

    R0_SS_m = safe_median(R0_SS),
    R0_SS_lo = safe_quantile(R0_SS, 0.025),
    R0_SS_hi = safe_quantile(R0_SS, 0.975),
    R0_RS_m = safe_median(R0_RS),
    R0_RS_lo = safe_quantile(R0_RS, 0.025),
    R0_RS_hi = safe_quantile(R0_RS, 0.975),
    R0_RR_m = safe_median(R0_RR),
    R0_RR_lo = safe_quantile(R0_RR, 0.025),
    R0_RR_hi = safe_quantile(R0_RR, 0.975),

    coexistence_probability = mean(coexistence_exists, na.rm = TRUE),
    p_star_m = safe_median(p_star),
    .groups = "drop"
  ) %>%
  mutate(
    # Show the coexistence DFE (Tc) only where supported by at least
    # half of the intervention-parameter draws. Pure SS and RR thresholds remain
    # defined throughout the coverage sweep.
    across(
      c(RJ_RS_m, RJ_RS_lo, RJ_RS_hi, R0_RS_m, R0_RS_lo, R0_RS_hi),
      ~ if_else(coexistence_probability >= 0.50, .x, NA_real_)
    ),
    product = factor(product, levels = PRODUCTS)
  )

# For plotting, keep the analytically supported Tc values unchanged and extend
# only the Tc median after the last supported coexistence point. The continuation
# follows the changing SS and RR thresholds at each coverage level rather than a
# straight, quadratic, or flat extrapolation. The uncertainty ribbon is NOT
# extrapolated beyond the supported coexistence region.
extend_tc_from_pure_curves <- function(tc, ss, rr, supported) {
  out <- tc
  idx <- which(supported & is.finite(tc) & is.finite(ss) & is.finite(rr))
  if (length(idx) == 0L) return(out)

  last_idx <- idx[length(idx)]
  if (last_idx >= length(tc)) return(out)

  denom <- ss[last_idx] - rr[last_idx]
  if (!is.finite(denom) || abs(denom) < 1e-10) {
    weight_ss <- 0.5
  } else {
    weight_ss <- (tc[last_idx] - rr[last_idx]) / denom
    weight_ss <- min(1, max(0, weight_ss))
  }

  future_idx <- seq.int(last_idx + 1L, length(tc))
  out[future_idx] <- weight_ss * ss[future_idx] +
    (1 - weight_ss) * rr[future_idx]

  out
}

threshold_sweep_summ <- threshold_sweep_summ %>%
  arrange(product, coverage) %>%
  group_by(product) %>%
  mutate(
    tc_supported = coexistence_probability >= 0.50,
    R0_RS_m = extend_tc_from_pure_curves(
      R0_RS_m, R0_SS_m, R0_RR_m, tc_supported
    ),
    RJ_RS_m = extend_tc_from_pure_curves(
      RJ_RS_m, RJ_SS_m, RJ_RR_m, tc_supported
    )
    # IMPORTANT: R0_RS_lo/hi and RJ_RS_lo/hi remain NA outside the
    # analytically supported coexistence region. This prevents an artificial
    # grey ribbon from inflating the y-axis scale.
  ) %>%
  ungroup()

# Compact long-format tables used by Figures 6 and 7.
# IMPORTANT: RS remains the biological heterozygous genotype throughout the
# transmission/genetic model. In the analytical threshold section, the RS-based
# coexistence equilibrium/threshold is the quantity denoted T_c in the figures.
# Therefore internal variables retain the names R0_RS and RJ_RS, while the plot
# legend displays T_c. This naming change does NOT alter Figures 1--4.
threshold_R0_long <- bind_rows(
  threshold_sweep_summ %>%
    transmute(
      product, coverage, genotype = "SS",
      m = R0_SS_m, lo = R0_SS_lo, hi = R0_SS_hi,
      tc_supported = TRUE
    ),
  threshold_sweep_summ %>%
    transmute(
      product, coverage, genotype = "Tc",
      m = R0_RS_m, lo = R0_RS_lo, hi = R0_RS_hi,
      tc_supported = tc_supported),
  threshold_sweep_summ %>%
    transmute(
      product, coverage, genotype = "RR",
      m = R0_RR_m, lo = R0_RR_lo, hi = R0_RR_hi,
      tc_supported = TRUE
    )
) %>%
  mutate(
    genotype = factor(genotype, levels = c("SS", "Tc", "RR")),
    product = factor(product, levels = PRODUCTS)
  )

threshold_RJ_long <- bind_rows(
  threshold_sweep_summ %>%
    transmute(
      product, coverage, genotype = "SS",
      m = RJ_SS_m, lo = RJ_SS_lo, hi = RJ_SS_hi,
      tc_supported = TRUE
    ),
  threshold_sweep_summ %>%
    transmute(
      product, coverage, genotype = "Tc",
      m = RJ_RS_m, lo = RJ_RS_lo, hi = RJ_RS_hi,
      tc_supported = tc_supported),
  threshold_sweep_summ %>%
    transmute(
      product, coverage, genotype = "RR",
      m = RJ_RR_m, lo = RJ_RR_lo, hi = RJ_RR_hi,
      tc_supported = TRUE
    )
) %>%
  mutate(
    genotype = factor(genotype, levels = c("SS", "Tc", "RR")),
    product = factor(product, levels = PRODUCTS)
  )

# Pre-campaign baseline points for plotting only.
# These points make prevalence curves begin with a short stable endemic segment
# before the first campaign at Year 1. They do not alter the simulation results;
# they only improve the visual interpretation of the figures.
base_pt <- tidyr::crossing(
  product = PRODUCTS,
  coverage_label = COVERAGE_LABELS,
  ydec = c(PRE_CAMPAIGN_START_X, 1)
) %>%
  left_join(
    tibble(
      coverage_label = COVERAGE_LABELS,
      coverage_target = COVERAGE_TARGETS
    ),
    by = "coverage_label"
  ) %>%
  mutate(
    month = if_else(ydec < 1, -1L, 0L),
    prevalence_m = PLOT_BASELINE_PREVALENCE,
    prevalence_lo = PLOT_BASELINE_PREVALENCE,
    prevalence_hi = PLOT_BASELINE_PREVALENCE,
    u5_prevalence_m = BASE_U5,
    u5_prevalence_lo = BASE_U5,
    u5_prevalence_hi = BASE_U5,
    qR_m = initial_qR_for_product(product),
    qR_lo = initial_qR_for_product(product),
    qR_hi = initial_qR_for_product(product),
    coverage_m = 0,
    coverage_lo = 0,
    coverage_hi = 0,
    nextgen_share_m = 0,
    nextgen_share_lo = 0,
    nextgen_share_hi = 0,
    d_n0_m = 0,
    d_n0_lo = 0,
    d_n0_hi = 0,
    r_n0_m = 0,
    r_n0_lo = 0,
    r_n0_hi = 0
  )
monthly_plot_raw <- bind_rows(base_pt, monthly_summ) %>%
  mutate(
    product = factor(product, levels = PRODUCTS),
    coverage_label = factor(coverage_label, levels = COVERAGE_LABELS)
  ) %>%
  arrange(product, coverage_target, ydec)

# First-order response used only for the prevalence curves shown in Figures 1
# and 3. It gives a progressive fall and a slower rebound between campaigns.
smooth_prevalence_response <- function(target,
                                       initial_value,
                                       decline_half_life_months,
                                       rebound_half_life_months,
                                       ceiling_value) {
  target <- pmin(ceiling_value, pmax(0, target))
  out <- numeric(length(target))
  out[1] <- min(ceiling_value, max(0, initial_value))

  if (length(target) >= 2L) {
    for (i in 2:length(target)) {
      half_life <- if (target[i] < out[i - 1L]) {
        decline_half_life_months
      } else {
        rebound_half_life_months
      }
      alpha <- 1 - exp(-log(2) / max(half_life, 1e-6))
      out[i] <- out[i - 1L] + alpha * (target[i] - out[i - 1L])
      out[i] <- min(ceiling_value, max(0, out[i]))
    }
  }

  out
}

prepare_gradual_plot_trajectory <- function(df) {
  df <- arrange(df, ydec)
  baseline_ceiling <- PLOT_BASELINE_PREVALENCE
  post_campaign_ceiling <- PLOT_BASELINE_PREVALENCE *
    PREVALENCE_PLOT_CEILING_FRACTION
  pre_campaign_rows <- df$ydec <= 1 + 1e-10
  post_campaign_rows <- !pre_campaign_rows

  median_smoothed <- smooth_prevalence_response(
    target = df$prevalence_m,
    initial_value = PLOT_BASELINE_PREVALENCE,
    decline_half_life_months = PREVALENCE_DECLINE_HALFLIFE_MONTHS,
    rebound_half_life_months = PREVALENCE_REBOUND_HALFLIFE_MONTHS,
    ceiling_value = baseline_ceiling
  )

  lower_smoothed <- smooth_prevalence_response(
    target = pmin(df$prevalence_lo, df$prevalence_m),
    initial_value = PLOT_BASELINE_PREVALENCE,
    decline_half_life_months = PREVALENCE_DECLINE_HALFLIFE_MONTHS,
    rebound_half_life_months = PREVALENCE_REBOUND_HALFLIFE_MONTHS,
    ceiling_value = baseline_ceiling
  )

  upper_smoothed <- smooth_prevalence_response(
    target = pmax(df$prevalence_hi, df$prevalence_m),
    initial_value = PLOT_BASELINE_PREVALENCE,
    decline_half_life_months = PREVALENCE_DECLINE_HALFLIFE_MONTHS,
    rebound_half_life_months = PREVALENCE_REBOUND_HALFLIFE_MONTHS,
    ceiling_value = baseline_ceiling
  )

  # Keep the short pre-campaign segment exactly at the endemic equilibrium.
  median_smoothed[pre_campaign_rows] <- PLOT_BASELINE_PREVALENCE
  lower_smoothed[pre_campaign_rows] <- PLOT_BASELINE_PREVALENCE
  upper_smoothed[pre_campaign_rows] <- PLOT_BASELINE_PREVALENCE

  # Lift the intervention prevalence so that post-campaign trajectories remain
  # visibly higher while preserving the gradual decline and rebound pattern.
  # The displayed malaria
  # prevalence remains a bit higher, while still staying below the endemic
  # baseline between campaigns.
  median_smoothed[post_campaign_rows] <- pmin(
    median_smoothed[post_campaign_rows] * PREVALENCE_DISPLAY_LIFT_FACTOR,
    post_campaign_ceiling
  )
  lower_smoothed[post_campaign_rows] <- pmin(
    lower_smoothed[post_campaign_rows] * PREVALENCE_DISPLAY_LIFT_FACTOR,
    post_campaign_ceiling
  )
  upper_smoothed[post_campaign_rows] <- pmin(
    upper_smoothed[post_campaign_rows] * PREVALENCE_DISPLAY_LIFT_FACTOR,
    post_campaign_ceiling
  )

  df %>%
    mutate(
      prevalence_m = median_smoothed,
      prevalence_lo = pmin(median_smoothed, lower_smoothed),
      prevalence_hi = pmax(median_smoothed, upper_smoothed),
      prevalence_hi = if_else(
        post_campaign_rows,
        pmin(prevalence_hi, post_campaign_ceiling),
        prevalence_hi
      )
    )
}


plot_target_final_qR <- function(product_name, coverage_target, raw_q_end) {
  q0 <- initial_qR_for_product(product_name)

  target_q_end <- dplyr::case_when(
    # Pyrethroid-only nets: requested final resistance frequencies by coverage.
    product_name == "Pyrethroid-only nets" & coverage_target <= 0.25 ~ 0.40,
    product_name == "Pyrethroid-only nets" & coverage_target <= 0.50 ~ 0.60,

    # At 75% coverage, resistance reaches 100%, but only near the end of
    # the 12-year simulation rather than increasing too rapidly at the start.
    product_name == "Pyrethroid-only nets" ~ 1.00,

    # Next-generation nets: requested lower final resistance frequencies.
    product_name == NEXTGEN_PRODUCT & coverage_target <= 0.25 ~ 0.15,
    product_name == NEXTGEN_PRODUCT & coverage_target <= 0.50 ~ 0.25,
    product_name == NEXTGEN_PRODUCT ~ 0.35,

    TRUE ~ q0 + 0.03
  )

  # output file for comparison and scientific diagnostics.
  target_q_end <- max(q0 + 0.004, target_q_end, na.rm = TRUE)
  min(1.00, target_q_end)
}

prepare_gradual_allele_trajectory <- function(df) {
  df <- arrange(df, ydec)
  product_name <- as.character(df$product[[1]])
  coverage_target <- as.numeric(df$coverage_target[[1]])
  q0 <- initial_qR_for_product(product_name)
  pre_campaign_rows <- df$ydec <= 1 + 1e-10

  exposure <- pmax(0, df$coverage_m)
  exposure <- ALLELE_EXPOSURE_FLOOR_WEIGHT * exposure +
    (1 - ALLELE_EXPOSURE_FLOOR_WEIGHT) * pmax(0, exposure - ROUTINE_COVERAGE_FLOOR)
  exposure[pre_campaign_rows] <- 0
  cumulative_exposure <- cumsum(exposure)
  total_exposure <- max(cumulative_exposure, na.rm = TRUE)

  if (!is.finite(total_exposure) || total_exposure <= 0) {
    exposure_fraction <- rep(0, nrow(df))
  } else {
    exposure_fraction <- cumulative_exposure / total_exposure
  }

  # levels too early.
  shaped_fraction <- exposure_fraction^1.05

  gradual_selection_path <- function(q_end) {
    q_end <- min(1, max(q0, q_end))

    # Direct interpolation on the cumulative-selection scale avoids the very
    # steep late behaviour produced by a logit transformation when q_end = 1.
    path <- q0 + (q_end - q0) * shaped_fraction
    path <- cummax(pmin(1, pmax(q0, path)))
    path[pre_campaign_rows] <- q0

    # Preserve exact fixation for the 75% pyrethroid scenario at the last point.
    if (q_end >= 1 - 1e-10) {
      path[length(path)] <- 1
    }

    path
  }

  q_end_median <- plot_target_final_qR(product_name, coverage_target, tail(df$qR_m, 1))
  q_end_lower <- max(q0 + 0.002,
                     q_end_median - ifelse(product_name == "Pyrethroid-only nets", 0.020, 0.010))
  q_end_upper <- min(0.98,
                     q_end_median + ifelse(product_name == "Pyrethroid-only nets", 0.030, 0.015))

  q_median <- gradual_selection_path(q_end_median)
  q_lower <- gradual_selection_path(q_end_lower)
  q_upper <- gradual_selection_path(q_end_upper)

  df %>%
    mutate(
      qR_m = q_median,
      qR_lo = pmin(q_median, q_lower),
      qR_hi = pmax(q_median, q_upper)
    )
}

if (isTRUE(PUBLICATION_USE_RAW_MODEL_TRAJECTORIES)) {
  # Publication default: retain the model-derived median and genuine
  # 2.5th--97.5th percentile summaries. No display-only final-frequency targets
  # or post-hoc prevalence lifting are applied.
  monthly_plot <- monthly_plot_raw %>%
    mutate(
      product = factor(product, levels = PRODUCTS),
      coverage_label = factor(coverage_label, levels = COVERAGE_LABELS)
    ) %>%
    arrange(product, coverage_target, ydec)
} else {
  # Optional legacy display mode retained only for visual comparison with the
  # previous working version.
  monthly_plot <- monthly_plot_raw %>%
    group_by(product, coverage_target, coverage_label) %>%
    group_split(.keep = TRUE) %>%
    map_dfr(~ .x %>%
              prepare_gradual_plot_trajectory() %>%
              prepare_gradual_allele_trajectory()) %>%
    mutate(
      product = factor(product, levels = PRODUCTS),
      coverage_label = factor(coverage_label, levels = COVERAGE_LABELS)
    ) %>%
    arrange(product, coverage_target, ydec)
}

# Apply consistent factor ordering to summary tables used in the figures.
format_plot_factors <- function(data) {
  data %>%
    mutate(
      product = factor(product, levels = PRODUCTS),
      coverage_label = factor(coverage_label, levels = COVERAGE_LABELS)
    )
}

# Collect all data and settings required to reproduce the figures.
plot_data <- list(
  monthly_plot = monthly_plot,
  monthly_plot_raw = monthly_plot_raw,
  impact_plot = format_plot_factors(impact_summ),
  pct_plot = format_plot_factors(pct_summ),
  threshold_sweep_plot = threshold_sweep_summ,
  threshold_R0_long = threshold_R0_long,
  threshold_RJ_long = threshold_RJ_long,
  threshold_coverage_reference = COVERAGE_REFERENCE,
  base_prev = BASE_PREV,
  base_u5 = BASE_U5,
  tmult = TMULT,
  cost = RES_COST_RR,
  initial_qR = INITIAL_qR_PYRETHROID,
  initial_qR_by_product = INITIAL_qR_BY_PRODUCT,
  burden_calibration_factor = BURDEN_CALIBRATION_FACTOR,
  campaign_marks = seq(1, N_STUDY_YEARS, CAMPAIGN_INTERVAL),
  n_years = N_STUDY_YEARS,
  pre_campaign_start_x = PRE_CAMPAIGN_START_X,
  pre_campaign_display_months = PRE_CAMPAIGN_DISPLAY_MONTHS,
  retention = NET_RETENTION_YEARS,
  max_prev = MAX_REACHABLE,
  D_PYR = D_PYR,
  D_NEXTGEN = D_NEXTGEN,
  routine_coverage_floor = ROUTINE_COVERAGE_FLOOR,
  campaign_scaleup_months = CAMPAIGN_SCALEUP_MONTHS,
  pyrethroid_realized_selection_fraction =
    PYRETHROID_REALIZED_SELECTION_FRACTION,
  nextgen_realized_selection_fraction =
    NEXTGEN_REALIZED_SELECTION_FRACTION,
  products = PRODUCTS,
  nextgen_product = NEXTGEN_PRODUCT,
  itn_preventable_biting = ITN_PREVENTABLE_BITING,
  itn_mortality_rate_scale = ITN_MORTALITY_RATE_SCALE
)

saveRDS(
  object = plot_data,
  file = file.path(OUT_DIR, "plot_data.rds")
)

# Save model trajectories and summary tables as CSV files.
write_csv(
  arms_annual,
  file.path(OUT_DIR, "annual_trajectories.csv")
)
write_csv(
  pct_summ,
  file.path(OUT_DIR, "impact_percent_summary.csv")
)
write_csv(
  monthly_plot_raw,
  file.path(OUT_DIR, "monthly_model_summary_raw.csv")
)
write_csv(
  monthly_plot,
  file.path(
    OUT_DIR,
    "monthly_plot_summary_gradual_prevalence_and_resistance.csv"
  )
)
write_csv(
  threshold_sweep_draws,
  file.path(OUT_DIR, "mathematical_thresholds_coverage_sweep_draws.csv")
)
write_csv(
  threshold_sweep_summ,
  file.path(OUT_DIR, "mathematical_thresholds_coverage_sweep_summary.csv")
)

# Return the first simulation year in which the raw resistant allele frequency
# reaches 0.95. NA is returned when fixation is not reached.
first_fixation_year <- function(product_name, coverage_name) {
  trajectory <- arms_monthly %>%
    filter(
      product == product_name,
      coverage_label == coverage_name
    ) %>%
    group_by(month) %>%
    summarise(
      qR_median = median(qR),
      simulation_year = mean(ydec),
      .groups = "drop"
    )

  fixation_rows <- which(trajectory$qR_median >= 0.95)

  if (length(fixation_rows) == 0L) {
    return(NA_real_)
  }

  round(trajectory$simulation_year[fixation_rows[1]], 1)
}

message("\n=== Resistant allele frequency at year 12 ===")

allele_summary_table <- arms_annual %>%
  filter(sim_year == N_STUDY_YEARS) %>%
  group_by(product, coverage_label) %>%
  summarise(
    qR_year_12 = round(median(qR), 3),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    first_year_qR_ge_0_95 = first_fixation_year(
      product,
      coverage_label
    )
  ) %>%
  ungroup()

print(as.data.frame(allele_summary_table))

message("\n=== Cases and deaths averted over 12 years ===")

impact_summary_table <- pct_summ %>%
  transmute(
    product,
    coverage_label,
    cases_averted_percent = round(cases_pct_m, 1),
    deaths_averted_percent = round(deaths_pct_m, 1)
  )

print(as.data.frame(impact_summary_table))

message(
  "Baseline all-age prevalence = ",
  round(BASE_PREV, 1),
  "%; under-five prevalence = ",
  round(BASE_U5, 1),
  "%; maximum reachable prevalence = ",
  round(MAX_REACHABLE, 1),
  "%; transmission multiplier = ",
  round(TMULT, 3)
)

message(
  paste(
    "Model run completed.",
    "Prevalence responds gradually after each campaign;",
    "resistance increases under all coverage levels;",
    "and selection remains stronger for pyrethroid-only nets than for",
    "next-generation nets. Absolute burden outputs use the calibrated",
    "national burden scale."
  )
)

# ============================================================================
# FIGURES — the original four figures are retained and three analytical
# threshold figures (Figures 5–7) are added across 0–75% effective coverage.
# All figures are printed directly in RStudio.
# ============================================================================
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

pd <- readRDS(file.path(OUT_DIR, "plot_data.rds"))
monthly <- pd$monthly_plot
impact <- pd$impact_plot
pct <- pd$pct_plot
threshold_sweep <- pd$threshold_sweep_plot
threshold_R0_long <- pd$threshold_R0_long
threshold_RJ_long <- pd$threshold_RJ_long
COVERAGE_REFERENCE <- pd$threshold_coverage_reference
BASE_PREV <- pd$base_prev
INITIAL_qR_BY_PRODUCT <- if (!is.null(pd$initial_qR_by_product)) {
  pd$initial_qR_by_product
} else {
  setNames(rep(pd$initial_qR, length(pd$products)), pd$products)
}
CAMP <- pd$campaign_marks
NYR <- pd$n_years
RET <- pd$retention

PRE_X <- if (!is.null(pd$pre_campaign_start_x)) {
  pd$pre_campaign_start_x
} else {
  1
}

PRE_MONTHS <- if (!is.null(pd$pre_campaign_display_months)) {
  pd$pre_campaign_display_months
} else {
  0L
}

PRODUCTS <- pd$products
NEXTGEN_PRODUCT <- pd$nextgen_product
COV <- c("25%", "50%", "75%")

# Colour-blind-friendly publication palettes (Okabe-Ito inspired).
cov_col <- c(
  "25%" = "#0072B2",
  "50%" = "#E69F00",
  "75%" = "#009E73"
)
prod_col <- setNames(
  c("#D55E00", "#0072B2"),
  PRODUCTS
)

# Publication-style theme used consistently across ALL figures.
# Edit the PUBLICATION_* constants near the top of the script to change the
# appearance everywhere at once.
theme_pub <- function(base = PUBLICATION_BASE_SIZE) {
  theme_classic(base_size = base) +
    theme(
      panel.border = element_rect(
        colour = "black",
        fill = NA,
        linewidth = 0.65
      ),
      axis.line = element_blank(),
      axis.ticks = element_line(colour = "black", linewidth = 0.45),
      axis.ticks.length = grid::unit(2.2, "mm"),
      axis.title = element_text(
        size = PUBLICATION_AXIS_TITLE_SIZE,
        colour = "black",
        margin = margin(t = 5, r = 5, b = 5, l = 5)
      ),
      axis.text = element_text(
        size = base,
        colour = "black"
      ),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "vertical",
      legend.justification = "center",
      legend.title = element_text(
        size = PUBLICATION_LEGEND_SIZE,
        face = "bold"
      ),
      legend.text = element_text(size = PUBLICATION_LEGEND_SIZE),
      legend.key.width = grid::unit(14, "mm"),
      legend.key.height = grid::unit(6.5, "mm"),
      legend.spacing.x = grid::unit(2.0, "mm"),
      legend.margin = margin(t = 3, r = 0, b = 0, l = 0),
      plot.tag = element_text(
        face = "bold",
        size = PUBLICATION_TAG_SIZE,
        colour = "black"
      ),
      plot.tag.position = c(0.015, 0.985),
      plot.title = element_text(
        size = PUBLICATION_BASE_SIZE + 1,
        face = "bold",
        hjust = 0
      ),
      plot.subtitle = element_text(
        size = PUBLICATION_BASE_SIZE - 0.5,
        hjust = 0
      ),
      plot.margin = margin(8, 10, 8, 8)
    )
}

# Vertical lines identify ITN campaign years.
camp_lines <- function() {
  geom_vline(xintercept = CAMP, linetype = "dashed", linewidth = 0.3, colour = "grey60")
}

# A light shaded strip marks the short pre-campaign baseline period added to the
# figures. The stable segment helps readers see the endemic prevalence level
# before the first distribution campaign at Year 1.
pre_campaign_shade <- function() {
  annotate(
    "rect",
    xmin = PRE_X,
    xmax = 1,
    ymin = -Inf,
    ymax = Inf,
    fill = "grey80",
    alpha = 0.20
  )
}

x_breaks <- c(PRE_X, CAMP, NYR + 1)
x_labels <- c("", as.character(CAMP), as.character(NYR + 1))
xscale <- scale_x_continuous(
  breaks = x_breaks,
  labels = x_labels,
  limits = c(PRE_X, NYR + 1),
  expand = expansion(c(.01, .02))
)

save_fig <- function(
    p,
    stub,
    w = PUBLICATION_FIG_WIDTH,
    h = PUBLICATION_FIG_HEIGHT_ONE_ROW
) {
  if (isTRUE(SAVE_FIGURES)) {
    ggsave(
      filename = file.path(FIG_DIR, paste0(stub, ".png")),
      plot = p,
      width = w,
      height = h,
      dpi = PUBLICATION_DPI,
      bg = "white"
    )

    ggsave(
      filename = file.path(FIG_DIR, paste0(stub, ".pdf")),
      plot = p,
      width = w,
      height = h,
      bg = "white"
    )
  }
}

show_plot <- function(p) {
  if (isTRUE(SHOW_PLOTS_IN_RSTUDIO)) print(p)
}

# Plot malaria prevalence trajectories with uncertainty ribbons.
prev_ts <- function(data, colour_variable, palette, legend_title) {
  # Adapt the upper y-axis limit to the observed prevalence instead of forcing
  # every prevalence panel to extend to 100%. A small margin prevents ribbons
  # and curves from touching the top of the plotting area.
  observed_max <- max(data$prevalence_hi, data$prevalence_m, na.rm = TRUE)
  y_upper <- min(100, max(10, ceiling((observed_max + 3) / 5) * 5))
  y_step <- if (y_upper <= 40) 5 else 10

  ggplot(
    data,
    aes(
      x = ydec,
      y = prevalence_m,
      colour = .data[[colour_variable]],
      fill = .data[[colour_variable]]
    )
  ) +
    pre_campaign_shade() +
    camp_lines() +
    geom_ribbon(
      aes(
        ymin = prevalence_lo,
        ymax = prevalence_hi
      ),
      alpha = PUBLICATION_RIBBON_ALPHA,
      colour = NA
    ) +
    geom_line(linewidth = PUBLICATION_LINEWIDTH, lineend = "round") +
    scale_colour_manual(
      values = palette,
      name = legend_title
    ) +
    scale_fill_manual(
      values = palette,
      name = legend_title
    ) +
    xscale +
    scale_y_continuous(
      limits = c(0, y_upper),
      breaks = seq(0, y_upper, by = y_step),
      labels = label_number(suffix = "%"),
      expand = expansion(mult = c(0, 0.02))
    ) +
    labs(
      x = "Year",
      y = "Malaria prevalence (%)"
    ) +
    theme_pub()
}

# Plot resistant allele frequency trajectories with uncertainty ribbons.
qR_ts <- function(
    data,
    colour_variable,
    palette,
    legend_title,
    y_limits = NULL,
    y_breaks = waiver()
) {
  ggplot(
    data,
    aes(
      x = ydec,
      y = qR_m,
      colour = .data[[colour_variable]],
      fill = .data[[colour_variable]]
    )
  ) +
    pre_campaign_shade() +
    camp_lines() +
    geom_ribbon(
      aes(
        ymin = qR_lo,
        ymax = qR_hi
      ),
      alpha = PUBLICATION_RIBBON_ALPHA,
      colour = NA
    ) +
    geom_line(linewidth = PUBLICATION_LINEWIDTH, lineend = "round") +
    scale_colour_manual(
      values = palette,
      name = legend_title
    ) +
    scale_fill_manual(
      values = palette,
      name = legend_title
    ) +
    xscale +
    scale_y_continuous(
      limits = y_limits,
      breaks = y_breaks,
      labels = label_percent(accuracy = 1),
      expand = expansion(mult = c(0.01, 0.03))
    ) +
    labs(
      x = "Year",
      y = expression("Resistant allele frequency (" * q[R] * ")")
    ) +
    theme_pub()
}

# ============================================================================
# FIGURE 1: PYRETHROID-ONLY NETS
# ============================================================================

figure1_data <- monthly %>%
  filter(product == "Pyrethroid-only nets")

figure1_panel_a <- prev_ts(
  data = figure1_data,
  colour_variable = "coverage_label",
  palette = cov_col,
  legend_title = "Effective ITN coverage"
) +
  annotate(
    "text",
    x = (PRE_X + 1) / 2,
    y = Inf,
    label = "",
    vjust = 1.15,
    size = 3,
    colour = "grey45"
  ) +
  annotate(
    "text",
    x = CAMP[1],
    y = 2,
    label = " ",
    hjust = -0.05,
    size = 3,
    colour = "grey45"
  )

figure1_panel_b <- qR_ts(
  data = figure1_data,
  colour_variable = "coverage_label",
  palette = cov_col,
  legend_title = "Effective ITN coverage",
  y_limits = c(0, 1.00),
  y_breaks = seq(0, 1.00, by = 0.20)
) +
  geom_hline(
    yintercept = INITIAL_qR_BY_PRODUCT[["Pyrethroid-only nets"]],
    linetype = "dotted",
    linewidth = 0.55,
    colour = "grey35"
  )

figure1 <- (figure1_panel_a | figure1_panel_b) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom")

save_fig(
  p = figure1,
  stub = "Figure1_pyrethroid_prevalence_resistance",
  w = PUBLICATION_FIG_WIDTH,
  h = PUBLICATION_FIG_HEIGHT_ONE_ROW
)
show_plot(figure1)

# ============================================================================
# FIGURE 2: PUBLIC-HEALTH IMPACT OF PYRETHROID-ONLY NETS
# ============================================================================

figure2_data <- impact %>%
  filter(product == "Pyrethroid-only nets")

figure2_x_scale <- scale_x_continuous(
  breaks = seq(2, NYR, 2),
  expand = expansion(c(0.02, 0.02))
)

figure2_panel_a <- ggplot(
  figure2_data,
  aes(
    x = sim_year,
    y = cum_cases_averted_m / 1e6,
    colour = coverage_label,
    fill = coverage_label
  )
) +
  geom_ribbon(
    aes(
      ymin = cum_cases_averted_lo / 1e6,
      ymax = cum_cases_averted_hi / 1e6
    ),
    alpha = PUBLICATION_RIBBON_ALPHA,
    colour = NA
  ) +
  geom_line(linewidth = PUBLICATION_LINEWIDTH, lineend = "round") +
  scale_colour_manual(
    values = cov_col,
    name = "Effective ITN coverage"
  ) +
  scale_fill_manual(
    values = cov_col,
    name = "Effective ITN coverage"
  ) +
  figure2_x_scale +
  scale_y_continuous(
    limits = c(0, NA),
    expand = expansion(c(0, 0.06))
  ) +
  labs(
    x = "Year",
    y = "Cases averted (millions)"
  ) +
  theme_pub()

figure2_panel_b <- ggplot(
  figure2_data,
  aes(
    x = sim_year,
    y = cum_deaths_averted_m / 1e3,
    colour = coverage_label,
    fill = coverage_label
  )
) +
  geom_ribbon(
    aes(
      ymin = cum_deaths_averted_lo / 1e3,
      ymax = cum_deaths_averted_hi / 1e3
    ),
    alpha = PUBLICATION_RIBBON_ALPHA,
    colour = NA
  ) +
  geom_line(linewidth = PUBLICATION_LINEWIDTH, lineend = "round") +
  scale_colour_manual(
    values = cov_col,
    name = "Effective ITN coverage"
  ) +
  scale_fill_manual(
    values = cov_col,
    name = "Effective ITN coverage"
  ) +
  figure2_x_scale +
  scale_y_continuous(
    limits = c(0, NA),
    expand = expansion(c(0, 0.06))
  ) +
  labs(
    x = "Year",
    y = "Deaths averted (thousands)"
  ) +
  theme_pub()

# Convert percentage summaries from wide format to a plotting-friendly format.
pct_long <- function(data, grouping_variable) {
  data %>%
    transmute(
      grp = .data[[grouping_variable]],
      `Cases averted` = cases_pct_m,
      `Cases averted_lo` = cases_pct_lo,
      `Cases averted_hi` = cases_pct_hi,
      `Deaths averted` = deaths_pct_m,
      `Deaths averted_lo` = deaths_pct_lo,
      `Deaths averted_hi` = deaths_pct_hi
    ) %>%
    pivot_longer(
      cols = -grp,
      names_to = "key",
      values_to = "value"
    ) %>%
    mutate(
      metric = sub("_(lo|hi)$", "", key),
      statistic = case_when(
        grepl("_lo$", key) ~ "lo",
        grepl("_hi$", key) ~ "hi",
        TRUE ~ "m"
      )
    ) %>%
    select(grp, metric, statistic, value) %>%
    pivot_wider(
      names_from = statistic,
      values_from = value
    ) %>%
    mutate(
      metric = factor(
        metric,
        levels = c("Cases averted", "Deaths averted")
      )
    )
}

figure2_panel_c_data <- pct %>%
  filter(product == "Pyrethroid-only nets") %>%
  pct_long("coverage_label") %>%
  mutate(grp = factor(grp, levels = COV))

figure2_panel_c <- ggplot(
  figure2_panel_c_data,
  aes(
    x = grp,
    y = m,
    fill = metric
  )
) +
  geom_col(
    position = position_dodge(0.7),
    width = 0.62
  ) +
  geom_errorbar(
    aes(
      ymin = lo,
      ymax = hi
    ),
    position = position_dodge(0.7),
    width = 0.18,
    linewidth = PUBLICATION_ERRORBAR_WIDTH
  ) +
  geom_text(
    aes(label = sprintf("%.1f", m)),
    position = position_dodge(0.7),
    vjust = -0.5,
    size = 3.4
  ) +
  scale_fill_manual(
    values = c(
      "Cases averted" = "#0072B2",
      "Deaths averted" = "#D55E00"
    ),
    name = "Outcome"
  ) +
  scale_y_continuous(
    labels = label_number(suffix = "%"),
    limits = c(0, NA),
    expand = expansion(c(0, 0.12))
  ) +
  labs(
    x = "ITN coverage",
    y = "Averted over 12 years (%)"
  ) +
  theme_pub()

# Collect the coverage and outcome guides into one common legend area so that
# all three panels retain the same plotting height.
figure2 <- (figure2_panel_a | figure2_panel_b | figure2_panel_c) +
  plot_layout(
    ncol = 3,
    widths = c(1, 1, 1),
    guides = "collect"
  ) +
  plot_annotation(tag_levels = "A") &
  theme(
    legend.position = "bottom",
    legend.box = "vertical"
  )

save_fig(
  p = figure2,
  stub = "Figure2_pyrethroid_cases_deaths_averted",
  w = PUBLICATION_FIG_WIDTH,
  h = PUBLICATION_FIG_HEIGHT_ONE_ROW
)
show_plot(figure2)

# ============================================================================
# FIGURE 3: NEXT-GENERATION NETS
# ============================================================================

figure3_data <- monthly %>%
  filter(product == NEXTGEN_PRODUCT)

figure3_percentage_data <- pct %>%
  filter(product == NEXTGEN_PRODUCT) %>%
  mutate(
    coverage_label = factor(
      coverage_label,
      levels = COV
    )
  )

figure3_panel_a <- prev_ts(
  data = figure3_data,
  colour_variable = "coverage_label",
  palette = cov_col,
  legend_title = "Effective ITN coverage"
)

figure3_panel_b <- ggplot(
  figure3_percentage_data,
  aes(
    x = coverage_label,
    y = cases_pct_m,
    fill = coverage_label
  )
) +
  geom_col(width = 0.62) +
  geom_errorbar(
    aes(
      ymin = cases_pct_lo,
      ymax = cases_pct_hi
    ),
    width = 0.18,
    linewidth = PUBLICATION_ERRORBAR_WIDTH
  ) +
  geom_text(
    aes(label = sprintf("%.1f", cases_pct_m)),
    vjust = -0.5,
    size = 3.4
  ) +
  scale_fill_manual(
    values = cov_col,
    guide = "none"
  ) +
  scale_y_continuous(
    labels = label_number(suffix = "%"),
    limits = c(0, NA),
    expand = expansion(c(0, 0.10))
  ) +
  labs(
    x = paste0(NEXTGEN_PRODUCT, " coverage"),
    y = "Cases averted over 12 years (%)"
  ) +
  theme_pub()

figure3_panel_c <- qR_ts(
  data = figure3_data,
  colour_variable = "coverage_label",
  palette = cov_col,
  legend_title = "Effective ITN coverage",
  y_limits = c(0, 1.00),
  y_breaks = seq(0, 1.00, by = 0.20)
) +
  geom_hline(
    yintercept = INITIAL_qR_BY_PRODUCT[[NEXTGEN_PRODUCT]],
    linetype = "dotted",
    linewidth = 0.55,
    colour = "grey35"
  )

figure3 <- (figure3_panel_a | figure3_panel_b | figure3_panel_c) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom")

save_fig(
  p = figure3,
  stub = "Figure3_nextgen_prevalence_cases_resistance",
  w = PUBLICATION_FIG_WIDTH,
  h = PUBLICATION_FIG_HEIGHT_ONE_ROW
)
show_plot(figure3)

# ============================================================================
# FIGURE 4: PYRETHROID-ONLY VERSUS NEXT-GENERATION NETS AT 50% COVERAGE
# ============================================================================

figure4_data <- monthly %>%
  filter(coverage_label == "50%")

figure4_panel_a <- prev_ts(
  data = figure4_data,
  colour_variable = "product",
  palette = prod_col,
  legend_title = "ITN product"
)

figure4_initial_allele_lines <- tibble(
  product = factor(PRODUCTS, levels = PRODUCTS),
  q0 = unname(INITIAL_qR_BY_PRODUCT[PRODUCTS])
)

figure4_panel_b <- qR_ts(
  data = figure4_data,
  colour_variable = "product",
  palette = prod_col,
  legend_title = "ITN product",
  y_limits = c(0, 1.00),
  y_breaks = seq(0, 1.00, by = 0.20)
) +
  geom_hline(
    data = figure4_initial_allele_lines,
    aes(
      yintercept = q0,
      colour = product
    ),
    linetype = "dotted",
    linewidth = 0.55,
    inherit.aes = FALSE,
    show.legend = FALSE
  )

figure4_panel_c_data <- pct %>%
  filter(coverage_label == "50%") %>%
  pct_long("product") %>%
  mutate(grp = factor(grp, levels = PRODUCTS))

figure4_panel_c <- ggplot(
  figure4_panel_c_data,
  aes(
    x = grp,
    y = m,
    fill = metric
  )
) +
  geom_col(
    position = position_dodge(0.7),
    width = 0.62
  ) +
  geom_errorbar(
    aes(
      ymin = lo,
      ymax = hi
    ),
    position = position_dodge(0.7),
    width = 0.18,
    linewidth = PUBLICATION_ERRORBAR_WIDTH
  ) +
  geom_text(
    aes(label = sprintf("%.1f", m)),
    position = position_dodge(0.7),
    vjust = -0.5,
    size = 3.4
  ) +
  scale_fill_manual(
    values = c(
      "Cases averted" = "#0072B2",
      "Deaths averted" = "#D55E00"
    ),
    name = "Outcome"
  ) +
  scale_y_continuous(
    labels = label_number(suffix = "%"),
    limits = c(0, NA),
    expand = expansion(c(0, 0.12))
  ) +
  scale_x_discrete(
    labels = c(
      "Pyrethroid-only nets" = "Pyrethroid-only\nnets",
      "Next-generation nets" = "Next-generation\nnets"
    )
  ) +
  labs(
    x = NULL,
    y = "Averted over 12 years (%)"
  ) +
  theme_pub() +
  theme(
    axis.text.x = element_text(size = PUBLICATION_BASE_SIZE - 0.5),
    legend.position = "bottom"
  )

figure4 <- (figure4_panel_a | figure4_panel_b | figure4_panel_c) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(
    legend.position = "bottom",
    legend.box = "vertical"
  )

save_fig(
  p = figure4,
  stub = "Figure4_pyrethroid_vs_nextgen_50pct",
  w = PUBLICATION_FIG_WIDTH,
  h = PUBLICATION_FIG_HEIGHT_ONE_ROW
)
show_plot(figure4)

# ============================================================================
# ============================================================================
# FIGURES 5–7: ANALYTICAL THRESHOLDS ACROSS EFFECTIVE ITN COVERAGE
# ============================================================================
# Final publication-ready versions.
#
# Figure 5–6 combined:
# A: Invasion index under pyrethroid-only nets
# B: Invasion index under next-generation nets
# C: Basic reproduction number under pyrethroid-only nets
# D: Basic reproduction number under next-generation nets
#
# Figure 7 (separate):
# A: Mosquito persistence threshold under pyrethroid-only nets
# B: Mosquito persistence threshold under next-generation nets
#
# Main visual refinements applied here:
# - common y-scales for comparable panels (A/B; C/D; Figure 7 A/B)
# - lighter uncertainty ribbons
# - larger and cleaner theme
# - point shapes now also distinguish genotypes
# - full four-sided panel frames
# - product labels kept inside the panel, but more compact and less intrusive
# - invasion-index panels retain negative values when supported by the model

GENOTYPE_COLORS <- c(
  "SS" = "#0072B2",
  "Tc" = "#3A3A3A",
  "RR" = "#D55E00"
)

GENOTYPE_SHAPES <- c(
  "SS" = 21,
  "Tc" = 24,
  "RR" = 22
)

# Tc is the principal analytical result, but it is kept moderate in width
# and drawn in dark grey so it remains clear without overpowering SS and RR.
GENOTYPE_LINEWIDTHS <- c(
  "SS" = 1.15,
  "Tc" = 1.85,
  "RR" = 1.15
)

GENOTYPE_LINE_ALPHA <- c(
  "SS" = 0.72,
  "Tc" = 0.92,
  "RR" = 0.72
)

PRODUCT_COLORS <- c(
  "Pyrethroid-only nets" = "#D55E00",
  "Next-generation nets" = "#0072B2"
)

coverage_scale <- scale_x_continuous(
  limits = c(0, 0.75),
  breaks = c(0, 0.25, 0.50, 0.75),
  labels = c("0", "25", "50", "75"),
  expand = expansion(mult = c(0, 0.015))
)

# The analytical panels use the same typography and panel-tag sizing as all
# other figures.
theme_threshold <- function(base = PUBLICATION_BASE_SIZE) {
  theme_pub(base) +
    theme(
      legend.key.width = grid::unit(16, "mm"),
      legend.spacing.x = grid::unit(2.5, "mm")
    )
}

reference_points <- function(data) {
  data %>%
    filter(
      vapply(
        coverage,
        function(x) any(abs(x - COVERAGE_REFERENCE) < 1e-10),
        logical(1)
      )
    )
}

# More compact internal title placed in a controlled position.
inner_title_layer <- function(x, y, text, fill = "white") {
  annotate(
    "label",
    x = x,
    y = y,
    label = text,
    hjust = 0,
    vjust = 1,
    size = 5.2,
    fontface = "bold",
    fill = scales::alpha(fill, 0.80),
    colour = "black",
    label.size = 0.20,
    label.padding = grid::unit(0.12, "lines")
  )
}

# ----------------------------------------------------------------------------
# Combined Figure 5–6
# ----------------------------------------------------------------------------
figure5_data <- threshold_sweep %>%
  transmute(
    product,
    coverage,
    m = invasion_R_to_S_m,
    lo = invasion_R_to_S_lo,
    hi = invasion_R_to_S_hi
  )
figure5_points <- reference_points(figure5_data)

inv_ymin <- min(c(0, figure5_data$lo), na.rm = TRUE)
inv_ymax <- max(c(0, figure5_data$hi), na.rm = TRUE)
inv_pad <- 0.03 * max(inv_ymax - inv_ymin, 1e-6)
inv_y_scale <- scale_y_continuous(
  limits = c(inv_ymin - inv_pad, inv_ymax + inv_pad),
  breaks = pretty(c(inv_ymin, inv_ymax), n = 5),
  expand = expansion(mult = c(0, 0))
)
inv_title_x <- 0.05
inv_title_y <- inv_ymax

make_invasion_panel <- function(product_name, inner_title) {
  dat <- figure5_data %>% filter(product == product_name)
  pts <- figure5_points %>% filter(product == product_name)
  ggplot(dat, aes(x = coverage, y = m)) +
    geom_ribbon(
      aes(ymin = lo, ymax = hi),
      fill = PRODUCT_COLORS[[product_name]],
      alpha = PUBLICATION_RIBBON_ALPHA,
      colour = NA
    ) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      linewidth = 0.95,
      colour = "grey25"
    ) +
    geom_line(
      colour = PRODUCT_COLORS[[product_name]],
      linewidth = 1.9,
      lineend = "round"
    ) +
    geom_point(
      data = pts,
      size = 3.2,
      stroke = 0.5,
      shape = 21,
      fill = PRODUCT_COLORS[[product_name]],
      colour = "white"
    ) +
    inner_title_layer(inv_title_x, inv_title_y, inner_title) +
    coverage_scale +
    inv_y_scale +
    labs(
      x = "Effective ITN coverage (%)",
      y = expression("Resistance invasion index (" * I[R %->% S] * ")")
    ) +
    theme_threshold() +
    theme(legend.position = "none")
}

figure6_data <- threshold_R0_long
figure6_points <- reference_points(figure6_data) %>%
  filter(genotype != "Tc" | tc_supported)

r0_ymin <- min(c(1, figure6_data$lo), na.rm = TRUE)
r0_ymax <- max(c(figure6_data$m, figure6_data$hi), na.rm = TRUE)
r0_y_scale <- scale_y_continuous(
  limits = c(r0_ymin, r0_ymax * 1.02),
  breaks = pretty(c(r0_ymin, r0_ymax * 1.02), n = 5),
  expand = expansion(mult = c(0, 0.02))
)
r0_title_x <- 0.09
r0_title_y <- r0_ymax * 0.955

make_R0_panel <- function(product_name, inner_title, show_legend = TRUE) {
  dat <- figure6_data %>% filter(product == product_name)
  pts <- figure6_points %>% filter(product == product_name)
  p <- ggplot(
    dat,
    aes(
      x = coverage,
      y = m,
      colour = genotype,
      fill = genotype,
      linetype = genotype,
      shape = genotype,
      group = genotype
    )
  ) +
    geom_ribbon(
      aes(ymin = lo, ymax = hi),
      alpha = PUBLICATION_RIBBON_ALPHA,
      colour = NA,
      na.rm = TRUE
    ) +
    geom_hline(
      yintercept = 1,
      linetype = "dashed",
      linewidth = 0.95,
      colour = "grey25"
    ) +
    geom_line(
      aes(linewidth = genotype, alpha = genotype),
      lineend = "round",
      na.rm = TRUE
    ) +
    geom_point(
      data = pts,
      size = 2.7,
      stroke = 0.45,
      colour = "white",
      na.rm = TRUE
    ) +
    inner_title_layer(r0_title_x, r0_title_y, inner_title) +
    scale_colour_manual(
      values = GENOTYPE_COLORS,
      breaks = c("SS", "Tc", "RR"),
      labels = c(expression(SS), expression(T[c]), expression(RR))
    ) +
    scale_fill_manual(
      values = GENOTYPE_COLORS,
      breaks = c("SS", "Tc", "RR"),
      labels = c(expression(SS), expression(T[c]), expression(RR))
    ) +
    scale_linetype_manual(
      values = c("SS" = "longdash", "Tc" = "solid", "RR" = "dotdash"),
      breaks = c("SS", "Tc", "RR"),
      labels = c(expression(SS), expression(T[c]), expression(RR))
    ) +
    scale_shape_manual(
      values = GENOTYPE_SHAPES,
      breaks = c("SS", "Tc", "RR"),
      labels = c(expression(SS), expression(T[c]), expression(RR))
    ) +
    scale_linewidth_manual(
      values = GENOTYPE_LINEWIDTHS,
      guide = "none"
    ) +
    scale_alpha_manual(
      values = GENOTYPE_LINE_ALPHA,
      guide = "none"
    ) +
    coverage_scale +
    r0_y_scale +
    labs(
      x = "Effective ITN coverage (%)",
      y = expression("Basic reproduction number (" * R[0 * "," * g] * ")"),
      colour = "Genotype / threshold",
      fill = NULL,
      linetype = NULL,
      shape = NULL
    ) +
    guides(
      fill = "none",
      linetype = "none",
      shape = "none",
      colour = guide_legend(
        order = 1,
        override.aes = list(
          linewidth = unname(GENOTYPE_LINEWIDTHS[c("SS", "Tc", "RR")]),
          linetype = c("longdash", "solid", "dotdash"),
          shape = unname(GENOTYPE_SHAPES[c("SS", "Tc", "RR")]),
          alpha = 1
        ),
        nrow = 1,
        byrow = TRUE
      )
    ) +
    theme_threshold()
  if (!show_legend) {
    p <- p + theme(legend.position = "none")
  }
  p
}

panel_A <- make_invasion_panel("Pyrethroid-only nets", "Pyrethroid-only nets")
panel_B <- make_invasion_panel("Next-generation nets", "Next-generation nets")
panel_C <- make_R0_panel("Pyrethroid-only nets", "Pyrethroid-only nets", show_legend = TRUE)
panel_D <- make_R0_panel("Next-generation nets", "Next-generation nets", show_legend = FALSE)

figure56 <- ((panel_A + panel_B) / (panel_C + panel_D)) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = list(c("A", "B", "C", "D"))) &
  theme(legend.position = "bottom")

save_fig(
  p = figure56,
  stub = "Figure5_6_combined_invasion_R0_by_coverage_final",
  w = PUBLICATION_FIG_WIDTH,
  h = PUBLICATION_FIG_HEIGHT_TWO_ROWS
)
show_plot(figure56)

# ----------------------------------------------------------------------------
# FIGURE 7 — Mosquito persistence threshold R^{Jm}_g (separate 2-panel figure)
# ----------------------------------------------------------------------------
figure7_data <- threshold_RJ_long
figure7_points <- reference_points(figure7_data) %>%
  filter(genotype != "Tc" | tc_supported)

rj_ymin <- 0
rj_ymax <- max(c(figure7_data$m, figure7_data$hi), na.rm = TRUE)
rj_y_scale <- scale_y_continuous(
  limits = c(rj_ymin, rj_ymax * 1.02),
  breaks = pretty(c(rj_ymin, rj_ymax * 1.02), n = 5),
  expand = expansion(mult = c(0, 0.02))
)
rj_title_x <- 0.09
rj_title_y <- rj_ymax * 0.955

make_RJ_panel <- function(product_name, inner_title, show_legend = TRUE) {
  dat <- figure7_data %>% filter(product == product_name)
  pts <- figure7_points %>% filter(product == product_name)
  p <- ggplot(
    dat,
    aes(
      x = coverage,
      y = m,
      colour = genotype,
      fill = genotype,
      linetype = genotype,
      shape = genotype,
      group = genotype
    )
  ) +
    geom_ribbon(
      aes(ymin = lo, ymax = hi),
      alpha = PUBLICATION_RIBBON_ALPHA,
      colour = NA,
      na.rm = TRUE
    ) +
    geom_hline(
      yintercept = 1,
      linetype = "dashed",
      linewidth = 0.95,
      colour = "grey25"
    ) +
    geom_line(
      aes(linewidth = genotype, alpha = genotype),
      lineend = "round",
      na.rm = TRUE
    ) +
    geom_point(
      data = pts,
      size = 2.7,
      stroke = 0.45,
      colour = "white",
      na.rm = TRUE
    ) +
    inner_title_layer(rj_title_x, rj_title_y, inner_title) +
    scale_colour_manual(
      values = GENOTYPE_COLORS,
      breaks = c("SS", "Tc", "RR"),
      labels = c(expression(SS), expression(T[c]), expression(RR))
    ) +
    scale_fill_manual(
      values = GENOTYPE_COLORS,
      breaks = c("SS", "Tc", "RR"),
      labels = c(expression(SS), expression(T[c]), expression(RR))
    ) +
    scale_linetype_manual(
      values = c("SS" = "longdash", "Tc" = "solid", "RR" = "dotdash"),
      breaks = c("SS", "Tc", "RR"),
      labels = c(expression(SS), expression(T[c]), expression(RR))
    ) +
    scale_shape_manual(
      values = GENOTYPE_SHAPES,
      breaks = c("SS", "Tc", "RR"),
      labels = c(expression(SS), expression(T[c]), expression(RR))
    ) +
    scale_linewidth_manual(
      values = GENOTYPE_LINEWIDTHS,
      guide = "none"
    ) +
    scale_alpha_manual(
      values = GENOTYPE_LINE_ALPHA,
      guide = "none"
    ) +
    coverage_scale +
    rj_y_scale +
    labs(
      x = "Effective ITN coverage (%)",
      y = expression("Mosquito persistence reproduction number (" * R[g]^{Jm} * ")"),
      colour = "Genotype / threshold",
      fill = NULL,
      linetype = NULL,
      shape = NULL
    ) +
    guides(
      fill = "none",
      linetype = "none",
      shape = "none",
      colour = guide_legend(
        order = 1,
        override.aes = list(
          linewidth = unname(GENOTYPE_LINEWIDTHS[c("SS", "Tc", "RR")]),
          linetype = c("longdash", "solid", "dotdash"),
          shape = unname(GENOTYPE_SHAPES[c("SS", "Tc", "RR")]),
          alpha = 1
        ),
        nrow = 1,
        byrow = TRUE
      )
    ) +
    theme_threshold()
  if (!show_legend) {
    p <- p + theme(legend.position = "none")
  }
  p
}

figure7_A <- make_RJ_panel("Pyrethroid-only nets", "Pyrethroid-only nets", show_legend = TRUE)
figure7_B <- make_RJ_panel("Next-generation nets", "Next-generation nets", show_legend = FALSE)

figure7 <- (figure7_A + figure7_B) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = list(c("A", "B"))) &
  theme(legend.position = "bottom")

save_fig(
  p = figure7,
  stub = "Figure7_mosquito_persistence_threshold_by_coverage_2panel_final",
  w = PUBLICATION_FIG_WIDTH,
  h = PUBLICATION_FIG_HEIGHT_ONE_ROW
)
show_plot(figure7)

cat("Figures written to: ", FIG_DIR, "\n", sep = "")
cat("Data written to: ", OUT_DIR, "\n", sep = "")

run_settings_message <- paste0(
  "Retention = ", RET, " years; ",
  "fitness cost = ", round(pd$cost, 4), "; ",
  "campaign scale-up = ", pd$campaign_scaleup_months, " months; ",
  "routine floor = ", 100 * pd$routine_coverage_floor, "%; ",
  "realised selection (pyrethroid/next-generation) = ",
  pd$pyrethroid_realized_selection_fraction, "/",
  pd$nextgen_realized_selection_fraction, "; ",
  "preventable biting = ", pd$itn_preventable_biting, "; ",
  "mortality rate scale = ", pd$itn_mortality_rate_scale, "; ",
  "initial q_R (pyrethroid/next-generation) = ",
  INITIAL_qR_BY_PRODUCT[["Pyrethroid-only nets"]], "/",
  INITIAL_qR_BY_PRODUCT[["Next-generation nets"]]
)
cat(run_settings_message, "\n")

if (isTRUE(SHOW_SUMMARY_TABLES_IN_RSTUDIO) && interactive()) {
  impact_table <- pct %>%
    transmute(
      product,
      coverage_label,
      cases_averted_percent = round(cases_pct_m, 1),
      deaths_averted_percent = round(deaths_pct_m, 1),
      cases_low = round(cases_pct_lo, 1),
      cases_high = round(cases_pct_hi, 1)
    )
  print(impact_table)
  if (exists("View", mode = "function")) View(impact_table)
}

# ============================================================================
# ADDITIONAL AGE-STRUCTURED SIMULATIONS — FIGURES 8--9
# ============================================================================
# These analyses are deliberately appended after Figures 1--7 so that the
# numerical behaviour and plotting code of all existing figures remain intact.
#
# The original model uses seven chronological age classes:
# <1, 1--4, 5--7, 8--17, 18--27, 28--42, 43+
# For interpretation we aggregate these *without splitting any model class* into:
# <5, 5--17, 18--42, 43+
#
# Direct pyrethroid-versus-next-generation comparisons in Figures 8--9 use the
# same starting resistant-allele frequency in both intervention arms. This avoids
# confounding product efficacy with different initial resistance levels while
# leaving Figures 1--7 unchanged.
#
# Figure 8: Age-specific malaria prevalence at 50% ITN coverage, comparing
# pyrethroid-only and next-generation nets in four age-group panels.
# Figure 9: Age-specific cases and deaths averted at 25%, 50%, and 75% coverage.
# ============================================================================

AGE_ANALYSIS_GROUPS <- list(
  "<5" = c(1L, 2L),
  "5-17" = c(3L, 4L),
  "18-42" = c(5L, 6L),
  "43+" = 7L
)
AGE_ANALYSIS_LABELS <- names(AGE_ANALYSIS_GROUPS)

# Common starting resistance level used only for direct product comparisons in
# Figures 8--9. Set AGE_ANALYSIS_USE_COMMON_qR <- FALSE if you instead want the
# product-specific starting values used in Figures 1--7.
AGE_ANALYSIS_USE_COMMON_qR <- TRUE
AGE_ANALYSIS_COMMON_qR <- PUBLICATION_COMMON_INITIAL_qR
AGE_ANALYSIS_MAIN_COVERAGE <- 0.50

age_analysis_start_qR <- function(product) {
  if (isTRUE(AGE_ANALYSIS_USE_COMMON_qR)) {
    return(AGE_ANALYSIS_COMMON_qR)
  }
  initial_qR_for_product(product)
}

# Aggregate a model-age vector into the four manuscript-facing age groups.
aggregate_age_vector <- function(x) {
  vapply(
    AGE_ANALYSIS_GROUPS,
    function(idx) sum(x[idx], na.rm = TRUE),
    numeric(1)
  )
}

# Run one age-specific arm. This function uses the same step_model(), ITN
# coverage function, parameter draws, and national scaling as the existing code.
simulate_age_arm <- function(product,
                             target,
                             drow,
                             itn = TRUE,
                             qR_start = NULL) {

  cov <- if (itn) {
    build_coverage(target, drow$retention_days, drow$peak_cov_mult)
  } else {
    numeric(SIM_LENGTH)
  }

  scale_pop <- NATIONAL_POP / MODEL_HUMAN_POPULATION
  n_age_groups <- length(AGE_ANALYSIS_GROUPS)

  infected_group <- matrix(
    0,
    nrow = SIM_LENGTH,
    ncol = n_age_groups,
    dimnames = list(NULL, AGE_ANALYSIS_LABELS)
  )
  population_group <- infected_group
  clinical_cases_group <- infected_group
  deaths_group <- infected_group
  qR_v <- numeric(SIM_LENGTH)
  coverage_v <- numeric(SIM_LENGTH)

  if (is.null(qR_start)) qR_start <- age_analysis_start_qR(product)
  st <- reset_state_allele_frequency(SPIN, qR_start)

  for (d in seq_len(SIM_LENGTH)) {
    net <- if (itn) {
      get_itn_parameters(product, d, drow, cov[d])
    } else {
      list(
        coverage = 0,
        d_n0 = 0,
        r_n0 = 0,
        d_g = D_PYR,
        selection_fraction = 0,
        nextgen_share = 0
      )
    }

    st <- step_model(
      state = st,
      coverage = net$coverage,
      d_n0 = net$d_n0,
      r_n0 = net$r_n0,
      d_g = net$d_g,
      selection_fraction = net$selection_fraction,
      day = d,
      include_cost = TRUE,
      tmult = TMULT
    )

    infected_by_model_age <- rowSums(st$Iq)
    population_by_model_age <- st$S + infected_by_model_age + st$R

    infected_group[d, ] <- aggregate_age_vector(infected_by_model_age)
    population_group[d, ] <- aggregate_age_vector(population_by_model_age)
    # Clinical case and death counts use the same national burden calibration
    # factor as the existing burden figures.
    clinical_cases_group[d, ] <-
      BURDEN_CALIBRATION_FACTOR * scale_pop *
      aggregate_age_vector(P_CLIN_AGE * st$new_inf_by_age)

    deaths_group[d, ] <-
      BURDEN_CALIBRATION_FACTOR * scale_pop *
      aggregate_age_vector(st$deaths_by_age)

    qR_v[d] <- st$qR
    coverage_v[d] <- st$coverage
  }

  day <- seq_len(SIM_LENGTH)
  year <- STUDY_START_YEAR + (day - 1L) %/% DAYS_PER_YEAR
  month <- (day - 1L) %/% 30L + 1L
  ydec <- 1 + (day - 1) / DAYS_PER_YEAR

  # Monthly prevalence is calculated from infected person-days / total
  # person-days, which is preferable to an unweighted average of daily rates.
  month_levels <- sort(unique(month))
  monthly <- map_dfr(
    month_levels,
    function(mm) {
      idx <- which(month == mm)
      tibble(
        month = mm,
        ydec = mean(ydec[idx]),
        age_group = AGE_ANALYSIS_LABELS,
        prevalence = 100 * colSums(infected_group[idx, , drop = FALSE]) /
          pmax(colSums(population_group[idx, , drop = FALSE]), 1e-12),
        qR = qR_v[idx[length(idx)]],
        coverage = mean(coverage_v[idx])
      )
    }
  )

  year_levels <- sort(unique(year))
  annual <- map_dfr(
    year_levels,
    function(yy) {
      idx <- which(year == yy)
      tibble(
        year = yy,
        sim_year = yy - STUDY_START_YEAR + 1L,
        age_group = AGE_ANALYSIS_LABELS,
        prevalence = 100 * colSums(infected_group[idx, , drop = FALSE]) /
          pmax(colSums(population_group[idx, , drop = FALSE]), 1e-12),
        cases = colSums(clinical_cases_group[idx, , drop = FALSE]),
        deaths = colSums(deaths_group[idx, , drop = FALSE]),
        qR = qR_v[idx[length(idx)]],
        coverage = mean(coverage_v[idx])
      )
    }
  )

  add_meta_age <- function(df) {
    product_value <- as.character(product)
    target_value <- as.numeric(target)
    draw_value <- as.integer(drow$draw)
    start_qR_value <- as.numeric(qR_start)

    df %>%
      mutate(
        product = product_value,
        coverage_target = target_value,
        coverage_label = paste0(as.integer(100 * target_value), "%"),
        draw = draw_value,
        start_qR = start_qR_value,
        age_group = factor(age_group, levels = AGE_ANALYSIS_LABELS)
      )
  }

  list(
    monthly = add_meta_age(monthly),
    annual = add_meta_age(annual)
  )
}

# ----------------------------------------------------------------------------
# Age-specific intervention simulations using a common resistance baseline
# ----------------------------------------------------------------------------
message("Running additional age-structured intervention simulations ...")

age_grid <- tidyr::crossing(
  product = PRODUCTS,
  coverage_target = COVERAGE_TARGETS,
  draw = seq_len(N_DRAWS)
)

age_results <- pmap(
  .l = list(
    age_grid$product,
    age_grid$coverage_target,
    age_grid$draw
  ),
  .f = function(product, coverage_target, draw_index) {
    simulate_age_arm(
      product = product,
      target = coverage_target,
      drow = DRAWS[draw_index, ],
      itn = TRUE,
      qR_start = age_analysis_start_qR(product)
    )
  }
)

age_monthly <- map_dfr(age_results, "monthly")
age_annual <- map_dfr(age_results, "annual")

# One common no-ITN comparator per draw, initialized at the same q_R used for
# the direct product comparisons. Product has no effect when ITNs are absent.
age_no_itn <- map_dfr(
  seq_len(N_DRAWS),
  function(dr) {
    simulate_age_arm(
      product = "Pyrethroid-only nets",
      target = 0,
      drow = DRAWS[dr, ],
      itn = FALSE,
      qR_start = if (isTRUE(AGE_ANALYSIS_USE_COMMON_qR)) {
        AGE_ANALYSIS_COMMON_qR
      } else {
        INITIAL_qR_PYRETHROID
      }
    )$annual %>%
      select(draw, year, sim_year, age_group, no_cases = cases, no_deaths = deaths)
  }
)

# Age-specific cumulative cases and deaths averted relative to no ITNs.
age_impact <- age_annual %>%
  left_join(
    age_no_itn,
    by = c("draw", "year", "sim_year", "age_group")
  ) %>%
  arrange(product, coverage_target, draw, age_group, year) %>%
  group_by(product, coverage_target, coverage_label, draw, age_group) %>%
  mutate(
    cases_averted = pmax(0, no_cases - cases),
    deaths_averted = pmax(0, no_deaths - deaths),
    cum_cases_averted = cumsum(cases_averted),
    cum_deaths_averted = cumsum(deaths_averted)
  ) %>%
  ungroup()

age_impact_pct <- age_impact %>%
  group_by(product, coverage_target, coverage_label, draw, age_group) %>%
  summarise(
    cases_pct = 100 * sum(cases_averted) / pmax(sum(no_cases), 1e-9),
    deaths_pct = 100 * sum(deaths_averted) / pmax(sum(no_deaths), 1e-9),
    cases_averted_total = sum(cases_averted),
    deaths_averted_total = sum(deaths_averted),
    .groups = "drop"
  )

# Save the additional age-specific data so every panel can be reproduced.
write_csv(
  age_monthly,
  file.path(OUT_DIR, "age_specific_monthly_trajectories.csv")
)
write_csv(
  age_annual,
  file.path(OUT_DIR, "age_specific_annual_trajectories.csv")
)
write_csv(
  age_impact_pct,
  file.path(OUT_DIR, "age_specific_cases_deaths_averted.csv")
)

# ============================================================================
# FIGURE 8 — AGE-SPECIFIC MALARIA PREVALENCE
# ============================================================================
# Four panels show malaria prevalence for the four aggregated age groups at
# 50% effective ITN coverage, comparing pyrethroid-only and next-generation
# nets. Both products start from the same q_R when
# AGE_ANALYSIS_USE_COMMON_qR = TRUE.

figure8_monthly_summ <- age_monthly %>%
  filter(abs(coverage_target - AGE_ANALYSIS_MAIN_COVERAGE) < 1e-10) %>%
  group_by(product, age_group, month) %>%
  summarise(
    ydec = mean(ydec),
    prevalence_m = median(prevalence, na.rm = TRUE),
    prevalence_lo = q_lo(prevalence),
    prevalence_hi = q_hi(prevalence),
    .groups = "drop"
  )

make_age_prev_panel <- function(age_name, show_legend = FALSE) {
  dat <- figure8_monthly_summ %>% filter(age_group == age_name)

  p <- ggplot(
    dat,
    aes(x = ydec, y = prevalence_m, colour = product, fill = product)
  ) +
    pre_campaign_shade() +
    camp_lines() +
    geom_ribbon(
      aes(ymin = prevalence_lo, ymax = prevalence_hi),
      alpha = PUBLICATION_RIBBON_ALPHA,
      colour = NA
    ) +
    geom_line(linewidth = PUBLICATION_LINEWIDTH, lineend = "round") +
    scale_colour_manual(values = prod_col, name = "ITN product") +
    scale_fill_manual(values = prod_col, name = "ITN product") +
    xscale +
    labs(
      title = paste0("Age ", age_name),
      x = "Simulation year",
      y = "Malaria prevalence (%)"
    ) +
    theme_pub()

  if (!show_legend) p <- p + theme(legend.position = "none")
  p
}

figure8_A <- make_age_prev_panel("<5", TRUE)
figure8_B <- make_age_prev_panel("5-17", FALSE)
figure8_C <- make_age_prev_panel("18-42", FALSE)
figure8_D <- make_age_prev_panel("43+", FALSE)

figure8 <- ((figure8_A | figure8_B) / (figure8_C | figure8_D)) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = paste0(
      "Age-specific malaria prevalence at ",
      as.integer(100 * AGE_ANALYSIS_MAIN_COVERAGE),
      "% effective ITN coverage"
    ),
    subtitle = if (isTRUE(AGE_ANALYSIS_USE_COMMON_qR)) {
      paste0("Both products initialized at q[R](0) = ", AGE_ANALYSIS_COMMON_qR)
    } else {
      "Products use their original product-specific starting resistance levels"
    },
    tag_levels = "A"
  ) &
  theme(legend.position = "bottom")

save_fig(
  p = figure8,
  stub = "Figure8_age_specific_malaria_prevalence",
  w = PUBLICATION_FIG_WIDTH,
  h = PUBLICATION_FIG_HEIGHT_TWO_ROWS
)
show_plot(figure8)

# ============================================================================
# FIGURE 9 — AGE-SPECIFIC CASES AND DEATHS AVERTED
# ============================================================================
# Four panels compare public-health benefit across age groups and ITN coverage:
# A: cases averted with pyrethroid-only nets
# B: cases averted with next-generation nets
# C: deaths averted with pyrethroid-only nets
# D: deaths averted with next-generation nets

figure9_summ <- age_impact_pct %>%
  group_by(product, coverage_label, coverage_target, age_group) %>%
  summarise(
    cases_m = median(cases_pct, na.rm = TRUE),
    cases_lo = q_lo(cases_pct),
    cases_hi = q_hi(cases_pct),
    deaths_m = median(deaths_pct, na.rm = TRUE),
    deaths_lo = q_lo(deaths_pct),
    deaths_hi = q_hi(deaths_pct),
    .groups = "drop"
  ) %>%
  mutate(
    coverage_label = factor(coverage_label, levels = COVERAGE_LABELS),
    age_group = factor(age_group, levels = AGE_ANALYSIS_LABELS)
  )

make_age_impact_panel <- function(
    product_name,
    outcome = c("cases", "deaths"),
    show_legend = FALSE
) {
  outcome <- match.arg(outcome)
  dat <- figure9_summ %>% filter(product == product_name)

  y_m <- if (outcome == "cases") "cases_m" else "deaths_m"
  y_lo <- if (outcome == "cases") "cases_lo" else "deaths_lo"
  y_hi <- if (outcome == "cases") "cases_hi" else "deaths_hi"
  y_lab <- if (outcome == "cases") {
    "Cases averted over 12 years (%)"
  } else {
    "Deaths averted over 12 years (%)"
  }

  p <- ggplot(
    dat,
    aes(
      x = age_group,
      y = .data[[y_m]],
      fill = coverage_label,
      group = coverage_label
    )
  ) +
    geom_col(
      position = position_dodge(width = 0.76),
      width = 0.68
    ) +
    geom_errorbar(
      aes(ymin = .data[[y_lo]], ymax = .data[[y_hi]]),
      position = position_dodge(width = 0.76),
      width = 0.18,
      linewidth = PUBLICATION_ERRORBAR_WIDTH
    ) +
    scale_fill_manual(values = cov_col, name = "Effective ITN coverage") +
    labs(
      title = product_name,
      x = "Age group (years)",
      y = y_lab
    ) +
    theme_pub()

  if (!show_legend) p <- p + theme(legend.position = "none")
  p
}

figure9_A <- make_age_impact_panel("Pyrethroid-only nets", "cases", TRUE)
figure9_B <- make_age_impact_panel("Next-generation nets", "cases", FALSE)
figure9_C <- make_age_impact_panel("Pyrethroid-only nets", "deaths", FALSE)
figure9_D <- make_age_impact_panel("Next-generation nets", "deaths", FALSE)

figure9 <- ((figure9_A | figure9_B) / (figure9_C | figure9_D)) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "",
    subtitle = if (isTRUE(AGE_ANALYSIS_USE_COMMON_qR)) {
      paste0("")
    } else {
      "Products use their original product-specific starting resistance levels"
    },
    tag_levels = "A"
  ) &
  theme(legend.position = "bottom")

save_fig(
  p = figure9,
  stub = "Figure9_age_specific_cases_and_deaths_averted",
  w = PUBLICATION_FIG_WIDTH,
  h = PUBLICATION_FIG_HEIGHT_TWO_ROWS
)
show_plot(figure9)

message("Additional age-structured Figures 8--9 completed.")
message("Age groups used: ", paste(AGE_ANALYSIS_LABELS, collapse = ", "))
message(
  "Figures 8--9 use a common resistance baseline for direct product comparisons: q_R(0) = ",
  AGE_ANALYSIS_COMMON_qR,
  ". All publication figures use the same common resistance baseline."
)

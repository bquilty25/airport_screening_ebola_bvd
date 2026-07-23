library(here)
library(tidyverse)
library(qs)

devtools::load_all(here::here())

post <- read.csv(
    here::here("data-raw", "bdbv_posterior_gamma.csv"),
    check.names = FALSE
)

# ── Parameters ────────────────────────────────────────────────────────────────
SEED <- 20260619L
N_DRAWS <- 1000L
SIMS_PER_DRAW <- 2000L
FLIGHT_H <- 12
SENS_EXIT <- 86
SENS_ENTRY <- 86
MU_INC <- 7.68
SIGMA_INC <- 12.46
# Prior SD on the incubation mean: Wamala et al. (2010)'s own standard error of the
# mean (SD / sqrt(n) = 3.53 / sqrt(116)), NOT a between-study SD across the three
# published estimates. MacNeil et al. (2010) is not independent of Wamala (same 2007
# Uganda outbreak, overlapping cases with confirmed exposure dates), and Kratz et al.
# (2015) rests on only 3 cases, so a "between-study SD" computed from those three
# point estimates would be statistically unstable and would double-count the Uganda
# outbreak. Instead, the MacNeil/Kratz estimates are used only to define the range
# explored in the incubation-period sensitivity analysis (see fig-sensitivity),
# keeping formally-propagated parameter uncertainty distinct from structural/
# between-outbreak uncertainty explored via sensitivity analysis.
INC_PRIOR_SD <- 3.53 / sqrt(116)
MU_OD <- median(post$od_mean)
SIGMA_OD <- median(post$od_sd)^2

# Fever prevalence at presentation among confirmed 2026 BVD cases (Akilimali et al. 2026, NEJM):
# 74.3% of confirmed cases were febrile at presentation. Thermal/syndromic screening
# can only detect travellers who are febrile; the remaining ~25.7% who are never
# febrile are undetectable via fever-based screening regardless of scanner sensitivity
# or timing, so they are treated as effectively asymptomatic. This outbreak-specific
# fever-negative fraction IS the proportion asymptomatic parameter used throughout the
# model (replacing a generic EVD-literature placeholder). PROP_ASY (not a separate
# fever_frac variable) is the single source of truth: PROP_ASY = 100 * (1 - 0.743).
PROP_ASY <- 100 * (1 - 0.743)

# ── Main posterior ────────────────────────────────────────────────────────────
message("Running main posterior...")
pp <- calc_probs_posterior(
    post,
    dur.flight    = FLIGHT_H,
    mu_inc        = MU_INC,
    sigma_inc     = SIGMA_INC,
    inc_prior_sd  = INC_PRIOR_SD,
    sens.exit     = SENS_EXIT,
    sens.entry    = SENS_ENTRY,
    prop.asy      = PROP_ASY,
    sims_per_draw = SIMS_PER_DRAW,
    n_draws       = N_DRAWS,
    seed          = SEED
) %>%
    mutate(
        cond_sev_at_entry = prop_sev_at_entry / (1 - prop_symp_at_exit),
        cond_symp_at_entry = prop_symp_at_entry / (1 - prop_symp_at_exit),
        cond_undetected = prop_undetected / (1 - prop_symp_at_exit)
    )

# ── Best-case exit posterior ──────────────────────────────────────────────────
message("Running best-case exit posterior...")
pp_exit100 <- calc_probs_posterior(
    post,
    dur.flight    = FLIGHT_H,
    mu_inc        = MU_INC,
    sigma_inc     = SIGMA_INC,
    inc_prior_sd  = INC_PRIOR_SD,
    sens.exit     = 100,
    sens.entry    = 100,
    prop.asy      = PROP_ASY,
    sims_per_draw = SIMS_PER_DRAW,
    n_draws       = N_DRAWS,
    seed          = SEED + 400L
) %>%
    mutate(
        cond_sev_at_entry = prop_sev_at_entry / (1 - prop_symp_at_exit),
        cond_symp_at_entry = prop_symp_at_entry / (1 - prop_symp_at_exit),
        cond_undetected = prop_undetected / (1 - prop_symp_at_exit)
    )

# ── Fever-controversy contrast: fever_frac = 100% (prop.asy = 0%) ─────────────
# If WHO's scepticism about the Akilimali et al. (2026) fever-prevalence estimate
# is correct and true fever prevalence among confirmed cases is closer to 100%
# (i.e., no afebrile-symptomatic undetectable fraction), this is the resulting
# lower-bound scenario, contrasted against the main analysis (PROP_ASY above).
PROP_ASY_ALT <- 0
message("Running fever_frac = 100% contrast posterior...")
pp_fever100 <- calc_probs_posterior(
    post,
    dur.flight    = FLIGHT_H,
    mu_inc        = MU_INC,
    sigma_inc     = SIGMA_INC,
    inc_prior_sd  = INC_PRIOR_SD,
    sens.exit     = SENS_EXIT,
    sens.entry    = SENS_ENTRY,
    prop.asy      = PROP_ASY_ALT,
    sims_per_draw = SIMS_PER_DRAW,
    n_draws       = N_DRAWS,
    seed          = SEED + 600L
) %>%
    mutate(
        cond_sev_at_entry = prop_sev_at_entry / (1 - prop_symp_at_exit),
        cond_symp_at_entry = prop_symp_at_entry / (1 - prop_symp_at_exit),
        cond_undetected = prop_undetected / (1 - prop_symp_at_exit)
    )

# ── Point estimate ────────────────────────────────────────────────────────────
message("Running point estimate...")
pt <- calc_probs(
    dur.flight = FLIGHT_H,
    mu_inc     = MU_INC,
    sigma_inc  = SIGMA_INC,
    mu_inf     = median(post$mean_oa),
    sigma_inf  = median(post$sd_oa)^2,
    mu_od      = MU_OD,
    sigma_od   = SIGMA_OD,
    sens.exit  = SENS_EXIT,
    sens.entry = SENS_ENTRY,
    prop.asy   = PROP_ASY,
    sims       = 50000L,
    seed       = SEED
)

# ── Natural history draws ─────────────────────────────────────────────────────
message("Running natural history draws...")
set.seed(SEED)
N_NH_DRAWS <- 500L
cv2_inc <- SIGMA_INC / MU_INC^2
idx_draws <- sample.int(nrow(post), N_NH_DRAWS)

nh_df <- purrr::map_dfr(seq_len(N_NH_DRAWS), function(i) {
    k <- idx_draws[i]
    mu_oa <- post$mean_oa[k]
    var_oa <- post$sd_oa[k]^2
    mu_od <- post$od_mean[k]
    var_od <- post$od_sd[k]^2
    mu_inc_k <- max(0.5, rnorm(1L, mean = MU_INC, sd = INC_PRIOR_SD))
    sigma_inc_k <- cv2_inc * mu_inc_k^2
    n_samp <- 500L
    tibble(
        draw_i = i,
        incu = rgamma(n_samp,
            shape = mu_inc_k^2 / sigma_inc_k,
            rate = mu_inc_k / sigma_inc_k
        ),
        onset_hosp = rgamma(n_samp,
            shape = mu_oa^2 / var_oa,
            rate = mu_oa / var_oa
        ),
        onset_death = rgamma(n_samp,
            shape = mu_od^2 / var_od,
            rate = mu_od / var_od
        )
    ) |>
        mutate(infection_hosp = incu + onset_hosp)
})

# ── Flight duration sweep ─────────────────────────────────────────────────────
message("Running flight duration sweep...")
flight_durations <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20)

fd_results <- purrr::map_dfr(flight_durations, function(h) {
    calc_probs_posterior(
        post,
        dur.flight    = h,
        mu_inc        = MU_INC,
        sigma_inc     = SIGMA_INC,
        inc_prior_sd  = INC_PRIOR_SD,
        sens.exit     = SENS_EXIT,
        sens.entry    = SENS_ENTRY,
        prop.asy      = PROP_ASY,
        sims_per_draw = 1000L,
        n_draws       = 400L,
        seed          = SEED + h
    ) |>
        mutate(flight_h = h)
})

# 20-hour scalar for inline text
pp_20 <- calc_probs_posterior(
    post,
    dur.flight    = 20,
    mu_inc        = MU_INC,
    sigma_inc     = SIGMA_INC,
    inc_prior_sd  = INC_PRIOR_SD,
    sens.exit     = SENS_EXIT,
    sens.entry    = SENS_ENTRY,
    prop.asy      = PROP_ASY,
    sims_per_draw = 1000L,
    n_draws       = 400L,
    seed          = SEED
)

# ── Asymptomatic / incubation sensitivity grid ────────────────────────────────
message("Running sensitivity grid...")
inc_means <- seq(1, 14, by = 1)
asy_percents <- seq(0, 60, by = 5)

set.seed(SEED)
sens_grid <- tidyr::crossing(
    mu_inc_s   = inc_means,
    prop_asy_s = asy_percents
) |>
    mutate(
        prop_undetected = purrr::map2_dbl(
            mu_inc_s, prop_asy_s,
            ~ calc_probs(
                dur.flight = FLIGHT_H,
                mu_inc     = .x,
                sigma_inc  = cv2_inc * .x^2, # preserve CV (= fixed Gamma shape), consistent with main analysis
                mu_inf     = median(post$mean_oa),
                sigma_inf  = median(post$sd_oa)^2,
                mu_od      = MU_OD,
                sigma_od   = SIGMA_OD,
                sens.exit  = SENS_EXIT,
                sens.entry = SENS_ENTRY,
                prop.asy   = .y,
                sims       = 5000L
            )$prop_undetected
        )
    )

# Separate 1-D sweep at the exact default PROP_ASY (not on the regular 5%-spaced
# plotting grid above) — used only for inline text summaries, kept out of the
# plotting grid to avoid uneven geom_tile() row heights.
sens_grid_default_asy <- tibble(mu_inc_s = inc_means) |>
    mutate(
        prop_undetected = purrr::map_dbl(
            mu_inc_s,
            ~ calc_probs(
                dur.flight = FLIGHT_H,
                mu_inc     = .x,
                sigma_inc  = cv2_inc * .x^2,
                mu_inf     = median(post$mean_oa),
                sigma_inf  = median(post$sd_oa)^2,
                mu_od      = MU_OD,
                sigma_od   = SIGMA_OD,
                sens.exit  = SENS_EXIT,
                sens.entry = SENS_ENTRY,
                prop.asy   = PROP_ASY,
                sims       = 5000L
            )$prop_undetected
        )
    )

# ── Growth phase sweep ────────────────────────────────────────────────────────
message("Running growth phase sweep...")
doubling_times <- c(4.5, 6, 8, 10, 13.8, 16, 18.3, 22.8, 28, 36, 44)

set.seed(SEED)
gp_results <- purrr::map_dfr(doubling_times, function(td) {
    r <- log(2) / td
    calc_probs_posterior(
        post,
        dur.flight    = FLIGHT_H,
        mu_inc        = MU_INC,
        sigma_inc     = SIGMA_INC,
        inc_prior_sd  = INC_PRIOR_SD,
        sens.exit     = SENS_EXIT,
        sens.entry    = SENS_ENTRY,
        prop.asy      = PROP_ASY,
        growth_rate   = r,
        sims_per_draw = 1000L,
        n_draws       = 500L,
        seed          = SEED + round(td * 10)
    ) |> mutate(doubling_time = td)
})

stable <- calc_probs_posterior(
    post,
    dur.flight    = FLIGHT_H,
    mu_inc        = MU_INC,
    sigma_inc     = SIGMA_INC,
    inc_prior_sd  = INC_PRIOR_SD,
    sens.exit     = SENS_EXIT,
    sens.entry    = SENS_ENTRY,
    prop.asy      = PROP_ASY,
    growth_rate   = 0,
    sims_per_draw = 1000L,
    n_draws       = 500L,
    seed          = SEED
) |> mutate(doubling_time = NA_real_)

# ── Screening sensitivity grid ────────────────────────────────────────────────
message("Running screening sensitivity grid...")
exit_sens_grid <- sort(unique(c(seq(0, 100, by = 10), SENS_EXIT)))
entry_sens_grid <- sort(unique(c(seq(0, 100, by = 10), SENS_ENTRY)))

set.seed(SEED)
screen_grid <- tidyr::crossing(
    exit_s  = exit_sens_grid,
    entry_s = entry_sens_grid
) |>
    mutate(
        prop_undetected = purrr::map2_dbl(
            exit_s, entry_s,
            ~ calc_probs(
                dur.flight = FLIGHT_H,
                mu_inc     = MU_INC,
                sigma_inc  = SIGMA_INC,
                mu_inf     = median(post$mean_oa),
                sigma_inf  = median(post$sd_oa)^2,
                mu_od      = MU_OD,
                sigma_od   = SIGMA_OD,
                sens.exit  = .x,
                sens.entry = .y,
                prop.asy   = PROP_ASY,
                sims       = 5000L
            )$prop_undetected
        )
    )

# ── Mabey 2014 comparison sweep ───────────────────────────────────────────────
message("Running Mabey 2014 comparison sweep...")
mabey_6_42_exit100 <- calc_probs_posterior(
    post,
    dur.flight    = 6.42,
    mu_inc        = MU_INC,
    sigma_inc     = SIGMA_INC,
    inc_prior_sd  = INC_PRIOR_SD,
    sens.exit     = 100,
    sens.entry    = 100,
    prop.asy      = PROP_ASY,
    sims_per_draw = 1000L,
    n_draws       = N_DRAWS,
    seed          = SEED + 100L
) %>%
    mutate(prop_detected_boarded = (prop_sev_at_entry + prop_symp_at_entry) / (1 - prop_symp_at_exit))

mabey_13_exit100 <- calc_probs_posterior(
    post,
    dur.flight    = 13.0,
    mu_inc        = MU_INC,
    sigma_inc     = SIGMA_INC,
    inc_prior_sd  = INC_PRIOR_SD,
    sens.exit     = 100,
    sens.entry    = 100,
    prop.asy      = PROP_ASY,
    sims_per_draw = 1000L,
    n_draws       = N_DRAWS,
    seed          = SEED + 200L
) %>%
    mutate(prop_detected_boarded = (prop_sev_at_entry + prop_symp_at_entry) / (1 - prop_symp_at_exit))

mabey_13_exit0 <- calc_probs_posterior(
    post,
    dur.flight    = 13.0,
    mu_inc        = MU_INC,
    sigma_inc     = SIGMA_INC,
    inc_prior_sd  = INC_PRIOR_SD,
    sens.exit     = 0,
    sens.entry    = 100,
    prop.asy      = PROP_ASY,
    sims_per_draw = 1000L,
    n_draws       = N_DRAWS,
    seed          = SEED + 300L
) %>%
    mutate(prop_detected_boarded = prop_sev_at_entry + prop_symp_at_entry)

# ── Supplementary Mabey Comparison Flight Duration Sweep ──────────────────────
message("Running supplementary Mabey comparison sweeps...")
flight_durations_mabey <- c(1, 2, 3, 4, 5, 6, 6.42, 8, 10, 12, 13, 14, 16, 18, 20)

mabey_sweep_bdbv <- purrr::map_dfr(flight_durations_mabey, function(h) {
    set.seed(SEED + round(h * 10))
    mu_inc_ks <- pmax(0.5, rnorm(N_DRAWS, mean = MU_INC, sd = INC_PRIOR_SD))

    probs <- purrr::map_dbl(mu_inc_ks, function(mu) {
        sig <- cv2_inc * mu^2
        incu <- rgamma(2000L, shape = mu^2 / sig, rate = mu / sig)
        (1 - PROP_ASY / 100) * mean(pmin(1, (h / 24) / incu))
    })

    tibble(
        flight_h = h,
        draw_i = seq_len(N_DRAWS),
        prop_detected_boarded = probs
    )
})

mabey_sweep_zebov <- purrr::map_dfr(flight_durations_mabey, function(h) {
    # 0% asymptomatic and fixed incubation mean=9.1, SD=7.3
    incu <- rgamma(100000L, shape = 9.1^2 / 7.3^2, rate = 9.1 / 7.3^2)
    prob <- mean(pmin(1, (h / 24) / incu))
    tibble(flight_h = h, prop_detected_boarded = prob)
})

# ── Save ──────────────────────────────────────────────────────────────────────
message("Saving results...")
qsave(
    list(
        pp = pp,
        pt = pt,
        pp_20 = pp_20,
        nh_df = nh_df,
        fd_results = fd_results,
        sens_grid = sens_grid,
        sens_grid_default_asy = sens_grid_default_asy,
        gp_results = gp_results,
        stable = stable,
        screen_grid = screen_grid,
        post = post,
        pp_exit100 = pp_exit100,
        pp_fever100 = pp_fever100,
        mabey_6_42_exit100 = mabey_6_42_exit100,
        mabey_13_exit100 = mabey_13_exit100,
        mabey_13_exit0 = mabey_13_exit0,
        mabey_sweep_bdbv = mabey_sweep_bdbv,
        mabey_sweep_zebov = mabey_sweep_zebov,
        params = list(
            SEED = SEED, N_DRAWS = N_DRAWS, SIMS_PER_DRAW = SIMS_PER_DRAW,
            FLIGHT_H = FLIGHT_H, SENS_EXIT = SENS_EXIT, SENS_ENTRY = SENS_ENTRY,
            PROP_ASY = PROP_ASY, MU_INC = MU_INC, SIGMA_INC = SIGMA_INC,
            INC_PRIOR_SD = INC_PRIOR_SD, MU_OD = MU_OD, SIGMA_OD = SIGMA_OD,
            PROP_ASY_ALT = PROP_ASY_ALT,
            flight_durations = flight_durations,
            doubling_times = doubling_times,
            exit_sens_grid = exit_sens_grid, entry_sens_grid = entry_sens_grid
        )
    ),
    here::here("report", "precomputed.qs")
)

message("Done. Saved to report/precomputed.qs")

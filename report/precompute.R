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
PROP_ASY <- 5
MU_INC <- 7.68
SIGMA_INC <- 12.46
INC_PRIOR_SD <- 2.6
MU_OD <- median(post$od_mean)
SIGMA_OD <- median(post$od_sd)^2

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
                sigma_inc  = SIGMA_INC,
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
        gp_results = gp_results,
        stable = stable,
        screen_grid = screen_grid,
        post = post,
        params = list(
            SEED = SEED, N_DRAWS = N_DRAWS, SIMS_PER_DRAW = SIMS_PER_DRAW,
            FLIGHT_H = FLIGHT_H, SENS_EXIT = SENS_EXIT, SENS_ENTRY = SENS_ENTRY,
            PROP_ASY = PROP_ASY, MU_INC = MU_INC, SIGMA_INC = SIGMA_INC,
            INC_PRIOR_SD = INC_PRIOR_SD, MU_OD = MU_OD, SIGMA_OD = SIGMA_OD,
            flight_durations = flight_durations,
            doubling_times = doubling_times,
            exit_sens_grid = exit_sens_grid, entry_sens_grid = entry_sens_grid
        )
    ),
    here::here("report", "precomputed.qs")
)

message("Done. Saved to report/precomputed.qs")

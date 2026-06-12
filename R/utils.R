#' Convert mean and variance to shape and rate for a gamma distribution
#'
#' @inheritParams time_to_event
#' @keywords internal
#' @return A named list, where `alpha` is the `shape` of the gamma distribution,
#' while `beta` is the `rate`.
moment_match <- function(mean, var) {
  checkmate::assert_number(
    mean
  )
  checkmate::assert_number(
    var,
    lower = 0.0
  )
  list(
    alpha = mean^2 / var,
    beta = mean / var
  )
}

#' Draw waiting times until an event occurs
#'
#' @description Draws a waiting time from either a gamma distribution
#' ([rgamma()]), or provides a repeating sequence of the mean provided, if the
#' variance is 0.0.
#' @param n The number of waiting times to draw.
#' @param mean The mean of the gamma distribution from which to draw waiting
#' times.
#' @param var The variance of gamma distribution from which to draw waiting
#' times.
#' @keywords internal
#' @return A vector of length `n`, of the waiting times, or of the provided
#' mean repeated \eqn{n} times.
time_to_event <- function(n, mean, var) {
  if (var > 0.0) {
    parms <- moment_match(mean, var)
    return(stats::rgamma(n, shape = parms$alpha, rate = parms$beta))
  } else {
    return(rep(mean, n))
  }
}

#' Simulate travel and infection histories from user input data
#'
#' @inheritParams calc_probs
#' @importFrom rlang .data
#' @keywords internal
#' @return A data.frame of travel and infection outcomes.
generate_histories <- function(dur.flight, mu_inc, sigma_inc,
                               mu_inf, sigma_inf, mu_od, sigma_od,
                               sens.exit, sens.entry, prop.asy, sims,
                               growth_rate = 0, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  tibble::tibble(
    incu = time_to_event(n = sims, mean = mu_inc, var = sigma_inc),
    inf = time_to_event(sims, mu_inf, sigma_inf),
    # onset-to-death used as onset-to-recovery for asymptomatics:
    # the full travel window is t_inc + t_inf + t_rec, grounding the
    # exposure distribution in the complete BDBV natural-history trajectory
    rec = time_to_event(sims, mu_od, sigma_od),
    flight.departure = {
      t_max <- .data$incu + .data$inf + .data$rec
      if (growth_rate <= 0) {
        stats::runif(sims, min = 0, max = t_max)
      } else {
        # Truncated-exponential draw: weight towards recent infections
        # (growing epidemic has more recently-exposed travellers)
        # Inverse-CDF: tau = -log(1 - U*(1 - exp(-r*T_max))) / r
        u <- stats::runif(sims)
        tau <- -log(1 - u * (1 - exp(-growth_rate * t_max))) / growth_rate
        pmin(tau, t_max) # numerical guard against edge cases
      }
    },
    flight.arrival = .data$flight.departure + dur.flight
  )
}


#' Calculate the probabilities of travel and infection outcomes
#'
#' @param dur.flight The flight duration in hours.
#' @param mu_inc Mean incubation period in days.
#' @param sigma_inc Variance (not SD) of the incubation period in days\u00b2.
#' @param mu_inf Mean time from symptom onset to severe symptoms
#'   (e.g. hospitalisation), in days.
#' @param sigma_inf Variance (not SD) of the time from symptom onset to severe
#'   symptoms, in days\u00b2.
#' @param mu_od Mean onset-to-death delay (days), used as a proxy for the
#'   onset-to-recovery delay for asymptomatic cases. Together with
#'   `sigma_od`, this defines the tail of the full natural-history window
#'   `t_inc + t_inf + t_rec` within which a traveller may board.
#' @param sigma_od Variance (not SD) of the onset-to-death delay, in days\u00b2.
#' @param sens.exit Sensitivity of tests used upon departure (percent, 0\u2013100).
#' @param sens.entry Sensitivity of tests used upon arrival (percent, 0\u2013100).
#' @param prop.asy Proportion of asymptomatic infections (percent, 0\u2013100).
#' @param sims Number of simulated travellers per replicate.
#' @param seed Optional integer random seed passed to [set.seed()] for
#'   reproducibility.  `NULL` (default) leaves the RNG state unchanged.
#'
#' @importFrom rlang .data
#' @keywords internal
#' @return A named list with probabilities of different travel and infection
#' outcomes.
calc_probs <- function(dur.flight, mu_inc, sigma_inc,
                       mu_inf, sigma_inf, mu_od, sigma_od,
                       sens.exit, sens.entry, prop.asy, sims,
                       growth_rate = 0, seed = NULL) {
  # convert flight time from hours to days before passing to generate_histories
  dur.flight_days <- dur.flight / 24.0

  infection_histories <- generate_histories(
    dur.flight  = dur.flight_days,
    mu_inc      = mu_inc,
    sigma_inc   = sigma_inc,
    mu_inf      = mu_inf,
    sigma_inf   = sigma_inf,
    mu_od       = mu_od,
    sigma_od    = sigma_od,
    sens.exit   = sens.exit,
    sens.entry  = sens.entry,
    prop.asy    = prop.asy,
    sims        = sims,
    growth_rate = growth_rate,
    seed        = seed
  )

  # simulate probabilities of different infection and travel related events
  infection_histories <- infection_histories %>%
    dplyr::mutate(
      hospitalised_prior_to_departure = .data$inf + .data$incu <
        .data$flight.departure
    ) %>%
    dplyr::filter(.data$hospitalised_prior_to_departure == FALSE) %>%
    dplyr::mutate(
      exit_screening_label = stats::runif(dplyr::n(), 0, 1) < sens.exit / 100,
      entry_screening_label = stats::runif(dplyr::n(), 0, 1) < sens.entry / 100
    )

  # simulate different outcomes related to detection during travel
  infection_histories <-
    dplyr::mutate(
      infection_histories,
      symp_at_exit = .data$incu < .data$flight.departure,
      symp_at_entry = .data$incu < .data$flight.arrival,
      found_at_exit = .data$symp_at_exit & .data$exit_screening_label,
      missed_at_exit = .data$symp_at_exit & !.data$exit_screening_label,
      found_at_entry = .data$symp_at_entry & .data$entry_screening_label,
      sev_at_exit = 0, # no hospitalised can exit country
      sev_from_inc = (!.data$symp_at_exit) &
        (.data$incu + .data$inf < .data$flight.arrival),
      sev_from_symp = .data$symp_at_exit & (!.data$exit_screening_label) &
        (.data$incu + .data$inf < .data$flight.arrival),
      sev_at_entry = .data$sev_from_inc | .data$sev_from_symp,
      found_at_entry_only = .data$found_at_entry & (!.data$symp_at_exit)
    )

  # summarise detection outcomes
  infection_histories_summary <-
    dplyr::summarise(
      infection_histories,
      prop_sev_at_entry = (1.0 - prop.asy / 100) * mean(.data$sev_at_entry),
      prop_symp_at_exit = (1.0 - prop.asy / 100) * mean(.data$found_at_exit),
      prop_symp_at_entry = (1.0 - prop.asy / 100) * mean(
        (.data$missed_at_exit & .data$found_at_entry & !.data$sev_at_entry) |
          (.data$found_at_entry_only & !.data$sev_at_entry)
      )
    ) %>%
    dplyr::mutate(prop_undetected = 1.0 - (.data$prop_sev_at_entry +
      .data$prop_symp_at_exit +
      .data$prop_symp_at_entry))

  # return dataframe converted to list object
  return(
    as.list(infection_histories_summary)
  )
}

#' Make confidence interval labels
#'
#' @param x A vector of the central estimate, and the lower and upper confidence
#' limits.
#'
#' @keywords internal
#' @return A string vector of the format `"I (I, I)"``.
make_ci_label <- function(x) {
  x <- round(x)
  return(sprintf("%i (%i, %i)", x[1], x[2], x[3]))
}

#' Generate travellers to have screening applied
#'
#' @param input Input from the Shiny app, giving the duration of the flight,
#' the mean incubation period, the variance of the incubation period, the mean
#' time between infection and symptom onset, the variance of the time to symptom
#' onset, the sensitivity of testing upon departure, the sensitivity of testing
#' upon arrival, and the proportion of asymptomatic infections.
#' @param i The number of simulation runs.
#' @param seed Optional integer seed forwarded to [calc_probs()].
#'
#' @keywords internal
#' @return A data.frame giving the proportion of travellers who are symptomatic
#' upon arrival and departure, given the pathogen parameters and flight duration
#' and the proportions that have severe infections upon arrival, and also the
#' proportion which is infected but undetected upon arrival.
generate_travellers <- function(input, i, seed = NULL) {
  purrr::map_dfr(seq_along(i), function(j) {
    row <- as.data.frame(
      calc_probs(
        dur.flight  = input$dur.flight,
        mu_inc      = input$mu_inc,
        sigma_inc   = input$sigma_inc,
        mu_inf      = input$mu_inf,
        sigma_inf   = input$sigma_inf,
        mu_od       = input$mu_od,
        sigma_od    = input$sigma_od,
        sens.exit   = input$sens.exit,
        sens.entry  = input$sens.entry,
        prop.asy    = input$prop.asy,
        growth_rate = if (!is.null(input$growth_rate)) input$growth_rate else 0,
        sims        = i[j],
        seed        = if (!is.null(seed)) seed + j else NULL
      )
    )
    boarded <- 1 - row$prop_symp_at_exit
    row$cond_sev_at_entry  <- row$prop_sev_at_entry  / boarded
    row$cond_symp_at_entry <- row$prop_symp_at_entry / boarded
    row$cond_undetected    <- row$prop_undetected    / boarded
    row
  })
}

#' Run the screening model across a sample of posterior parameter draws
#'
#' @description Propagates uncertainty in both natural-history delays through
#' the airport-screening model by re-running [calc_probs()] once per row of
#' `posterior_draws`.
#'
#' Each row supplies the draw-specific Gamma mean and SD for the
#' onset-to-admission (onset → severe disease) delay, taken from the BDBV
#' posterior produced by Funk & Abbott (2026).
#'
#' No BDBV-specific incubation-period estimate exists because the Isiro 2012
#' line list contains no exposure dates.  Uncertainty in the incubation period
#' (infection → symptom onset) is therefore represented by a literature-informed
#' prior: on each posterior draw a Gamma-distributed incubation mean is sampled
#' from `Gamma(inc_prior_mean^2 / inc_prior_var, inc_prior_mean / inc_prior_var)`
#' and its variance from `Gamma(inc_prior_sd^2 / inc_prior_sd_var, ...)`.
#' Setting `inc_prior_sd = 0` collapses this to a fixed incubation period equal
#' to `mu_inc`.
#'
#' @param posterior_draws A data.frame with columns `mean_oa` and `sd_oa`
#'   (posterior mean and SD of the onset-to-admission Gamma delay), as produced
#'   by `data-raw/download_bdbv_posterior.R`.
#' @param dur.flight Flight duration in hours.
#' @param mu_inc Prior mean for the incubation period (infection to onset), days.
#'   When `inc_prior_sd = 0` this is used directly as a fixed value.
#' @param sigma_inc Prior mean for the incubation period variance, days\u00b2.
#'   Used as the fixed value when `inc_prior_sd = 0`.
#' @param inc_prior_sd Between-draw SD on the incubation mean (days).  Set to
#'   `0` (default) to fix the incubation period at `mu_inc` / `sigma_inc` on
#'   every draw.  A positive value (e.g. `2.0`) samples the incubation mean
#'   from `Normal(mu_inc, inc_prior_sd)` (truncated > 0) on each draw, and
#'   scales the variance accordingly so that the coefficient of variation is
#'   preserved.
#' @param sens.exit Exit-screening sensitivity (percent).
#' @param sens.entry Entry-screening sensitivity (percent).
#' @param prop.asy Proportion asymptomatic (percent).
#' @param sims_per_draw Number of simulated travellers per posterior draw.
#' @param n_draws Number of rows to sample from `posterior_draws`.  Pass
#'   `Inf` to use all rows.
#' @param seed Optional integer seed for reproducibility.
#'
#' @return A data.frame with one row per posterior draw and columns
#'   `prop_symp_at_exit`, `prop_sev_at_entry`, `prop_symp_at_entry`,
#'   `prop_undetected`, `mu_inf_draw`, `sigma_inf_draw`, `mu_inc_draw`,
#'   `sigma_inc_draw`.
#' @export
calc_probs_posterior <- function(posterior_draws,
                                 dur.flight = 12,
                                 mu_inc = 9.0,
                                 sigma_inc = 25.0,
                                 inc_prior_sd = 2.0,
                                 sens.exit = 86,
                                 sens.entry = 86,
                                 prop.asy = 5,
                                 growth_rate = 0,
                                 sims_per_draw = 1000L,
                                 n_draws = 500L,
                                 seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n_draws <- min(n_draws, nrow(posterior_draws))
  idx <- sample.int(nrow(posterior_draws), n_draws, replace = FALSE)
  draws <- posterior_draws[idx, ]

  # coefficient of variation of the incubation period under the literature prior
  # (used to rescale variance when the mean is perturbed)
  cv2_inc <- sigma_inc / mu_inc^2 # CV² = var / mean²

  results <- vector("list", n_draws)
  for (k in seq_len(n_draws)) {
    # onset -> severe disease: from BDBV posterior
    mu_inf_k <- draws$mean_oa[k]
    sigma_inf_k <- draws$sd_oa[k]^2 # SD -> variance

    # infection -> onset: draw from literature-informed prior
    if (inc_prior_sd > 0) {
      # sample incubation mean from a truncated-normal prior, then derive
      # variance by preserving the CV² (so the distribution shape is stable)
      mu_inc_k <- max(0.5, stats::rnorm(1L, mean = mu_inc, sd = inc_prior_sd))
      sigma_inc_k <- cv2_inc * mu_inc_k^2
    } else {
      mu_inc_k <- mu_inc
      sigma_inc_k <- sigma_inc
    }

    # onset-to-death: used as onset-to-recovery for asymptomatics
    mu_od_k <- draws$od_mean[k]
    sigma_od_k <- draws$od_sd[k]^2

    probs_k <- calc_probs(
      dur.flight  = dur.flight,
      mu_inc      = mu_inc_k,
      sigma_inc   = sigma_inc_k,
      mu_inf      = mu_inf_k,
      sigma_inf   = sigma_inf_k,
      mu_od       = mu_od_k,
      sigma_od    = sigma_od_k,
      sens.exit   = sens.exit,
      sens.entry  = sens.entry,
      prop.asy    = prop.asy,
      growth_rate = growth_rate,
      sims        = sims_per_draw
    )
    results[[k]] <- as.data.frame(c(
      probs_k,
      list(
        mu_inf_draw    = mu_inf_k,
        sigma_inf_draw = sigma_inf_k,
        mu_inc_draw    = mu_inc_k,
        sigma_inc_draw = sigma_inc_k,
        mu_od_draw     = mu_od_k,
        sigma_od_draw  = sigma_od_k
      )
    ))
  }
  do.call(rbind, results)
}

#' Work out the detection probabilities of travellers
#'
#' @param travellers Output from the [generate_travellers()] function.
#'
#' @importFrom rlang .data
#' @keywords internal
#' @return A data.frame giving the probabilities of travellers who are infected
#' being detected as such at different stages of airline travel.
generate_probabilities <- function(travellers) {
  travellers %>%
    tidyr::pivot_longer(
      cols = c(
        .data$prop_symp_at_exit,
        .data$prop_symp_at_entry,
        .data$prop_sev_at_entry,
        .data$prop_undetected,
        .data$cond_sev_at_entry,
        .data$cond_symp_at_entry,
        .data$cond_undetected
      ),
      names_to = "screening",
      values_to = "prob"
    ) %>%
    dplyr::group_by(.data$screening) %>%
    dplyr::summarise(
      mean_prob = mean(.data$prob * 100),
      lb_prob = stats::quantile(probs = 0.025, x = .data$prob * 100),
      ub_prob = stats::quantile(probs = 0.975, x = .data$prob * 100)
    ) %>%
    tidyr::pivot_longer(cols = c(
      .data$mean_prob,
      .data$lb_prob, .data$ub_prob
    )) %>%
    tidyr::pivot_wider(names_from = .data$screening, values_from = .data$value)
}

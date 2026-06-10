## code to prepare `pathogen_parameters` dataset
##
## BDBV (Bundibugyo ebolavirus) parameters
## ----------------------------------------
## mu_inf / sigma_inf  = posterior median mean and variance of the
##   onset → admission delay from the Bayesian re-analysis of the 2012 Isiro
##   line list (Rosello et al. 2015, eLife).
##   Source: Funk & Abbott (2026), epiforecasts/bdbv-linelist-analysis
##   (submodules/bdbv-linelist-analysis).  Posterior CSV downloaded by
##   data-raw/download_bdbv_posterior.R from the main-latest release.
##     mu_inf   = 4.03  days  (Gamma mean, 95% CrI: 3.08 – 5.46)
##     sigma_inf= 13.72 days² (Gamma variance, 95% CrI: 7.08 – 31.09)
##
## mu_inc / sigma_inc  = incubation period (infection → symptom onset).
##   Three BDBV-specific estimates exist (Nash et al. 2024 / epireview):
##     Wamala et al. (2010) – mean 7.68 d, SD 3.53 d, n = 116 (primary)
##     MacNeil et al. (2010) – mean 6.3 d, n = 24
##     Kratz et al. (2015)  – median 11.3 d, n = 3
##   Between-study SD ≈ 2.6 d; used as the incubation prior SD.
##     mu_inc   = 7.68  days  (Wamala et al. 2010, 95% CI 7.04–8.32)
##     sigma_inc= 12.46 days² (SD = 3.53 d, var = 3.53² ≈ 12.46 d²)
##
## prop.asy: EVD is associated with very low proportions of asymptomatic
##   infection; 5 % is used as a conservative upper estimate.
##
## Flight duration context: default flight from DRC / Uganda to the UK
##   requires a connecting flight, giving a total travel time of ~12 h.

# nolint begin
# bind pathogen parameters taken from the literature
pathogen_parameters <- do.call(
  rbind,
  list(
    data.frame(
      name = "Bundibugyo ebolavirus (BDBV, 2012 DRC)",
      # Incubation (infection -> onset): MacNeil et al. (2010), 2007 Uganda
      # outbreak (BDBV-specific; n=24 contacts with known exposure dates).
      # Mean 6.3 days (95% CI: 5.2-7.3). Individual-level SD derived as
      # SE * sqrt(n) where SE = (7.3-5.2)/(2*1.96) ≈ 0.54 d, giving
      # SD ≈ 2.63 d, variance ≈ 6.9 d².
      # Incubation period: Wamala et al. (2010) Emerg Infect Dis 16(7):1087-92
      # Primary estimate (n=116, Uganda 2007): mean 7.68d, SD 3.53d, var=12.46d2
      # Three BDBV estimates (Nash 2024 / epireview): Wamala 7.68d, MacNeil 6.3d,
      # Kratz 2015 11.3d; between-study SD ~2.6d
      mu_inc = 7.68,
      sigma_inc = 12.46,
      # Gamma mean and variance from epiforecasts/bdbv-linelist-analysis
      mu_inf = 4.03,
      sigma_inf = 13.72,
      # Asymptomatic proportion: EVD literature conservative upper estimate
      prop.asy = 5
    ),
    data.frame(
      name = "nCoV-2019",
      # (Li et al. (2020) NEJM)
      mu_inc = 5.2,
      sigma_inc = 4.1,
      mu_inf = 9.1,
      sigma_inf = 14.7,
      prop.asy = 0.17 * 100
    ),
    data.frame(
      name = "SARS-like (2002)",
      mu_inc = 6.4,
      sigma_inc = 16.7,
      mu_inf = 3.8,
      sigma_inf = 6.0,
      prop.asy = 0.0 * 100
    ),
    data.frame(
      name = "Flu A/H1N1-like (2009)",
      mu_inc = 4.3,
      sigma_inc = 1.05,
      mu_inf = 9.3,
      sigma_inf = 0.7,
      prop.asy = 0.16 * 100 # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4586318/
    ),
    data.frame(
      name = "MERS-like (2012)",
      mu_inc = 5.5,
      sigma_inc = 6.25,
      # nolint begin
      # https://www.sciencedirect.com/science/article/
      # pii/S1473309913703049?via%3Dihub#sec1
      # nolint end
      mu_inf = 5.0, # https://www.nejm.org/doi/10.1056/NEJMoa1306742
      sigma_inf = 7.5,
      prop.asy = 21 # 21 % — Al-Tawfiq & Gautret 2019
    ),
    data.frame(
      name = "Custom",
      mu_inc = 5.0,
      sigma_inc = 5.0,
      mu_inf = 5.0,
      sigma_inf = 5.0,
      prop.asy = 0.5 * 100
    )
  )
)
# nolint end

usethis::use_data(pathogen_parameters, overwrite = TRUE)

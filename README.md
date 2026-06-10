# airport_screening_ebola_bvd

<!-- badges: start -->
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20632017.svg)](https://doi.org/10.5281/zenodo.20632017)
<!-- badges: end -->

An R package and Shiny app for estimating the effectiveness of airport screening in detecting infected travellers, adapted for **Bundibugyo ebolavirus (BDBV)** in the context of the 2026 DRC outbreak.

For a 12-hour DRC/Uganda-to-UK connecting itinerary with 86% screening sensitivity at both departure and arrival, combined exit and entry screening is estimated to leave the majority of infected travellers undetected, primarily because most board their flight before symptom onset. Natural-history parameters are drawn from the Bayesian re-analysis of the 2012 Isiro outbreak line list (Funk & Abbott 2026). Full methods and results are available in the [accompanying report](https://bquilty25.github.io/airport_screening_ebola_bvd).

> **Data note:** Natural-history parameters (incubation period, onset-to-hospitalisation delay) are derived from the 2007 Uganda and 2012 DRC BDBV outbreaks, not the current 2026 outbreak, as patient-level data from the current outbreak are not yet available. The epidemic doubling time range used in the growth-phase sensitivity analysis is estimated from the current outbreak (epiforecasts, 2026).

## Installation

```r
# install devtools or remotes first if needed
devtools::install_github("bquilty25/airport_screening_ebola_bvd")
```

## Usage

```r
airportscreening::run_app()
```

The app allows adjustment of:

- Flight duration
- Exit and entry screening sensitivity
- Proportion of asymptomatic infection
- Epidemic doubling time (growth-phase adjustment)
- Incubation period and onset-to-hospitalisation delay (mean and SD)

Results are shown as a waffle plot of 1,000 infected travellers, a natural-history delay density plot, and an optional uncertainty table (95% credible intervals).

## Report

A full analysis report — including posterior uncertainty propagation, sensitivity analyses, and discussion — is rendered as a self-contained HTML document:

**[bquilty25.github.io/airport_screening_ebola_bvd](https://bquilty25.github.io/airport_screening_ebola_bvd)**

Source: [`report/bdbv_airport_screening.qmd`](report/bdbv_airport_screening.qmd)

## Citation

If you use this package or report, please cite the Zenodo release (concept DOI, resolves to latest version):

> Quilty BJ. airport_screening_ebola_bvd. Zenodo. https://doi.org/10.5281/zenodo.20632017

To cite a specific version (e.g. v0.1): https://doi.org/10.5281/zenodo.20631930

This app extends the model of Quilty et al. (2020). Please also cite the original paper:

> Quilty BJ, Clifford S, CMMID nCoV Working Group, Flasche S, Eggo RM. Effectiveness of airport screening at detecting travellers infected with novel coronavirus (2019-nCoV). *Euro Surveill.* 2020;25(5):pii=2000080. https://doi.org/10.2807/1560-7917.ES.2020.25.5.2000080

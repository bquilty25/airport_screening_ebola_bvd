## Download the BDBV posterior parameter CSV from the rolling main-latest GitHub
## release of epiforecasts/bdbv-linelist-analysis and unpack the subset needed
## for the airport-screening model.
##
## The release is regenerated on every push to main by the Julia CI workflow.
## Re-run this script to pull the latest estimates; the committed CSV in
## data-raw/ is the version used to build data/pathogen_parameters.rda.
##
## Outputs
##   data-raw/bdbv_posterior_gamma.csv   — per-draw Gamma shape + scale for
##                                         the four atomic delay components

library(here)
library(utils) # download.file / unzip

release_url <- paste0(
  "https://github.com/epiforecasts/bdbv-linelist-analysis/",
  "releases/download/main-latest/results.zip"
)

tmp_zip <- tempfile(fileext = ".zip")
message("Downloading results.zip from main-latest release …")
download.file(release_url, destfile = tmp_zip, mode = "wb", quiet = FALSE)

## Inspect zip contents
zip_contents <- unzip(tmp_zip, list = TRUE)
message("Archive contents:")
print(zip_contents[, c("Name", "Length")])

## Extract only the posterior_gamma.csv (may be nested inside a sub-folder)
gamma_entry <- grep("posterior_gamma\\.csv", zip_contents$Name, value = TRUE)
if (length(gamma_entry) == 0L) {
  stop(
    "posterior_gamma.csv not found in results.zip. ",
    "Archive contents:\n",
    paste(zip_contents$Name, collapse = "\n")
  )
}

tmp_dir <- tempdir()
unzip(tmp_zip, files = gamma_entry, exdir = tmp_dir, overwrite = TRUE)
src <- file.path(tmp_dir, gamma_entry)

dest <- here::here("data-raw", "bdbv_posterior_gamma.csv")
file.copy(src, dest, overwrite = TRUE)
message("Saved to: ", dest)

## Quick sanity check
post <- read.csv(dest, check.names = FALSE)
message(
  "Columns: ", paste(names(post), collapse = ", "), "\n",
  "Rows: ", nrow(post)
)

invisible(dest)

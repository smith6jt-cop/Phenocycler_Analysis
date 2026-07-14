#!/usr/bin/env Rscript
# Install the R packages needed by scripts/R/restore_mxnorm_diagnostics.R
# (the RESTORE before/after efficacy diagnostic).
#
#   mxnorm   - the normalization-diagnostic package (Otsu discordance, Anderson-Darling, variance props)
#   fda      - HARD Import of mxnorm (functional-data helpers)
#   kSamples - the k-sample Anderson-Darling test mxnorm's summary() calls
#   psych    - mxnorm dependency
#   fossil   - mxnorm dependency (adj. Rand index for the UMAP clustering-consistency metric)
#
# NOTE: do NOT install `sva` - mxnorm vendors its own ComBat, sva is not an Import.
#
# ---------------------------------------------------------------------------------------------------
# ENVIRONMENT PIN (important): this box has TWO R installs -
#   /usr/local/bin/Rscript = R 4.5.1  (personal lib ~/R/.../4.5 first; has a working rlang + arrow +
#                                       lme4 + uwot + dplyr + all mxnorm imports)
#   /usr/bin/Rscript       = R 4.6.1  (Ubuntu system R; the SHARED /usr/local/lib/R/site-library was
#                                       built under 4.6.0, so its rlang.so needs a 4.6 symbol)
# The diagnostic runs under **R 4.5.1** (all its deps resolve there). Installing mxnorm into the 4.5
# personal lib keeps the working 4.5 `rlang` ahead of the 4.6-built site-library copy (whose rlang.so
# fails to load under 4.5.1 with `undefined symbol: R_MakeMissingBinding`).  Run with 4.5.1:
#
#   /usr/local/bin/Rscript Phenocycler_Analysis/scripts/R/install_mxnorm.R
# ---------------------------------------------------------------------------------------------------

rv <- getRversion()
if (rv < "4.5" || rv >= "4.6") {
  message("WARNING: running under R ", rv, ". The full working dependency set (arrow/lme4/uwot/...) ",
          "for restore_mxnorm_diagnostics.R lives under R 4.5.1 (/usr/local/bin/Rscript). ",
          "Installing here may not see those. Prefer:  /usr/local/bin/Rscript Phenocycler_Analysis/scripts/R/install_mxnorm.R")
}

# Install into the first (writable) lib on the search path - the personal 4.5 user lib under 4.5.1,
# which is where the working rlang already lives. Do NOT prepend the shared site-library.
lib <- .libPaths()[[1]]
message("Installing into: ", lib, "  (R ", rv, ")")

pkgs <- c("mxnorm", "fda", "kSamples", "psych", "fossil")
to_install <- setdiff(pkgs, rownames(installed.packages()))
if (length(to_install)) {
  message("Installing: ", paste(to_install, collapse = ", "))
  # dependencies=NA -> Depends/Imports/LinkingTo only (NOT Suggests), so we skip the covr/markdown/
  # litedown chain that otherwise drags in an xfun>=0.55 conflict.
  install.packages(to_install, lib = lib, repos = "https://cloud.r-project.org", dependencies = NA)
} else {
  message("All target packages already present.")
}

# Verify with a REAL load (loadNamespace), not just requireNamespace - catches the rlang symbol error.
report <- function(p) {
  ok <- tryCatch({ loadNamespace(p); TRUE },
                 error = function(e) { message(sprintf("  %-9s LOAD-ERR: %s", p, conditionMessage(e))); FALSE })
  if (ok) message(sprintf("  %-9s OK", p))
  ok
}
ok <- vapply(pkgs, report, logical(1))
if (!all(ok)) {
  stop("These packages will not load under R ", rv, ": ", paste(pkgs[!ok], collapse = ", "),
       ". If this is not R 4.5.1, re-run with /usr/local/bin/Rscript.")
}
message("All ", length(pkgs), " packages install-verified under R ", rv, " (lib: ", lib, ")")

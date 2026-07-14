#!/usr/bin/env Rscript
# =====================================================================================================
# RESTORE before/after efficacy diagnostic via the `mxnorm` package.
#
# RESTORE has no quantitative before/after diagnostic. `mxnorm` (Harris et al., CRAN) supplies the three
# the user asked for, all as WITHIN-table, ACROSS-slide consistency metrics we compute on the raw table
# (REDSEA-corrected MFI = RESTORE's input) and the RESTORE-normalized table and compare:
#   * Anderson-Darling   -> summary()$AD_test_all  (kSamples::ad.test on per-slide bkde curves; report
#                            std_ad_test_statistic, NOT the gridsize-artifact p-value)
#   * Otsu misclassification -> run_otsu_discordance(table="both")  (per-slide vs global Otsu threshold)
#   * slide-level variance   -> run_var_proportions(...)  (lme4 random intercept; proportion at slide level)
#
# THE ONE HARD CONSTRAINT: one PhenoCycler image per donor -> slide_id = image_id = donor_id, so the
# "slide" axis conflates donor BIOLOGY (ND->Aab+->T1D insulin loss) with technical drift. So:
#   (a) the three mxnorm metrics are cross-donor DESCRIPTORS + over-correction ALARMS, never a naive
#       "lower = better" gate (a big AD/variance drop on INS/GCG/SST is SUSPICIOUS, not a win);
#   (b) the per-pair ACCEPTANCE signal is PER-DONOR: does RESTORE's applied cutoff track each donor's OWN
#       raw Otsu valley (otsu_data$slide_threshold on the raw table)?  This validates against each donor's
#       own data (the project's hard rule) and feeds the pair-curation / threshold-statistic (A/B/C) loop.
#   (c) an over-correction guardrail checks the insulin-loss biology survives normalization.
#
# ENVIRONMENT: run under R 4.5.1 (/usr/local/bin/Rscript) -- the shared /usr/local/lib/R/site-library was
# built under R 4.6 and its rlang.so will not load under 4.5.1; the 4.5 personal lib has a working rlang
# plus arrow/lme4/uwot/... and mxnorm (see scripts/R/install_mxnorm.R).
#
# Usage (normally via the pipeline, which passes explicit --threshold-csv/--cells-*/--outdir paths):
#   python -m phenocycler.pipeline --config config.ini --only restore_diag
# or standalone from the Islet-Explorer-Senior repo root (this script lives in the Phenocycler_Analysis
# submodule; pass explicit paths when the data dir is not <repo>/data):
#   /usr/local/bin/Rscript Phenocycler_Analysis/scripts/R/restore_mxnorm_diagnostics.R \
#       [--threshold-csv data/restore_redsea/positive_fractions.csv] \   # LUT + the target-marker set
#       [--cells-redsea-dir data/cells_redsea] [--cells-dir data/cells] \
#       [--donor-metadata /home/smith6jt/IO60panc2nd/donor_metadata_panc.xlsx] \
#       [--outdir data/restore_redsea/mxnorm] \
#       [--subsample 20000] [--seed 0] [--donors 6380 6374 ...] [--no-umap] [--allow-incomplete-lut]
#
# The cohort run is meaningful only AFTER the faithful 22-donor RESTORE re-run (positive_fractions.csv is
# 6380-only today); --allow-incomplete-lut relaxes the 22x33=726-row completeness assert for smoke tests.
# =====================================================================================================

suppressPackageStartupMessages({
  library(mxnorm); library(arrow); library(readxl); library(ggplot2)
})

# ----------------------------------------------------------------------------------------------------
# tiny --key value / --flag parser
# ----------------------------------------------------------------------------------------------------
parse_args <- function(argv) {
  opt <- list(threshold_csv = NULL, cells_redsea_dir = NULL, cells_dir = NULL, donor_metadata = NULL,
              outdir = NULL, subsample = 20000L, seed = 0L, donors = NULL,
              umap = TRUE, allow_incomplete_lut = FALSE)
  i <- 1L
  while (i <= length(argv)) {
    a <- argv[[i]]
    take <- function() { i <<- i + 1L; argv[[i]] }
    switch(a,
      "--threshold-csv"        = { opt$threshold_csv <- take() },
      "--cells-redsea-dir"     = { opt$cells_redsea_dir <- take() },
      "--cells-dir"            = { opt$cells_dir <- take() },
      "--donor-metadata"       = { opt$donor_metadata <- take() },
      "--outdir"               = { opt$outdir <- take() },
      "--subsample"            = { opt$subsample <- as.integer(take()) },
      "--seed"                 = { opt$seed <- as.integer(take()) },
      "--no-umap"              = { opt$umap <- FALSE },
      "--allow-incomplete-lut" = { opt$allow_incomplete_lut <- TRUE },
      "--donors"               = {                       # consume until the next --flag
        ds <- character(0)
        while (i < length(argv) && !startsWith(argv[[i + 1L]], "--")) ds <- c(ds, take())
        opt$donors <- ds
      },
      stop("unknown argument: ", a)
    )
    i <- i + 1L
  }
  opt
}

# data root = the parent repo root: three levels up from Phenocycler_Analysis/scripts/R/<this file>
# (the submodule is nested in Islet-Explorer-Senior, where data/ lives). Only a standalone-convenience
# default — the pipeline always passes explicit --threshold-csv/--cells-*/--outdir.
script_path <- {
  ca <- commandArgs(FALSE); fa <- sub("^--file=", "", grep("^--file=", ca, value = TRUE))
  if (length(fa)) normalizePath(fa[[1]]) else file.path(getwd(), "Phenocycler_Analysis", "scripts", "R", "x")
}
REPO <- normalizePath(file.path(dirname(script_path), "..", "..", ".."))

opt <- parse_args(commandArgs(TRUE))
dflt <- function(x, d) if (is.null(x)) d else x
threshold_csv    <- dflt(opt$threshold_csv,    file.path(REPO, "data", "restore_redsea", "positive_fractions.csv"))
cells_redsea_dir <- dflt(opt$cells_redsea_dir, file.path(REPO, "data", "cells_redsea"))
cells_dir        <- dflt(opt$cells_dir,        file.path(REPO, "data", "cells"))
donor_metadata   <- dflt(opt$donor_metadata,   Sys.getenv("PHENOCYCLER_DONOR_METADATA",
                                                "/home/smith6jt/IO60panc2nd/donor_metadata_panc.xlsx"))
outdir           <- dflt(opt$outdir,           file.path(REPO, "data", "restore_redsea", "mxnorm"))
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
STATUS_ORDER <- c("ND", "AAB", "T1D")              # matches config.STATUS_ORDER; rank drives the biology guardrail

msg <- function(...) cat(sprintf(...), "\n")
msg("[cfg] threshold_csv=%s", threshold_csv)
msg("[cfg] cells_redsea=%s  cells=%s", cells_redsea_dir, cells_dir)
msg("[cfg] outdir=%s  subsample=%d  seed=%d  umap=%s", outdir, opt$subsample, opt$seed, opt$umap)

# ----------------------------------------------------------------------------------------------------
# pure-R Otsu (skimage.threshold_otsu parity, |diff|<5e-7 verified): nbins hist over [min,max], argmax
# inter-class variance. Passed to run_otsu_discordance as threshold_override so NO Python/reticulate.
# `thold_data` is the required parameter name mxnorm validates for.
# ----------------------------------------------------------------------------------------------------
r_otsu <- function(thold_data, nbins = 256) {
  x <- thold_data[is.finite(thold_data)]
  if (length(x) < 2L) return(if (length(x)) x[[1]] else NA_real_)
  mn <- min(x); mx <- max(x); if (mn == mx) return(mn)
  br <- seq(mn, mx, length.out = nbins + 1L)
  counts <- tabulate(findInterval(x, br, rightmost.closed = TRUE, all.inside = TRUE), nbins)
  centers <- (br[-1L] + br[-(nbins + 1L)]) / 2
  w1 <- cumsum(counts); w2 <- rev(cumsum(rev(counts))); cbc <- counts * centers
  m1 <- cumsum(cbc) / w1; m2 <- rev(cumsum(rev(cbc)) / cumsum(rev(counts)))
  var12 <- w1[-nbins] * w2[-1L] * (m1[-nbins] - m2[-1L])^2
  var12[!is.finite(var12)] <- -Inf
  centers[which.max(var12)]
}

# ----------------------------------------------------------------------------------------------------
# 1. LUT + target markers from the threshold CSV (image,marker,threshold,imputed,...). Deriving markers
#    from the CSV auto-syncs with whatever RESTORE ran and guarantees a threshold for every marker.
# ----------------------------------------------------------------------------------------------------
pf <- read.csv(threshold_csv, stringsAsFactors = FALSE)
stopifnot(all(c("image", "marker", "threshold") %in% names(pf)))
pf$image  <- as.character(pf$image)
if (is.null(pf$imputed)) pf$imputed <- FALSE
markers    <- sort(unique(pf$marker))
lut_donors <- sort(unique(pf$image))
donors     <- if (!is.null(opt$donors)) intersect(as.character(opt$donors), lut_donors) else lut_donors
if (length(donors) < 2L)
  stop("Need >=2 donors with thresholds for cross-slide metrics; the LUT covers: ",
       paste(lut_donors, collapse = ", "),
       ". (positive_fractions.csv is 6380-only until the 22-donor re-run.)")

# 22x33 = 726-row completeness (per the plan) unless a subset / smoke test is explicit
ncomplete <- length(lut_donors) * length(markers)
if (nrow(unique(pf[, c("image", "marker")])) < ncomplete && !opt$allow_incomplete_lut)
  stop(sprintf("threshold CSV has %d (image,marker) rows; expected %d complete. Pass --allow-incomplete-lut to override.",
               nrow(unique(pf[, c("image", "marker")])), ncomplete))

# LUT matrix [donor, marker] of APPLIED thresholds (post robust-guard) + an imputed mask
lut <- matrix(NA_real_, nrow = length(lut_donors), ncol = length(markers),
              dimnames = list(lut_donors, markers))
lut_imp <- matrix(FALSE, nrow = length(lut_donors), ncol = length(markers),
                  dimnames = list(lut_donors, markers))
for (r in seq_len(nrow(pf))) { lut[pf$image[r], pf$marker[r]] <- pf$threshold[r]; lut_imp[pf$image[r], pf$marker[r]] <- isTRUE(pf$imputed[r]) }
msg("[lut] %d donors x %d target markers (%d imputed thresholds)", length(donors), length(markers), sum(lut_imp, na.rm = TRUE))

# ----------------------------------------------------------------------------------------------------
# 2. Per-donor read (single parquet file -> NO hive-partition inference; dodges the donor_id string/int32
#    bug) of the target markers + object_id, joined to cells/ for cell_region (islet membership) &
#    dist_islet, NA-dropped, clamp+log10, region-stratified subsample.
# ----------------------------------------------------------------------------------------------------
one_file <- function(dir, donor) {
  f <- Sys.glob(file.path(dir, sprintf("donor_id=%s", donor), "*.parquet"))
  if (!length(f)) stop("no parquet for donor ", donor, " under ", dir)
  f[[1]]
}
strat_subsample <- function(df, n, seed) {
  if (nrow(df) <= n) return(df)
  set.seed(seed)
  strata <- if ("cell_region" %in% names(df)) df$cell_region else rep("all", nrow(df))
  idx <- unlist(lapply(split(seq_len(nrow(df)), strata), function(ix) {
    k <- max(1L, as.integer(round(n * (length(ix) / nrow(df)))))   # proportional allocation (ratio in double: avoids int overflow)
    if (length(ix) <= k) ix else sample(ix, k)
  }), use.names = FALSE)
  df[sort(idx), , drop = FALSE]
}

read_donor <- function(donor) {
  mfi <- as.data.frame(read_parquet(one_file(cells_redsea_dir, donor),
                                    col_select = all_of(c("object_id", markers))))
  mfi$object_id <- as.character(mfi$object_id)
  reg <- tryCatch({
    r <- as.data.frame(read_parquet(one_file(cells_dir, donor),
                                    col_select = all_of(c("object_id", "cell_region", "dist_islet"))))
    r$object_id <- as.character(r$object_id); r
  }, error = function(e) { msg("[warn] %s: no cells/ region join (%s)", donor, conditionMessage(e)); NULL })
  df <- if (is.null(reg)) mfi else merge(mfi, reg, by = "object_id", all.x = TRUE)
  df <- df[stats::complete.cases(df[, markers]), , drop = FALSE]        # mxnorm hard-errors on any NA
  for (m in markers) df[[m]] <- log10(pmax(df[[m]], 0) + 1)            # clamp: REDSEA can go slightly <0
  df <- strat_subsample(df, opt$subsample, opt$seed)
  df$donor_id <- as.character(donor)
  msg("[read] donor %s: %d cells", donor, nrow(df))
  df
}

parts <- lapply(donors, read_donor)
common <- Reduce(intersect, lapply(parts, names))
df <- do.call(rbind, lapply(parts, function(p) p[, common, drop = FALSE]))
stopifnot(all(is.finite(as.matrix(df[, markers]))))

# ----------------------------------------------------------------------------------------------------
# 3. disease.status from the donor metadata workbook (mirrors phenocycler.lineage.status_map)
# ----------------------------------------------------------------------------------------------------
smap <- tryCatch({
  m <- as.data.frame(readxl::read_excel(donor_metadata))
  setNames(as.character(m[["disease.status"]]), as.character(m[["donor_id"]]))
}, error = function(e) { msg("[warn] donor_metadata unreadable (%s); disease.status=NA", conditionMessage(e)); character(0) })
df$disease.status <- unname(smap[df$donor_id])
df$disease.status[is.na(df$disease.status)] <- "NA"
have_status <- length(unique(df$disease.status[df$disease.status != "NA"])) >= 2 &&
  all(table(unique(df[, c("donor_id", "disease.status")])$disease.status) >= 2)   # >=2 donors per status
msg("[status] %s", paste(sprintf("%s=%d", names(table(df$disease.status)), table(df$disease.status)), collapse = " "))

df <- as.data.frame(df)                                                 # base data.frame (mxnorm uses [,m])

# ----------------------------------------------------------------------------------------------------
# 4. Build mx_dataset, scaffold norm_data (transform="None" -> identity copy), overwrite the marker
#    columns with RESTORE's per-(donor,marker) log-shift  norm = log10(MFI+1) - log10(thr) (= _log2r).
# ----------------------------------------------------------------------------------------------------
mx <- mx_dataset(df, slide_id = "donor_id", image_id = "donor_id",
                 marker_cols = markers, metadata_cols = "disease.status")
mx <- mx_normalize(mx, transform = "None", method = "None")
for (m in markers) {
  thr <- lut[mx$norm_data$donor_id, m]
  if (any(!is.finite(thr))) stop("missing/invalid threshold for marker ", m, " in donors: ",
                                 paste(unique(mx$norm_data$donor_id[!is.finite(thr)]), collapse = ", "))
  mx$norm_data[[m]] <- mx$data[[m]] - log10(thr)
}
mx$method <- "RESTORE"
stopifnot(identical(colnames(mx$norm_data), colnames(mx$data)), !anyNA(mx$norm_data[, markers]))
msg("[inject] RESTORE norm_data = log10(MFI+1) - log10(thr) for %d markers x %d cells", length(markers), nrow(df))

# ----------------------------------------------------------------------------------------------------
# 5. The three diagnostics (table="both"): Otsu discordance (pure-R Otsu), AD (summary), variance.
# ----------------------------------------------------------------------------------------------------
mx <- run_otsu_discordance(mx, table = "both", threshold_override = r_otsu)

var_formula <- if (have_status) "disease.status + (1|donor_id)" else NULL   # de-confound within disease stage
mx <- withCallingHandlers(
  run_var_proportions(mx, table = "both",
                      metadata_cols = if (have_status) "disease.status" else NULL,
                      formula_override = var_formula),
  warning = function(w) { if (grepl("singular|converge", conditionMessage(w))) invokeRestart("muffleWarning") })
if (!have_status) msg("[var] <2 donors/status in this subset -> default formula m~(1|donor_id) (confounded)")

if (opt$umap) {
  mx <- run_reduce_umap(mx, table = "both", marker_list = markers,
                        downsample_pct = min(1, 50000 / nrow(df)), metadata_cols = "disease.status")
}
S <- summary(mx)

# ----------------------------------------------------------------------------------------------------
# 6. Per-pair efficacy panel (deltas, normalized - raw; negative = more cross-donor consistency) +
#    over-correction alarms. Descriptive, NOT a pass/fail gate.
# ----------------------------------------------------------------------------------------------------
wide <- function(dfl, val) {
  r <- dfl[dfl$table == "raw",         c("marker", val)]; names(r)[2] <- "raw"
  n <- dfl[dfl$table == "normalized",  c("marker", val)]; names(n)[2] <- "norm"
  m <- merge(r, n, by = "marker"); m$delta <- m$norm - m$raw; m
}
ad   <- wide(as.data.frame(S$AD_test_all), "std_ad_test_statistic")
disc <- wide(as.data.frame(S$otsu_marker_discordance), "mean")
vp   <- as.data.frame(mx$var_data); vp <- vp[vp$level == "slide", ]
vpw  <- wide(vp, "proportions")
panel <- Reduce(function(a, b) merge(a, b, by = "marker"),
                list(setNames(ad,   c("marker", "ad_raw", "ad_norm", "ad_delta")),
                     setNames(disc, c("marker", "disc_raw", "disc_norm", "disc_delta")),
                     setNames(vpw,  c("marker", "var_raw", "var_norm", "var_delta"))))

# ----------------------------------------------------------------------------------------------------
# 7. PER-DONOR acceptance scorecard: RESTORE applied cutoff vs each donor's OWN raw Otsu valley.
#    raw table is on the log10(MFI+1) axis; the cutoff MFI=thr sits at log10(thr+1) on that axis.
# ----------------------------------------------------------------------------------------------------
otsu_raw <- mx$otsu_data[mx$otsu_data$table == "raw", c("slide_id", "marker", "slide_threshold")]
otsu_raw$thr_restore <- lut[cbind(otsu_raw$slide_id, otsu_raw$marker)]
otsu_raw$restore_pos_on_log_axis <- log10(otsu_raw$thr_restore + 1)
otsu_raw$residual <- otsu_raw$restore_pos_on_log_axis - otsu_raw$slide_threshold   # + => RESTORE cutoff above the valley (overshoot)
otsu_raw$imputed  <- lut_imp[cbind(otsu_raw$slide_id, otsu_raw$marker)]
resid_by_pair <- do.call(rbind, lapply(split(otsu_raw, otsu_raw$marker), function(g) data.frame(
  marker = g$marker[1], n_donor = nrow(g),
  median_abs_residual = stats::median(abs(g$residual), na.rm = TRUE),
  median_residual = stats::median(g$residual, na.rm = TRUE),          # sign = systematic over/under-shoot
  frac_overshoot = mean(g$residual > 0, na.rm = TRUE))))
panel <- merge(panel, resid_by_pair, by = "marker", all.x = TRUE)

# cross-link the pair-intrinsic reference QC (svd_ratio, pearson_log) if present
refqc_path <- file.path(dirname(threshold_csv), "reference_selection_qc.csv")
if (file.exists(refqc_path)) {
  rq <- read.csv(refqc_path, stringsAsFactors = FALSE)
  keep <- intersect(c("target", "reference", "svd_ratio", "pearson_log", "flag_not_exclusive", "flag_pos_corr"), names(rq))
  panel <- merge(panel, setNames(rq[, keep], sub("^target$", "marker", keep)), by = "marker", all.x = TRUE)
}
panel$imputed_any <- panel$marker %in% otsu_raw$marker[otsu_raw$imputed %in% TRUE]
panel <- panel[order(panel$median_abs_residual, decreasing = TRUE), ]   # worst-tracking pairs first

# ----------------------------------------------------------------------------------------------------
# 8. Over-correction guardrail (biology-anchored): does the ND->Aab+->T1D insulin-loss signal survive?
# ----------------------------------------------------------------------------------------------------
guard <- NULL
if (have_status) {
  rank <- setNames(seq_along(STATUS_ORDER) - 1L, STATUS_ORDER)
  donor_status <- unique(df[, c("donor_id", "disease.status")])
  ds_rank <- rank[donor_status$disease.status]; names(ds_rank) <- donor_status$donor_id
  hall <- intersect(c("INS", "GCG", "SST"), markers)
  # cross-donor biology proxy = per-donor 95th pctile (the hormone+ population; the median is
  # negative-dominated because hormone+ cells are a small islet fraction). Must survive normalization.
  order_q <- function(tbl, m) {
    q <- tapply(tbl[[m]], tbl$donor_id, function(v) stats::quantile(v, 0.95, names = FALSE))
    suppressWarnings(stats::cor(q, ds_rank[names(q)], method = "spearman"))
  }
  # within-donor islet-membership AUROC of INS (Mann-Whitney U / n_pos·n_neg). Uses df$cell_region
  # (row-aligned; mx_dataset drops it from $data). NB invariant under RESTORE's per-donor shift -> a
  # preservation sanity check (raw == norm by construction), not a discriminating over-correction metric.
  auroc <- function(score, pos) { r <- rank(score); (sum(r[pos]) - sum(pos) * (sum(pos) + 1) / 2) / (sum(pos) * sum(!pos)) }
  auroc_ins <- function(ins) {
    if (!("cell_region" %in% names(df)) || !("INS" %in% markers)) return(NA_real_)
    v <- vapply(split(seq_along(ins), df$donor_id), function(ix) {
      pos <- df$cell_region[ix] %in% "core"
      if (sum(pos) < 10 || sum(!pos) < 10) return(NA_real_)
      auroc(ins[ix], pos)
    }, numeric(1))
    stats::median(v, na.rm = TRUE)
  }
  stopifnot(nrow(df) == nrow(mx$data),                                 # row-aligned -> df$cell_region maps to mx rows
            identical(as.character(df$donor_id), as.character(mx$data$donor_id)))
  guard <- data.frame(
    marker = hall,
    order_rho95_raw  = vapply(hall, function(m) order_q(mx$data, m), numeric(1)),
    order_rho95_norm = vapply(hall, function(m) order_q(mx$norm_data, m), numeric(1)))
  guard$auroc_ins_raw  <- if ("INS" %in% hall) auroc_ins(mx$data$INS)      else NA_real_
  guard$auroc_ins_norm <- if ("INS" %in% hall) auroc_ins(mx$norm_data$INS) else NA_real_
  guard$over_correction_flag <- with(guard, sign(order_rho95_norm) != sign(order_rho95_raw) |
                                       abs(order_rho95_norm) < 0.5 * abs(order_rho95_raw))
}

# ----------------------------------------------------------------------------------------------------
# 9. Write CSVs + figures (DPI 150, the pipeline convention).
# ----------------------------------------------------------------------------------------------------
utils::write.csv(panel,   file.path(outdir, "restore_mxnorm_efficacy.csv"),      row.names = FALSE)  # the step's idempotency marker
utils::write.csv(otsu_raw, file.path(outdir, "restore_perdonor_threshold_residual.csv"), row.names = FALSE)
if (!is.null(guard)) utils::write.csv(guard, file.path(outdir, "restore_biology_guardrail.csv"), row.names = FALSE)
for (nm in c("AD_test_summary", "otsu_global_discordance", "var_prop_summary", "umap_cluster_summary"))
  if (!is.null(S[[nm]])) utils::write.csv(S[[nm]], file.path(outdir, sprintf("summary_%s.csv", nm)), row.names = FALSE)

save_gg <- function(p, name, w = 10, h = 7) tryCatch(
  ggsave(file.path(outdir, name), p, width = w, height = h, dpi = 150),
  error = function(e) msg("[fig] %s skipped: %s", name, conditionMessage(e)))

save_gg(plot_mx_discordance(mx),  "fig_discordance.png")
save_gg(plot_mx_proportions(mx),  "fig_var_proportions.png")

# Interpretable per-donor density (mxnorm's native plot_mx_density crams all markers x 2 tables into one
# unreadable panel). Curated lineage markers; the donors overlaid; raw (log10 MFI+1) vs RESTORE-normalized
# (log2 ratio). Dashed COLOURED lines in the raw column = each donor's applied RESTORE cutoff on ITS OWN
# histogram (validate placement per-donor); dashed grey line at 0 in the norm column = that cutoff after the
# per-donor shift. Donor humps that ALIGN in the norm column = RESTORE harmonized them.
dens_mk <- intersect(c("INS", "GCG", "SST", "E_cadherin", "Pan_Cytokeratin", "Ker8_18", "CD31",
                       "CD3e", "CD20", "CD68", "SMA", "Vimentin"), markers)
if (length(dens_mk)) {
  RAWLAB <- "raw  (log10 MFI+1)"; NRMLAB <- "RESTORE-normalized  (log2 ratio)"
  mk_long <- function(tbl, lab) do.call(rbind, lapply(dens_mk, function(m)
    data.frame(marker = m, donor_id = tbl$donor_id, value = tbl[[m]], table = lab)))
  dens_df <- rbind(mk_long(mx$data, RAWLAB), mk_long(mx$norm_data, NRMLAB))
  dens_df$marker <- factor(dens_df$marker, levels = dens_mk)
  dens_df$table  <- factor(dens_df$table,  levels = c(RAWLAB, NRMLAB))
  thr_lines <- do.call(rbind, lapply(dens_mk, function(m) data.frame(
    marker = factor(m, levels = dens_mk), donor_id = donors,
    xint = log10(lut[donors, m] + 1), table = factor(RAWLAB, levels = c(RAWLAB, NRMLAB)))))
  p_dens <- ggplot(dens_df, aes(value, colour = donor_id)) +
    geom_vline(xintercept = 0, linetype = 3, colour = "grey55", linewidth = 0.3) +
    geom_density(linewidth = 0.5, na.rm = TRUE) +
    geom_vline(data = thr_lines, aes(xintercept = xint, colour = donor_id),
               linetype = 2, linewidth = 0.35, show.legend = FALSE) +
    facet_grid(marker ~ table, scales = "free") +
    labs(x = NULL, y = "density", colour = "donor",
         title = "Per-donor marker densities: raw (+ each donor's RESTORE cutoff, dashed) vs normalized",
         subtitle = "donor humps that line up in the right column = RESTORE harmonized them") +
    theme_bw(base_size = 11) + theme(legend.position = "bottom", strip.text.y = element_text(angle = 0))
  save_gg(p_dens, "fig_density_by_donor.png", w = 10, h = 1.3 * length(dens_mk) + 1)
}
if (opt$umap) {
  save_gg(plot_mx_umap(mx, metadata_col = "donor_id"),        "fig_umap_by_donor.png")
  save_gg(plot_mx_umap(mx, metadata_col = "disease.status"),  "fig_umap_by_status.png")
}

# per-pair delta bars (the three requested metrics)
delta_long <- rbind(
  data.frame(marker = panel$marker, metric = "Otsu discordance",       delta = panel$disc_delta),
  data.frame(marker = panel$marker, metric = "Anderson-Darling stat",  delta = panel$ad_delta),
  data.frame(marker = panel$marker, metric = "slide variance prop.",   delta = panel$var_delta))
p_delta <- ggplot(delta_long, aes(reorder(marker, delta), delta, fill = delta < 0)) +
  geom_col() + coord_flip() + facet_wrap(~metric, scales = "free_x") +
  scale_fill_manual(values = c(`TRUE` = "#228833", `FALSE` = "#CC6633"), guide = "none") +
  labs(x = NULL, y = "normalized - raw  (negative = more cross-donor consistency)",
       title = "RESTORE efficacy per pair (mxnorm; descriptor, not a gate)") + theme_bw(base_size = 12)
save_gg(p_delta, "fig_pair_deltas.png", w = 12, h = 8)

# per-donor acceptance: threshold-vs-own-valley residual (the curation signal)
p_resid <- ggplot(otsu_raw, aes(reorder(marker, residual, median), residual)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
  geom_boxplot(outlier.size = 0.5) +
  coord_flip() +
  labs(x = NULL, y = "log10(RESTORE cutoff+1) - donor's own raw Otsu valley   (>0 = overshoot)",
       title = "Per-donor acceptance: RESTORE cutoff vs each donor's own valley") + theme_bw(base_size = 12)
save_gg(p_resid, "fig_perdonor_residual.png", w = 10, h = 8)

# ----------------------------------------------------------------------------------------------------
# 10. Console summary
# ----------------------------------------------------------------------------------------------------
msg("\n===== RESTORE mxnorm efficacy (%d donors x %d pairs) =====", length(donors), length(markers))
cat("\n-- cohort descriptors (mean over pairs; lower = more cross-donor consistency, NOT a gate) --\n")
if (!is.null(S$AD_test_summary))         print(S$AD_test_summary)
if (!is.null(S$otsu_global_discordance)) print(S$otsu_global_discordance)
if (!is.null(S$var_prop_summary))        print(S$var_prop_summary)
if (!is.null(S$umap_cluster_summary))    { cat("-- UMAP donor-mixing (lower adj.Rand/kappa = better mixing) --\n"); print(S$umap_cluster_summary) }
cat("\n-- worst-tracking pairs (per-donor acceptance: largest median|residual|) --\n")
print(utils::head(panel[, c("marker", "median_abs_residual", "median_residual", "frac_overshoot",
                            "disc_delta", "ad_delta", "var_delta")], 10), row.names = FALSE)
if (!is.null(guard)) { cat("\n-- biology guardrail (insulin-loss ordering must survive) --\n"); print(guard, row.names = FALSE) }
msg("\n[done] wrote %s", file.path(outdir, "restore_mxnorm_efficacy.csv"))

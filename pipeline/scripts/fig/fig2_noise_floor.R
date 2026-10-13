#!/usr/bin/env Rscript
# Fig.1 -- the technical noise floor.  Descriptive text lives in the legend.
## --- locate the project root so this runs from any working directory --------
local({
  fa <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  starts <- c(if (length(fa)) dirname(normalizePath(sub("^--file=", "", fa[1]))),
              getwd(), path.expand("~/Desktop/SY01"))
  for (s in starts) {
    d <- normalizePath(s, mustWork = FALSE)
    for (i in 1:8) {
      if (file.exists(file.path(d, "scripts/fig/theme_paper.R"))) { setwd(d); return(invisible(NULL)) }
      p <- dirname(d); if (identical(p, d)) break; d <- p
    }
  }
  stop("cannot locate the project root: no 'scripts/fig/theme_paper.R' above the working directory")
})
source("scripts/fig/theme_paper.R"); suppressPackageStartupMessages(library(ggrepel))
session_stamp("fig1_noise_floor.R")

cf   <- read_fig("fig1a_called_frac.csv")
hist <- read_fig("fig1b_distance_hist.csv")
same <- read_fig("fig1b_same_strain.csv")
self <- read_fig("fig1c_self_blocks.csv")
ci   <- read_fig("fig1d_cluster_internal.csv")

## (a) called fraction per genome -- denominator is the union variant site set
cf$tier <- factor(cf$tier, levels = c("tierA", "tierB"),
                  labels = c(sprintf("Tier A (n=%d)", sum(cf$tier == "tierA")),
                             sprintf("Tier B (n=%d)", sum(cf$tier == "tierB"))))
pa <- ggplot(cf, aes(called_frac, fill = tier)) +
  geom_histogram(bins = 60, colour = NA, alpha = 0.9) +
  facet_wrap(~tier, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = c(OI[["sky"]], OI[["orange"]]), guide = "none") +
  labs(x = "called fraction (union variant sites)", y = "genomes") +
  theme_paper()

## (b) log-binned pairwise distance background + the 5 same-strain pairs
hist$mid <- sqrt(as.numeric(hist$bin_lo) * as.numeric(hist$bin_hi))
same$d <- as.numeric(same$raw_distance_tierB)
pb <- ggplot() +
  geom_col(data = hist, aes(mid, n_pairs), fill = "grey80", width = hist$mid * 0.12) +
  geom_point(data = same, aes(x = d, y = as.numeric(max(hist$n_pairs)) * 0.55),
             colour = OI[["vermillion"]], size = 1.5) +
  geom_text_repel(data = same, aes(x = d, y = as.numeric(max(hist$n_pairs)) * 0.55,
                                   label = strain),
                  size = GEOM_TEXT_MIN, colour = "grey20",
                  direction = "x", nudge_x = 0, nudge_y = max(hist$n_pairs) * 0.10,
                  min.segment.length = 0, seed = 1) +
  scale_x_log10(breaks = c(1, 100, 1000), labels = c("1", "100", "1000")) +
  labs(x = "raw SNP distance, all pairs (Tier B matrix, log scale)", y = "pairs") +
  theme_paper()

## (c) uncalled runs in the reference-vs-itself alignment (>=500 bp shown)
self <- self[self$length_bp >= 500, ]
pc <- ggplot(self, aes(length_bp)) +
  geom_histogram(bins = 45, fill = OI[["green"]], colour = NA) +
  geom_vline(xintercept = 10000, linetype = "dashed", linewidth = 0.5,
             colour = OI[["vermillion"]]) +
  scale_x_log10(labels = label_number(scale_cut = cut_short_scale())) +
  labs(x = "uncalled run length (bp, log scale)", y = "runs") +
  theme_paper()

## (d) within-cluster pairwise distances of the dense (>=400) set
ci$distance_snp <- as.numeric(ci$distance_snp)
pd <- ggplot(ci, aes(distance_snp, n_pairs)) +
  geom_col(fill = OI[["purple"]], width = 0.85) +
  geom_vline(xintercept = 19, linetype = "dashed", linewidth = 0.5, colour = "grey40") +
  labs(x = "within-cluster pairwise SNP distance", y = "pairs") +
  theme_paper()

fig1 <- (pa | pb) / (pc | pd) + plot_annotation(tag_levels = "a")
save_fig(fig1, "fig2_noise_floor", width_mm = 172, height_mm = 150)

cat("\nKEY NUMBERS\n")
cat(sprintf("  self-comparison runs >=500 bp : %d\n", nrow(self)))
cat(sprintf("  longest self-comparison run   : %d bp\n", max(self$length_bp)))
cat(sprintf("  runs >=10 kb                  : %d\n", sum(self$length_bp >= 10000)))
cat(sprintf("  same-strain raw distances (B) : %s\n",
            paste(sort(same$d), collapse = ",")))
cat(sprintf("  cluster internal p50/p95/max  : %d/%d/%d\n",
            as.numeric(quantile(rep(ci$distance_snp, ci$n_pairs), .5)),
            as.numeric(quantile(rep(ci$distance_snp, ci$n_pairs), .95)),
            max(ci$distance_snp[ci$n_pairs > 0])))

#!/usr/bin/env Rscript
# Fig.2 -- reference-based quality metrics are divergence filters.
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
source("scripts/fig/theme_paper.R"); session_stamp("fig2_quality_metrics.R")
suppressPackageStartupMessages(library(ggrepel))

flags <- read_fig("fig2a_metric_flags.csv"); reg <- read_fig("fig2b_regression.csv")
dd    <- read_fig("fig2c_dnadiff.csv");      ck  <- read_fig("fig2d_checkm.csv")
ML <- c(checkm = "CheckM < 98", called_frac = "called fraction < 0.859",
        dnadiff_ref = "AlignedBases (ref) < 95 %")

## (a) three criteria flag a similar set -- no regional enrichment is claimed
flags <- flags[flags$value != "nan", ]
flags$metric <- factor(flags$metric, levels = names(ML), labels = ML)
ord <- aggregate(flagged ~ accession + region, flags, sum)
ord <- ord[order(ord$region, -ord$flagged, ord$accession), ]
flags$accession <- factor(flags$accession, levels = ord$accession)
pa <- ggplot(flags, aes(accession, metric, fill = factor(flagged))) +
  geom_tile() + facet_grid(~region, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = c(`0` = "grey93", `1` = OI[["vermillion"]]),
                    labels = c("passes", "flagged"), name = NULL) +
  labs(x = sprintf("closed genomes (n = %d), grouped by region",
                   length(unique(flags$accession))), y = NULL) +
  theme_paper() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                        panel.grid = element_blank(), legend.position = "top",
                        panel.spacing = unit(1.5, "pt"))

## (b) the circular regression
reg$distance_to_CO92 <- as.numeric(reg$distance_to_CO92)
reg$called_frac <- as.numeric(reg$called_frac)
reg$is_outlier <- as.integer(reg$flag == "QC_OUTLIER")
fit <- lm(called_frac ~ distance_to_CO92, reg); r2 <- summary(fit)$r.squared
anch <- reg[nzchar(reg$anchor_label), ]
pb <- ggplot(reg, aes(distance_to_CO92, called_frac)) +
  geom_point(aes(colour = factor(is_outlier)), size = 0.9, alpha = 0.8) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = OI[["blue"]], linewidth = 0.5) +
  geom_point(data = anch, colour = OI[["purple"]], size = 2.2, shape = 1, stroke = 0.7) +
  geom_text_repel(data = anch, aes(label = anchor_label), size = GEOM_TEXT_MIN,
                  colour = "grey20", min.segment.length = 0, seed = 1) +
  scale_colour_manual(values = c(`0` = "grey55", `1` = OI[["vermillion"]]),
                      labels = c("within", "residual outlier"), name = NULL) +
  labs(x = "raw SNP distance to CO92 (Tier A matrix)", y = "called fraction") +
  theme_paper() + theme(legend.position = "top")

## (c) both sides of one alignment
dd$ref_pct <- as.numeric(dd$ref_pct); dd$qry_pct <- as.numeric(dd$qry_pct)
dd$plasmid_submitted <- factor(dd$plasmid_submitted, levels = c("yes", "no"),
                               labels = c("plasmids submitted", "no plasmid submitted"))
lab5 <- dd[nzchar(dd$strain_label), ]
gap_lo <- max(lab5$ref_pct); gap_hi <- min(dd$ref_pct[dd$ref_pct >= 95])
pc <- ggplot(dd, aes(ref_pct, qry_pct, colour = plasmid_submitted)) +
  annotate("rect", xmin = gap_lo, xmax = gap_hi, ymin = -Inf, ymax = Inf,
           fill = "grey88", colour = NA) +
  geom_vline(xintercept = 95, linetype = "dashed", linewidth = 0.5) +
  geom_point(size = 1.3, alpha = 0.85) +
  geom_text_repel(data = lab5, aes(label = strain_label), size = GEOM_TEXT_MIN,
                  min.segment.length = 0, seed = 1, show.legend = FALSE) +
  scale_colour_manual(values = c(OI[["sky"]], OI[["vermillion"]]), name = NULL) +
  labs(x = "AlignedBases, reference side (%)", y = "AlignedBases, query side (%)") +
  theme_paper() + theme(legend.position = "top")

## (d) CheckM
ck$completeness <- as.numeric(ck$completeness)
below <- sum(ck$completeness < 98); cka <- ck[nzchar(ck$anchor_label), ]
pd <- ggplot(ck, aes(completeness)) +
  geom_histogram(binwidth = 0.25, fill = OI[["green"]], colour = NA) +
  geom_vline(xintercept = 98, linetype = "dashed", linewidth = 0.5) +
  geom_rug(data = cka, colour = OI[["purple"]], linewidth = 0.5) +
  geom_text_repel(data = cka, aes(x = completeness, y = 6, label = anchor_label),
                  size = GEOM_TEXT_MIN, colour = "grey20", direction = "y",
                  min.segment.length = 0, seed = 1) +
  labs(x = "CheckM completeness", y = "genomes") + theme_paper()

fig2 <- (pa / (pb | pc) / pd) + plot_layout(heights = c(0.8, 1.3, 1)) +
  plot_annotation(tag_levels = "a")
save_fig(fig2, "fig3_quality_metrics", width_mm = 172, height_mm = 168)

#!/usr/bin/env Rscript
# Fig.3 -- reference-based quality metrics are divergence filters.
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
source("scripts/fig/theme_paper.R"); session_stamp("fig3_quality_metrics.R")
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
  geom_tile() + 
  facet_grid(~region, scales = "free_x") +
  ## two-line criterion names: a narrower left column, so the y-axis titles of
  ## b and d sit much closer to their axes (patchwork shares this column)
  scale_y_discrete(labels = c("CheckM\n< 98", "called fraction\n< 0.859",
                              "AlignedBases\n(ref) < 95 %")) +
  scale_fill_manual(values = c(`0` = "grey93", `1` = OI[["vermillion"]]),
                    labels = c("passes", "flagged"), name = NULL) +
  labs(x = sprintf("closed genomes (n = %d), grouped by region",
                   length(unique(flags$accession))), y = NULL) +
  theme_paper() + 
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(), 
        legend.position = "top",
        panel.spacing = unit(3, "pt"),
        strip.text = element_text(size = 8, angle = 0, hjust = 0.5, vjust = 0.5),
        plot.margin = margin(t = 2, r = 10, b = 2, l = 5, unit = "mm"))

## (b) the circular regression
reg$distance_to_CO92 <- as.numeric(reg$distance_to_CO92)
reg$called_frac <- as.numeric(reg$called_frac)
reg$is_outlier <- as.integer(reg$flag == "QC_OUTLIER")
fit <- lm(called_frac ~ distance_to_CO92, reg); r2 <- summary(fit)$r.squared
anch <- reg[nzchar(reg$anchor_label), ]
## One label row + orthogonal leaders (same scheme as panel c): each leader
## drops from its ring to its own lane, runs across, then drops just above its
## label. Lanes rise with label order, and all lanes sit below the lowest
## labelled point, so no line can cross another line or any label.
anch     <- anch[order(anch$distance_to_CO92), ]
anch$lx  <- seq(25, 1000, length.out = nrow(anch))
anch$ly  <- 0.790
anch$lane <- 0.800 + 0.005 * (seq_len(nrow(anch)) - 1)
## the drop stops just left of the label (small gap, no overlap)
lead <- do.call(rbind, lapply(seq_len(nrow(anch)), function(i) {
  r <- anch[i, ]
  data.frame(x = c(r$distance_to_CO92, r$distance_to_CO92, r$lx - 5, r$lx - 5),
             y = c(r$called_frac, r$lane, r$lane, 0.797), grp = i)
}))
pb <- ggplot(reg, aes(distance_to_CO92, called_frac)) +
  geom_point(aes(colour = factor(is_outlier)), size = 0.35, alpha = 0.8) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = OI[["blue"]], linewidth = 0.5) +
  ## hollow ring = whitelist anchor; mapped to the shape scale only so it gets
  ## its own legend key (the colour scale stays within / residual outlier)
  geom_point(data = anch, aes(shape = "anchor"), colour = OI[["purple"]],
             size = 2.2, stroke = 0.7, fill = NA) +
  ## orthogonal leaders from each ring to its label
  geom_path(data = lead, aes(x, y, group = grp), colour = "grey45",
            linewidth = 0.3, lineend = "round", linejoin = "round") +
  geom_text(data = anch, aes(x = lx, y = ly, label = anchor_label), hjust = 0,
            size = GEOM_TEXT_MIN, colour = "grey20", fontface = "plain") +
  scale_colour_manual(values = c(`0` = "grey55", `1` = OI[["vermillion"]]),
                      labels = c("within", "residual outlier"), name = NULL) +
  scale_shape_manual(values = c(anchor = 1), name = NULL,
                     labels = c(anchor = "whitelist anchor (ringed)")) +
  labs(x = "raw SNP distance to CO92 (Tier A matrix)", y = "called fraction") +
  theme_paper() + theme(legend.position = "top", legend.box = "horizontal",
                        legend.spacing.x = unit(6, "pt"),
                        ## free(type="space") frees the axis-space area, so reserve
                        ## room for b's own x-axis title explicitly
                        plot.margin = margin(t = 2, r = 2, b = 9, l = 2, unit = "mm"))

## (c) both sides of one alignment
dd$ref_pct <- as.numeric(dd$ref_pct); dd$qry_pct <- as.numeric(dd$qry_pct)
dd$plasmid_submitted <- factor(dd$plasmid_submitted, levels = c("yes", "no"),
                               labels = c("plasmids submitted", "no plasmid submitted"))
lab5 <- dd[nzchar(dd$strain_label), ]
## The five reference-side flags are a tight cluster (ref 93.0-94.2, qry 99.5-
## 100.0), so their labels fan out below-left. They are ordered left-to-right by
## the POINTS' x and laid out on a monotone grid, which guarantees that no two
## leader lines can cross (repel cannot guarantee this).
lab5    <- lab5[order(lab5$ref_pct), ]
## all five short labels in ONE row, and each leader routed orthogonally:
## point -> down to its own lane (above the label row) -> across -> down to the
## label. The lanes rise with label order, which provably means no line crosses
## another line or any label box.
## x positions chosen so no label straddles the x = 95 % reference rule
## (KM 567 ends before 95, C-830 starts after it); assumes the five flags
stopifnot(nrow(lab5) == 5L)
lab5$lx   <- min(dd$ref_pct) + c(0.05, 1.05, 2.35, 3.60, 4.85)
lab5$ly   <- 98.60
## text spans ly +/- ~0.07 units, so lanes and drops stay at >= 98.70
lab5$lane <- 98.74 + 0.05 * (seq_len(nrow(lab5)) - 1)
## the vertical drop ends ~0.06 units LEFT of the text start so the line can
## never touch a glyph even allowing for font metric rounding
lead5 <- do.call(rbind, lapply(seq_len(nrow(lab5)), function(i) {
  r  <- lab5[i, ]
  xe <- r$lx - 0.06
  data.frame(x = c(r$ref_pct, r$ref_pct, xe, xe),
             y = c(r$qry_pct, r$lane, r$lane, 98.70), grp = i)
}))
gap_lo <- max(lab5$ref_pct); gap_hi <- min(dd$ref_pct[dd$ref_pct >= 95])

pc <- ggplot(dd, aes(ref_pct, qry_pct, colour = plasmid_submitted)) +
  annotate("rect", xmin = gap_lo, xmax = gap_hi, ymin = -Inf, ymax = Inf,
           fill = "grey88", colour = NA) +
  geom_vline(xintercept = 95, linetype = "dashed", linewidth = 0.5) +
  geom_point(size = 0.7, alpha = 0.85) +
  ## orthogonal leader (see lead5 above)
  geom_path(data = lead5, aes(x, y, group = grp), colour = "grey45",
            linewidth = 0.3, lineend = "round", linejoin = "round") +
  geom_text(data = lab5, aes(x = lx, y = ly, label = strain_label), hjust = 0,
            size = GEOM_TEXT_MIN, colour = "grey20", fontface = "plain") +
  scale_colour_manual(values = c(OI[["sky"]], OI[["vermillion"]]), name = NULL) +
  labs(x = "AlignedBases, reference side (%)", y = "AlignedBases, query side (%)") +
  theme_paper() + theme(legend.position = "top",
                        ## free(type="space") frees the axis-space area, so reserve
                        ## room for c's own x-axis title explicitly
                        plot.margin = margin(t = 2, r = 2, b = 9, l = 2, unit = "mm"))

## (d) CheckM
ck$completeness <- as.numeric(ck$completeness)
below <- sum(ck$completeness < 98); cka <- ck[nzchar(ck$anchor_label), ]

## anchors are marked with coloured asterisks + a legend instead of text, so no
## label has to sit on the x = 98 % rule
ANC_COL <- c(Angola = OI[["blue"]], Antiqua = OI[["purple"]],
             `Pestoides F` = OI[["vermillion"]])
stopifnot(all(cka$anchor_label %in% names(ANC_COL)))
ck_brk  <- seq(floor(min(ck$completeness)), ceiling(max(ck$completeness)) + 0.25, by = 0.25)
## keep clear headroom so the asterisks sit well above the tallest bar
ymark   <- max(table(cut(ck$completeness, ck_brk, right = FALSE))) * 1.16
pd <- ggplot(ck, aes(completeness)) +
  geom_histogram(binwidth = 0.25, fill = OI[["green"]], colour = NA) +
  geom_vline(xintercept = 98, linetype = "dashed", linewidth = 0.5) +
  geom_point(data = cka, aes(x = completeness, y = ymark, colour = anchor_label),
             shape = 8, size = 2.6, stroke = 0.9) +
  scale_colour_manual(values = ANC_COL, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.10))) +
  labs(x = "CheckM completeness", y = "genomes") + theme_paper() +
  theme(legend.position = "top", legend.margin = margin(0, 0, 0, 0))

## panel b gets more width than c so its 256 points spread out
## b gets its own full-width row (its 7 accessions need ~130 mm of plot width to
## sit in one non-overlapping label row); c and d share the third row.
## free(type="space") releases b/d from the shared axis-space column so their
## y-axis titles sit right next to their own tick labels.
fig2 <- (pa / free(pb, type = "space") / (free(pc, type = "space") | pd)) +
  plot_layout(heights = c(0.5, 1.5, 1.28)) +
  plot_annotation(tag_levels = "a")

# 采用和满意图一致的宽扁横版比例（宽 220 mm，高 115 mm）
save_fig(fig2, "fig3_quality_metrics", width_mm = 180, height_mm = 152)

## explicit output paths (results/figures -> ../figures)
for (f in c("fig3_quality_metrics.pdf", "fig3_quality_metrics.png")) {
  p <- file.path(FIGDIR, f)
  cat(sprintf("  output: %s\n", normalizePath(p, mustWork = FALSE)))
}

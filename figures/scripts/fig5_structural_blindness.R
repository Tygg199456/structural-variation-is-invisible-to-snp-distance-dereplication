#!/usr/bin/env Rscript
# Fig.5 -- dereplication is blind to structural variation (headline).
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
session_stamp("fig5_structural_blindness.R")

trk  <- read_fig("fig4b_co92_tracks.csv")
con  <- read_fig("fig4c_constraints.csv")
sw   <- read_fig("fig4d_sweep.csv")

XL <- 4653728 / 1000            # chromosome length in kb

## (b) the CO92 pair drawn on the reference coordinate
XL <- 4653728 / 1000            # chromosome length, kb
trk$start <- as.numeric(trk$start); trk$end <- as.numeric(trk$end)
trk$y <- ifelse(grepl("LANL", trk$label), 2, 1)
pgm <- trk[trk$length_bp > 100000, ]
pb <- ggplot() +
  annotate("rect", xmin = 0, xmax = XL, ymin = 1.75, ymax = 2.25,
           fill = "grey93", colour = "grey55", linewidth = 0.5) +
  geom_segment(data = trk, aes(x = start / 1000, xend = end / 1000,
                               y = y, yend = y),
               linewidth = 2.4, colour = OI[["vermillion"]]) +
  ## short tick labels: these sit in the shared left column, so long strings
  ## would push the y-axis titles of panels b and c far from their axes
  scale_y_continuous(breaks = c(1, 2),
                     labels = c("CO92 reference",
                                "CO92 (LANL,\npgm\u2212/pPCP1\u2212)")) +
  coord_cartesian(ylim = c(0.6, 2.75)) +
  labs(x = "CO92 chromosome position (kb)", y = NULL) + theme_paper() +
  theme(panel.grid = element_blank())

## (c) the four constraints
con$constraint <- factor(con$constraint, levels = con$constraint)
pc <- ggplot(con, aes(constraint, retained)) +
  geom_col(fill = OI[["green"]], width = 0.6) +
  geom_text(aes(label = retained), vjust = -0.35, size = GEOM_TEXT_MIN) +
  scale_y_continuous(limits = c(0, max(con$retained) * 1.12)) +
  labs(x = NULL, y = "genomes retained (Tier A)") + theme_paper() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

## (d) the three conditions of the threshold-selection rule, one panel each
sw$threshold <- as.numeric(sw$threshold)
sw$dup <- as.numeric(sw$dup_pairs_merged)
sw$max_size <- as.numeric(sw$max_size)
sw$max_diameter <- as.numeric(sw$max_diameter)
t_sel <- 20
mk <- function(y, yint, ylab) {
  ggplot(sw, aes(threshold, .data[[y]])) +
    geom_line(colour = OI[["blue"]], linewidth = 0.5) +
    geom_point(colour = OI[["blue"]], size = 1.2) +
    geom_hline(yintercept = yint, linetype = "dashed", linewidth = 0.5, colour = "grey40") +
    geom_vline(xintercept = t_sel, linetype = "dotted", linewidth = 0.5,
               colour = OI[["vermillion"]]) +
    scale_x_continuous(breaks = c(0, 10, 20, 30, 50, 75, 100)) +
    labs(x = "SNP threshold t", y = ylab) + theme_paper()
}
d1 <- mk("dup", 5, "same-strain pairs\nmerged (of 5)")
d2 <- mk("max_diameter", 40, "max cluster\ndiameter")
d3 <- mk("max_size", 25, "max cluster\nsize")
pd <- (free(d1, type = "space") | d2 | d3) + plot_layout(axes = "collect")

## free(type="space") keeps panel b's and c's y-axis titles next to their own
## axes instead of at the far left of the shared column
fig4 <- pb / free(pc, type = "space") / pd +
  plot_layout(heights = c(1.0, 0.75, 0.85)) +
  plot_annotation(tag_levels = "a")
save_fig(fig4, "fig5_structural_blindness", width_mm = 172, height_mm = 168)

cat(sprintf("\nKEY NUMBERS\n  constraints: %s\n",
            paste(sprintf("%s=%d", con$constraint, con$retained), collapse = " ")))
cat(sprintf("  CO92 tracks rows (>=10 kb runs): %d\n", nrow(trk)))

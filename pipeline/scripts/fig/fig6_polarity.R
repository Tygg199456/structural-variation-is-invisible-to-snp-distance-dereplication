#!/usr/bin/env Rscript
# Fig.5 -- polarity and the background baseline for polarity calls.
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
source("scripts/fig/theme_paper.R"); session_stamp("fig5_polarity.R")
suppressPackageStartupMessages(library(ggrepel))

reg <- read_fig("fig5a_region_polarity.csv")
bg  <- read_fig("fig5b_background.csv")
bd  <- read_fig("fig5b_block_delta.csv")

POL <- c(derived_LOSS = "derived loss", CO92_lineage_GAIN = "CO92-lineage gain",
         AMBIGUOUS = "ambiguous", OUTGROUP_UNINFORMATIVE = "outgroup uninformative")
PCOL <- c(`derived loss` = OI[["blue"]], `CO92-lineage gain` = OI[["vermillion"]],
          ambiguous = OI[["yellow"]], `outgroup uninformative` = OI[["grey"]])

## (a) the four-way rule, in the plane of the two outgroups' callable fractions
reg$og1_callable <- as.numeric(reg$og1_callable)
reg$og2_callable <- as.numeric(reg$og2_callable)
reg$n_genomes <- as.numeric(reg$n_genomes)
reg$polarity <- factor(POL[reg$polarity], levels = POL)
hl1 <- reg[reg$highlight == "primary", ]; hl2 <- reg[reg$highlight == "secondary", ]
pa <- ggplot(reg, aes(og1_callable, og2_callable)) +
  annotate("rect", xmin = 0.70, xmax = 1.02, ymin = -0.02, ymax = 1.02,
           fill = OI[["blue"]], alpha = 0.06) +
  annotate("rect", xmin = -0.02, xmax = 1.02, ymin = 0.70, ymax = 1.02,
           fill = OI[["blue"]], alpha = 0.06) +
  annotate("rect", xmin = -0.02, xmax = 0.20, ymin = -0.02, ymax = 0.20,
           fill = OI[["vermillion"]], alpha = 0.08) +
  geom_vline(xintercept = c(0.20, 0.70), linetype = "dotted", linewidth = 0.5) +
  geom_hline(yintercept = c(0.20, 0.70), linetype = "dotted", linewidth = 0.5) +
  geom_point(aes(size = n_genomes, fill = polarity), shape = 21, colour = "grey25",
             stroke = 0.3, alpha = 0.9) +
  geom_point(data = hl1, shape = 21, size = 6, colour = OI[["vermillion"]],
             fill = NA, stroke = 0.8) +
  geom_text_repel(data = rbind(hl1, hl2),
                  aes(label = sprintf("%s +%d kb\n%d genomes", format(start, big.mark = ","),
                                      round(length_bp / 1000), n_genomes)),
                  size = GEOM_TEXT_MIN, min.segment.length = 0, seed = 1, lineheight = 0.95) +
  scale_fill_manual(values = PCOL, name = NULL) +
  scale_size_area(max_size = 7, name = "genomes") +
  coord_cartesian(xlim = c(-0.02, 1.02), ylim = c(-0.02, 1.02)) +
  labs(x = "callable fraction, outgroup 1 (IP 32953)",
       y = "callable fraction, outgroup 2 (IP 31758)") +
  theme_paper() + theme(legend.position = "right", legend.box = "vertical",
                        legend.spacing.y = unit(1, "pt"))

## (b) measure the background first
bg$delta <- as.numeric(bg$delta); bd$delta <- as.numeric(bd$delta)
q95 <- quantile(bg$delta, 0.95)
bd$within <- factor(bd$within_background, levels = c(1, 0),
                    labels = c("within background", "exceeds p95"))
pb <- ggplot(bg, aes(delta)) +
  geom_histogram(bins = 60, fill = "grey82", colour = NA) +
  geom_vline(xintercept = q95, linetype = "dashed", linewidth = 0.5,
             colour = OI[["vermillion"]]) +
  geom_rug(data = bd, aes(x = delta, colour = within), linewidth = 0.6,
           length = unit(5, "pt")) +
  scale_colour_manual(values = c(`within background` = "grey30",
                                 `exceeds p95` = OI[["vermillion"]]), name = NULL) +
  scale_x_sqrt() +
  labs(x = "|callable(OG1) - callable(OG2)| over 25 kb windows (sqrt scale)", y = "windows") +
  theme_paper() + theme(legend.position = "top")

fig5 <- (pa | pb) + plot_layout(widths = c(1.15, 1)) + plot_annotation(tag_levels = "a")
save_fig(fig5, "fig6_polarity", width_mm = 172, height_mm = 85)

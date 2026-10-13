#!/usr/bin/env Rscript
# Fig.3 -- spatial morphology distinguishes deletion from mapping failure.
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
source("scripts/fig/theme_paper.R"); session_stamp("fig3_spatial_morphology.R")
suppressPackageStartupMessages(library(ggrepel))

bl <- read_fig("fig3a_block_lengths.csv"); rec <- read_fig("fig3b_recurrence.csv")
lv <- read_fig("fig3b_level_frac.csv");    pa6 <- read_fig("fig3c_pgm_align.csv")
gn <- read_fig("fig3c_pgm_genes.csv")

## (a) block-length distribution -- IS-scale scatter vs a few large blocks
bl$length_bp <- as.numeric(bl$length_bp)
pa <- ggplot(bl, aes(length_bp)) +
  geom_histogram(bins = 50, fill = OI[["sky"]], colour = NA) +
  geom_vline(xintercept = c(3000, 10000), linetype = c("dotted", "dashed"),
             linewidth = 0.5, colour = c("grey40", OI[["vermillion"]])) +
  scale_x_log10(labels = label_number(scale_cut = cut_short_scale())) +
  labs(x = "uncalled run length (bp, log scale)", y = "runs") + theme_paper()

## (b) recurrence spectrum + level breakdown inset
rec <- rec[order(-rec$n_genomes), ]; rec$rank <- seq_len(nrow(rec))
n_once <- sum(rec$n_genomes == 1); n_rec10 <- sum(rec$n_genomes >= 10)
pb_main <- ggplot(rec, aes(rank, n_genomes)) +
  geom_col(fill = OI[["blue"]], width = 0.9) +
  geom_hline(yintercept = 10, linetype = "dashed", linewidth = 0.5, colour = "grey40") +
  scale_y_log10() +
  labs(x = "distinct >=10 kb coordinate, ranked", y = "genomes sharing it") +
  theme_paper()
lv$frac <- as.numeric(lv$frac)
lv$level <- factor(lv$level, levels = c("complete", "chromosome", "scaffold", "contig"))
pb <- pb_main
## (c) pgm positive control: assemblies aligned at the pgm locus + gene track
pa6$length_bp <- as.numeric(pa6$length_bp)
pa6 <- pa6[order(-pa6$length_bp), ]
pa6$strain <- factor(pa6$strain, levels = rev(unique(pa6$strain)))
x0 <- min(as.numeric(pa6$start))
pc_top <- ggplot(pa6) +
  geom_segment(aes(x = (as.numeric(start) - x0) / 1000, xend = (as.numeric(end) - x0) / 1000,
                   y = strain, yend = strain), linewidth = 3, colour = OI[["vermillion"]]) +
  labs(x = NULL, y = NULL) + theme_paper() +
  theme(axis.text = element_blank(), axis.ticks = element_blank())
gn$group <- factor(gn$group, levels = c("ybt", "irp", "fyuA", "pga"))
pat <- "(Ybt[A-Z]|Irp[0-9]|Pga[A-D]|HMWP[0-9]|FyuA|Psn)"
gn$lab <- ifelse(grepl(pat, gn$gene), sub(paste0("^.*?(", pat, ").*$"), "\\1", gn$gene),
                 ifelse(grepl("salycil-AMP", gn$gene), "YbtE",
                        ifelse(grepl("thioesterase", gn$gene), "YbtT",
                               ifelse(grepl("non-ribosomal", gn$gene), "HMWP2",
                                      ifelse(grepl("salicylate synthase", gn$gene), "Irp9",
                                             ifelse(grepl("signal transducer", gn$gene), "YbtX",
                                                    ifelse(grepl("N-deacetylase", gn$gene), "PgaB",
                                                           ifelse(grepl("porin", gn$gene), "PgaA",
                                                                  ifelse(grepl("biosynthesis protein", gn$gene), "PgaD", "")))))))))
pc_bot <- ggplot(gn) +
  geom_rect(aes(xmin = (as.numeric(start) - x0) / 1000, xmax = (as.numeric(end) - x0) / 1000,
                ymin = 0, ymax = 1, fill = group), colour = NA) +
  geom_text_repel(aes(x = (as.numeric(start) + as.numeric(end)) / 2 / 1000 - x0 / 1000,
                      y = 1, label = lab), size = 1.9, direction = "x", angle = 90,
                  vjust = 0, segment.size = 0.15, max.overlaps = 30, seed = 1) +
  scale_fill_manual(values = unname(OI[c("blue", "green", "purple", "orange")]), name = NULL) +
  scale_y_continuous(limits = c(0, 3.4), breaks = NULL) +
  labs(x = sprintf("position relative to %s (kb)", format(x0, big.mark = ",")), y = NULL) +
  theme_paper() + theme(legend.position = "bottom", panel.grid = element_blank())
pc <- pc_top / pc_bot + plot_layout(heights = c(1, 1.1))

fig3 <- (pa | pb) / pc + plot_layout(heights = c(1, 1.25)) +
  plot_annotation(tag_levels = "a")
save_fig(fig3, "fig4_spatial_morphology", width_mm = 172, height_mm = 165)

#!/usr/bin/env Rscript
# Fig.4 -- spatial morphology distinguishes deletion from mapping failure.
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
source("scripts/fig/theme_paper.R"); session_stamp("fig4_spatial_morphology.R")
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
## Label only named strains, and only when at least 3 rows clear of the last
## one -- the 33 numbered isolates would be unreadable at this panel height.
## The explanatory line goes in the subtitle so it cannot sit on the labels.
lev      <- levels(pa6$strain)
labs_vec <- ifelse(grepl("[A-Za-z]", lev), lev, "")
keep <- rep(FALSE, length(labs_vec)); last <- -99
for (i in seq_along(labs_vec)) if (nzchar(labs_vec[i]) && i - last >= 3) { keep[i] <- TRUE; last <- i }
labs_vec[!keep] <- ""
pc_top <- ggplot(pa6) +
  geom_segment(aes(x = (as.numeric(start) - x0) / 1000, xend = (as.numeric(end) - x0) / 1000,
                   y = strain, yend = strain), linewidth = 3, colour = OI[["vermillion"]]) +
  scale_y_discrete(labels = labs_vec) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.01))) +
  labs(x = "position relative to 2,135,548 (kb)", y = NULL,
       subtitle = sprintf(paste0("%d assemblies, each carrying the same 102.5 kb block ",
                                 "(block starts within %d bp, lengths within %d bp)"),
                          nrow(pa6), diff(range(as.numeric(pa6$start))),
                          diff(range(as.numeric(pa6$length_bp))))) +
  theme_paper() +
  theme(axis.ticks.y = element_blank(), panel.grid = element_blank(),
        plot.subtitle = element_text(size = BASE_SIZE - 1.5, colour = "grey25"))
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
## Gene labels: each sits directly above its own gene, no repel. In the dense
## ybt run (FyuA/Irp2/YbtP/YbtQ) repel had shoved labels sideways and drawn
## diagonal leader lines across the text; instead alternate the labels between
## two rows so nothing needs to move and no connector is drawn.
gn$mid  <- (as.numeric(gn$start) + as.numeric(gn$end)) / 2 / 1000 - x0 / 1000
gn      <- gn[order(gn$mid), ]
## stagger by VISIBLE label (not by gene row): the unlabelled rows would
## otherwise flip the parity, e.g. pushing PgaD up while PgaA sits low
vis      <- nzchar(gn$lab)
gn$ypos  <- NA_real_
gn$ypos[vis] <- 1.08 + ((cumsum(vis)[vis] - 1) %% 2) * 0.62
pc_bot <- ggplot(gn) +
  geom_rect(aes(xmin = (as.numeric(start) - x0) / 1000, xmax = (as.numeric(end) - x0) / 1000,
                ymin = 0, ymax = 1, fill = group), colour = NA) +
  ## hjust = 0 + angle = 90 -> rotated text grows upward from ypos, clear of the bars
  geom_text(data = gn[vis, ], aes(x = mid, y = ypos, label = lab),
            size = 1.9, angle = 90, hjust = 0) +
  scale_fill_manual(values = unname(OI[c("blue", "green", "purple", "orange")]), name = NULL) +
  scale_y_continuous(limits = c(0, 3.4), breaks = NULL) +
  labs(x = sprintf("position relative to %s (kb)", format(x0, big.mark = ",")), y = NULL) +
  theme_paper() + theme(legend.position = "bottom", panel.grid = element_blank())
pc <- pc_top / pc_bot + plot_layout(heights = c(1.35, 1))

fig3 <- (pa | pb) / pc + plot_layout(heights = c(1, 1.25)) +
  plot_annotation(tag_levels = "a")
save_fig(fig3, "fig4_spatial_morphology", width_mm = 172, height_mm = 165)

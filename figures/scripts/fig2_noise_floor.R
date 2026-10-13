#!/usr/bin/env Rscript
# Fig.2 -- the technical noise floor.  Descriptive text lives in the legend.
#
# Changes in this version:
#   - Wrapped long x-axis label in panel (b) into two lines to prevent truncation.
#   - Maintained enlarged fonts and vertically staggered strain labels.

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
session_stamp("fig2_noise_floor.R")

cf   <- read_fig("fig1a_called_frac.csv")
hist <- read_fig("fig1b_distance_hist.csv")
same <- read_fig("fig1b_same_strain.csv")
self <- read_fig("fig1c_self_blocks.csv")
ci   <- read_fig("fig1d_cluster_internal.csv")

## --- 定义全局字体放大主题扩展（可应用至所有子图）
theme_larger <- theme(
  axis.title = element_text(size = 11),   # 放大坐标轴标题
  axis.text  = element_text(size = 9),    # 放大坐标轴刻度数字
  strip.text = element_text(size = 10)    # 放大分面标签标题
)

## ---------------------------------------------------------------- (a)
cf$tier <- factor(cf$tier, levels = c("tierA", "tierB"),
                  labels = c(sprintf("Tier A (n=%d)", sum(cf$tier == "tierA")),
                             sprintf("Tier B (n=%d)", sum(cf$tier == "tierB"))))
pa <- ggplot(cf, aes(called_frac, fill = tier)) +
  geom_histogram(bins = 60, colour = NA, alpha = 0.9) +
  facet_wrap(~tier, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = c(OI[["sky"]], OI[["orange"]]), guide = "none") +
  labs(x = "called fraction (union variant sites)", y = "genomes") +
  theme_paper() + theme_larger

## ---------------------------------------------------------------- (b)
hist$bin_lo      <- as.numeric(hist$bin_lo)
hist$bin_hi      <- as.numeric(hist$bin_hi)
hist$n_pairs     <- as.numeric(hist$n_pairs)
hist$bin_lo_plot <- pmax(hist$bin_lo, 0.5)      # 0 cannot sit on a log axis
same$d           <- as.numeric(same$raw_distance_tierB)

ymax <- max(hist$n_pairs, na.rm = TRUE)
xhi  <- max(hist$bin_hi, na.rm = TRUE)

pb <- ggplot() +
  geom_rect(data = hist,
            aes(xmin = bin_lo_plot, xmax = bin_hi, ymin = 0, ymax = n_pairs),
            fill = "grey80", colour = NA) +
  # 底部短刻度线（rug ticks）
  geom_segment(data = same,
               aes(x = d, xend = d, y = 0, yend = ymax * 0.04),
               colour = OI[["vermillion"]], linewidth = 0.8) +
  # 顶部标签放大并优化排布，确保上下错开、绝不重叠
  geom_text_repel(data = same,
                  aes(x = d, y = ymax * 0.04, label = strain),
                  size = 3.3,                     # 适当放大字体
                  colour = OI[["vermillion"]],
                  fontface = "bold",              # 加粗字体
                  direction = "both",             # 允许横向与纵向自由排布
                  force = 8,                      # 强排斥力，让标签充分散开
                  box.padding = 0.5,              # 增加文本框留白
                  point.padding = 0.3,
                  max.overlaps = Inf,             # 允许所有标签完整显示
                  min.segment.length = 0, segment.size = 0.35,
                  segment.colour = OI[["vermillion"]], seed = 1) +
  scale_x_log10(breaks = c(1, 10, 100, 1000),
                labels = c("1", "10", "100", "1000")) +
  scale_y_continuous(labels = label_number(scale_cut = cut_short_scale()),
                     expand = expansion(mult = c(0, 0.38))) + # 扩大顶部留白，容纳多层错开的标签
  coord_cartesian(xlim = c(0.5, xhi), clip = "off") +
  # 标题已通过 \n 换行，防止被截断
  labs(x = "raw SNP distance, all pairs\n(Tier B matrix, log scale)",
       y = "pairs") +
  theme_paper() + theme_larger

## ---------------------------------------------------------------- (c)
self$length_bp <- as.numeric(self$length_bp)
self <- self[self$length_bp >= 500, ]
pc <- ggplot(self, aes(length_bp)) +
  geom_histogram(bins = 45, fill = OI[["green"]], colour = NA) +
  geom_vline(xintercept = 10000, linetype = "dashed", linewidth = 0.5,
             colour = OI[["vermillion"]]) +
  scale_x_log10(labels = label_number(scale_cut = cut_short_scale())) +
  labs(x = "uncalled run length (bp, log scale)", y = "runs") +
  theme_paper() + theme_larger

## ---------------------------------------------------------------- (d)
ci$distance_snp <- as.numeric(ci$distance_snp)
ci$n_pairs      <- as.numeric(ci$n_pairs)
pd <- ggplot(ci, aes(distance_snp, n_pairs)) +
  geom_col(fill = OI[["purple"]], width = 0.85) +
  geom_vline(xintercept = 19, linetype = "dashed", linewidth = 0.5,
             colour = "grey40") +
  labs(x = "within-cluster pairwise SNP distance", y = "pairs") +
  theme_paper() + theme_larger

fig2 <- (pa | pb) / (pc | pd) + plot_annotation(tag_levels = "a")
save_fig(fig2, "fig2_noise_floor", width_mm = 172, height_mm = 150)

## ---------------------------------------------------------------- checks
n_zero   <- sum(hist$n_pairs[hist$bin_lo <= 0 & hist$bin_hi > 0], na.rm = TRUE)
n_500    <- nrow(self)
n_10k    <- sum(self$length_bp >= 10000)
ci_long  <- rep(ci$distance_snp, ci$n_pairs)
ci_max   <- max(ci$distance_snp[ci$n_pairs > 0])

cat("\nKEY NUMBERS\n")
n_500_chr <- nrow(self[self$replicon == "NC_003143", ])
cat(sprintf("  self-comparison runs >=500 bp : %d all replicons, %d on the chromosome    <- text says 165 / 157\n", n_500, n_500_chr))
cat(sprintf("  longest self-comparison run   : %s bp\n",
            format(max(self$length_bp), big.mark = ",")))
cat(sprintf("  runs >=10 kb                  : %d    <- must be 0\n", n_10k))
cat(sprintf("  same-strain raw distances (B) : %s\n",
            paste(sort(same$d), collapse = ", ")))
cat(sprintf("  pairs at distance 0           : %s    <- quote this in the legend\n",
            format(n_zero, big.mark = ",")))
cat(sprintf("  cluster internal p50/p95/max  : %d / %d / %d    <- text says 6 / 11 / 19\n",
            as.integer(quantile(ci_long, .50)),
            as.integer(quantile(ci_long, .95)),
            ci_max))

if (n_10k != 0)
  stop("self-comparison produced a run of >=10 kb; Results 3.1 must be rechecked")
if (ci_max > 19)
  warning(sprintf("largest within-cluster distance is %d, but the text and the dashed line in (d) say 19", ci_max))
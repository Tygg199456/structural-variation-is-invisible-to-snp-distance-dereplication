#!/usr/bin/env Rscript
# Shared theme + helpers for paper 1 figures.
# R version pinned: /workdir/conda/envs/r-4.4.3/bin/Rscript
#
# Conventions
#   - every figure reads only CSVs from results/figures/ (written by
#     scripts/fig/prep_figure_data.py); no .npy/.npz/GFF access from R
#   - panel letters come ONLY from plot_annotation(tag_levels = "a"); panel
#     titles are NOT used (theme_paper clears them) -- descriptive text belongs
#     in the figure legend
#   - output: results/figures/figN_<name>.pdf (cairo_pdf, vector) and .png
#     (ragg, 600 dpi)
#   - colourblind-safe palette (Okabe-Ito); full width 180 mm, half width 85 mm

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(data.table)
  library(scales)
})

FIGDIR <- "results/figures"
BASE_SIZE <- 8
BASE_FAMILY <- "Arial"
GEOM_TEXT_MIN <- 2.2          # pt (~6 pt); never go below this

# Okabe-Ito, colourblind safe
OI <- c(orange = "#E69F00", sky = "#56B4E9", green = "#009E73",
        yellow = "#F0E442", blue = "#0072B2", vermillion = "#D55E00",
        purple = "#CC79A7", grey = "#999999")

theme_paper <- function(base_size = BASE_SIZE) {
  theme_bw(base_size = base_size, base_family = BASE_FAMILY) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      strip.background = element_rect(fill = "grey95", colour = NA),
      strip.text = element_text(face = "bold", size = base_size - 0.5),
      axis.text = element_text(colour = "grey20"),
      axis.title = element_text(size = base_size),
      legend.key.size = unit(0.35, "cm"),
      legend.background = element_blank(),
      plot.tag = element_text(face = "bold", size = base_size + 1),
      plot.title = element_blank()      # panel titles are a convention, not a choice
    )
}

read_fig <- function(name) fread(file.path(FIGDIR, name), data.table = FALSE)

save_fig <- function(plot, name, width_mm = 172, height_mm = 120) {
  w <- width_mm / 25.4
  h <- height_mm / 25.4
  pdf_f <- file.path(FIGDIR, paste0(name, ".pdf"))
  png_f <- file.path(FIGDIR, paste0(name, ".png"))
  ggsave(pdf_f, plot, width = w, height = h, device = cairo_pdf)
  ggsave(png_f, plot, width = w, height = h, dpi = 600, device = ragg::agg_png)
  if (!file.exists(pdf_f) || file.info(pdf_f)$size < 10240)
    stop(sprintf("save_fig: %s missing or suspiciously small", pdf_f))
  cat(sprintf("  wrote %s  (%d x %d mm)\n", pdf_f, width_mm, height_mm))
  cat(sprintf("  wrote %s\n", png_f))
  ## show in the RStudio Plots pane when run interactively (Source / Run)
  if (interactive()) print(plot)
}

session_stamp <- function(script) {
  cat(sprintf("== %s ==\n", script))
  cat(sprintf("   R %s.%s ; data from %s/\n",
              R.version$major, R.version$minor, FIGDIR))
}

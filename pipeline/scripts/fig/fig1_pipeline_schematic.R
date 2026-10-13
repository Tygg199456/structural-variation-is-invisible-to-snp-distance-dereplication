#!/usr/bin/env Rscript
# Fig.1 -- how the distance is constructed, and what it cannot see.
#
# Replaces two monospace blocks previously in Methods 2.1/2.2: character-aligned
# text cannot survive journal typesetting, and at revision stage a figure is a
# separate file.
#
# CONVENTIONS (two deliberate relaxations, documented so they are not read as
# drift): (1) ggrepel is NOT used -- in a schematic a label's position carries
# meaning and must stay inside its box; (2) panels contain text, but no
# explanatory text -- every string is an element NAME or a value COMPUTED from
# the CSVs. All explanatory sentences live in the legend.
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
source("scripts/fig/theme_paper.R")
session_stamp("fig1_pipeline_schematic.R")

read_fig2 <- function(name) fread(file.path(FIGDIR, name), sep = ";", data.table = FALSE)

nodes <- read_fig2("fig_m1_nodes.csv")
edges <- read_fig2("fig_m1_edges.csv")
elabs <- read_fig2("fig_m1_edge_labels.csv")
sites <- read_fig("fig_m1_sites.csv")

FILL <- c(plain = "white",
          blue = scales::alpha(OI[["sky"]], 0.18),
          orange = scales::alpha(OI[["vermillion"]], 0.14))
BORDER <- c(plain = "grey30", blue = OI[["blue"]], orange = OI[["vermillion"]])

## =============================================================== PANEL A
nodes$sublabel <- gsub("\\\\n", "\n", nodes$sublabel)
nodes$cx <- (nodes$xmin + nodes$xmax) / 2
nodes$n_sub <- ifelse(nzchar(nodes$sublabel), lengths(strsplit(nodes$sublabel, "\n")), 0)
nodes$y_lab <- ifelse(nodes$n_sub == 0, (nodes$ymin + nodes$ymax) / 2, nodes$ymax - 4.2)
nodes$y_sub <- nodes$ymax - 4.2 - 3.6

ed <- edges[order(edges$edge_id, edges$seq), ]
ed_arrow <- ed[ed$arrow_end == 1 & ed$dashed == 0, ]
ed_dash  <- ed[ed$dashed == 1, ]
ed_plain <- ed[ed$arrow_end == 0 & ed$dashed == 0, ]
AR <- arrow(length = unit(1.5, "mm"), type = "closed")

pa <- ggplot() +
  geom_path(data = ed_plain, aes(x, y, group = edge_id), colour = "grey30", linewidth = 0.5) +
  geom_path(data = ed_arrow, aes(x, y, group = edge_id), colour = "grey30",
            linewidth = 0.5, arrow = AR, lineend = "butt", linejoin = "mitre") +
  geom_path(data = ed_dash, aes(x, y, group = edge_id), colour = "grey55",
            linewidth = 0.5, linetype = "22", arrow = AR) +
  geom_rect(data = nodes, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                              fill = style, colour = style), linewidth = 0.5) +
  geom_text(data = nodes, aes(cx, y_lab, label = label), size = GEOM_TEXT_MIN + 0.45,
            fontface = "bold", colour = "grey10", lineheight = 0.92) +
  geom_text(data = nodes[nodes$n_sub > 0, ], aes(cx, y_sub, label = sublabel),
            size = GEOM_TEXT_MIN, colour = "grey30", lineheight = 1.0, vjust = 1) +
  geom_text(data = elabs, aes(x, y, label = label, hjust = hjust), size = GEOM_TEXT_MIN,
            colour = "grey35") +
  scale_fill_manual(values = FILL, guide = "none") +
  scale_colour_manual(values = BORDER, guide = "none") +
  coord_cartesian(xlim = c(-1, 101), ylim = c(0, 108), expand = FALSE, clip = "off") +
  labs(x = NULL, y = NULL) + theme_paper() +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        panel.grid = element_blank(), panel.border = element_blank(),
        plot.margin = margin(2, 2, 2, 2))

## =============================================================== PANEL B
N <- max(sites$site_index)
PITCH <- 5.25; W <- 4.75; X0 <- 14
sites$xmin <- X0 + (sites$site_index - 1) * PITCH
sites$xmax <- sites$xmin + W
sites$cx <- (sites$xmin + sites$xmax) / 2
YY <- c(A = 78, B = 62); H <- 10
sites$ymin <- as.numeric(YY[as.character(sites$genome)])
sites$ymax <- as.numeric(sites$ymin) + H

SFILL <- c(called_match = scales::alpha(OI[["sky"]], 0.30),
           called_diff = OI[["vermillion"]],
           uncalled = "white")
SLAB <- c(called_match = "called, matches reference",
          called_diff = "called, differs",
          uncalled = "uncalled (N, - or masked)")
sites$state <- factor(sites$state, levels = names(SFILL), labels = SLAB)

a <- as.character(sites$state[sites$genome == "A"])
b <- as.character(sites$state[sites$genome == "B"])
both_called <- a != SLAB[["uncalled"]] & b != SLAB[["uncalled"]]
differs <- both_called & (a == SLAB[["called_diff"]] | b == SLAB[["called_diff"]])
C_n <- sum(both_called)
d_n <- sum(differs)
wide <- data.frame(
  site_index = as.integer(sites$site_index[sites$genome == "A"]),
  enters = ifelse(!both_called, "\u2013", ifelse(differs, "C+d", "C")))
wide$cx <- X0 + (wide$site_index - 1) * PITCH + W / 2

gap <- sites[sites$genome == "B" & sites$state == SLAB[["uncalled"]], ]
gap_lo <- min(gap$xmin); gap_hi <- max(gap$xmax); gap_n <- nrow(gap)
brk <- data.frame(x = c(gap_lo, gap_lo, gap_hi, gap_hi), y = c(57.5, 54.5, 54.5, 57.5))
gap_lab <- data.frame(x = (gap_lo + gap_hi) / 2, y = 50,
                      label = sprintf("deletion in genome B (%d sites)", gap_n))
diag_ln <- gap
tally <- data.frame(x = 17, y = c(29, 23.5, 17),
                    label = c(sprintf("jointly called sites C = %d", C_n),
                              sprintf("differing sites d = %d", d_n),
                              sprintf("raw distance = %d", d_n)),
                    bold = c(FALSE, FALSE, TRUE))

pb <- ggplot() +
  geom_rect(data = sites, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                              fill = state), colour = "grey35", linewidth = 0.5) +
  geom_segment(data = diag_ln, aes(x = xmin + 0.6, xend = xmax - 0.6,
                                   y = ymin + 0.6, yend = ymax - 0.6),
               colour = "grey65", linewidth = 0.5) +
  geom_path(data = brk, aes(x, y), colour = "grey55", linewidth = 0.5) +
  geom_text(data = gap_lab, aes(x, y, label = label), size = GEOM_TEXT_MIN, colour = "grey35") +
  geom_text(data = wide, aes(cx, 40, label = as.character(enters)), size = GEOM_TEXT_MIN, colour = "grey25") +
  geom_rect(aes(xmin = X0, xmax = X0 + N * PITCH - (PITCH - W), ymin = 11, ymax = 34),
            fill = NA, colour = "grey35", linewidth = 0.5) +
  geom_text(data = tally, aes(as.numeric(x), as.numeric(y), label = label,
                              fontface = ifelse(bold, "bold", "plain")),
            size = GEOM_TEXT_MIN + 0.3, colour = "grey10", hjust = 0) +
  scale_fill_manual(values = setNames(unname(SFILL), SLAB), name = NULL, breaks = SLAB) +
  scale_y_continuous(breaks = c(83, 67, 40), labels = c("genome A", "genome B", "enters"),
                     limits = c(8, 95), expand = expansion(0)) +
  scale_x_continuous(limits = c(0, 101), expand = expansion(0)) +
  labs(x = NULL, y = NULL) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(colour = "grey35"))) +
  theme_paper() +
  theme(axis.text.x = element_blank(), axis.ticks = element_blank(),
        axis.text.y = element_text(size = BASE_SIZE - 1.5, colour = "grey25"),
        panel.grid = element_blank(), panel.border = element_blank(),
        legend.position = "bottom", legend.margin = margin(0, 0, 0, 0),
        plot.margin = margin(2, 2, 2, 2))

fig1 <- (pa | pb) + plot_layout(widths = c(1, 1.18)) + plot_annotation(tag_levels = "a")
save_fig(fig1, "fig1_pipeline_schematic", width_mm = 172, height_mm = 96)

cat("\nKEY NUMBERS (counted from fig_m1_sites.csv)\n")
cat(sprintf("  sites drawn          : %d\n", N))
cat(sprintf("  jointly called C     : %d\n", C_n))
cat(sprintf("  differing d          : %d\n", d_n))
cat(sprintf("  uncalled in genome B : %d\n", gap_n))
cat(sprintf("  raw distance         : %d\n", d_n))
stopifnot(C_n + gap_n == N)

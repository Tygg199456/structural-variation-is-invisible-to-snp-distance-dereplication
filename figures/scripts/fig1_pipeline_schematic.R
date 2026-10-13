#!/usr/bin/env Rscript
# Fig.1 -- how the distance is constructed, and what it cannot see.
#
# Panel (a) is drawn in the CheckM Fig.2 idiom (Parks et al. 2015):
#   * two shaded stage regions ("Pass 1" / "Pass 2") with a title in the headroom
#   * node SHAPE encodes role, and the key decodes it with terms that are true
#     of EVERY node carrying that shape:
#         trapezoid = step
#         cylinder  = snippy output file       (snps.vcf, snps.aligned.fa)
#         rectangle = other data
#   * ONE node is accented in colour -- `raw distance`, the quantity the
#     dereplication threshold is applied to. That is the figure's point.
#   * arrows are orthogonal; every label is wrapped to the usable width of its
#     own shape, and every text size is >= GEOM_TEXT_MIN.
# Geometry (regions/nodes/edges) is DESIGN, not data, and is therefore kept
# in-script where it is easy to tune; panel (b) still reads its CSV.
#
# CONVENTIONS: ggrepel is NOT used (in a schematic a label's position carries
# meaning and must stay inside its box); panels contain element names and counts
# only -- explanatory sentences live in the figure legend.
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

## ============================================================ helpers ======
wrap_txt <- function(x, width) {
  if (length(width) == 1) width <- rep(width, length(x))
  vapply(seq_along(x), function(i) {
    if (!nzchar(x[i])) return("")
    if (grepl("\n", x[i], fixed = TRUE)) return(gsub(" *\n *", "\n", x[i]))
    paste(strwrap(x[i], width = width[i]), collapse = "\n")
  }, character(1))
}

## text sizes: the project floor is a FLOOR, never a target to undercut
LAB_SIZE   <- GEOM_TEXT_MIN + 0.2      # node names
SUB_SIZE   <- GEOM_TEXT_MIN            # node sub-labels (counts, dimensions)
TITLE_SIZE <- GEOM_TEXT_MIN + 0.6      # stage titles and "Key"

## shape vocabulary
SH_STEP <- "step"; SH_FILE <- "file"; SH_DATA <- "data"
SLANT_F <- 0.22                        # trapezoid slant, as a fraction of h
RY_F    <- 0.12                        # cylinder cap depth, as a fraction of h

NODE_FILL   <- "white"
NODE_BORDER <- "grey25"
ACC_BORDER  <- OI[["blue"]]
ACC_FILL    <- scales::alpha(OI[["blue"]], 0.13)
REG_FILL    <- "grey91"
REG_BORDER  <- "grey72"

## --- node geometry ---------------------------------------------------------
## x,y = centre; w,h = box; wrap = characters per line (set from the USABLE
## width of the shape at the text's own height, not from w)
nodes <- data.frame(
  id    = c("assembly", "snippy",  "vcf",    "pass1",  "aligned", "union",
            "pass2",    "gmatrix", "raw",    "scaled"),
  x     = c(12,         38,        63,       88,       62,        88,
            61,         87,        87,       87),
  y     = c(86,         86,        86,       86,       68,        68,
            38,         38,        22,       8),
  w     = c(20,         22,        22,       20,       22,        20,
            24,         22,        22,       22),
  h     = c(13.5,       13.5,      13.5,     13.5,     16,        16,
            14,         14,        11,       11),
  shape = c(SH_DATA,    SH_STEP,   SH_FILE,  SH_STEP,  SH_FILE,   SH_DATA,
            SH_STEP,    SH_DATA,   SH_DATA,  SH_DATA),
  label = c("assembly", "snippy 4.6.0", "snps.vcf", "pass 1", "snps.aligned.fa",
            "union variant sites", "pass 2", "genotype matrix",
            "raw distance", "scaled distance"),
  sub   = c("one of 3,441", "CO92 reference", "differences only",
            "reads snps.vcf", "reference length,\nthree states",
            "Tier A 12,926\nTier B 67,711",
            "reads union sites\nACGT\u21921\u20134, else 0",
            "int8\n256 \u00d7 12,926\n3,441 \u00d7 67,711",
            "threshold applied here", "raw \u00d7 n_sites / C\ndiagnostic only"),
  wrap  = c(15,         16,        16,       15,       16,        15,
            18,         17,        16,       17),
  stringsAsFactors = FALSE)

## a cylinder loses its top band to the cap, so its text sits lower
nodes$ylab <- ifelse(nodes$shape == SH_FILE, 0.14, 0.22)
nodes$ysub <- ifelse(nodes$shape == SH_FILE, -0.22, -0.20)

## the one accented node, named once so nothing else hard-codes it
ACCENT_ID    <- "raw"
nodes$fill   <- ifelse(nodes$id == ACCENT_ID, ACC_FILL,   NODE_FILL)
nodes$border <- ifelse(nodes$id == ACCENT_ID, ACC_BORDER, NODE_BORDER)

nodes$label <- wrap_txt(nodes$label, nodes$wrap)
nodes$sub   <- wrap_txt(nodes$sub,   nodes$wrap)

## --- stage regions ---------------------------------------------------------
regions <- data.frame(
  xmin = c(1, 1), xmax = c(99, 99), ymin = c(57, 1), ymax = c(100, 54),
  title = c("Pass 1 \u2014 variant calling with snippy",
            "Pass 2 \u2014 matrix and distance"),
  stringsAsFactors = FALSE)

## --- edges -----------------------------------------------------------------
## H = horizontal, V = vertical, L = down then right into the target's left
## edge, E = down, across a lane, then down into the target's top edge.
## `lane` is the y of an E cross-bar; `tx` overrides the entry x.
edges <- data.frame(
  from = c("assembly", "snippy", "snippy",  "vcf",   "pass1", "aligned",
           "union",    "pass2",  "gmatrix", "raw"),
  to   = c("snippy",   "vcf",    "aligned", "pass1", "union", "pass2",
           "pass2",    "gmatrix", "raw",     "scaled"),
  kind = c("H",        "H",      "L",       "H",     "V",     "V",
           "E",        "H",      "V",       "V"),
  lane = c(NA,         NA,       NA,        NA,      NA,      NA,
           55.5,       NA,       NA,        NA),
  tx   = c(NA,         NA,       NA,        NA,      NA,      62,
           68,         NA,       NA,        NA),
  stringsAsFactors = FALSE)

edge_path <- function(e) {
  s <- nodes[nodes$id == e$from, ]; t <- nodes[nodes$id == e$to, ]
  if (e$kind == "H") {
    ## force a single y so a 'horizontal' edge is never a shallow diagonal
    yh <- s$y
    data.frame(x = c(s$x + s$w / 2, t$x - t$w / 2), y = c(yh, yh))
  } else if (e$kind == "V") {
    xv <- if (is.na(e$tx)) s$x else e$tx
    data.frame(x = c(xv, xv), y = c(s$y - s$h / 2, t$y + t$h / 2))
  } else if (e$kind == "L") {
    xv <- if (is.na(e$tx)) s$x else e$tx
    data.frame(x = c(xv, xv, t$x - t$w / 2), y = c(s$y - s$h / 2, t$y, t$y))
  } else {
    xe <- if (is.na(e$tx)) t$x else e$tx
    data.frame(x = c(s$x, s$x, xe, xe),
               y = c(s$y - s$h / 2, e$lane, e$lane, t$y + t$h / 2))
  }
}
epath <- do.call(rbind, lapply(seq_len(nrow(edges)), function(i) {
  d <- edge_path(edges[i, ]); d$edge_id <- i; d
}))

## --- shape polygons --------------------------------------------------------
cyl_pts <- function(cx, cy, w, h, ry, lower_only = FALSE, n = 60) {
  a <- seq(0, 2 * pi, length.out = n)
  if (lower_only) a <- a[a >= pi]
  cyc <- if (lower_only) cy - h / 2 + ry else cy + h / 2 - ry
  data.frame(x = cx + (w / 2) * cos(a), y = cyc + ry * sin(a))
}
trap_pts <- function(cx, cy, w, h, slant = SLANT_F * h) {
  data.frame(x = c(cx - w / 2 + slant, cx + w / 2 - slant, cx + w / 2, cx - w / 2),
             y = c(cy + h / 2, cy + h / 2, cy - h / 2, cy - h / 2))
}
## carries fill/border through so the accent survives into the polygon layer
shape_polys <- function(nd) {
  out <- list(); k <- 0
  for (i in seq_len(nrow(nd))) {
    r <- nd[i, ]
    p <- if (r$shape == SH_STEP) trap_pts(r$x, r$y, r$w, r$h)
    else if (r$shape == SH_FILE) cyl_pts(r$x, r$y, r$w, r$h, RY_F * r$h)
    else NULL
    if (is.null(p)) next
    k <- k + 1
    p$grp <- k; p$fill <- r$fill; p$border <- r$border
    out[[length(out) + 1]] <- p
  }
  if (!length(out)) return(NULL)
  do.call(rbind, out)
}

polys    <- shape_polys(nodes)
boxes    <- nodes[nodes$shape == SH_DATA, ]
cyl_body <- nodes[nodes$shape == SH_FILE, ]; cyl_body$ry <- RY_F * cyl_body$h
cyl_low  <- do.call(rbind, lapply(seq_len(nrow(cyl_body)), function(i) {
  r <- cyl_body[i, ]
  p <- cyl_pts(r$x, r$y, r$w, r$h, r$ry, lower_only = TRUE)
  p$grp <- i; p$border <- r$border; p
}))

## --- key (shape legend), compact, in the dead lower-left of region 2 -------
KEY <- list(xmin = 4, xmax = 48, ymin = 3, ymax = 31)
key_nodes <- data.frame(
  id    = c("k_step", "k_file", "k_data"),
  x     = 12,
  y     = c(23, 15.5, 8),
  w     = 14,
  h     = 5.5,
  shape = c(SH_STEP, SH_FILE, SH_DATA),
  label = c("step", "snippy output file", "other data"),
  fill   = NODE_FILL,
  border = NODE_BORDER,
  stringsAsFactors = FALSE)
key_polys <- shape_polys(key_nodes)
key_box   <- key_nodes[key_nodes$shape == SH_DATA, ]
key_cyl   <- key_nodes[key_nodes$shape == SH_FILE, ]; key_cyl$ry <- RY_F * key_cyl$h
key_low   <- cyl_pts(key_cyl$x, key_cyl$y, key_cyl$w, key_cyl$h, key_cyl$ry,
                     lower_only = TRUE)

## ============================================================ panel A =====
pa <- ggplot() +
  ## stage regions
  geom_rect(data = regions, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = REG_FILL, colour = REG_BORDER, linewidth = 0.4) +
  geom_text(data = regions, aes(x = xmin + 2, y = ymax - 3.5, label = title),
            hjust = 0, size = TITLE_SIZE, fontface = "bold", colour = "grey15") +
  ## arrows, behind the shapes
  geom_path(data = epath, aes(x, y, group = edge_id), colour = "grey30",
            linewidth = 0.5, lineend = "round", linejoin = "round",
            arrow = arrow(length = unit(1.4, "mm"), type = "closed")) +
  ## cylinder bodies and sides, before their caps
  geom_rect(data = cyl_body, aes(xmin = x - w / 2, xmax = x + w / 2,
                                 ymin = y - h / 2 + ry, ymax = y + h / 2 - ry,
                                 fill = fill), colour = NA) +
  geom_segment(data = cyl_body, aes(x = x - w / 2, xend = x - w / 2,
                                    y = y - h / 2 + ry, yend = y + h / 2 - ry,
                                    colour = border), linewidth = 0.5) +
  geom_segment(data = cyl_body, aes(x = x + w / 2, xend = x + w / 2,
                                    y = y - h / 2 + ry, yend = y + h / 2 - ry,
                                    colour = border), linewidth = 0.5) +
  geom_path(data = cyl_low, aes(x, y, group = grp, colour = border), linewidth = 0.5) +
  ## trapezoids and cylinder caps
  geom_polygon(data = polys, aes(x, y, group = grp, fill = fill, colour = border),
               linewidth = 0.5) +
  ## plain data rectangles
  geom_rect(data = boxes, aes(xmin = x - w / 2, xmax = x + w / 2,
                              ymin = y - h / 2, ymax = y + h / 2,
                              fill = fill, colour = border), linewidth = 0.5) +
  ## node text
  geom_text(data = nodes, aes(x, y + h * ylab, label = label),
            size = LAB_SIZE, fontface = "bold", colour = "grey10", lineheight = 0.92) +
  geom_text(data = nodes[nzchar(nodes$sub), ], aes(x, y + h * ysub, label = sub),
            size = SUB_SIZE, colour = "grey30", lineheight = 1.05, vjust = 0.5) +
  ## --- key ---------------------------------------------------------------
geom_rect(data = as.data.frame(KEY),
          aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
          fill = "white", colour = "grey60", linewidth = 0.4) +
  geom_text(data = data.frame(x = KEY$xmin + 2, y = KEY$ymax - 2.5, label = "Key"),
            aes(x, y, label = label), hjust = 0,
            size = TITLE_SIZE, fontface = "bold", colour = "grey15") +
  geom_rect(data = key_cyl, aes(xmin = x - w / 2, xmax = x + w / 2,
                                ymin = y - h / 2 + ry, ymax = y + h / 2 - ry),
            fill = NODE_FILL, colour = NA) +
  geom_segment(data = key_cyl, aes(x = x - w / 2, xend = x - w / 2,
                                   y = y - h / 2 + ry, yend = y + h / 2 - ry),
               colour = NODE_BORDER, linewidth = 0.5) +
  geom_segment(data = key_cyl, aes(x = x + w / 2, xend = x + w / 2,
                                   y = y - h / 2 + ry, yend = y + h / 2 - ry),
               colour = NODE_BORDER, linewidth = 0.5) +
  geom_path(data = key_low, aes(x, y), colour = NODE_BORDER, linewidth = 0.5) +
  geom_polygon(data = key_polys, aes(x, y, group = grp), fill = NODE_FILL,
               colour = NODE_BORDER, linewidth = 0.5) +
  geom_rect(data = key_box, aes(xmin = x - w / 2, xmax = x + w / 2,
                                ymin = y - h / 2, ymax = y + h / 2),
            fill = NODE_FILL, colour = NODE_BORDER, linewidth = 0.5) +
  geom_text(data = key_nodes, aes(x + w / 2 + 2.5, y, label = label), hjust = 0,
            size = SUB_SIZE, colour = "grey25") +
  scale_fill_identity() + scale_colour_identity() +
  coord_cartesian(xlim = c(0, 100), ylim = c(0, 100), expand = FALSE, clip = "off") +
  labs(x = NULL, y = NULL) + theme_paper() +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        panel.grid = element_blank(), panel.border = element_blank(),
        plot.margin = margin(2, 2, 2, 2))

## ============================================================ panel B =====
sites <- read_fig("fig_m1_sites.csv")
## genome A and genome B are read off as two parallel vectors below, so the
## row order has to be made explicit rather than trusted
sites <- sites[order(as.character(sites$genome), sites$site_index), ]
N <- max(sites$site_index)
stopifnot(identical(sort(unique(as.character(sites$genome))), c("A", "B")),
          nrow(sites) == 2L * N,
          !anyDuplicated(paste(sites$genome, sites$site_index)))

## pitch is derived from N so the row always fills the panel
X0 <- 7; XR <- 99
PITCH <- (XR - X0) / N; W <- PITCH * 0.86
sites$xmin <- X0 + (sites$site_index - 1) * PITCH
sites$xmax <- sites$xmin + W
sites$cx   <- (sites$xmin + sites$xmax) / 2
YY <- c(A = 78, B = 62); H <- 10
sites$ymin <- as.numeric(YY[as.character(sites$genome)])
sites$ymax <- sites$ymin + H

SFILL <- c(called_match = scales::alpha(OI[["sky"]], 0.30),
           called_diff  = OI[["vermillion"]],
           uncalled     = "white")
SLAB <- c(called_match = "called, matches reference",
          called_diff  = "called, differs from ref.",
          uncalled     = "uncalled (N, -, or masked)")
stopifnot(all(as.character(sites$state) %in% names(SFILL)))
sites$state <- factor(sites$state, levels = names(SFILL), labels = SLAB)

a <- as.character(sites$state[sites$genome == "A"])
b <- as.character(sites$state[sites$genome == "B"])
both_called <- a != SLAB[["uncalled"]] & b != SLAB[["uncalled"]]

## A three-state cartoon cannot say whether two genomes that BOTH differ from
## the reference differ from each other -- that depends on the allele. So the
## illustration must not contain such a site; if it does, say so loudly rather
## than silently counting it as a difference.
both_diff <- sum(a == SLAB[["called_diff"]] & b == SLAB[["called_diff"]])
if (both_diff > 0)
  stop(sprintf(paste0("fig_m1_sites.csv has %d site(s) where BOTH genomes are 'called_diff'. ",
                      "With three states the figure cannot say whether they differ from each ",
                      "other (same alternative allele -> no pairwise difference). Set one of ",
                      "the two to 'called_match', or add an allele column and redefine ",
                      "`differs` on alleles."), both_diff))
differs <- both_called & (a == SLAB[["called_diff"]] | b == SLAB[["called_diff"]])
C_n <- sum(both_called)
d_n <- sum(differs)

wide <- data.frame(
  site_index = as.integer(sites$site_index[sites$genome == "A"]),
  enters = ifelse(!both_called, "\u2013", ifelse(differs, "C+d", "C")),
  stringsAsFactors = FALSE)
wide$cx <- X0 + (wide$site_index - 1) * PITCH + W / 2

## the deletion bracket assumes one contiguous uncalled run in genome B
gap <- sites[sites$genome == "B" & sites$state == SLAB[["uncalled"]], ]
gap_n <- nrow(gap)
gap_idx <- sort(as.integer(gap$site_index))
stopifnot(gap_n > 0, identical(gap_idx, seq(min(gap_idx), max(gap_idx))))
gap_lo <- min(gap$xmin); gap_hi <- max(gap$xmax)
brk <- data.frame(x = c(gap_lo, gap_lo, gap_hi, gap_hi), y = c(57.5, 54.5, 54.5, 57.5))
gap_lab <- data.frame(x = (gap_lo + gap_hi) / 2, y = 50,
                      label = sprintf("deletion in genome B (%d sites)", gap_n))

## hatch EVERY uncalled cell, in either genome, so the hatch means one thing
hatch <- sites[sites$state == SLAB[["uncalled"]], ]
row_box <- data.frame(xmin = X0, xmax = X0 + N * PITCH - (PITCH - W),
                      ymin = 11, ymax = 34)
tally <- data.frame(
  x = X0 + 3, y = c(29, 23.5, 17),
  label = c(sprintf("jointly called sites C = %d", C_n),
            sprintf("differing sites d = %d", d_n),
            sprintf("raw distance = %d", d_n)),
  face = c("plain", "plain", "bold"),
  stringsAsFactors = FALSE)

pb <- ggplot() +
  geom_rect(data = sites, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                              fill = state), colour = "grey35", linewidth = 0.5) +
  geom_segment(data = hatch, aes(x = xmin + 0.6, xend = xmax - 0.6,
                                 y = ymin + 0.6, yend = ymax - 0.6),
               colour = "grey65", linewidth = 0.5) +
  geom_path(data = brk, aes(x, y), colour = "grey55", linewidth = 0.5) +
  geom_text(data = gap_lab, aes(x, y, label = label), size = GEOM_TEXT_MIN,
            colour = "grey35") +
  geom_text(data = wide, aes(cx, 40, label = enters), size = GEOM_TEXT_MIN,
            colour = "grey25") +
  geom_rect(data = row_box, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = NA, colour = "grey35", linewidth = 0.5) +
  geom_text(data = tally, aes(x, y, label = label), fontface = tally$face,
            size = GEOM_TEXT_MIN + 0.3, colour = "grey10", hjust = 0) +
  ## breaks/labels unnamed: a named `breaks` vector makes ggplot echo the names
  ## (the internal tokens) instead of the descriptive labels
  scale_fill_manual(values = setNames(unname(SFILL), SLAB), name = NULL,
                    breaks = unname(SLAB), labels = unname(SLAB)) +
  scale_y_continuous(breaks = c(83, 67, 40),
                     labels = c("genome A", "genome B", "enters"),
                     limits = c(8, 95), expand = expansion(0)) +
  scale_x_continuous(limits = c(0, 100), expand = expansion(0)) +
  labs(x = NULL, y = NULL) +
  guides(fill = guide_legend(nrow = 3, byrow = TRUE,
                             override.aes = list(colour = "grey35"))) +
  theme_paper() +
  theme(axis.text.x = element_blank(), axis.ticks = element_blank(),
        axis.text.y = element_text(size = BASE_SIZE - 1.5, colour = "grey25"),
        panel.grid = element_blank(), panel.border = element_blank(),
        ## the legend lives in the empty right half of the tally box
        ## (box spans x 7-98, y 11-34; the three counts occupy x 10-45)
        legend.position = c(0.70, 0.17), legend.justification = "center",
        legend.background = element_rect(fill = "white", colour = "grey85",
                                         linewidth = 0.3),
        legend.margin = margin(1.5, 2.5, 1.5, 2.5),
        legend.key.size = unit(0.26, "cm"),
        legend.text = element_text(size = BASE_SIZE - 2),
        plot.margin = margin(2, 2, 2, 2))

## ================================================= checks, THEN write =====
## Assertions run before save_fig so a failure cannot leave a wrong PDF behind.
uncalled_either <- sum(a == SLAB[["uncalled"]] | b == SLAB[["uncalled"]])
cat("\nKEY NUMBERS (counted from fig_m1_sites.csv)\n")
cat(sprintf("  sites drawn             : %d\n", N))
cat(sprintf("  jointly called C        : %d\n", C_n))
cat(sprintf("  differing d             : %d\n", d_n))
cat(sprintf("  uncalled in genome B    : %d\n", gap_n))
cat(sprintf("  uncalled in A or B      : %d\n", uncalled_either))
cat(sprintf("  raw distance            : %d\n", d_n))
## general identity: a site is either jointly called or uncalled in at least one
## genome. C + gap_n == N only holds while genome A has no uncalled site.
stopifnot(C_n + uncalled_either == N, d_n <= C_n)
if (uncalled_either != gap_n)
  message("note: genome A also has uncalled sites; the bracket labels genome B only")

fig1 <- (pa | pb) + plot_layout(widths = c(1.32, 1)) + plot_annotation(tag_levels = "a")
save_fig(fig1, "fig1_pipeline_schematic", width_mm = 172, height_mm = 108)

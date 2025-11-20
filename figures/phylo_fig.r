# =========================
# Libraries
# =========================
library(ggplot2)
library(ape)
library(ggtree)
library(dplyr)
library(readr)
library(ggtreeExtra)
library(RColorBrewer)
library(ggnewscale)
library(scales)
library(cowplot)
library(stringr)
library(beeswarm)

# =========================
# 0) Load tree + metadata
# =========================
tree <- read.tree("../data/full_prune.tree")
Taxa <- read.csv("../data/full_phenodat.csv")

# =========================
# 1) Clean metadata
# =========================
Taxa <- Taxa %>%
  mutate(
    genus_species = gsub(" ", "_", genus_species),
    mating_system = tolower(trimws(mating_system)),
    Haploid.Chrom.. = suppressWarnings(as.numeric(Haploid.Chrom..))
  ) %>%
  filter(!is.na(genus_species))

# =========================
# 2) Match tree and metadata
# =========================
overlap <- intersect(tree$tip.label, Taxa$genus_species)

tree_pruned <- drop.tip(tree, setdiff(tree$tip.label, overlap))
meta_pruned <- Taxa %>% filter(genus_species %in% overlap)

tip_metadata <- tibble(label = tree_pruned$tip.label) %>%
  left_join(meta_pruned, by = c("label" = "genus_species")) %>%
  mutate(
    # ✅ Relabel mating system
    mating_system = recode_factor(
      mating_system,
      "sc" = "Self-compatible",
      "si" = "Self-incompatible"
    ),
    Haploid.Chrom.. = suppressWarnings(as.numeric(Haploid.Chrom..))
  )

cat("Mating system counts:\n")
print(table(tip_metadata$mating_system, useNA = "ifany"))

# =========================
# 3) Add normalized chromosome number
# =========================
tip_metadata <- tip_metadata %>%
  mutate(chrom_norm = Haploid.Chrom.. / max(Haploid.Chrom.., na.rm = TRUE))

# =========================
# 4) Group tree by family
# =========================
grouped_tree <- ladderize(tree_pruned, right = TRUE)
fam2tips <- split(as.character(tip_metadata$label),
                  as.character(tip_metadata$Family))
grouped_tree_col <- ggtree::groupOTU(grouped_tree, fam2tips,
                                     group_name = "branch_family")

families    <- sort(unique(as.character(tip_metadata$Family)))
family_cols <- setNames(hue_pal(l = 70, c = 60)(length(families)), families)

# =========================
# 5) Plot setup
# =========================
# ✅ Updated palette with full names
mating_cols <- c(
  "Self-compatible" = "#2c7fb8",
  "Self-incompatible" = "#99d8c9"
)

rm_axis  <- list(axis = "x", text = NA, line = NA, ticks = NA)

# Ring geometry
off_ms <- 0.02   # distance from tree
pw_ms  <- 0.20   # mating system band thickness

off_bar <- off_ms   # bars close to mating system
pw_bar  <- 0.45     # taller chromosome bars

pad_out <- off_bar + pw_bar + 0.05

tight_tree_theme <- ggtree::theme_tree() +
  theme(plot.margin = margin(2,2,2,2,"pt"),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank())

radius <- 3

# =========================
# 6) Base tree
# =========================
p0 <- (ggtree(grouped_tree_col, layout = "circular",
              branch.length = "none", open.angle = 0,
              size = 0) %<+% tip_metadata) +
  xlim_tree(radius) +
  tight_tree_theme +
  scale_x_continuous(expand = expansion(mult = 0, add = pad_out)) +
  scale_y_continuous(expand = expansion(mult = 0, add = 0))

p <- p0 +
  geom_tree(
    aes(color = branch_family),
    linewidth = 0.12,
    lineend = "butt",
    linejoin = "bevel",
    show.legend = TRUE
  ) +
  scale_color_manual(
    values = family_cols,
    na.translate = FALSE,
    name = "Family",
    guide = guide_legend(
      override.aes = list(linewidth = 3),
      keywidth = unit(28, "pt"),
      keyheight = unit(12, "pt")
    )
  ) +
  theme(
    legend.key.width  = unit(28, "pt"),
    legend.key.height = unit(12, "pt")
  )

# =========================
# 7) Rings
# =========================
# Ring 1: mating system band
p1 <- p + 
  ggnewscale::new_scale_fill() +
  geom_fruit(
    geom = geom_col,
    aes(y = label, x = 2, fill = mating_system),
    orientation = "y",
    pwidth = pw_ms, offset = off_ms,
    color = NA, linewidth = 0,
    axis.params = rm_axis
  ) +
  scale_fill_manual(
    values = mating_cols,
    na.translate = FALSE,
    name = "Mating System",
    guide = guide_legend(
      override.aes = list(linewidth = 3),
      keywidth = unit(28, "pt"),
      keyheight = unit(12, "pt")
    )
  ) +
  theme(
    legend.key.height = unit(12, "pt"),
    legend.key.width  = unit(28, "pt"),
    legend.text       = element_text(size = 12),
    legend.title      = element_text(size = 14, face = "bold"),
    legend.key        = element_rect(fill = "white", colour = NA)
  )

# Ring 2: chromosome numbers (normalized)
p_final <- p1 + ggnewscale::new_scale_fill() +
  geom_fruit(
    geom = geom_col,
    aes(y = label, x = chrom_norm),
    orientation = "y",
    pwidth = pw_bar, offset = off_bar,
    fill = "gray40", color = NA,
    axis.params = rm_axis
  )

# =========================
# 8) Export
# =========================
ggsave("tree_final.pdf", p_final, width = 8, height = 8, bg = "white")
svglite::svglite("tree_final.svg", width = 8, height = 8, bg = "white")
print(p_final); dev.off()

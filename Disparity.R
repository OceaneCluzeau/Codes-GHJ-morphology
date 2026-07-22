# ============================================================
# MORPHOLOGICAL DISPARITY ANALYSIS BY LOCOMOTOR GROUP
# ============================================================

# ============================================================
# DISPARITY ANALYSIS USING YOUR CSV STRUCTURE
# ============================================================

library(geomorph)

# ------------------------------------------------------------
# 1. Read metadata CSV
# ------------------------------------------------------------
meta <- read.csv(file = "./Datas/id_loco_humerus.csv", stringsAsFactors = FALSE)
meta <- read.csv(file = "./Datas/id_loco_scapula.csv", stringsAsFactors = FALSE)

# Convert grouping variables to factors
meta$species    <- factor(meta$species)
meta$family     <- factor(meta$family)
meta$locomotion <- factor(meta$locomotion)

# ------------------------------------------------------------
# 2. Load SLIDED landmark object
# ------------------------------------------------------------
load("./Datas/landmarks_slided/landmarks_humerus_surf_slided.Rdata")
load("./Datas/landmarks_slided/landmarks_scapula_surf_slided.Rdata")
# Extract slid landmarks
coords <- landmarks_humerus_slided$dataslide
coords <- landmarks_scapula_slided$dataslide
# ------------------------------------------------------------
# 3. Match landmark order to metadata (CRITICAL)
# ------------------------------------------------------------
# Check specimen names stored in the landmark array
dimnames(coords)[[3]]

# ------------------------------------------------------------
# 4. Generalized Procrustes Analysis
# ------------------------------------------------------------
gpa <- gpagen(coords, print.progress = FALSE)

# ------------------------------------------------------------
# 5. Morphological disparity by locomotor group
# ------------------------------------------------------------
group <- meta$locomotion

disparity <- morphol.disparity(
  gpa$coords ~ group,
  groups = group,
  iter = 999,
  boot = TRUE
)

summary(disparity)

raw_disp <- morphol.disparity(
  gpa$coords ~ group,
  groups = group,
  iter = 999,
  boot = FALSE
)
cbind(
  Raw = raw_disp$Procrustes.var,
  Bootstrapped = disparity$Procrustes.var
)

# Extract pairwise p-values
pvals <- disparity$PV.dist.Pval

# Correction FDR
pvals_fdr <- p.adjust(pvals[lower.tri(pvals)], method = "fdr")

# Reconstruire la matrice corrigée
pvals_fdr_mat <- matrix(NA, nrow = 5, ncol = 5,
                        dimnames = dimnames(pvals))
pvals_fdr_mat[lower.tri(pvals_fdr_mat)] <- pvals_fdr
pvals_fdr_mat[upper.tri(pvals_fdr_mat)] <- t(pvals_fdr_mat)[upper.tri(pvals_fdr_mat)]

# Vérifier
pvals_fdr_mat

library(reshape2)
library(ggplot2)
# Heatmap mise à jour
pvals_long <- melt(pvals_fdr_mat)
cols <- colorRampPalette(c("blue", "red"))

ggplot(pvals_long, aes(Var2, Var1, fill = value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.3f", value)), size = 4) +
  scale_fill_gradientn(
    colours = cols(100),   # generate colors
    name = "P-value"
  ) +
  labs(
    title = "Pairwise Morphological Disparity (P-values)",
    x = "",
    y = ""
  ) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 12),
    plot.title = element_text(hjust = 0.5)
  )
# ------------------------------------------------------------
# 6. Plot disparity
# ------------------------------------------------------------
# Procrustes-aligned coordinates
Y <- gpa$coords

groups_levels <- levels(group)

group_means <- lapply(groups_levels, function(g) {
  mshape(Y[, , group == g])
})

names(group_means) <- groups_levels

# Procrustes distance to group mean
dist_to_mean <- sapply(1:dim(Y)[3], function(i) {
  sqrt(sum(
    (Y[, , i] - group_means[[ as.character(group[i]) ]])^2
  ))
})

# Boxplot
boxplot(dist_to_mean ~ group,
        ylab = "Procrustes distance to group mean",
        xlab = "Locomotor group",
        las = 2,
        #col = "darkorange")# scapula
        col= "steelblue") # humerus

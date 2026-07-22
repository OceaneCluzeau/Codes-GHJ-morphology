# ── Scapula shape visualisation on PC axes ─────────────────────────────────
# Uses: pPCA_scap (gm.prcomp on allometry-corrected species means)
#       scap_corrected_array (species mean landmark array, p x k x n)
#       mesh_scap  : your scapula .ply mesh
#       landmark_scap : reference landmarks on the mesh (p x k matrix)

library(Morpho)
library(rgl)
library(geomorph)
library(RANN)
library(Rvcg)
library(ape)

phy <- read.nexus(file = "./2_Block_analysis/data/mcctree.nexus")

load(file = "./outputs/Yhum.Rdata")
load(file = "./outputs/Yscap.Rdata")
load(file = "./outputs/Y_sb.Rdata")

load(file = "./outputs/data_all_humerus.Rdata")
load(file = "./outputs/data_all_scapula.Rdata")
load(file = "./outputs/data_all_sb_nocsize.Rdata")

###############Humerus######################
PCA_hum <- gm.prcomp(Yhum, phy = phy)

# ── 1. Extract PC scores ────────────────────────────────────────────────────
pc1_scores <- data.frame(species = rownames(PCA_hum$x),
                         scores  = PCA_hum$x[, 1])
pc2_scores <- data.frame(species = rownames(PCA_hum$x),
                         scores  = PCA_hum$x[, 2])

# ── 2. Identify species at extremes ────────────────────────────────────────
sp_pc1_min <- pc1_scores$species[which.min(pc1_scores$scores)]
sp_pc1_max <- pc1_scores$species[which.max(pc1_scores$scores)]
sp_pc2_min <- pc2_scores$species[which.min(pc2_scores$scores)]
sp_pc2_max <- pc2_scores$species[which.max(pc2_scores$scores)]

cat("PC1 min:", sp_pc1_min, "\n")
cat("PC1 max:", sp_pc1_max, "\n")
cat("PC2 min:", sp_pc2_min, "\n")
cat("PC2 max:", sp_pc2_max, "\n")

# ── 3. Extract corresponding species mean landmark configurations ───────────
lma_pc1_min <- Yhum[,, sp_pc1_min]
lma_pc1_max <- Yhum[,, sp_pc1_max]
lma_pc2_min <- Yhum[,, sp_pc2_min]
lma_pc2_max <- Yhum[,, sp_pc2_max]

# ── Load template landmark and mesh ────────────────────────────────
template_lm <- as.matrix(read.table("./2_Block_analysis/data_blocks/Humerus/all_pts/Nasalis_larvatus_m_FMNH_68684_humerus_L.pts"))
mesh_template <- vcgImport("./Datas/Resampling/all_pts/all_ply/Nasalis_larvatus_m_FMNH_68684_humerus_L.ply")

# ── Sanity check ───────────────────────────────────────────────────
stopifnot(nrow(template_lm) == dim(Yhum)[1])   # p landmarks must match
stopifnot(ncol(template_lm) == dim(Yhum)[2])   # k dimensions must match

# ── Then use exactly as before ─────────────────────────────────────
mean_shape  <- mshape(Yhum)
warp_mean   <- tps3d(mesh_template, template_lm, mean_shape)

open3d()
shade3d(warp_mean, col = "#B3B3B3", alpha = 1, specular = 1)


# ── 5. Set and save viewing angles interactively ────────────────────────────
# Rotate to desired view in the rgl window, then run:
lateral_view <- list(zoom       = par3d()$zoom,
                     userMatrix = par3d()$userMatrix,
                     windowRect = par3d()$windowRect)

dorsal_view  <- list(zoom       = par3d()$zoom,
                     userMatrix = par3d()$userMatrix,
                     windowRect = par3d()$windowRect)

ventral_view <- list(zoom       = par3d()$zoom,
                     userMatrix = par3d()$userMatrix,
                     windowRect = par3d()$windowRect)

shapes_views <- list(lateral  = lateral_view,
                     dorsal   = dorsal_view,
                     ventral = ventral_view)

save(shapes_views, file = "./outputs/shapes_views_humerus.Rdata")  # .Rdata not .R
load(file = "./outputs/shapes_views_humerus.Rdata")              # reload if needed
mesh_hum <- mesh_template
landmark_hum <- template_lm
# ── 6. Helper function to overlay two warped meshes ────────────────────────
plot_shape_comparison <- function(lma_plus, lma_minus, view, snapshot_file) {
  open3d(zoom       = view$zoom,
         userMatrix = view$userMatrix,
         windowRect = view$windowRect)
  # Positive extreme — red
  warp_plus <- tps3d(mesh_hum, landmark_hum, lma_plus)
  shade3d(warp_plus,  col = "#fa0000", alpha = 0.5, specular = 1)
  # Negative extreme — blue
  warp_minus <- tps3d(mesh_hum, landmark_hum, lma_minus)
  shade3d(warp_minus, col = "#0022fa", alpha = 0.5, specular = 1)
  rgl.snapshot(snapshot_file, fmt = "png")
}

# ── 7. PC1 visualisation ────────────────────────────────────────────────────
plot_shape_comparison(
  lma_plus      = lma_pc1_max,
  lma_minus     = lma_pc1_min,
  view          = shapes_views$lateral,
  snapshot_file = "humerus_shapes_pc1_lateral.png"
)
plot_shape_comparison(
  lma_plus      = lma_pc1_max,
  lma_minus     = lma_pc1_min,
  view          = shapes_views$dorsal,
  snapshot_file = "humerus_shapes_pc1_dorsal.png"
)
plot_shape_comparison(
  lma_plus      = lma_pc1_max,
  lma_minus     = lma_pc1_min,
  view          = shapes_views$ventral,
  snapshot_file = "humerus_shapes_pc1_ventral.png"
)

# ── 8. PC2 visualisation ────────────────────────────────────────────────────
plot_shape_comparison(
  lma_plus      = lma_pc2_max,
  lma_minus     = lma_pc2_min,
  view          = shapes_views$lateral,
  snapshot_file = "humerus_shapes_pc2_lateral.png"
)
plot_shape_comparison(
  lma_plus      = lma_pc2_max,
  lma_minus     = lma_pc2_min,
  view          = shapes_views$dorsal,
  snapshot_file = "humerus_shapes_pc2_dorsal.png"
)
plot_shape_comparison(
  lma_plus      = lma_pc2_max,
  lma_minus     = lma_pc2_min,
  view          = shapes_views$ventral,
  snapshot_file = "humerus_shapes_pc2_ventral.png"
)

# ── one extreme vs template (mean shape) ───────────────────────────
plot_shape_vs_template <- function(lma_extreme, color, view, snapshot_file) {
  open3d(zoom       = view$zoom,
         userMatrix = view$userMatrix,
         windowRect = c(0, 0, 2000, 2000))  # bigger = higher-res PNG
  # Template (mean shape) — grey, solid
  warp_mean <- tps3d(mesh_hum, landmark_hum, mean_shape)
  shade3d(warp_mean, col = "#B3B3B3", alpha = 0.4, specular = 1)
  
  # Extreme shape — species colour, semi-transparent
  warp_extreme <- tps3d(mesh_hum, landmark_hum, lma_extreme)
  shade3d(warp_extreme, col = color, alpha = 0.5, specular = 1)
  
  rgl.snapshot(snapshot_file, fmt = "png")
}

# ── Define one colour per extreme species ──────────────────────────────────
col_pc1_min <- "#9895DB"   # purple = QA
col_pc1_max <- "#93D175"   # green = S
col_pc2_min <- "#93D175"   # green = s
col_pc2_max <- "#ED3BAC"   # pink = QT

# ── Plot all 4 extremes vs template ────────────────────────────────────────
##PC1-
plot_shape_vs_template(lma_pc1_min, col_pc1_min,
                       shapes_views$lateral,
                       "humerus_pc1_min_vs_mean.png")
plot_shape_vs_template(lma_pc1_min, col_pc1_min,
                       shapes_views$dorsal,
                       "humerus_pc1_min_vs_mean2.png")
##PC1+
plot_shape_vs_template(lma_pc1_max, col_pc1_max,
                       shapes_views$lateral,
                       "humerus_pc1_max_vs_mean.png")
plot_shape_vs_template(lma_pc1_max, col_pc1_max,
                       shapes_views$dorsal,
                       "humerus_pc1_max_vs_mean2.png")
##PC2-
plot_shape_vs_template(lma_pc2_min, col_pc2_min,
                       shapes_views$lateral,
                       "humerus_pc2_min_vs_mean.png")
plot_shape_vs_template(lma_pc2_min, col_pc2_min,
                       shapes_views$dorsal,
                       "humerus_pc2_min_vs_mean2.png")
##PC2+
plot_shape_vs_template(lma_pc2_max, col_pc2_max,
                       shapes_views$lateral,
                       "humerus_pc2_max_vs_mean.png")
plot_shape_vs_template(lma_pc2_max, col_pc2_max,
                       shapes_views$dorsal,
                       "humerus_pc2_max_vs_mean2.png")

###############Scapula######################
pPCA_scap <- gm.prcomp(scap_corrected_array, phy = phy)

# ── 1. Extract PC scores ────────────────────────────────────────────────────
pc1_scores <- data.frame(species = rownames(pPCA_scap$x),
                         scores  = pPCA_scap$x[, 1])
pc2_scores <- data.frame(species = rownames(pPCA_scap$x),
                         scores  = pPCA_scap$x[, 2])

# ── 2. Identify species at extremes ────────────────────────────────────────
sp_pc1_min <- pc1_scores$species[which.min(pc1_scores$scores)]
sp_pc1_max <- pc1_scores$species[which.max(pc1_scores$scores)]
sp_pc2_min <- pc2_scores$species[which.min(pc2_scores$scores)]
sp_pc2_max <- pc2_scores$species[which.max(pc2_scores$scores)]

cat("PC1 min:", sp_pc1_min, "\n")
cat("PC1 max:", sp_pc1_max, "\n")
cat("PC2 min:", sp_pc2_min, "\n")
cat("PC2 max:", sp_pc2_max, "\n")

# ── 3. Extract corresponding species mean landmark configurations ───────────
lma_pc1_min <- scap_corrected_array[,, sp_pc1_min]
lma_pc1_max <- scap_corrected_array[,, sp_pc1_max]
lma_pc2_min <- scap_corrected_array[,, sp_pc2_min]
lma_pc2_max <- scap_corrected_array[,, sp_pc2_max]

# ── Load template landmark and mesh ────────────────────────────────
template_lm <- as.matrix(read.table("./2_Block_analysis/data_blocks/Scapula/all_pts/Colobus_guereza_f_AMNH_52241_scapula_L.pts"))
mesh_template <- vcgImport("./Datas/Resampling/all_pts/all_ply/Colobus_guereza_f_AMNH_52241_scapula_L.ply")

# ── Sanity check ───────────────────────────────────────────────────
stopifnot(nrow(template_lm) == dim(Yscap)[1])   # p landmarks must match
stopifnot(ncol(template_lm) == dim(Yscap)[2])   # k dimensions must match

# ── Then use exactly as before ─────────────────────────────────────
mean_shape  <- mshape(Yscap)
warp_mean   <- tps3d(mesh_template, template_lm, mean_shape)

open3d()
shade3d(warp_mean, col = "#B3B3B3", alpha = 1, specular = 1)


# ── 5. Set and save viewing angles interactively ────────────────────────────
# Rotate to desired view in the rgl window, then run:
lateral_view <- list(zoom       = par3d()$zoom,
                     userMatrix = par3d()$userMatrix,
                     windowRect = par3d()$windowRect)

dorsal_view  <- list(zoom       = par3d()$zoom,
                     userMatrix = par3d()$userMatrix,
                     windowRect = par3d()$windowRect)

caudal_view <- list(zoom       = par3d()$zoom,
                     userMatrix = par3d()$userMatrix,
                     windowRect = par3d()$windowRect)
ventral_view <- list(zoom       = par3d()$zoom,
                    userMatrix = par3d()$userMatrix,
                    windowRect = par3d()$windowRect)

shapes_views <- list(lateral  = lateral_view,
                     dorsal   = dorsal_view,
                     caudal  = caudal_view,
                     ventral = ventral_view)

save(shapes_views, file = "./outputs/shapes_views_scapula.Rdata")  # .Rdata not .R
load(file = "./outputs/shapes_views_scapula.Rdata")              # reload if needed
mesh_scap <- mesh_template
landmark_scap <- template_lm
# ── 6. Helper function to overlay two warped meshes ────────────────────────
plot_shape_comparison <- function(lma_plus, lma_minus, view, snapshot_file) {
  open3d(zoom       = view$zoom,
         userMatrix = view$userMatrix,
         windowRect = view$windowRect)
  # Positive extreme — red
  warp_plus <- tps3d(mesh_scap, landmark_scap, lma_plus)
  shade3d(warp_plus,  col = "#fa0000", alpha = 0.5, specular = 1)
  # Negative extreme — blue
  warp_minus <- tps3d(mesh_scap, landmark_scap, lma_minus)
  shade3d(warp_minus, col = "#0022fa", alpha = 0.5, specular = 1)
  rgl.snapshot(snapshot_file, fmt = "png")
}

# ── 7. PC1 visualisation ────────────────────────────────────────────────────
plot_shape_comparison(
  lma_plus      = lma_pc1_max,
  lma_minus     = lma_pc1_min,
  view          = shapes_views$dorsal,
  snapshot_file = "scapula_shapes_pc1_dorsal.png"
)
plot_shape_comparison(
  lma_plus      = lma_pc1_max,
  lma_minus     = lma_pc1_min,
  view          = shapes_views$caudal,
  snapshot_file = "scapula_shapes_pc1_caudal.png"
)


# ── 8. PC2 visualisation ────────────────────────────────────────────────────
plot_shape_comparison(
  lma_plus      = lma_pc2_max,
  lma_minus     = lma_pc2_min,
  view          = shapes_views$dorsal,
  snapshot_file = "scapula_shapes_pc2_dorsal.png"
)
plot_shape_comparison(
  lma_plus      = lma_pc2_max,
  lma_minus     = lma_pc2_min,
  view          = shapes_views$caudal,
  snapshot_file = "scapula_shapes_pc2_caudal.png"
)

# ── Helper: one extreme vs template (mean shape) ───────────────────────────
plot_shape_vs_template <- function(lma_extreme, color, view, snapshot_file) {
  open3d(zoom       = view$zoom,
         userMatrix = view$userMatrix,
         windowRect = view$windowRect)
  
  # Template (mean shape) — grey, solid
  warp_mean <- tps3d(mesh_scap, landmark_scap, mean_shape)
  shade3d(warp_mean, col = "#B3B3B3", alpha = 0.4, specular = 1)
  
  # Extreme shape — species colour, semi-transparent
  warp_extreme <- tps3d(mesh_scap, landmark_scap, lma_extreme)
  shade3d(warp_extreme, col = color, alpha = 0.5, specular = 1)
  
  rgl.snapshot(snapshot_file, fmt = "png")
}

# ── Define one colour per extreme species ──────────────────────────────────
col_pc1_min <- "#9895DB"   # purple = QA
col_pc1_max <- "#ED3BAC"   # pink = QT
col_pc2_min <- "#93D175"   # green = s
col_pc2_max <- "#299157"   # dark green = B

# ── Plot all 4 extremes vs template ────────────────────────────────────────
##PC1-
plot_shape_vs_template(lma_pc1_min, col_pc1_min,
                       shapes_views$caudal,
                       "humerus_pc1_min_vs_mean.png")
plot_shape_vs_template(lma_pc1_min, col_pc1_min,
                       shapes_views$dorsal,
                       "humerus_pc1_min_vs_mean2.png")
##PC1+
plot_shape_vs_template(lma_pc1_max, col_pc1_max,
                       shapes_views$caudal,
                       "humerus_pc1_max_vs_mean.png")
plot_shape_vs_template(lma_pc1_max, col_pc1_max,
                       shapes_views$dorsal,
                       "humerus_pc1_max_vs_mean2.png")
##PC2-
plot_shape_vs_template(lma_pc2_min, col_pc2_min,
                       shapes_views$caudal,
                       "humerus_pc2_min_vs_mean.png")
plot_shape_vs_template(lma_pc2_min, col_pc2_min,
                       shapes_views$dorsal,
                       "humerus_pc2_min_vs_mean2.png")
##PC2+
plot_shape_vs_template(lma_pc2_max, col_pc2_max,
                       shapes_views$caudal,
                       "humerus_pc2_max_vs_mean.png")
plot_shape_vs_template(lma_pc2_max, col_pc2_max,
                       shapes_views$dorsal,
                       "humerus_pc2_max_vs_mean2.png")

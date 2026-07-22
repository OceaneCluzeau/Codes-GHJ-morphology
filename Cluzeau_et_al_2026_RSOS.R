### Library
library(devtools)
#install.packages("remotes")
#devtools::install_github("aharmer/morphoBlocks", build_vignettes = TRUE, dependencies = TRUE)
library(rgl)
library(morphoBlocks)
library(FactoMineR)
library(abind)
library(geomorph)
####################################################
#### Morphoblock code # adapted from Thomas et al. 2023
####################################################

for (d in c(
  "./2_Block_analysis_TEST/blocks/blocks_without_surf",
  "./2_Block_analysis_TEST/blocks/blocks_with_surf",
  "./2_Block_analysis_TEST/mean_data",
  "./2_Block_analysis_TEST/modelsevol",
  "./2_Block_analysis_TEST/data_blocks/Humerus/all_pts",
  "./2_Block_analysis_TEST/data_blocks/Scapula/all_pts",
  "./2_Block_analysis_TEST/data",
  "./outputs_TEST",
  "./results_TEST"
)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ── Shared functions (used by both pipelines) ─────────────────────────────────
remove_template <- function(block) {
  block@gpa.3D   <- block@gpa.3D[,,  -1]
  block@raw      <- block@raw[,,    -1]
  block@gpa.2D   <- block@gpa.2D[   -1, ]
  block@centroid <- block@centroid[  -1]
  block@n        <- dim(block@gpa.3D)[3]
  return(block)
}

rename_block <- function(block_obj, spec_names) {
  if (length(spec_names) != block_obj@n) stop("Specimen count mismatch.")
  dimnames(block_obj@gpa.3D)[[3]] <- spec_names
  dimnames(block_obj@raw)[[3]]    <- spec_names
  rownames(block_obj@gpa.2D)      <- spec_names
  return(block_obj)
}

rename_blocklist <- function(blocklist_obj, spec_names) {
  if (length(spec_names) != blocklist_obj@n[1]) stop("Specimen count mismatch.")
  for (i in seq_along(blocklist_obj@block.list)) {
    rownames(blocklist_obj@block.list[[i]]) <- spec_names
  }
  return(blocklist_obj)
}

# Create output directories
dir.create("./2_Block_analysis_TEST/blocks/blocks_without_surf", recursive = TRUE, showWarnings = FALSE)
dir.create("./2_Block_analysis_TEST/blocks/blocks_with_surf",    recursive = TRUE, showWarnings = FALSE)

####------------------------------------------------
#### WITHOUT surface semi landmarks
####------------------------------------------------
dirpath1 <- "./2_Block_analysis/new_pts_blocks/humerus_prep"
dirpath2 <- "./2_Block_analysis/new_pts_blocks/scapula_prep"

# Read data (no GPA yet)
humerus <- readPts(dirpath1, gpa = FALSE)
scapula <- readPts(dirpath2, gpa = FALSE)

# GPA + sliding of curve semilandmarks
block1 <- formatBlock(humerus@raw, curves = humerus@curves, k = 87, gpa = TRUE)
block2 <- formatBlock(scapula@raw, curves = scapula@curves, k = 87, gpa = TRUE)

# Remove duplicate template specimen
block1 <- remove_template(block1)
block2 <- remove_template(block2)

# Combine into blocklist
blocklist <- combineBlocks(blocks = c(block1, block2))

# Rename specimens
spec_data  <- read.csv("./2_Block_analysis/data/id_loco_humerus.csv", sep = ";")
block1     <- rename_block(block1, spec_data$filename)

spec_data  <- read.csv("./2_Block_analysis/data/id_loco_scapula.csv", sep = ",")
block2     <- rename_block(block2, spec_data$filename)

spec_data  <- read.csv("./2_Block_analysis/data/id_loco_general.csv", sep = ";")
blocklist  <- rename_blocklist(blocklist, spec_data$filename)

# Sanity checks
cat("Block1 specimen 1:", dimnames(block1@gpa.3D)[[3]][1], "\n")
cat("Block2 specimen 1:", dimnames(block2@gpa.3D)[[3]][1], "\n")
cat("Blocklist row 1:", rownames(blocklist@block.list$block_A)[1], "\n")

# ✅ Save WITHOUT surface
save(block1,    file = "./2_Block_analysis_TEST/blocks/blocks_without_surf/block1.Rdata")
save(block2,    file = "./2_Block_analysis_TEST/blocks/blocks_without_surf/block2.Rdata")
save(blocklist, file = "./2_Block_analysis_TEST/blocks/blocks_without_surf/blocklist.Rdata")
cat("✔ Saved blocks WITHOUT surface\n")

####------------------------------------------------
#### WITH surface semi landmarks
####------------------------------------------------

# Load already-slided surface coordinates
load(file = "./2_Block_analysis/data_blocks/Humerus/landmarks_humerus_surf_slided.Rdata")
load(file = "./2_Block_analysis/data_blocks/Scapula/landmarks_scapula_surf_slided.Rdata")

hum_coords  <- landmarks_humerus_slided$dataslide
scap_coords <- landmarks_scapula_slided$dataslide

# ── Write .pts files from dataslide arrays ────────────────────────────────────
write_pts_from_dataslide <- function(coords, outdir) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  n          <- dim(coords)[3]
  spec_names <- dimnames(coords)[[3]]
  if (is.null(spec_names)) spec_names <- paste0("specimen_", seq_len(n))
  for (i in seq_len(n)) {
    write.table(coords[,,i],
                file      = file.path(outdir, paste0(spec_names[i], ".pts")),
                row.names = FALSE,
                col.names = FALSE)
  }
  cat("✔ Exported", n, ".pts files to", outdir, "\n")
}

hum_outdir  <- "./2_Block_analysis_TEST/data_blocks/Humerus/all_pts"
scap_outdir <- "./2_Block_analysis_TEST/data_blocks/Scapula/all_pts"

write_pts_from_dataslide(hum_coords,  hum_outdir)
write_pts_from_dataslide(scap_coords, scap_outdir)

# ── Add version header + S numbering to .pts files ───────────────────────────
add_pts_header <- function(input_dir) {
  output_dir <- file.path(input_dir, "numbered_pts")
  dir.create(output_dir, showWarnings = FALSE)
  pts_files  <- list.files(input_dir, pattern = "\\.pts$", full.names = TRUE)
  
  for (file in pts_files) {
    coords   <- read.table(file, header = FALSE)
    n_lm     <- nrow(coords)
    labels   <- sprintf("S%03d", 0:(n_lm - 1))
    new_data <- cbind(labels, coords)
    header   <- c("Version 1.0", as.character(n_lm))
    out_file <- file.path(output_dir, basename(file))
    writeLines(header, out_file)
    write.table(new_data, file = out_file, append = TRUE,
                row.names = FALSE, col.names = FALSE, quote = FALSE)
  }
  cat("✔ Numbered .pts files written to", output_dir, "\n")
  return(output_dir)
}

hum_numbered  <- add_pts_header(hum_outdir)
scap_numbered <- add_pts_header(scap_outdir)

# ── Read numbered .pts files back into arrays ─────────────────────────────────
readPts_no_sl <- function(folder) {
  files       <- sort(list.files(folder, pattern = "\\.pts$", full.names = TRUE))
  coords_list <- lapply(files, function(f) {
    dat <- read.table(f, skip = 2)
    as.matrix(dat[, 2:4])
  })
  simplify2array(coords_list)
}

humerus <- readPts_no_sl(hum_numbered)
scapula <- readPts_no_sl(scap_numbered)

# GPA (no sliding needed — already slided)
block1 <- formatBlock(humerus, k = 86, gpa = TRUE)
block2 <- formatBlock(scapula, k = 86, gpa = TRUE)


# Combine into blocklist
blocklist <- combineBlocks(blocks = c(block1, block2))

# Rename specimens
spec_data <- read.csv("./2_Block_analysis/data/id_loco_humerus.csv", sep = ";")
block1    <- rename_block(block1, spec_data$filename)

spec_data <- read.csv("./2_Block_analysis/data/id_loco_scapula.csv", sep = ",")
block2    <- rename_block(block2, spec_data$filename)

spec_data <- read.csv("./2_Block_analysis/data/id_loco_general.csv", sep = ";")
blocklist <- rename_blocklist(blocklist, spec_data$filename)

# Sanity checks
cat("Block1 specimen 1:", dimnames(block1@gpa.3D)[[3]][1], "\n")
cat("Block2 specimen 1:", dimnames(block2@gpa.3D)[[3]][1], "\n")
cat("Blocklist row 1:",   rownames(blocklist@block.list$block_A)[1], "\n")

# ✅ Save WITH surface
save(block1,    file = "./2_Block_analysis_TEST/blocks/blocks_with_surf/block1.Rdata")
save(block2,    file = "./2_Block_analysis_TEST/blocks/blocks_with_surf/block2.Rdata")
save(blocklist, file = "./2_Block_analysis_TEST/blocks/blocks_with_surf/blocklist.Rdata")
cat("✔ Saved blocks WITH surface\n")
#############Analysis###############################################################
library(geomorph)
library(ape)

phy <- read.nexus(file = "./2_Block_analysis/data/mcctree.nexus")

### Metadata (shared by both pipelines)
meta          <- read.csv("./2_Block_analysis/data/plot.csv", header = TRUE, sep = ";")
meta$species  <- as.factor(meta$species)
stopifnot(nrow(meta) == 86)

# ── Shared functions ──────────────────────────────────────────────────────────
bootstrap_disparity <- function(coords_3D, species, n_boot = 999, min_n = NULL) {
  species   <- droplevels(as.factor(species))
  sp_levels <- levels(species)
  sp_counts <- table(species)
  if (is.null(min_n)) min_n <- max(2, min(sp_counts))
  keep_sp <- names(sp_counts)[sp_counts >= min_n]
  if (length(keep_sp) < length(sp_levels)) {
    message("Species excluded (n < min_n=", min_n, "): ",
            paste(setdiff(sp_levels, keep_sp), collapse = ", "))
  }
  proc_var <- function(mat) {
    mn <- rowMeans(mat); mean(colSums((mat - mn)^2))
  }
  flat <- function(arr) {
    d <- dim(arr); matrix(arr, nrow = d[1]*d[2], ncol = d[3])
  }
  coords_flat <- flat(coords_3D)
  boot_mat <- matrix(NA, nrow = n_boot, ncol = length(keep_sp),
                     dimnames = list(NULL, keep_sp))
  for (b in seq_len(n_boot))
    for (sp in keep_sp) {
      idx <- which(species == sp)
      boot_mat[b, sp] <- proc_var(coords_flat[, sample(idx, min_n), drop = FALSE])
    }
  list(mean   = colMeans(boot_mat),
       sd     = apply(boot_mat, 2, sd),
       ci_low = apply(boot_mat, 2, quantile, 0.025),
       ci_up  = apply(boot_mat, 2, quantile, 0.975),
       boot   = boot_mat, min_n = min_n, n_boot = n_boot)
}

plot_boot_disparity <- function(boot_res, title, col_bars = "steelblue",
                                col_ci = "grey30") {
  mn <- boot_res$mean; ci_low <- boot_res$ci_low; ci_up <- boot_res$ci_up
  bp <- barplot(mn, main = title,
                ylab = paste0("Procrustes Variance\n(bootstrapped, n=",
                              boot_res$min_n, "/spp, ", boot_res$n_boot, " iters)"),
                las = 2, col = col_bars, border = NA,
                ylim = c(0, max(ci_up) * 1.15), cex.names = 0.75)
  arrows(x0 = bp, y0 = ci_low, x1 = bp, y1 = ci_up,
         angle = 90, code = 3, length = 0.05, col = col_ci, lwd = 1.5)
}

compute_dimorphism <- function(coords_3D, species, sex) {
  species <- droplevels(as.factor(species))
  sex     <- factor(toupper(sex))
  valid   <- !is.na(sex)
  coords_3D <- coords_3D[,,valid]; species <- species[valid]; sex <- sex[valid]
  cat("Specimens after removing NA sex:", sum(valid), "\n")
  print(table(species, sex))
  sp_list <- levels(species)
  dimorphism <- sapply(sp_list, function(sp) {
    idx <- species == sp
    if (sum(idx & sex == "M") > 0 && sum(idx & sex == "F") > 0) {
      sqrt(sum((mshape(coords_3D[,,idx & sex == "M", drop = FALSE]) -
                  mshape(coords_3D[,,idx & sex == "F", drop = FALSE]))^2))
    } else NA
  })
  setNames(dimorphism, sp_list)
}

reshape_resid <- function(resid_mat, ref_array) {
  array(t(resid_mat), dim = dim(ref_array))
}

plot_dimorphism_disparity <- function(dimorphism, boot_sexcor, title,
                                      col_pts, sink_file) {
  shared <- intersect(names(dimorphism), names(boot_sexcor$mean))
  df <- na.omit(data.frame(
    species    = shared,
    dimorphism = dimorphism[shared],
    disparity  = boot_sexcor$mean[shared],
    disp_low   = boot_sexcor$ci_low[shared],
    disp_up    = boot_sexcor$ci_up[shared]))
  cat("\nSpecies in analysis:", nrow(df), "\n"); print(df)
  model <- lm(disparity ~ dimorphism, data = df)
  r2 <- round(summary(model)$r.squared, 3)
  pv <- round(summary(model)$coefficients[2, 4], 3)
  par(mar = c(5, 5, 3, 1))
  plot(df$dimorphism, df$disparity, pch = 19, col = col_pts,
       xlab = "Sexual dimorphism (Procrustes distance M vs F)",
       ylab = "Within-species shape disparity\n(sex-corrected, bootstrapped)",
       main = title,
       ylim = range(c(df$disp_low, df$disp_up)) * c(0.95, 1.05))
  arrows(x0 = df$dimorphism, y0 = df$disp_low,
         x1 = df$dimorphism, y1 = df$disp_up,
         angle = 90, code = 3, length = 0.04, col = "grey50", lwd = 1.5)
  text(df$dimorphism, df$disparity, labels = df$species, pos = 2, cex = 0.7)
  abline(model, col = "red", lwd = 2, lty = 2)
  legend("topleft", legend = paste0("R² = ", r2, "  p = ", pv), bty = "n", cex = 0.85)
  sink(sink_file); print(summary(model)); sink()
}

compute_species_means <- function(coords_3D, centroid, species_levels, meta_species) {
  n_sp  <- length(species_levels)
  n_lm  <- dim(coords_3D)[1]
  n_dim <- dim(coords_3D)[2]
  mean_shape <- array(NA, dim = c(n_lm, n_dim, n_sp),
                      dimnames = list(NULL, NULL, species_levels))
  mean_cs    <- setNames(numeric(n_sp), species_levels)
  for (sp in species_levels) {
    idx <- which(meta_species == sp)
    mean_shape[,,sp] <- apply(coords_3D[,,idx, drop = FALSE], c(1,2), mean)
    mean_cs[sp]      <- mean(centroid[idx])
  }
  list(shape = mean_shape, cs = mean_cs)
}

####--------------------------------------------------------------------
#### WITHOUT SURFACE
####--------------------------------------------------------------------
load(file = "./2_Block_analysis_TEST/blocks/blocks_without_surf/block1.Rdata")
load(file = "./2_Block_analysis_TEST/blocks/blocks_without_surf/block2.Rdata")
load(file = "./2_Block_analysis_TEST/blocks/blocks_without_surf/blocklist.Rdata")

# ── Allometry ────────────────────────────────────────────────────────────────
fit_allo_hum  <- procD.lm(block1@gpa.3D ~ log(block1@centroid))
fit_allo_scap <- procD.lm(block2@gpa.3D ~ log(block2@centroid))
summary(fit_allo_hum)
summary(fit_allo_scap)

# ── Raw bootstrapped disparity ────────────────────────────────────────────────
set.seed(42)
boot1_nosurf <- bootstrap_disparity(block1@gpa.3D, meta$species, n_boot = 999)
boot2_nosurf <- bootstrap_disparity(block2@gpa.3D, meta$species, n_boot = 999)

cat("Humerus (no surf) — subsampled to n =", boot1_nosurf$min_n, "per species\n")
cat("Scapula (no surf) — subsampled to n =", boot2_nosurf$min_n, "per species\n")

par(mar = c(8, 5, 3, 1))
plot_boot_disparity(boot1_nosurf, "Humerus Disparity — no surface (bootstrapped)")
plot_boot_disparity(boot2_nosurf, "Scapula Disparity — no surface (bootstrapped)")

# ── Sexual dimorphism ─────────────────────────────────────────────────────────
dimorphism_hum_nosurf  <- compute_dimorphism(block1@gpa.3D, meta$species, meta$sex)
dimorphism_scap_nosurf <- compute_dimorphism(block2@gpa.3D, meta$species, meta$sex)

# ── Sex-corrected bootstrapped disparity ──────────────────────────────────────
fit_sex_hum_nosurf  <- procD.lm(block1@gpa.3D ~ meta$sex)
fit_sex_scap_nosurf <- procD.lm(block2@gpa.3D ~ meta$sex)

set.seed(42)
boot1_nosurf_sexcor <- bootstrap_disparity(
  reshape_resid(fit_sex_hum_nosurf$residuals,  block1@gpa.3D), meta$species, n_boot = 999)
boot2_nosurf_sexcor <- bootstrap_disparity(
  reshape_resid(fit_sex_scap_nosurf$residuals, block2@gpa.3D), meta$species, n_boot = 999)

par(mar = c(8, 5, 3, 1))
plot_boot_disparity(boot1_nosurf_sexcor,
                    "Humerus Disparity — no surface (sex-corrected, bootstrapped)",
                    col_bars = "steelblue")
plot_boot_disparity(boot2_nosurf_sexcor,
                    "Scapula Disparity — no surface (sex-corrected, bootstrapped)",
                    col_bars = "darkorange")

# ── Dimorphism vs sex-corrected disparity ─────────────────────────────────────
plot_dimorphism_disparity(dimorphism_hum_nosurf,  boot1_nosurf_sexcor,
                          "Humerus: Dimorphism vs Disparity — no surface (sex-corrected)",
                          col_pts = "steelblue",
                          sink_file = "./outputs_TEST/disparity_dimorphism_nosurf_hum.txt")

plot_dimorphism_disparity(dimorphism_scap_nosurf, boot2_nosurf_sexcor,
                          "Scapula: Dimorphism vs Disparity — no surface (sex-corrected)",
                          col_pts = "darkorange",
                          sink_file = "./outputs_TEST/disparity_dimorphism_nosurf_scap.txt")

# ── Phylogenetic allometry on species means ───────────────────────────────────
species_levels <- levels(meta$species)

means_hum_nosurf  <- compute_species_means(block1@gpa.3D, block1@centroid,
                                           species_levels, meta$species)
means_scap_nosurf <- compute_species_means(block2@gpa.3D, block2@centroid,
                                           species_levels, meta$species)

cat("In data not in tree:", setdiff(species_levels, phy$tip.label), "\n")
cat("In tree not in data:", setdiff(phy$tip.label,  species_levels), "\n")

phy_pruned <- keep.tip(phy, species_levels)
tip_order  <- phy_pruned$tip.label

fit_allo_pgls_hum_nosurf <- procD.pgls(
  means_hum_nosurf$shape[,,tip_order] ~ log(means_hum_nosurf$cs[tip_order]),
  phy = phy_pruned, iter = 999)
summary(fit_allo_pgls_hum_nosurf)

fit_allo_pgls_scap_nosurf <- procD.pgls(
  means_scap_nosurf$shape[,,tip_order] ~ log(means_scap_nosurf$cs[tip_order]),
  phy = phy_pruned, iter = 999)
summary(fit_allo_pgls_scap_nosurf)

####--------------------------------------------------------------------
#### WITH SURFACE
####--------------------------------------------------------------------
load(file = "./2_Block_analysis_TEST/blocks/blocks_with_surf/block1.Rdata")
load(file = "./2_Block_analysis_TEST/blocks/blocks_with_surf/block2.Rdata")
load(file = "./2_Block_analysis_TEST/blocks/blocks_with_surf/blocklist.Rdata")

# ── Allometry ────────────────────────────────────────────────────────────────
fit_allo_hum  <- procD.lm(block1@gpa.3D ~ log(block1@centroid))
fit_allo_scap <- procD.lm(block2@gpa.3D ~ log(block2@centroid))
summary(fit_allo_hum)
summary(fit_allo_scap)

# ── Raw bootstrapped disparity ────────────────────────────────────────────────
set.seed(42)
boot1_surf <- bootstrap_disparity(block1@gpa.3D, meta$species, n_boot = 999)
boot2_surf <- bootstrap_disparity(block2@gpa.3D, meta$species, n_boot = 999)

cat("Humerus (surf) — subsampled to n =", boot1_surf$min_n, "per species\n")
cat("Scapula (surf) — subsampled to n =", boot2_surf$min_n, "per species\n")

par(mar = c(8, 5, 3, 1))
plot_boot_disparity(boot1_surf, "Humerus Disparity — with surface (bootstrapped)")
plot_boot_disparity(boot2_surf, "Scapula Disparity — with surface (bootstrapped)")

# ── Sexual dimorphism ─────────────────────────────────────────────────────────
dimorphism_hum_surf  <- compute_dimorphism(block1@gpa.3D, meta$species, meta$sex)
dimorphism_scap_surf <- compute_dimorphism(block2@gpa.3D, meta$species, meta$sex)

# ── Sex-corrected bootstrapped disparity ──────────────────────────────────────
fit_sex_hum_surf  <- procD.lm(block1@gpa.3D ~ meta$sex)
fit_sex_scap_surf <- procD.lm(block2@gpa.3D ~ meta$sex)

set.seed(42)
boot1_surf_sexcor <- bootstrap_disparity(
  reshape_resid(fit_sex_hum_surf$residuals,  block1@gpa.3D), meta$species, n_boot = 999)
boot2_surf_sexcor <- bootstrap_disparity(
  reshape_resid(fit_sex_scap_surf$residuals, block2@gpa.3D), meta$species, n_boot = 999)

par(mar = c(8, 5, 3, 1))
plot_boot_disparity(boot1_surf_sexcor,
                    "Humerus Disparity — with surface (sex-corrected, bootstrapped)",
                    col_bars = "steelblue")
plot_boot_disparity(boot2_surf_sexcor,
                    "Scapula Disparity — with surface (sex-corrected, bootstrapped)",
                    col_bars = "darkorange")

# ── Dimorphism vs sex-corrected disparity ─────────────────────────────────────
plot_dimorphism_disparity(dimorphism_hum_surf,  boot1_surf_sexcor,
                          "Humerus: Dimorphism vs Disparity — with surface (sex-corrected)",
                          col_pts = "steelblue",
                          sink_file = "./outputs_TEST/disparity_dimorphism_surf_hum.txt")

plot_dimorphism_disparity(dimorphism_scap_surf, boot2_surf_sexcor,
                          "Scapula: Dimorphism vs Disparity — with surface (sex-corrected)",
                          col_pts = "darkorange",
                          sink_file = "./outputs_TEST/disparity_dimorphism_surf_scap.txt")

# ── Phylogenetic allometry on species means ───────────────────────────────────
means_hum_surf  <- compute_species_means(block1@gpa.3D, block1@centroid,
                                         species_levels, meta$species)
means_scap_surf <- compute_species_means(block2@gpa.3D, block2@centroid,
                                         species_levels, meta$species)

fit_allo_pgls_hum_surf <- procD.pgls(
  means_hum_surf$shape[,,tip_order] ~ log(means_hum_surf$cs[tip_order]),
  phy = phy_pruned, iter = 999)
summary(fit_allo_pgls_hum_surf)

fit_allo_pgls_scap_surf <- procD.pgls(
  means_scap_surf$shape[,,tip_order] ~ log(means_scap_surf$cs[tip_order]),
  phy = phy_pruned, iter = 999)
summary(fit_allo_pgls_scap_surf)
########### RCPCA #################
plot <- read.csv("./2_Block_analysis/data/plot.csv", sep=";")
plot <- as.data.frame(plot)

library(RColorBrewer)

Lococolo <- as.factor(plot$locomotion)
palette <- brewer.pal(nlevels(Lococolo), "Dark2")
names(palette)<- levels(Lococolo)

color_loco <-palette[Lococolo]

result_RCPCA <- analyseBlocks(blocklist, option = "rcpca",ncomp = 5) 
result_RCPCA$result$AVE # variance for each PC axis (interest = superblock here)

n_comp <- 5
block_contrib <- data.frame(
  Component = paste0("PC", 1:n_comp),
  Block_A    = sapply(1:n_comp, function(i)
    cor(result_RCPCA$result$Y$block_A[, i],
        result_RCPCA$result$Y$superblock[, i])),
  Block_B    = sapply(1:n_comp, function(i)
    cor(result_RCPCA$result$Y$block_B[, i],
        result_RCPCA$result$Y$superblock[, i]))
)

# Absolute contributions (sign is arbitrary in RGCCA)
block_contrib$Block_A_abs <- abs(block_contrib$Block_A)
block_contrib$Block_B_abs <- abs(block_contrib$Block_B)

# Relative contribution (% share between the two blocks)
block_contrib$Block_A_pct <- block_contrib$Block_A_abs / 
  (block_contrib$Block_A_abs + block_contrib$Block_B_abs) * 100
block_contrib$Block_B_pct <- block_contrib$Block_B_abs / 
  (block_contrib$Block_A_abs + block_contrib$Block_B_abs) * 100

print(block_contrib)

scoresPlot(result_RCPCA,pcol = color_loco)#plabels = names$Species)# plot per block
#writexl::write_xlsx(as.data.frame(result_PCA$result$x), "RCPCA_scores.xlsx")
legend("topright",
       legend = levels(Lococolo),
       pt.bg = palette,
       pch = 21,
       pt.cex = 1,
       bty = "n")

scores <- result_RCPCA$result$Y$superblock
dim(scores)
df_scores <- data.frame(
  Comp1 = scores[,1],
  Comp2 = scores[,2],
  Comp3 = scores[,3],
  species = meta$species
)
boxplot(Comp1 ~ species,
        data = df_scores,
        las = 2,
        main = "RCPCA Axis 1 per species")
aggregate(Comp1 ~ species,
          data = df_scores,
          var)
#### PCA
result_PCA <- analyseBlocks(blocklist, option = "pca",ncomp = 5) ## Superblock
scoresPlot(result_PCA,pcol = color_loco)#plabels = names$Species)# plot per block

####################################################
####--------Data preparation ----------
##extract new data cordinates
#--humerus--
Ycoords_hum <- block1@gpa.3D
csize_hum <- block1@centroid
#--scapula--
Ycoords_scap <- block2@gpa.3D
csize_scap <- block2@centroid
#--superblock--
Ycoords_sb <- blocklist@block.list$superblock

###Create mean##############################
build_moyenne_csv <- function(folder_path, out_file) {
  files <- list.files(folder_path, pattern = "\\.pts$", full.names = FALSE)
  filename <- tools::file_path_sans_ext(files)
  species <- sapply(strsplit(filename, "_"), function(x) paste(x[1:2], collapse = "_"))
  df <- data.frame(filename = filename, species = species, stringsAsFactors = FALSE)
  write.csv(df, out_file, row.names = FALSE)
  df
}
build_moyenne_csv("./2_Block_analysis_TEST/data_blocks/Humerus/all_pts/numbered_pts",
                  "./2_Block_analysis_TEST/data/moyenne_humerus.csv")
build_moyenne_csv("./2_Block_analysis_TEST/data_blocks/Scapula/all_pts/numbered_pts",
                  "./2_Block_analysis_TEST/data/moyenne_scapula.csv")

identH <- read.csv("./2_Block_analysis_TEST/data/moyenne_humerus.csv",
                   sep = ",", header = TRUE, stringsAsFactors = FALSE)
identS <- read.csv("./2_Block_analysis_TEST/data/moyenne_scapula.csv",
                   sep = ",", header = TRUE, stringsAsFactors = FALSE)
identSB <- read.csv("./2_Block_analysis/data/moyenne_superblock.csv",
                    sep = ";", header = TRUE, stringsAsFactors = FALSE)

spnamesH  <- identH[,2]
spnamesS  <- identS[,2]
spnamesSB <- identSB[,2]

# transform data into 2d array
Yhum <- two.d.array(Ycoords_hum)
Yscap <- two.d.array(Ycoords_scap)
Ysuperblock <- (Ycoords_sb)

##### creating mean SHAPE coordinates per species (for arrayspecs) #####
#-----humerus-------
means_hum <- rowsum(Yhum, spnamesH)/as.vector(table(spnamesH))
means_csize_hum <- rowsum(csize_hum, spnamesH) / as.vector(table(spnamesH))
means_csize_df_hum <- data.frame(
  species = rownames(means_csize_hum),
  csize   = means_csize_hum[, 1],
  row.names = NULL
)
#-----scapula-------
means_scap <- rowsum(Yscap, spnamesS)/as.vector(table(spnamesS))
means_csize_scap <- rowsum(csize_scap, spnamesS) / as.vector(table(spnamesS))
means_csize_df_scap <- data.frame(
  species = rownames(means_csize_scap),
  csize   = means_csize_scap[, 1],
  row.names = NULL
)
#-----superblock-------
means_sb <- rowsum(Ysuperblock, spnamesSB)/as.vector(table(spnamesSB))

# import identification file
ident2 <- read.csv(file = "./2_Block_analysis/data/tab2.csv", header = FALSE, sep = ",")
library(ape)
library(geiger)
phy <- read.nexus(file = "./2_Block_analysis/data/mcctree.nexus")
plot(phy)

setdiff(rownames(means_sb), phy$tip.label)
setdiff(phy$tip.label, rownames(means_sb))

td <- treedata(phy, means_sb)
phy <- td$phy

rownames(ident2) <- ident2$V1
means_hum  <- means_hum[match(phy$tip.label, rownames(means_hum)),]
means_scap <- means_scap[match(phy$tip.label, rownames(means_scap)),]
means_sb   <- means_sb[match(phy$tip.label, rownames(means_sb)),]
ident2     <- ident2[match(phy$tip.label, rownames(ident2)),]

write.table(means_csize_df_hum,
            file = "./2_Block_analysis_TEST/mean_data/csize_mean_humerus.csv", sep = ";")
write.table(means_csize_df_scap,
            file = "./2_Block_analysis_TEST/mean_data/csize_mean_scapula.csv", sep = ";")

# create landmark array (directly from the SHAPE means computed above -
# do NOT reload from the csize CSV, that file only has species+csize,
# not landmark coordinates)
Yhum        <- arrayspecs(means_hum, 131, 3)  #humerus
Yscap       <- arrayspecs(means_scap, 172, 3) #scapula
Ysuperblock <- arrayspecs(means_sb, 303, 3)   #superblock

save(Yhum,        file = "./outputs_TEST/Yhum.Rdata")  
save(Yscap,       file = "./outputs_TEST/Yscap.Rdata")  
save(Ysuperblock, file = "./outputs_TEST/Y_sb.Rdata")

##### Prepare data for analyses ====================================================================

# import superimposed coordonates
load(file = "./outputs_TEST/Yhum.Rdata")
load(file = "./outputs_TEST/Yscap.Rdata")
load(file = "./outputs_TEST/Y_sb.Rdata")
# identify outliers
library(geomorph)
plotOutliers(Yhum)
plotOutliers(Yscap)
plotOutliers(Ysuperblock)

# loading identifications
ident <- read.csv("./2_Block_analysis/data/id_loco_humerus.csv", sep= ";")
ident_hum <- ident

ident <- read.csv("./2_Block_analysis/data/id_loco_scapula.csv", sep =",")
ident_scap <- ident

ident <- read.csv("./2_Block_analysis//data/id_loco_general.csv", sep =";")
ident_general <- ident

# import phylogenetic tree
library(ape)
phy <- read.nexus(file = "./2_Block_analysis/data/mcctree.nexus")
plot(phy,
     show.tip = TRUE, 
     cex = 1,
     type = "fan")

# ordering identification to match coordinates order
library(geomorph)
Y_2d_hum  <- two.d.array(Yhum)
Y_2d_scap <- two.d.array(Yscap)
Y_2d_sb   <- two.d.array(Ysuperblock)# transform 3D array into 2D matrix #choose the right Y
library(geiger)
td <- treedata(phy, Y_2d_sb) # match the phylogeny with dataset - superblock used as reference set of species
phy = td$phy # extract phylogy
plot(phy,
     show.tip = TRUE, 
     cex = 1,
     type = "fan")

ident_ordered_hum <- ident_hum[match(phy$tip.label, ident_hum$species),] # order identifications
ident_ordered_scap <- ident_scap[match(phy$tip.label, ident_scap$species),]
ident_ordered_sb <- ident_general[match(phy$tip.label, ident_general$species),]

# # saving table
write.table(ident_ordered_hum,
            file = paste0("./2_Block_analysis_TEST/mean_data/ident_ordered_hum.csv"),
            sep = ";")
write.table(ident_ordered_scap,
            file = paste0("./2_Block_analysis_TEST/mean_data/ident_ordered_scap.csv"),
            sep = ";")
write.table(ident_ordered_sb,
            file = paste0("./2_Block_analysis_TEST/mean_data/ident_ordered_sb.csv"),
            sep = ";")


####################################################
#### Multivariate analysis
####################################################

# load ident ordered
ident_ordered_hum  <- read.csv("./2_Block_analysis_TEST/mean_data/ident_ordered_hum.csv", sep=";")
ident_ordered_scap <- read.csv("./2_Block_analysis_TEST/mean_data/ident_ordered_scap.csv", sep=";")
ident_ordered_sb   <- read.csv("./2_Block_analysis_TEST/mean_data/ident_ordered_sb.csv", sep=";")

# load sizes - choose the bone 
csize_hum  <- read.csv("./2_Block_analysis_TEST/mean_data/csize_mean_humerus.csv", sep=";")
csize_scap <- read.csv("./2_Block_analysis_TEST/mean_data/csize_mean_scapula.csv", sep=";")


# load the phylogeny
library(ape)
phy <- read.nexus(file = "./2_Block_analysis/data/mcctree.nexus")

# load landmarks
load(file = "./outputs_TEST/Yhum.Rdata")
load(file = "./outputs_TEST/Yscap.Rdata")
load(file = "./outputs_TEST/Y_sb.Rdata")

library(geomorph)
Y_2d_hum  <- two.d.array(Yhum) # transform 3D array into 2D matrix
Y_2d_scap <- two.d.array(Yscap)
Y_2d_sb   <- two.d.array(Ysuperblock)
### arrange data
# species name order (humerus)
sp_names_hum <- ident_ordered_hum$species
family_names_hum <- ident_ordered_hum$family
names(family_names_hum) <- sp_names_hum
loco_hum <- ident_ordered_hum$locomotion
names(loco_hum) <- sp_names_hum

# species name order (scapula)
sp_names_scap <- ident_ordered_scap$species
family_names_scap <- ident_ordered_scap$family
names(family_names_scap) <- sp_names_scap
loco_scap <- ident_ordered_scap$locomotion
names(loco_scap) <- sp_names_scap

# species name order (superblock)
sp_names_sb <- ident_ordered_sb$species
family_names_sb <- ident_ordered_sb$family
names(family_names_sb) <- sp_names_sb
loco_sb <- ident_ordered_sb$locomotion
names(loco_sb) <- sp_names_sb

# centroid size (humerus)
colnames(csize_hum) <- c("sp", "csize")
library(dplyr)
data_csize_hum <- csize_hum %>% arrange(match(sp, sp_names_hum))
data_csize_hum <- data_csize_hum[,-1]
names(data_csize_hum) <- sp_names_hum

# centroid size (scapula)
colnames(csize_scap) <- c("sp", "csize")
data_csize_scap <- csize_scap %>% arrange(match(sp, sp_names_scap))
data_csize_scap <- data_csize_scap[,-1]
names(data_csize_scap) <- sp_names_scap

# create list for analyses
data_all_humerus <- list(Y = Y_2d_hum,
                         sp_names = sp_names_hum,
                         family_names = family_names_hum,
                         loco = loco_hum,
                         csize = data_csize_hum)
save(data_all_humerus,
     file = "./outputs_TEST/data_all_humerus.Rdata")

data_all_scapula <- list(Y = Y_2d_scap,
                         sp_names = sp_names_scap,
                         family_names = family_names_scap,
                         loco = loco_scap,
                         csize = data_csize_scap)
save(data_all_scapula,
     file = "./outputs_TEST/data_all_scapula.Rdata")

data_all_superblock <- list(Y = Y_2d_sb,
                            sp_names = sp_names_sb,
                            family_names = family_names_sb,
                            loco = loco_sb)
save(data_all_superblock,
     file = "./outputs_TEST/data_all_sb_nocsize.Rdata")

load(file = "./outputs_TEST/data_all_humerus.Rdata")
load(file = "./outputs_TEST/data_all_scapula.Rdata")
load(file = "./outputs_TEST/data_all_sb_nocsize.Rdata")

##### fitting evolutionary models ==================================================================
# Installing mvMORPH
library(devtools)
#install_github("JClavel/mvMORPH", build_vignettes = FALSE)
library(mvMORPH)
phy <- read.nexus(file = "./2_Block_analysis/data/mcctree.nexus")
#-----------------humerus---------------------
### identify the best model
# Brownian motion
# This is like mvBM(…,model="BM1")*
fit_bm <- mvgls(Y ~ loco + csize,
                data = data_all_humerus,
                tree = phy,
                model = "BM",
                method = "PL-LOOCV")
save(fit_bm,
     file = "./2_Block_analysis_TEST/modelsevol/bm_hum.Rdata")
# Ornstein-Uhlenbeck 
# This is like mvOU(…,model="OU1",param=list(decomp="equaldiagonal", vcv="fixedRoot")*
fit_ou <- mvgls(Y ~ loco + csize,
                data = data_all_humerus,
                tree = phy,
                model = "OU",
                method = "PL-LOOCV")
save(fit_ou,
     file = "./2_Block_analysis_TEST/modelsevol/ou_hum.Rdata")
# early burst 
# The rates decay is estimated jointly for all the traits. 
# When the “upper” limit is set >0, this corresponds to the ACDC model 
# (rates can increase). This is like mvEB(…)*
fit_eb <- mvgls(Y ~ loco + csize,
                data = data_all_humerus,
                tree = phy,
                model = "EB",
                method = "PL-LOOCV")
save(fit_eb,
     file = "./2_Block_analysis_TEST/modelsevol/eb_hum.Rdata")
# Pagel's lambda - a measure of “phylogenetic signal
fit_pl <- mvgls(Y ~ loco + csize,
                data = data_all_humerus,
                tree = phy,
                model = "lambda",
                method = "PL-LOOCV")
save(fit_pl,
     file = "./2_Block_analysis_TEST/modelsevol/pl_hum.Rdata")

load(file= "./2_Block_analysis_TEST/modelsevol/bm_hum.Rdata")
load(file= "./2_Block_analysis_TEST/modelsevol/ou_hum.Rdata")
load(file= "./2_Block_analysis_TEST/modelsevol/eb_hum.Rdata")
load(file= "./2_Block_analysis_TEST/modelsevol/pl_hum.Rdata")

logLik(fit_bm)
logLik(fit_ou)
logLik(fit_eb)
logLik(fit_pl)

## model comparison
summ_fit_bm <- summary(fit_bm) 
capture.output(summ_fit_bm, 
               file = "./outputs_TEST/summary_fit_bm.txt")
summ_fit_ou <- summary(fit_ou)
capture.output(summ_fit_ou, 
               file = "./outputs_TEST/summary_fit_ou.txt")
summ_fit_eb <- summary(fit_eb)
capture.output(summ_fit_eb, 
               file = "./outputs_TEST/summary_fit_eb.txt")
summ_fit_pl <- summary(fit_pl)
capture.output(summ_fit_pl, 
               file = "./outputs_TEST/summary_fit_pl.txt")

# Compare the models using the GIC criterion 
GIC(fit_bm)
GIC(fit_ou)
GIC(fit_eb)
GIC(fit_pl)

EIC(fit_ou)
EIC(fit_pl)

#-------------Scapula---------------------------------
### identify the best model
# Brownian motion
fit_bm <- mvgls(Y ~ loco + csize,
                data = data_all_scapula,
                tree = phy,
                model = "BM",
                method = "PL-LOOCV")
save(fit_bm,
     file = "./2_Block_analysis_TEST/modelsevol/bm_scap.Rdata")
# Ornstein-Uhlenbeck 
fit_ou <- mvgls(Y ~ loco + csize,
                data = data_all_scapula,
                tree = phy,
                model = "OU",
                method = "PL-LOOCV")
save(fit_ou,
     file = "./2_Block_analysis_TEST/modelsevol/ou_scap.Rdata")
# early burst 
fit_eb <- mvgls(Y ~ loco + csize,
                data = data_all_scapula,
                tree = phy,
                model = "EB",
                method = "PL-LOOCV")
save(fit_eb,
     file = "./2_Block_analysis_TEST/modelsevol/eb_scap.Rdata")
# Pagel's lambda - a measure of “phylogenetic signal
fit_pl <- mvgls(Y ~ loco + csize,
                data = data_all_scapula,
                tree = phy,
                model = "lambda",
                method = "PL-LOOCV")
save(fit_pl,
     file = "./2_Block_analysis_TEST/modelsevol/pl_scap.Rdata")

load(file= "./2_Block_analysis_TEST/modelsevol/bm_scap.Rdata")
load(file= "./2_Block_analysis_TEST/modelsevol/ou_scap.Rdata")
load(file= "./2_Block_analysis_TEST/modelsevol/eb_scap.Rdata")
load(file= "./2_Block_analysis_TEST/modelsevol/pl_scap.Rdata")

logLik(fit_bm)
logLik(fit_ou)
logLik(fit_eb)
logLik(fit_pl)

## model comparison
summ_fit_bm <- summary(fit_bm) 
capture.output(summ_fit_bm, 
               file = "./outputs_TEST/summary_fit_bm_scap.txt")
summ_fit_ou <- summary(fit_ou)
capture.output(summ_fit_ou, 
               file = "./outputs_TEST/summary_fit_ou_scap.txt")
summ_fit_eb <- summary(fit_eb)
capture.output(summ_fit_eb, 
               file = "./outputs_TEST/summary_fit_eb_scap.txt")
summ_fit_pl <- summary(fit_pl)
capture.output(summ_fit_pl, 
               file = "./outputs_TEST/summary_fit_pl_scap.txt")
# Compare the models using the GIC criterion 
GIC(fit_bm)
GIC(fit_ou)
GIC(fit_eb)
GIC(fit_pl)

EIC(fit_ou)
EIC(fit_pl)

#-------------superblock---------------------------------
### identify the best model
# Brownian motion
fit_bm <- mvgls(Y ~ loco,
                data = data_all_superblock,
                tree = phy,
                model = "BM",
                method = "PL-LOOCV")
save(fit_bm,
     file = "./2_Block_analysis_TEST/modelsevol/bm_sb.Rdata")
# Ornstein-Uhlenbeck 
fit_ou <- mvgls(Y ~ loco,
                data = data_all_superblock,
                tree = phy,
                model = "OU",
                method = "PL-LOOCV")
save(fit_ou,
     file = "./2_Block_analysis_TEST/modelsevol/ou_sb.Rdata")
# early burst 
fit_eb <- mvgls(Y ~ loco,
                data = data_all_superblock,
                tree = phy,
                model = "EB",
                method = "PL-LOOCV")
save(fit_eb,
     file = "./2_Block_analysis_TEST/modelsevol/eb_sb.Rdata")
# Pagel's lambda - a measure of “phylogenetic signal
fit_pl <- mvgls(Y ~ loco,
                data = data_all_superblock,
                tree = phy,
                model = "lambda",
                method = "PL-LOOCV")
save(fit_pl,
     file = "./2_Block_analysis_TEST/modelsevol/pl_sb.Rdata")

load("./outputs_TEST/data_all_sb_nocsize.Rdata")

load(file= "./2_Block_analysis_TEST/modelsevol/bm_sb.Rdata")
load(file= "./2_Block_analysis_TEST/modelsevol/ou_sb.Rdata")
load(file= "./2_Block_analysis_TEST/modelsevol/eb_sb.Rdata")
load(file= "./2_Block_analysis_TEST/modelsevol/pl_sb.Rdata")

logLik(fit_bm)
logLik(fit_ou)
logLik(fit_eb)
logLik(fit_pl)

## model comparison
summ_fit_bm <- summary(fit_bm) 
capture.output(summ_fit_bm, 
               file = "./outputs_TEST/summary_fit_bm_sb.txt")
summ_fit_ou <- summary(fit_ou)
capture.output(summ_fit_ou, 
               file = "./outputs_TEST/summary_fit_ou_sb.txt")
summ_fit_eb <- summary(fit_eb)
capture.output(summ_fit_eb, 
               file = "./outputs_TEST/summary_fit_eb_sb.txt")
summ_fit_pl <- summary(fit_pl)
capture.output(summ_fit_pl, 
               file = "./outputs_TEST/summary_fit_pl_sb.txt")
GIC(fit_bm)
GIC(fit_ou)
GIC(fit_eb)
GIC(fit_pl)

EIC(fit_bm)
EIC(fit_ou)
EIC(fit_eb)
EIC(fit_pl)

#######simulation###### for illustrations
library(mvMORPH)
###humerus####
pca <- prcomp(data_all_humerus$Y)

# keep PCs explaining ~95% variance
summary(pca)

Y_pc <- pca$x[,1:15]
fit_BM <- mvBM(phy, Y_pc)

sim_BM <- mvSIM(phy, nsim = 10, model = "BM", param = fit_BM)
sim_mat <- sapply(sim_BM, function(x) x[,1])

library(phytools)

h <- nodeHeights(phy)
time <- max(h) - h[,2]

matplot(sim_mat,
  type = "l",
  lwd = 2,
  lty = 1,
  xlab = "evolutionary distance",
  ylab = "Humeral Trait value",
  main = "Brownian Motion simulations"
)
### scapula####
pca <- prcomp(data_all_scapula$Y)
summary(pca)
var <- cumsum(pca$sdev^2 / sum(pca$sdev^2))
n_pc <- which(var >= 0.95)[1]

Y_pc <- pca$x[,1:n_pc]
fit_OU <- mvOU(phy, Y_pc)

# simulate trait evolution
sim_OU <- mvSIM(phy, nsim = 10, model = "OU", param = fit_OU)
sim_mat <- sapply(sim_OU, function(x) x[,1])

library(phytools)

h <- nodeHeights(phy)
time <- max(h) - h[,2]

matplot(sim_mat,
        type = "l",
        lwd = 2,
        lty = 1,
        xlab = "evolutionary distance",
        ylab = "Scapular Trait value",
        main = "Orstein Ulhenbeck simulations"
)

################ MANOVAs ##############
#---------humerus---------------------
load(file = "./2_Block_analysis_TEST/modelsevol/pl_hum.Rdata")
results <- manova.gls(
  fit_pl,
  test = "Wilks",
  nperm = 999,
  nbcores = 4L, # to speed up the calculation
  type = "II")

capture.output(results, 
               file = "./results_TEST/manova_all_hum_surf.txt")
save(results,
     file = "./results_TEST/manova_humerus_surf.Rdata")

#---------scapula---------------------
load(file = "./2_Block_analysis_TEST/modelsevol/ou_scap.Rdata")
results <- manova.gls(
  fit_ou,
  test = "Wilks",
  nperm = 999,
  nbcores = 4L, # to speed up the calculation
  type = "II")

capture.output(results, 
               file = "./results_TEST/manova_all_scap_surf.txt")
save(results,
     file = "./results_TEST/manova_scapula_surf.Rdata")

#---------superblock---------------------
load(file = "./2_Block_analysis_TEST/modelsevol/pl_sb.Rdata")
results <- manova.gls(
  fit_pl,
  test = "Wilks",
  nperm = 999,
  nbcores = 4L, # to speed up the calculation
  type = "II")

capture.output(results, 
               file = "./results_TEST/manova_all_sb_surf.txt")
save(results,
     file = "./results_TEST/manova_superblock_surf.Rdata")

################ ACPs ####################
##chargement des données 
library(ape)
library(dplyr)
#load phylogeny
phy <- read.nexus(file = "./2_Block_analysis/data/mcctree.nexus")

# load ident ordered
ident_ordered <- read.csv("./2_Block_analysis/data/ident_ordered.csv")

# load landmarks
load(file = "./outputs_TEST/Yhum.Rdata")
load(file = "./outputs_TEST/Yscap.Rdata")
load(file = "./outputs_TEST/Y_sb.Rdata")

library(geomorph)
Y_2d_hum  <- two.d.array(Yhum)
Y_2d_scap <- two.d.array(Yscap)
Y_2d_sb   <- two.d.array(Ysuperblock)


# Phylogenetic allometry test
fit_allo_pgls <- procD.pgls(data_all_humerus$Y ~ log(data_all_humerus$csize), 
                            phy = phy, iter = 999)
summary(fit_allo_pgls)

fit_allo_scap_simple <- procD.pgls(data_all_scapula$Y ~ log(data_all_scapula$csize), 
                            phy = phy, iter = 999)
summary(fit_allo_scap_simple)

#######correcting allometry for scapula ##########################
# Step 1 — correct mean shape (2D matrix, no mshape needed)
fit_allo_scap_simple <- procD.lm(
  data_all_scapula$Y ~ log(data_all_scapula$csize), 
  iter = 999
)
mean_shape_scap <- colMeans(data_all_scapula$Y)  # 1 x (p*k) mean

scap_corrected <- fit_allo_scap_simple$residuals + 
  matrix(mean_shape_scap,
         nrow = nrow(fit_allo_scap_simple$residuals),
         ncol = ncol(fit_allo_scap_simple$residuals),
         byrow = TRUE)

# Step 2 — check dimensions before arrayspecs
n_spec <- nrow(scap_corrected)          # number of species (rows)
n_vars <- ncol(scap_corrected)          # total variables = p * k

# What are p and k for your scapula?
# Since block2 is 3D landmarks in X,Y,Z: k = 3, p = n_vars / 3
k <- 3
p <- n_vars / k
cat("n =", n_spec, "| p =", p, "| k =", k, "| p*k =", p*k, "\n")

# Step 3 — convert to array with correct dims
scap_corrected_array <- arrayspecs(scap_corrected, p = p, k = k)

# Step 4 — add species names (critical for gm.prcomp + phylogeny matching)
dimnames(scap_corrected_array)[[3]] <- rownames(scap_corrected)

# Step 5 — re-run pPCA
PCA_scap_corrected <- gm.prcomp(scap_corrected_array, phy = phy)

#######pPCA

#install.packages("RColorBrewer")  # only if not installed
library(RColorBrewer)

library(phytools)

lococo <- as.factor(ident_ordered$Locomotion)
# Number of locomotion categories
n_loco <- nlevels(lococo)

# Generate colors
color_loco <- brewer.pal(n_loco, "Dark2")

#pPCA sur locomotion - humerus
PCA<-gm.prcomp(
  Y_2d_hum,
  phy = phy,
)


plot(PCA, main = "PCA - Humerus")
d<-PCA$x
PCA_df <- as.data.frame(d)
x<-PCA_df$Comp1
y<-PCA_df$Comp2

#humerus
phylomorphospace(tree = phy,
                 X = d[, 1:2],
                 A = NULL,
                 xlab = "Principal component 1 (52.95%)",
                 ylab = "Principal component 2 (24.3%)",
                 node.size = 0,
                 label = "NULL", # label of each point
                 lwd = 0.5, # line widths
)
points(x = x,
       y = y,
       col = "black",
       pch =21, 
       bg = color_loco[lococo],
       cex = 1) 

legend("topright",
       legend = levels(lococo),
       pt.bg = color_loco,
       pch = 21,
       pt.cex = 1,
       bty = "n",
       xjust = 1,
       yjust = 1)


#scapula
plot(PCA_scap_corrected, main = "PCA - Scapula")
d<-PCA_scap_corrected$x
PCA_df <- as.data.frame(d)
x<-PCA_df$Comp1
y<-PCA_df$Comp2

phylomorphospace(tree = phy,
                 X = d[, 1:2],
                 A = NULL,
                 xlab = "Principal component 1 (42.18)",
                 ylab = "Principal component 2 (19.26)",
                 node.size = 0,
                 label = "NULL",        
                 lwd = 0.5)

points(x = x,
       y = y,
       col = "black",
       pch =21, 
       bg = color_loco[lococo],
       cex = 1) 

legend("topright",
       legend = levels(lococo),
       pt.bg = color_loco,
       pch = 21,
       pt.cex = 1,
       bty = "n",
       xjust = 1,
       yjust = 1)

#superblock
PCA_sb <- gm.prcomp(Y_2d_sb, phy = phy)
plot(PCA_sb, main = "PCA - Superblock")
d <- PCA_sb$x
PCA_df <- as.data.frame(d)
x <- PCA_df$Comp1
y <- PCA_df$Comp2

phylomorphospace(tree = phy,
                 X = d[, 1:2],
                 A = NULL,
                 xlab = "Principal component 1 (56.82%)",
                 ylab = "Principal component 2 (12.78%)",
                 node.size = 0,
                 label = "NULL", # label of each point
                 lwd = 0.5, # line widths
)
points(x = x,
       y = y,
       col = "black",
       pch =21, 
       bg = color_loco[lococo],
       cex = 1) 

legend("topright",
       legend = levels(lococo),
       pt.bg = color_loco,
       pch = 21,
       pt.cex = 1,
       bty = "n",
       xjust = 1,
       yjust = 1)


#### defomation ####
#-------humerus------
PCA <- gm.prcomp(Yhum, phy = phy)
p <- dim(Yhum)[1]
k <- dim(Yhum)[2]

mean_shape <- mshape(Yhum)

pc1 <- PCA$rotation[,1]
sd_pc1 <- sd(PCA$x[,1])
shape_pc1_plus  <- mean_shape + matrix(pc1 * sd_pc1 * 2, p, k, byrow = TRUE)
shape_pc1_minus <- mean_shape - matrix(pc1 * sd_pc1 * 2, p, k, byrow = TRUE)

#PC1-
plotRefToTarget(
  mean_shape,
  shape_pc1_minus,
  method = "TPS",
  mag = 1
)
#PC1+
plotRefToTarget(
  mean_shape,
  shape_pc1_plus,
  method = "TPS",
  mag = 1
)

pc2 <- PCA$rotation[,2]
sd_pc2 <- sd(PCA$x[,2])
shape_pc2_plus  <- mean_shape + matrix(pc2 * sd_pc2 * 2, p, k, byrow = TRUE)
shape_pc2_minus <- mean_shape - matrix(pc2 * sd_pc2 * 2, p, k, byrow = TRUE)
#PC2-
plotRefToTarget(
  mean_shape,
  shape_pc2_minus,
  method = "TPS",
  mag = 1
)
#PC2+
plotRefToTarget(
  mean_shape,
  shape_pc2_plus,
  method = "TPS",
  mag = 1
)

#-------scapula------
# ── Scapula pPCA on allometry-CORRECTED coordinates ──────────────
pPCA_scap <- gm.prcomp(scap_corrected_array, phy = phy)
p <- dim(scap_corrected_array)[1]   # number of landmarks
k <- dim(scap_corrected_array)[2]   # dimensions (3 for 3D)

mean_shape <- mshape(scap_corrected_array)   # mean shape of corrected coords
pc1 <- pPCA_scap$rotation[, 1]      # PC1 eigenvector (p*k vector)
sd_pc1 <- sd(pPCA_scap$x[, 1])      # SD of scores along PC1
pc1_matrix <- matrix(pc1, p, k)

shape_pc1_plus  <- mean_shape + pc1_matrix * sd_pc1 * 2   # +2SD
shape_pc1_minus <- mean_shape - pc1_matrix * sd_pc1 * 2   # -2SD

#PC1-
plotRefToTarget(
  mean_shape,
  shape_pc1_minus,
  method = "TPS",
  mag = 1
)
#PC1+
plotRefToTarget(
  mean_shape,
  shape_pc1_plus,
  method = "TPS",
  mag = 1
)

pc2 <- pPCA_scap$rotation[,2]
sd_pc2 <- sd(pPCA_scap$x[, 2])      # SD of scores along PC1
pc2_matrix <- matrix(pc2, p, k)

shape_pc2_plus  <- mean_shape + pc2_matrix * sd_pc2 * 2   # +2SD
shape_pc2_minus <- mean_shape - pc2_matrix * sd_pc2 * 2   # -2SD
#PC2-
plotRefToTarget(
  mean_shape,
  shape_pc2_minus,
  method = "TPS",
  mag = 1
)
#PC2+
plotRefToTarget(
  mean_shape,
  shape_pc2_plus,
  method = "TPS",
  mag = 1
)

#### Loading of each landmarks ####
load(file = "./outputs_TEST/Yhum.Rdata")
load(file = "./outputs_TEST/Yscap.Rdata")
load(file = "./outputs_TEST/Y_sb.Rdata")
#------ humerus surf -------------
phy <- read.nexus(file = "./2_Block_analysis/data/mcctree.nexus")

PCA <- gm.prcomp(Yhum, phy = phy)
k <- dim(Yhum)[2]

######PC1########
rot <- PCA$rotation[,1]
# reshape to landmark × dimension
rot_mat <- matrix(rot, ncol = k, byrow = TRUE)
# compute magnitude
loading <- sqrt(rowSums(rot_mat^2))

loading_scaled <- (loading - min(loading)) /
  (max(loading) - min(loading))
cols <- colorRampPalette(c("#56BF47", "purple"))(100)
landmark_colors <- cols[round(loading_scaled * 99) + 1]

library(Morpho)
library(rgl)
library(geomorph)
library(RANN)

template_lm <- as.matrix(read.table("./2_Block_analysis_TEST/data_blocks/Humerus/all_pts/Nasalis_larvatus_m_FMNH_68684_humerus_L.pts"))
mean_shape <- mshape(Yhum)
mesh_template <- vcgImport("./Datas/Resampling/all_pts/all_ply/Nasalis_larvatus_m_FMNH_68684_humerus_L.ply")

mesh_mean <- tps3d(mesh_template,
                   refmat = template_lm,
                   tarmat = mean_shape)
inflation <- 1.02
mean_mesh_inflated <- mean_shape * inflation

open3d()
shade3d(mesh_mean, color="beige", alpha = 0.8)
points3d(mean_mesh_inflated,
         col = landmark_colors,
         size = 10)
rgl.snapshot("./outputs_TEST/humerus_PC1_loadings.png", fmt="png")

barplot(loading,
        xlab = "Landmarks",
        ylab = "Loading magnitude",
        main = "Landmark contributions to PC1")
######PC2####################
rot <- PCA$rotation[,2]
# reshape to landmark × dimension
rot_mat <- matrix(rot, ncol = k, byrow = TRUE)
# compute magnitude
loading <- sqrt(rowSums(rot_mat^2))

loading_scaled <- (loading - min(loading)) /
  (max(loading) - min(loading))
cols <- colorRampPalette(c("#56BF47", "purple"))(100)
landmark_colors <- cols[round(loading_scaled * 99) + 1]

template_lm <- as.matrix(read.table("./2_Block_analysis_TEST/data_blocks/Humerus/all_pts/Nasalis_larvatus_m_FMNH_68684_humerus_L.pts"))
mean_shape <- mshape(Yhum)
mesh_template <- vcgImport("./Datas/Resampling/all_pts/all_ply/Nasalis_larvatus_m_FMNH_68684_humerus_L.ply")

mesh_mean <- tps3d(mesh_template,
                   refmat = template_lm,
                   tarmat = mean_shape)
inflation <- 1.02
mean_mesh_inflated <- mean_shape * inflation

open3d()
shade3d(mesh_mean, color="beige", alpha = 0.8)
points3d(mean_mesh_inflated,
         col = landmark_colors,
         size = 10)
rgl.snapshot("./outputs_TEST/humerus_PC2_loadings_ventral.png", fmt="png")

barplot(loading,
       xlab = "Landmarks",
       ylab = "Loading magnitude",
       main = "Landmark contributions to PC2")
#------ scapula surf  -------------
phy <- read.nexus(file = "./2_Block_analysis/data/mcctree.nexus")

PCA <- gm.prcomp(Yscap, phy = phy)
k <- dim(Yscap)[2]

######PC1########
rot <- PCA$rotation[,1]
# reshape to landmark × dimension
rot_mat <- matrix(rot, ncol = k, byrow = TRUE)
# compute magnitude
loading <- sqrt(rowSums(rot_mat^2))

loading_scaled <- (loading - min(loading)) /
  (max(loading) - min(loading))
cols <- colorRampPalette(c("#56BF47", "purple"))(100)
landmark_colors <- cols[round(loading_scaled * 99) + 1]

library(Morpho)
library(rgl)
library(geomorph)
library(RANN)

template_lm <- as.matrix(read.table("./2_Block_analysis_TEST/data_blocks/Scapula/all_pts/Colobus_guereza_f_AMNH_52241_scapula_L.pts"))
mean_shape <- mshape(Yscap)
mesh_template <- vcgImport("./Datas/Resampling/all_pts/all_ply/Colobus_guereza_f_AMNH_52241_scapula_L.ply")

mesh_mean <- tps3d(mesh_template,
                   refmat = template_lm,
                   tarmat = mean_shape)

open3d()
shade3d(mesh_mean, color="beige", alpha = 0.8)
points3d(mean_shape,
         col = landmark_colors,
         size = 13)
rgl.snapshot("./outputs_TEST/scapula_PC1_loadings_ventral.png", fmt="png")

barplot(loading,
        xlab = "Landmarks",
        ylab = "Loading magnitude",
        main = "Landmark contributions to PC1")
######PC2####################
rot <- PCA$rotation[,2]
# reshape to landmark × dimension
rot_mat <- matrix(rot, ncol = k, byrow = TRUE)
# compute magnitude
loading <- sqrt(rowSums(rot_mat^2))

loading_scaled <- (loading - min(loading)) /
  (max(loading) - min(loading))
cols <- colorRampPalette(c("#56BF47", "purple"))(100)
landmark_colors <- cols[round(loading_scaled * 99) + 1]

library(Morpho)
library(rgl)
library(geomorph)
library(RANN)

template_lm <- as.matrix(read.table("./2_Block_analysis_TEST/data_blocks/Scapula/all_pts/Colobus_guereza_f_AMNH_52241_scapula_L.pts"))
mean_shape <- mshape(Yscap)
mesh_template <- vcgImport("./Datas/Resampling/all_pts/all_ply/Colobus_guereza_f_AMNH_52241_scapula_L.ply")

mesh_mean <- tps3d(mesh_template,
                   refmat = template_lm,
                   tarmat = mean_shape)

open3d()
shade3d(mesh_mean, color="beige", alpha = 0.8)
points3d(mean_shape,
         col = landmark_colors,
         size = 10)
rgl.snapshot("./outputs_TEST/scapula_PC2_loadings_ventral.png", fmt="png")

barplot(loading,
        xlab = "Landmarks",
        ylab = "Loading magnitude",
        main = "Landmark contributions to PC2")

########## creating figure paper ############
legend_cols <- colorRampPalette(c("#56BF47", "purple"))(100)

png("./outputs_TEST/legend.png", width=400, height=120)

par(mar=c(2,2,2,2))
image(
  z = matrix(seq(0,1,length=100), ncol=1),
  col = legend_cols,
  axes = FALSE
)

axis(1, at=c(0,1), labels=c("Low loading","High loading"))

dev.off()

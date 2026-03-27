air_rearrangement <- read.delim(file = "VDJ/airr_rearrangement.tsv")
all_contig <- read.delim("VDJ/all_contig_annotations.csv", sep = ",")
concensus <- read.csv("VDJ/consensus_annotations.csv")


##  Setting directories
lib3 <- "Lib3"
objectsDir <- file.path("outputObjects", lib3)
figuresDir <- file.path("outputFigures", lib3)
dir.create(objectsDir, recursive = TRUE, showWarnings = FALSE)
dir.create(figuresDir, recursive = TRUE, showWarnings = FALSE)
getwd()
library(Seurat)
##  Importing 10x data      --------------------------------------
rawData <- Read10X("C:/Users/ktari/Desktop/Analysis scRNAseq/Lib3/filtered_feature_bc_matrix/")
library(dplyr)
dataDir <- "C:/Users/ktari/Desktop/Analysis scRNAseq/Lib3"
vdjData <- read.delim(
  file.path(dataDir, "VDJ", "airr_rearrangement.tsv"),
  header = TRUE,
  sep = "\t"
) %>%
  dplyr::select(-c(sequence, sequence_aa, sequence_alignment, germline_alignment))
vdjDataSupp <- read.csv(
  file.path(dataDir, "VDJ", "all_contig_annotations.csv"),
  header = TRUE
)
barcodes <- intersect(colnames(rawData$`Gene Expression`),
                      colnames(rawData$`Antibody Capture`))
# Extraire RNA et HTO à partir de rawData et garder les mêmes barcodes
umi <- rawData$`Gene Expression`[, barcodes]
hto <- rawData$`Antibody Capture`[, barcodes]
##  Creating seurat object & normalizing RNA data (log normalizing)   -------------
seuratObj <- CreateSeuratObject(counts = umi)
seuratObj <- NormalizeData(seuratObj)

##  Find and scale variable features        ---------------------------------
seuratObj <- FindVariableFeatures(seuratObj, selection.method = "mean.var.plot", nfeatures = 2000)
seuratObj <- ScaleData(seuratObj, features = VariableFeatures(seuratObj))

##  Add HTO data as an independent assay & normalize      --------------------------
seuratObj[["HTO"]] <- CreateAssayObject(counts = hto)
seuratObj <- NormalizeData(seuratObj, assay = "HTO", normalization.method = "CLR")

## QC for HTODemux ---------------------- 
library(ggplot2)
library(dplyr)
# créer fontion colpalettegenerator 
colPaletteGenerator <- function(n) { 
  grDevices::colorRampPalette(c(
    "#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02", "#a6761d", "#666666" ))(n) }




# Créer le sous-dossier QC
qcDir <- file.path(figuresDir, "QC")
dir.create(qcDir, recursive = TRUE, showWarnings = FALSE)

pdf(
  file = file.path(qcDir, paste0(lib3, "_QC_HTODemux.pdf")),
  width = 8,
  height = 6,
  title = "QC Demultiplexing"
)

for (i in c(0.95, 0.96, 0.97, 0.98, 0.99, 0.995)) {
  
  seuratTmp <- HTODemux(
    seuratObj,
    assay = "HTO",
    positive.quantile = i
  )
  
  # Plot 1: répartition hash.ID
  p1 <- seuratTmp@meta.data %>%
    dplyr::count(hash.ID) %>%
    ggplot() +
    geom_bar(
      aes(x = hash.ID, y = n, fill = hash.ID),
      width = 0.7,
      stat = "identity",
      position = "stack"
    ) +
    scale_fill_manual(values = colPaletteGenerator(4)) +
    labs(
      x = "",
      y = "",
      title = paste0("Positive.quantile = ", i),
      fill = "HTO"
    ) +
    theme_bw() +
    theme(
      axis.text = element_text(size = 11, face = "bold"),
      legend.title = element_text(size = 12, face = "bold", colour = "darkred"),
      plot.title = element_text(size = 15, face = "bold", colour = "darkred")
    )
  
  print(p1)
  
  # Plot 2: distribution UMI par hash.ID
  p2 <- VlnPlot(
    seuratTmp,
    features = "nCount_RNA",
    pt.size = 0.1,
    log = TRUE,
    group.by = "hash.ID"
  ) +
    labs(
      title = paste0("UMI distributions at ", i),
      x = ""
    ) +
    scale_fill_manual(values = colPaletteGenerator(4))
  
  print(p2)
}

dev.off()

##  Selecting a proper threshold ------------------
seuratObj <- HTODemux(seuratObj, assay = "HTO", positive.quantile = 0.99)
table(seuratObj$HTO_classification.global)
table(seuratObj$hash.ID)  
colnames(seuratObj@meta.data)
Assays(seuratObj)
seuratObj[["HTO"]]
seuratObj <- HTODemux(
  object = seuratObj,
  assay = "HTO",
  positive.quantile = 0.99
)
colnames(seuratObj@meta.data)



##  Filtering : keep singlet   ------------------------
Idents(seuratObj) <- seuratObj$HTO_classification.global
seuratObj <- subset(seuratObj, idents = "Singlet", invert = FALSE)

##  QC for MT genes expression     -------------------------------
seuratObj <- PercentageFeatureSet(seuratObj,
                                  pattern = "^MT-",
                                  col.name = "percent.mt",
                                  assay = "RNA")

##  Visualizing MT pct genes & Importing as PDF    ---------------------------
# Dossier QC dans outputFigures/Lib3
qcDir <- file.path(figuresDir, "QC")
dir.create(qcDir, recursive = TRUE, showWarnings = FALSE)

# Export PDF dans Lib3/QC
pdf(
  file = file.path(qcDir, paste0(lib3, "_MT_pct_checkout.pdf")),
  width = 8,
  height = 6,
  title = "Quality Control"
)

print(
  VlnPlot(
    object = seuratObj,
    features = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
    pt.size = 0.01,
    log = TRUE
  )
)

print(
  ggplot(seuratObj@meta.data, aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt)) +
    geom_point(size = 0.8) +
    scale_y_log10(breaks = c(90, 100, 300, 1000, 3000)) +
    scale_color_gradientn(
      colours = c("green", "gold", "orange", "darkred"),
      values = c(0, 0.2, 0.5, 1),
      limits = c(0, 50)
    ) +
    ggtitle("QC plot", "Number of detected genes in function of number of UMI") +
    labs(x = "Number of UMI by cell", y = "Number of detected genes by cell") +
    theme_bw()
)

dev.off()

##  Filtering cells and visualization as PDF    -----------------------------
mtThreshold <- 15

qcDir <- file.path(figuresDir, "QC")
dir.create(qcDir, recursive = TRUE, showWarnings = FALSE)

seuratObj <- subset(
  seuratObj,
  subset = percent.mt < mtThreshold & nFeature_RNA < 3500 & nFeature_RNA >200 & nCount_RNA > 500
)

pdf(
  file = file.path(qcDir, paste0(lib3, "_MT_pct_filtered.pdf")),
  width = 8,
  height = 6,
  title = "Quality Control - filtered"
)

print(
  VlnPlot(
    object = seuratObj,
    features = c("nCount_RNA","nFeature_RNA","percent.mt"),
    pt.size = 0.01,
    log = TRUE
  )
)

print(
  ggplot(seuratObj@meta.data, aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt)) +
    geom_point(size = 0.8) +
    scale_y_log10(breaks = c(90,100,300,1000,3000)) +
    scale_color_gradientn(
      colours = c("green","gold","orange","darkred"),
      values = c(0,0.2,0.5,1),
      limits = c(0,50)
    ) +
    ggtitle("QC plot", "Number of detected genes in function of number of UMI") +
    labs(
      x = "Number of UMI by cell",
      y = "Number of detected genes by cell",
      caption = paste0("Chosen threshold: ", mtThreshold)
    ) +
    theme_bw()
)

dev.off()

colnames(vdjData)

##  Add vdj data to seurat obj  
# 1) Match VDJ rows to Seurat cells
vdj_use <- vdjData %>%
  filter(cell_id %in% colnames(seuratObj))

cat("Matched VDJ rows  :", nrow(vdj_use), "\n")
cat("Matched VDJ cells :", n_distinct(vdj_use$cell_id), "\n")

# 2) Productive flag
vdj_use <- vdj_use %>%
  mutate(productive_flag = trimws(tolower(as.character(productive))) %in% c("true", "t", "1", "yes"))

# 3) Simple summary per cell
vdj_meta <- vdj_use %>%
  group_by(cell_id) %>%
  summarise(
    n_vdj_contigs = n(),
    n_productive_contigs = sum(productive_flag, na.rm = TRUE),
    clone_id_productive = paste(unique(clone_id[productive_flag & !is.na(clone_id) & clone_id != ""]), collapse = ";"),
    .groups = "drop"
  ) %>%
  mutate(clone_id_productive = ifelse(clone_id_productive == "", NA, clone_id_productive)) %>%
  distinct(cell_id, .keep_all = TRUE)

# 4) Convert to data.frame BEFORE rownames
vdj_meta_df <- as.data.frame(vdj_meta)
rownames(vdj_meta_df) <- vdj_meta_df$cell_id
vdj_meta_df$cell_id <- NULL

# 5) Align to Seurat cells
vdj_meta_aligned <- vdj_meta_df[colnames(seuratObj), , drop = FALSE]

cat("Exact rowname match:", identical(rownames(vdj_meta_aligned), colnames(seuratObj)), "\n")

# 6) Add to Seurat
seuratObj <- AddMetaData(seuratObj, metadata = vdj_meta_aligned)

# 7) Check
cat("Seurat cells total      :", ncol(seuratObj), "\n")
cat("Cells with VDJ metadata :", sum(!is.na(seuratObj$n_vdj_contigs)), "\n")
cat("Cells without VDJ       :", sum(is.na(seuratObj$n_vdj_contigs)), "\n")
table(is.na(seuratObj$n_vdj_contigs))


##  SCtransforme the data     -------------------------------------------------
seuratObj <- SCTransform(
  seuratObj,
  vars.to.regress = "percent.mt",
  verbose = TRUE,
  ncells = min(5000, ncol(seuratObj)),  
  variable.features.n = 3000
)

DefaultAssay(seuratObj) <- "SCT"

##   Clustering        ----------------------------------------------------------
##  Linear dim reduction      
DefaultAssay(seuratObj) <- "SCT"

seuratObj <- RunPCA(
  seuratObj,
  npcs = 50,
  ndims.print = 1:2,
  nfeatures.print = 10,
  verbose = TRUE
)

##  Elbow plot for number of PCs selection  ------------------------------------
qcDir <- file.path(figuresDir, "QC")
dir.create(qcDir, recursive = TRUE, showWarnings = FALSE)

pdf(
  file = file.path(qcDir, paste0(lib3, "_elbowPlot.pdf")),
  width = 8,
  height = 6,
  title = "Elbow plot"
)

print(ElbowPlot(seuratObj, ndims = 50))

dev.off()

##  non-linear dim reduction (UMAP)  -------------------------------------------
ndims <- 30

seuratObj <- FindNeighbors(seuratObj, dims = 1:ndims)
seuratObj <- FindClusters(seuratObj, resolution = 0.5)
seuratObj <- RunUMAP(seuratObj, dims = 1:ndims)

#enregistrer UMAP
qcDir <- file.path(figuresDir, "QC")
dir.create(qcDir, recursive = TRUE, showWarnings = FALSE)

pdf(file = file.path(qcDir, paste0(lib3, "_UMAP_plots.pdf")), width = 8, height = 6)

print(DimPlot(seuratObj, reduction = "umap", group.by = "seurat_clusters", label = TRUE))
print(DimPlot(seuratObj, reduction = "umap", group.by = "hash.ID"))

seuratObj$VDJ_present <- ifelse(is.na(seuratObj$n_vdj_contigs), "No VDJ", "VDJ")
print(DimPlot(seuratObj, reduction = "umap", group.by = "VDJ_present"))

dev.off()


##  Clustering      -------------------------------------------------------------
seuratObj <- FindNeighbors(seuratObj, dims = 1:ndims, verbose = FALSE)

seuratObj <- FindClusters(
  seuratObj,
  verbose = FALSE,
  resolution = c(seq(0.1, 0.6, by = 0.1), seq(0.7, 1.5, by = 0.2))
)

## ================== CLEAN OBJECT REBUILD ==================

# Extraire les counts RNA
counts <- GetAssayData(seuratObj, assay = "RNA", layer = "counts")
# Recréer un objet propre
seuratObj_clean <- CreateSeuratObject(counts = counts)

# Ajouter metadata (IMPORTANT : aligné automatiquement par noms)
seuratObj_clean <- AddMetaData(seuratObj_clean, metadata = seuratObj@meta.data)

## ================== REFAIRE PIPELINE ==================

# SCTransform
seuratObj_clean <- SCTransform(
  seuratObj_clean,
  vars.to.regress = "percent.mt",
  verbose = FALSE
)

# PCA
seuratObj_clean <- RunPCA(seuratObj_clean, npcs = 50, verbose = FALSE)

# Choisir dimensions raisonnables
ndims <- 30

# Neighbors + clustering
seuratObj_clean <- FindNeighbors(seuratObj_clean, dims = 1:ndims)
seuratObj_clean <- FindClusters(seuratObj_clean, resolution = 0.5)

# UMAP
seuratObj_clean <- RunUMAP(seuratObj_clean, dims = 1:ndims, verbose = FALSE)

## ================== CHECK ==================

cat("Nombre de cellules :", ncol(seuratObj_clean), "\n")
cat("PCA cellules :", ncol(seuratObj_clean[["pca"]]), "\n")
cat("UMAP cellules :", ncol(seuratObj_clean[["umap"]]), "\n")

## ================== PLOTS ==================

DimPlot(seuratObj_clean, reduction = "umap", group.by = "seurat_clusters")
DimPlot(seuratObj_clean, reduction = "umap", group.by = "hash.ID")
DimPlot(seuratObj_clean, reduction = "umap", group.by = "VDJ_present")
##  Cell Cycle      -------------------------------------------------------
# Faire le scoring sur RNA (plus stable pour les signatures de cycle cellulaire)
DefaultAssay(seuratObj) <- "RNA"

seuratObj <- CellCycleScoring(
  seuratObj,
  s.features = Seurat::cc.genes.updated.2019$s.genes,
  g2m.features = Seurat::cc.genes.updated.2019$g2m.genes
)

# Revenir sur SCT
DefaultAssay(seuratObj) <- "SCT"

table(seuratObj$Phase)
head(seuratObj@meta.data[, c("S.Score", "G2M.Score", "Phase")])
# enregistrer UMAP
pdf(file = file.path(qcDir, paste0(lib3, "_UMAP_plots.pdf")), width = 8, height = 6)

print(DimPlot(seuratObj, reduction = "umap", group.by = "seurat_clusters", label = TRUE))
print(DimPlot(seuratObj, reduction = "umap", group.by = "hash.ID"))

seuratObj$VDJ_present <- ifelse(is.na(seuratObj$n_vdj_contigs), "No VDJ", "VDJ")
print(DimPlot(seuratObj, reduction = "umap", group.by = "VDJ_present"))

print(DimPlot(seuratObj, reduction = "umap", group.by = "Phase"))
print(FeaturePlot(seuratObj, features = c("S.Score", "G2M.Score"), reduction = "umap"))

dev.off()

#  Recréer l'assay SCT -----------------------
if (!"percent.mt" %in% colnames(seuratObj@meta.data)) {
  seuratObj <- PercentageFeatureSet(seuratObj, pattern = "^MT-", col.name = "percent.mt", assay = "RNA")
}

seuratObj <- SCTransform(
  seuratObj,
  assay = "RNA",
  new.assay.name = "SCT",
  vars.to.regress = "percent.mt",
  verbose = TRUE,
  ncells = min(5000, ncol(seuratObj)),
  variable.features.n = 3000
)

DefaultAssay(seuratObj) <- "SCT"

# PCA + Neighbors + Clustering (multi-résolutions) + UMAP
seuratObj <- RunPCA(
  seuratObj,
  npcs = 50,
  verbose = FALSE
)

ndims <- 30

seuratObj <- FindNeighbors(
  seuratObj,
  dims = 1:ndims,
  verbose = FALSE
)

# Clustering principal + multi-résolutions
seuratObj <- FindClusters(seuratObj, resolution = 0.5, verbose = FALSE)

seuratObj <- FindClusters(
  seuratObj,
  resolution = c(seq(0.1, 0.6, by = 0.1), seq(0.7, 1.5, by = 0.2)),
  verbose = FALSE
)

seuratObj <- RunUMAP(
  seuratObj,
  dims = 1:ndims,
  verbose = FALSE
)
# Vérifier les colonnes de clustering créées
print(Assays(seuratObj))
print(colnames(seuratObj@meta.data)[grep("snn_res|seurat_clusters", colnames(seuratObj@meta.data))])
print(table(Idents(seuratObj)))
colnames(seuratObj@meta.data)
table(seuratObj$hash.ID)
##  Enrich metadata : Condition + Clonotype   -------------------------------------------


# 5) Ton label clustering personnalisé 
seuratObj$Condition <- seuratObj$hash.ID
seuratObj@meta.data <- seuratObj@meta.data %>%
  mutate(
    Clust_res0.1 = as.factor(paste0(SCT_snn_res.0.1, "_", Condition)),
    Clust_res0.2 = as.factor(paste0(SCT_snn_res.0.2, "_", Condition)),
    Clust_res0.3 = as.factor(paste0(SCT_snn_res.0.3, "_", Condition)),
    Clust_res0.4 = as.factor(paste0(SCT_snn_res.0.4, "_", Condition)),
    Clust_res0.5 = as.factor(paste0(SCT_snn_res.0.5, "_", Condition)),
    Clust_res0.6 = as.factor(paste0(SCT_snn_res.0.6, "_", Condition)),
    Clust_res0.7 = as.factor(paste0(SCT_snn_res.0.7, "_", Condition)),
    Clust_res0.9 = as.factor(paste0(SCT_snn_res.0.9, "_", Condition)),
    Clust_res1.1 = as.factor(paste0(SCT_snn_res.1.1, "_", Condition)),
    Clust_res1.3 = as.factor(paste0(SCT_snn_res.1.3, "_", Condition)),
    Clust_res1.5 = as.factor(paste0(SCT_snn_res.1.5, "_", Condition))
  )
#  label à partir du clustering actif
seuratObj@meta.data <- seuratObj@meta.data %>%
  mutate(
    Clust_current = as.factor(paste0(seurat_clusters, substr(Condition,1,1)))
  )

# Vérifs
table(seuratObj$Condition, useNA = "ifany")

# Repartir de vdj_use ----------------

vdj_use <- vdjData %>%
  filter(cell_id %in% colnames(seuratObj)) %>%
  mutate(productive_flag = trimws(tolower(as.character(productive))) %in% c("true","t","1","yes"))

vdj_meta2 <- vdj_use %>%
  group_by(cell_id) %>%
  summarise(
    v_call_productive = paste(unique(v_call[productive_flag & !is.na(v_call) & v_call != ""]), collapse = ";"),
    .groups = "drop"
  ) %>%
  mutate(v_call_productive = ifelse(v_call_productive == "", NA, v_call_productive))

vdj_meta2_df <- as.data.frame(vdj_meta2)
rownames(vdj_meta2_df) <- vdj_meta2_df$cell_id
vdj_meta2_df$cell_id <- NULL

vdj_meta2_aligned <- vdj_meta2_df[colnames(seuratObj), , drop = FALSE]

seuratObj <- AddMetaData(seuratObj, metadata = vdj_meta2_aligned)

seuratObj$TCR <- seuratObj$v_call_productive

# Garder les versions originales
seuratObj$Clonotype_raw <- seuratObj$clone_id_productive
seuratObj$TCR_raw <- seuratObj$TCR

# Fonction de recodage
group_by_frequency <- function(x) {
  x_chr <- as.character(x)
  tab <- table(x_chr, useNA = "no")
  
  out <- x_chr
  out[is.na(x)] <- "missing"
  
  if (length(tab) > 0) {
    unique_vals <- names(tab)[tab == 1]
    low_vals    <- names(tab)[tab >= 2 & tab < 10]
    
    out[!is.na(x) & x_chr %in% unique_vals] <- "unique"
    out[!is.na(x) & x_chr %in% low_vals]    <- "low_count"
  }
  
  out
}

# Créer nouvelles colonnes groupées
seuratObj$Clonotype_grouped <- group_by_frequency(seuratObj$Clonotype_raw)
seuratObj$TCR_grouped       <- group_by_frequency(seuratObj$TCR_raw)

# Vérifs
table(seuratObj$Clonotype_grouped, useNA = "ifany")
table(seuratObj$TCR_grouped, useNA = "ifany")


## Adjusting clonotypes and TCR annotations
seuratObj@meta.data$TCR <- seuratObj@meta.data$v_call_productive

# ── ÉTAPE 1 : Extraire V-J depuis vdjData ──────────────────────
vdj_tcr <- vdjData %>%
  filter(trimws(tolower(as.character(productive))) %in% c("true", "t", "1", "yes")) %>%
  group_by(cell_id) %>%
  reframe(
    TCR = paste(
      paste(unique(na.omit(v_call)), unique(na.omit(j_call)), sep = "_"),
      collapse = ";"
    )
  ) %>%
  distinct(cell_id, .keep_all = TRUE) %>%  # garder 1 ligne par cellule
  as.data.frame()

rownames(vdj_tcr) <- vdj_tcr$cell_id
vdj_tcr$cell_id <- NULL
# ── ÉTAPE 2 : Ajouter TCR au Seurat object ─────────────────────
seuratObj <- AddMetaData(seuratObj, metadata = vdj_tcr)

# ── ÉTAPE 3 : Créer Clonotype ──────────────────────────────────
seuratObj$Clonotype <- seuratObj$clone_id_productive

# ── ÉTAPE 4 : Vérifier avant la boucle ────────────────────────
head(seuratObj@meta.data[, c("Clonotype", "TCR")])
table(is.na(seuratObj$TCR))

# ── BOUCLE FOR (ton code existant) ────────────────────────────
seuratObj$Clonotype_raw <- seuratObj$Clonotype
seuratObj$TCR_raw       <- seuratObj$TCR

for (c in c("Clonotype", "TCR")) {
  
  tmp <- seuratObj@meta.data %>% 
    count(.data[[c]], name = "n") %>% 
    mutate(
      x = case_when(
        n == 1 ~ "unique",
        n >= 2 & n < 10 ~ "few",
        TRUE ~ as.character(.data[[c]])
      )
    )
  
  unique_cl <- tmp[[c]][tmp$x == "unique"]
  low_cl    <- tmp[[c]][tmp$x == "few"]
  
  seuratObj@meta.data <- seuratObj@meta.data %>% 
    dplyr::mutate(
      !!sym(c) := case_when(
        is.na(.data[[c]]) ~ "missing",
        .data[[c]] %in% unique_cl ~ "unique",
        .data[[c]] %in% low_cl ~ "low_count",
        TRUE ~ as.character(.data[[c]])
      )
    )
}

# Vérifications
table(seuratObj@meta.data$Condition, seuratObj@meta.data$Clonotype, useNA = "ifany")
table(seuratObj@meta.data$Condition, seuratObj@meta.data$TCR, useNA = "ifany")


qcDir <- file.path(figuresDir, "QC")
dir.create(qcDir, recursive = TRUE, showWarnings = FALSE)

pdf(file = file.path(qcDir, paste0(lib3, "_UMAP_clonotype_TCR_summary.pdf")),
    width = 10, height = 8)

# UMAPs principaux
print(DimPlot(seuratObj, reduction = "umap", group.by = "Condition") + ggtitle("Condition"))
print(DimPlot(seuratObj, reduction = "umap", group.by = "Clonotype") + ggtitle("Clonotype (grouped)"))
print(DimPlot(seuratObj, reduction = "umap", group.by = "TCR") + ggtitle("TCR (grouped)"))

# Split par condition
print(DimPlot(seuratObj, reduction = "umap", split.by = "Condition", group.by = "Clonotype") + ggtitle("Clonotype by Condition"))
print(DimPlot(seuratObj, reduction = "umap", split.by = "Condition", group.by = "TCR") + ggtitle("TCR by Condition"))

if ("Clonal_status2" %in% colnames(seuratObj@meta.data)) {
  print(DimPlot(seuratObj, reduction = "umap", group.by = "Clonal_status2") + ggtitle("Clonal_status2"))
  print(DimPlot(seuratObj, reduction = "umap", split.by = "Condition", group.by = "Clonal_status2") + ggtitle("Clonal_status2 by Condition"))
}

dev.off()


#creation de la fonction checkMissingClonoCells -----------------------

checkMissingClonoCells_simple <- function(seuratObj, vdjDataSupp) {
  
  if (!"Clonotype" %in% colnames(seuratObj@meta.data)) {
    stop("Clonotype absent de seuratObj@meta.data")
  }
  if (!"barcode" %in% colnames(vdjDataSupp)) {
    stop("barcode absent de vdjDataSupp (all_contig_annotations.csv)")
  }
  
  missing_cells <- seuratObj@meta.data %>%
    tibble::rownames_to_column("seurat_barcode") %>%
    mutate(cell_id_use = seurat_barcode) %>%
    filter(Clonotype == "missing") %>%
    pull(cell_id_use) %>%
    unique()
  
  vdj_sub <- vdjDataSupp %>%
    filter(barcode %in% missing_cells)
  
  has_chain <- "chain" %in% colnames(vdj_sub)
  has_prod  <- "productive" %in% colnames(vdj_sub)
  
  if (nrow(vdj_sub) == 0) {
    return(data.frame(cell_id = character(0), reason = character(0)))
  }
  
  out <- vdj_sub %>%
    group_by(barcode) %>%
    summarise(
      n_contigs = n(),
      any_productive = if (has_prod) any(tolower(as.character(productive)) %in% c("true","t","1","yes"), na.rm = TRUE) else NA,
      chains_seen = if (has_chain) paste(sort(unique(chain)), collapse = ";") else NA_character_,
      .groups = "drop"
    ) %>%
    mutate(
      reason = case_when(
        !is.na(any_productive) & any_productive ~ "contig_present_productive_but_no_clonotype",
        TRUE ~ "contig_present_no_productive_clonotype"
      )
    ) %>%
    rename(cell_id = barcode)
  
  as.data.frame(out)
}
# check the cells with missing clonotypes --------------------------

missingClonotypesInfo <- checkMissingClonoCells_simple(seuratObj, vdjDataSupp)

write.csv(
  missingClonotypesInfo,
  file = file.path(objectsDir, "missingClonoInfo.csv"),
  row.names = FALSE
)

head(missingClonotypesInfo)
nrow(missingClonotypesInfo)


## annotation  ----------
# Créer TCR_chain avant le mutate
seuratObj$TCR_chain <- case_when(
  grepl("TRAV", seuratObj$TCR) & grepl("TRBV", seuratObj$TCR) ~ "AB",
  grepl("TRBV", seuratObj$TCR) ~ "B",
  grepl("TRAV", seuratObj$TCR) ~ "A",
  TRUE ~ "unknown"
)

# Mutate corrigé
seuratObj@meta.data <- seuratObj@meta.data %>%
  mutate(
    Clonotype = case_when(
      rownames(seuratObj@meta.data) %in% missingClonotypesInfo$cell_id ~ "immature_TCR",
      TRUE ~ Clonotype
    ),
    Clonal_status = case_when(
      Clonotype == "unique"                          ~ "no_blast",
      Clonotype %in% c("missing", "immature_TCR")   ~ "missing_clonotype",
      Clonotype == "low_count"                       ~ "tiny_blast",
      TRUE                                           ~ "blast"
    ),
    VDJ_like_class1 = case_when(
      TCR_chain == "AB" & n_productive_contigs > 0  ~ "late_post_B_checkpoint",
      TCR_chain == "B"  & n_productive_contigs > 0  ~ "early_post_B_checkpoint",
      TCR_chain == "A"  & n_productive_contigs > 0  ~ "only_A_chain",
      TRUE                                           ~ "uncertain"
    ),
    VDJ_like_class2 = case_when(
      Clonal_status == "tiny_blast"                                                         ~ "tiny_blast",
      Clonal_status == "no_blast"                                                           ~ "no_blast",
      Clonal_status == "missing_clonotype" & rownames(seuratObj@meta.data) %in% missingClonotypesInfo$cell_id ~ "pre_B_checkpoint_candidate",
      TRUE                                                                                  ~ VDJ_like_class1
    )
  )

#plots metadata ---------------------

## Plot clusters---------------
# Compatibilité avec ton ancien style
lib <- lib3

# Assurer que le dossier existe
dir.create(figuresDir, recursive = TRUE, showWarnings = FALSE)

# Palette (si pas déjà définie)
if (!exists("colPaletteGenerator")) {
  colPaletteGenerator <- function(n) {
    grDevices::colorRampPalette(c(
      "#1b9e77", "#d95f02", "#7570b3", "#e7298a",
      "#66a61e", "#e6ab02", "#a6761d", "#666666"
    ))(n)
  }
}

# Thème UMAP (si pas déjà défini)
if (!exists("My_umap_theme")) {
  My_umap_theme <- function() {
    theme_bw() +
      theme(
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "bottom",
        strip.background = element_rect(fill = "grey95"),
        strip.text = element_text(face = "bold")
      )
  }
}

# Vérifier UMAP
if (!"umap" %in% names(seuratObj@reductions)) {
  stop("UMAP absent dans seuratObj. Lance RunUMAP() avant ce bloc.")
}

# Ajouter umap_1 / umap_2 à meta.data si absents
if (!all(c("umap_1", "umap_2") %in% colnames(seuratObj@meta.data))) {
  umap_df <- as.data.frame(Embeddings(seuratObj, "umap"))
  colnames(umap_df) <- c("umap_1", "umap_2")
  
  seuratObj@meta.data <- seuratObj@meta.data %>%
    tibble::rownames_to_column("cell_barcode_tmp") %>%
    left_join(
      umap_df %>% tibble::rownames_to_column("cell_barcode_tmp"),
      by = "cell_barcode_tmp"
    ) %>%
    tibble::column_to_rownames("cell_barcode_tmp")
}

# Vérifs colonnes nécessaires
if (!"Condition" %in% colnames(seuratObj@meta.data)) {
  stop("La colonne 'Condition' est absente de seuratObj@meta.data.")
}

# Colonnes de clustering SCT
sctRes <- grep("^SCT_snn_res\\.", colnames(seuratObj@meta.data), value = TRUE)

if (length(sctRes) == 0) {
  warning("Aucune colonne SCT_snn_res.* trouvée. Le PDF contiendra seulement le plot Condition.")
}


## Export PDF

pdf(
  file = file.path(figuresDir, paste0(lib, "_clusters.pdf")),
  width = 9,
  height = 6,
  title = "Clusters - no.filt"
)

# 1) UMAP global par Condition
n_cond <- length(unique(seuratObj@meta.data$Condition))

print(
  DimPlot(
    seuratObj,
    reduction = "umap",
    group.by = "Condition",
    cols = colPaletteGenerator(n_cond),
    pt.size = 0.7,
    alpha = 0.8
  ) +
    labs(
      title = "General conditions",
      colour = "Conditions"
    ) +
    My_umap_theme()
)

# 2) UMAP facetté pour chaque résolution SCT
for (i in sctRes) {
  
  ncolors <- length(unique(seuratObj@meta.data[[i]]))
  
  p <- seuratObj@meta.data %>%
    ggplot(aes(x = umap_1, y = umap_2)) +
    geom_point(
      aes(colour = .data[[i]]),
      size = 0.7,
      alpha = 0.8
    ) +
    facet_grid(~ Condition) +
    My_umap_theme() +
    ggtitle(i) +
    scale_colour_manual(values = colPaletteGenerator(ncolors)) +
    guides(
      colour = guide_legend(
        position = "bottom",
        title = "",
        nrow = ifelse(ncolors <= 16, 1, 2),
        override.aes = list(size = 3)
      )
    ) +
    theme(
      plot.margin = margin(t = 0.05, b = 0.05, r = 0.05, l = 0.05, unit = "in"),
      plot.title = element_text(margin = margin(b = 0.1, unit = "in"))
    )
  
  print(p)
}

dev.off()

cat("PDF généré :", file.path(figuresDir, paste0(lib, "_clusters.pdf")), "\n")

# resolution 0.5------
Idents(seuratObj) <- "SCT_snn_res.0.5"

# mettre aussi seurat_clusters = res 0.5
seuratObj$seurat_clusters <- seuratObj$SCT_snn_res.0.5

# Vérif
table(Idents(seuratObj))


# plot cell cycle ----------

## 1) Calculer le Cell Cycle Scoring (si Phase absente)
if (!"Phase" %in% colnames(seuratObj@meta.data)) {
  cat("Phase absente -> lancement de CellCycleScoring...\n")
  
  # Le scoring se fait sur RNA
  DefaultAssay(seuratObj) <- "RNA"
  
  seuratObj <- CellCycleScoring(
    seuratObj,
    s.features = Seurat::cc.genes.updated.2019$s.genes,
    g2m.features = Seurat::cc.genes.updated.2019$g2m.genes,
    set.ident = FALSE
  )
  
  # Revenir à SCT pour la suite
  if ("SCT" %in% Assays(seuratObj)) {
    DefaultAssay(seuratObj) <- "SCT"
  } else {
    DefaultAssay(seuratObj) <- "RNA"
    warning("Assay SCT absent. Je reste sur RNA.")
  }
} else {
  cat("Phase déjà présente, pas de recalcul.\n")
}

# Vérif
cat("Répartition des phases:\n")
print(table(seuratObj$Phase, useNA = "ifany"))
print(head(seuratObj@meta.data[, c("S.Score", "G2M.Score", "Phase")]))

## 2) Préparation plotting


# Compatibilité avec ton style
lib <- lib3

# Dossier de sortie
dir.create(figuresDir, recursive = TRUE, showWarnings = FALSE)

# Palette (si pas déjà définie)
if (!exists("colPaletteGenerator")) {
  colPaletteGenerator <- function(n) {
    grDevices::colorRampPalette(c(
      "#1b9e77", "#d95f02", "#7570b3", "#e7298a",
      "#66a61e", "#e6ab02", "#a6761d", "#666666"
    ))(n)
  }
}

# Thème UMAP (si pas déjà défini)
if (!exists("My_umap_theme")) {
  My_umap_theme <- function() {
    theme_bw() +
      theme(
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "bottom",
        strip.background = element_rect(fill = "grey95"),
        strip.text = element_text(face = "bold")
      )
  }
}

# Vérifs
if (!"umap" %in% names(seuratObj@reductions)) {
  stop("UMAP absent dans seuratObj. Lance RunUMAP() avant ce bloc.")
}
if (!"Condition" %in% colnames(seuratObj@meta.data)) {
  stop("La colonne 'Condition' est absente de seuratObj@meta.data.")
}
if (!"Phase" %in% colnames(seuratObj@meta.data)) {
  stop("La colonne 'Phase' est toujours absente après CellCycleScoring().")
}

# Ajouter coordonnées UMAP à meta.data si absentes
if (!all(c("UMAP_1_plot", "UMAP_2_plot") %in% colnames(seuratObj@meta.data))) {
  umap_df <- as.data.frame(Embeddings(seuratObj, "umap"))
  colnames(umap_df) <- c("UMAP_1_plot", "UMAP_2_plot")
  
  seuratObj@meta.data <- seuratObj@meta.data %>%
    tibble::rownames_to_column("cell_barcode_tmp") %>%
    left_join(
      umap_df %>% tibble::rownames_to_column("cell_barcode_tmp"),
      by = "cell_barcode_tmp"
    ) %>%
    tibble::column_to_rownames("cell_barcode_tmp")
}

##Export PDF Cell Cycle


outfile <- file.path(figuresDir, paste0(lib, "_cellCycle.pdf"))

# Si un ancien PDF vide existe, on le remplace
if (file.exists(outfile)) {
  try(file.remove(outfile), silent = TRUE)
}

pdf(
  file = outfile,
  width = 9,
  height = 6,
  title = "Cell Cycle - no.filt"
)

# Plot 1 : UMAP global par phase 
n_phase <- length(unique(as.character(seuratObj@meta.data$Phase)))

p1 <- DimPlot(
  seuratObj,
  reduction = "umap",
  group.by = "Phase",
  cols = colPaletteGenerator(n_phase),
  pt.size = 0.7,
  alpha = 0.8
) +
  labs(
    title = "Cell Cycle exploration",
    colour = "Phases"
  ) +
  My_umap_theme()

print(p1)

# Plot 2 : UMAP facetté par Condition 
ncolors <- length(unique(as.character(seuratObj@meta.data$Phase)))

p2 <- seuratObj@meta.data %>%
  ggplot(aes(x = UMAP_1_plot, y = UMAP_2_plot)) +
  geom_point(
    aes(colour = Phase),
    size = 0.7,
    alpha = 0.8
  ) +
  facet_grid(~ Condition) +
  My_umap_theme() +
  ggtitle("Cell cycle phases") +
  scale_colour_manual(values = colPaletteGenerator(ncolors)) +
  guides(
    colour = guide_legend(
      position = "bottom",
      title = "",
      nrow = ifelse(ncolors <= 16, 1, 2),
      override.aes = list(size = 3)
    )
  ) +
  theme(
    plot.margin = margin(t = 0.05, b = 0.05, r = 0.05, l = 0.05, unit = "in"),
    plot.title = element_text(margin = margin(b = 0.1, unit = "in"))
  )

print(p2)

dev.off()


# plot proportions cell cycle ----------------

# Dossiers de sortie
qcDir <- file.path(figuresDir, "QC")
dir.create(qcDir, recursive = TRUE, showWarnings = FALSE)

tabDir <- file.path(objectsDir, "tables")
dir.create(tabDir, recursive = TRUE, showWarnings = FALSE)

# Vérifs
if (!"Phase" %in% colnames(seuratObj@meta.data)) {
  stop("La colonne 'Phase' est absente. Lance CellCycleScoring() avant ce bloc.")
}
if (!"Condition" %in% colnames(seuratObj@meta.data)) {
  stop("La colonne 'Condition' est absente de seuratObj@meta.data.")
}

# 1) Tables

tab_phase <- table(seuratObj$Phase, useNA = "ifany")
tab_cond_phase <- table(seuratObj$Condition, seuratObj$Phase, useNA = "ifany")
prop_cond_phase <- prop.table(tab_cond_phase, margin = 1)

# Affichage console
cat("\n=== Table Phase ===\n")
print(tab_phase)

cat("\n=== Table Condition x Phase ===\n")
print(tab_cond_phase)

cat("\n=== Proportions Phase par Condition ===\n")
print(round(prop_cond_phase, 3))


# 2) Sauvegarde CSV
write.csv(
  as.data.frame(tab_phase),
  file = file.path(tabDir, paste0(lib3, "_table_Phase.csv")),
  row.names = FALSE
)

write.csv(
  as.data.frame.matrix(tab_cond_phase),
  file = file.path(tabDir, paste0(lib3, "_table_Condition_by_Phase_counts.csv"))
)

write.csv(
  round(as.data.frame.matrix(prop_cond_phase), 4),
  file = file.path(tabDir, paste0(lib3, "_table_Condition_by_Phase_proportions.csv"))
)

# 3) Préparer dataframes pour plots
df_phase <- as.data.frame(tab_phase)
colnames(df_phase) <- c("Phase", "Count")

df_cond_phase <- as.data.frame(tab_cond_phase)
colnames(df_cond_phase) <- c("Condition", "Phase", "Count")

df_prop_cond_phase <- as.data.frame(prop_cond_phase)
colnames(df_prop_cond_phase) <- c("Condition", "Phase", "Proportion")

# Optionnel : pourcentages
df_prop_cond_phase <- df_prop_cond_phase %>%
  mutate(Percent = round(Proportion * 100, 1))

# 4) Sauvegarde des figures en PDF
pdf(
  file = file.path(qcDir, paste0(lib3, "_Phase_Condition_tables_plots.pdf")),
  width = 9,
  height = 6
)

# Plot A : distribution globale des phases
print(
  ggplot(df_phase, aes(x = Phase, y = Count, fill = Phase)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = Count), vjust = -0.2, fontface = "bold") +
    labs(
      title = "Distribution globale des phases du cycle cellulaire",
      x = "Phase",
      y = "Nombre de cellules"
    ) +
    theme_bw() +
    theme(legend.position = "none")
)

# Plot B : counts par condition (barres côte à côte)
print(
  ggplot(df_cond_phase, aes(x = Condition, y = Count, fill = Phase)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    geom_text(
      aes(label = Count),
      position = position_dodge(width = 0.8),
      vjust = -0.2,
      size = 3
    ) +
    labs(
      title = "Nombre de cellules par phase et par condition",
      x = "Condition",
      y = "Nombre de cellules"
    ) +
    theme_bw()
)

# Plot C : proportions par condition (stacked 100%)
print(
  ggplot(df_prop_cond_phase, aes(x = Condition, y = Proportion, fill = Phase)) +
    geom_col(position = "fill", width = 0.7) +
    geom_text(
      aes(label = paste0(Percent, "%")),
      position = position_fill(vjust = 0.5),
      color = "white",
      fontface = "bold",
      size = 3
    ) +
    scale_y_continuous(labels = function(x) paste0(round(x * 100), "%")) +
    labs(
      title = "Proportions des phases par condition",
      x = "Condition",
      y = "Proportion"
    ) +
    theme_bw()
)

dev.off()



## Plot clonotypes -----------------

# Compatibilité style ancien code
lib <- lib3

# Assurer dossier de sortie
dir.create(figuresDir, recursive = TRUE, showWarnings = FALSE)

# Palette (si pas déjà définie)
if (!exists("colPaletteGenerator")) {
  colPaletteGenerator <- function(n) {
    grDevices::colorRampPalette(c(
      "#1b9e77", "#d95f02", "#7570b3", "#e7298a",
      "#66a61e", "#e6ab02", "#a6761d", "#666666"
    ))(n)
  }
}

# Thème UMAP (si pas déjà défini)
if (!exists("My_umap_theme")) {
  My_umap_theme <- function() {
    theme_bw() +
      theme(
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "bottom",
        strip.background = element_rect(fill = "grey95"),
        strip.text = element_text(face = "bold")
      )
  }
}

# Thème barplot (si pas déjà défini)
if (!exists("My_bplot_theme")) {
  My_bplot_theme <- function() {
    theme_bw() +
      theme(
        panel.grid.major.x = element_blank(),
        legend.position = "bottom",
        plot.title = element_text(face = "bold"),
        axis.text.x = element_text(face = "bold")
      )
  }
}

# RotatedAxis (si pas déjà défini)
if (!exists("RotatedAxis")) {
  RotatedAxis <- function() {
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }
}

# Vérifs essentielles
needed_cols <- c("Clonotype", "Condition")
missing_needed <- setdiff(needed_cols, colnames(seuratObj@meta.data))
if (length(missing_needed) > 0) {
  stop("Colonnes manquantes dans meta.data : ", paste(missing_needed, collapse = ", "))
}

if (!"umap" %in% names(seuratObj@reductions)) {
  stop("UMAP absent dans seuratObj. Lance RunUMAP() avant ce bloc.")
}

# Ajouter coordonnées UMAP de plot si absentes
if (!all(c("UMAP_1_plot", "UMAP_2_plot") %in% colnames(seuratObj@meta.data))) {
  umap_df <- as.data.frame(Embeddings(seuratObj, "umap"))
  colnames(umap_df) <- c("UMAP_1_plot", "UMAP_2_plot")
  
  seuratObj@meta.data <- seuratObj@meta.data %>%
    rownames_to_column("cell_barcode_tmp") %>%
    left_join(
      umap_df %>% rownames_to_column("cell_barcode_tmp"),
      by = "cell_barcode_tmp"
    ) %>%
    column_to_rownames("cell_barcode_tmp")
}

# Harmoniser Clonotype en character (évite certains bugs de facteur)
seuratObj@meta.data$Clonotype <- as.character(seuratObj@meta.data$Clonotype)
seuratObj@meta.data$Condition <- as.character(seuratObj@meta.data$Condition)

# Nombre de couleurs
ncolors <- length(unique(seuratObj@meta.data$Clonotype))

# Fichier de sortie
outfile <- file.path(figuresDir, paste0(lib, "_clonotypes.pdf"))
if (file.exists(outfile)) {
  try(file.remove(outfile), silent = TRUE)
}

pdf(
  file = outfile,
  width = 9,
  height = 6,
  title = "Clonotypes - no.filt"
)


## 01) UMAP global par Clonotype

p1 <- DimPlot(
  seuratObj,
  reduction = "umap",
  group.by = "Clonotype",
  cols = colPaletteGenerator(ncolors),
  pt.size = 0.7,
  alpha = 0.8
) +
  labs(title = "Clonotypes", colour = "") +
  My_umap_theme()

print(p1)

## 02) Clone size mapping (log10)
p2 <- seuratObj@meta.data %>%
  group_by(Clonotype) %>%
  mutate(clone_size = n()) %>%
  ungroup() %>%
  ggplot(aes(x = UMAP_1_plot, y = UMAP_2_plot)) +
  geom_point(aes(colour = log10(clone_size)), size = 0.7, alpha = 0.8) +
  labs(
    title = "Clones size mapping in log10",
    colour = ""
  ) +
  My_umap_theme() +
  theme(
    plot.margin = margin(t = 0.05, b = 0.05, r = 0.05, l = 0.05, unit = "in"),
    plot.title = element_text(margin = margin(b = 0.1, unit = "in"))
  ) +
  scale_colour_gradient(low = "#bfd", high = "#064")

print(p2)

## 03) UMAP Clonotype facetté par Condition

p3 <- seuratObj@meta.data %>%
  ggplot(aes(x = UMAP_1_plot, y = UMAP_2_plot)) +
  geom_point(
    aes(colour = Clonotype),
    size = 0.7,
    alpha = 0.8
  ) +
  facet_grid(~ Condition) +
  My_umap_theme() +
  ggtitle("Clonotypes") +
  scale_colour_manual(values = colPaletteGenerator(ncolors)) +
  guides(
    colour = guide_legend(
      position = "bottom",
      title = "",
      nrow = ifelse(ncolors <= 16, 1, 2),
      override.aes = list(size = 3)
    )
  ) +
  theme(
    plot.margin = margin(t = 0.05, b = 0.05, r = 0.05, l = 0.05, unit = "in"),
    plot.title = element_text(margin = margin(b = 0.1, unit = "in"))
  )

print(p3)


## 04) Distribution des clonotypes par Condition (stacked bar)
t <- 150

p4 <- seuratObj@meta.data %>%
  count(Condition, Clonotype) %>%
  ggplot(aes(x = Condition, y = n, fill = Clonotype)) +
  geom_bar(
    position = "stack",
    stat = "identity",
    width = 0.8
  ) +
  geom_text(
    aes(label = ifelse(n >= t, n, "")),
    position = position_stack(vjust = 0.5, reverse = FALSE),
    colour = "white",
    size = 4,
    fontface = "bold"
  ) +
  scale_fill_manual(values = colPaletteGenerator(ncolors)) +
  labs(
    title = "Clonotypes distribution",
    y = "",
    x = "",
    caption = paste0("Label count min threshold : ", t)
  ) +
  My_bplot_theme()

print(p4)


## 05) Distribution clonotypes / VDJ_like_class2 (si disponible)
if ("VDJ_like_class2" %in% colnames(seuratObj@meta.data)) {
  
  t <- 150
  
  p5 <- seuratObj@meta.data %>%
    count(VDJ_like_class2, Clonotype, Condition) %>%
    ggplot(aes(x = VDJ_like_class2, y = n, fill = Clonotype)) +
    geom_bar(
      position = "stack",
      stat = "identity",
      width = 0.8
    ) +
    geom_text(
      aes(label = ifelse(n >= t, n, "")),
      position = position_stack(vjust = 0.5, reverse = FALSE),
      colour = "white",
      size = 3,
      fontface = "bold"
    ) +
    scale_fill_manual(values = colPaletteGenerator(ncolors)) +
    labs(
      title = "Clonotypes distribution / VDJ-like classification",
      y = "",
      x = "",
      caption = paste0("Label count min threshold : ", t)
    ) +
    My_bplot_theme() +
    RotatedAxis() +
    facet_grid(~ Condition)
  
  print(p5)
  
} else {
  message("VDJ_like_class2 absent : plot 05 ignoré.")
}

dev.off()


## Plot TCR (V gene) ------------------

# Compatibilité avec ton ancien style
lib <- lib3

# Assurer dossier de sortie
dir.create(figuresDir, recursive = TRUE, showWarnings = FALSE)

# Palette (si pas déjà définie)
if (!exists("colPaletteGenerator")) {
  colPaletteGenerator <- function(n) {
    grDevices::colorRampPalette(c(
      "#1b9e77", "#d95f02", "#7570b3", "#e7298a",
      "#66a61e", "#e6ab02", "#a6761d", "#666666"
    ))(n)
  }
}

# Thème UMAP (si pas déjà défini)
if (!exists("My_umap_theme")) {
  My_umap_theme <- function() {
    theme_bw() +
      theme(
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "bottom",
        strip.background = element_rect(fill = "grey95"),
        strip.text = element_text(face = "bold")
      )
  }
}

# Thème barplot (si pas déjà défini)
if (!exists("My_bplot_theme")) {
  My_bplot_theme <- function() {
    theme_bw() +
      theme(
        panel.grid.major.x = element_blank(),
        legend.position = "bottom",
        plot.title = element_text(face = "bold"),
        axis.text.x = element_text(face = "bold")
      )
  }
}

# RotatedAxis (si pas déjà défini)
if (!exists("RotatedAxis")) {
  RotatedAxis <- function() {
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }
}

# Vérifs essentielles
needed_cols <- c("TCR", "Condition", "Clonotype")
missing_needed <- setdiff(needed_cols, colnames(seuratObj@meta.data))
if (length(missing_needed) > 0) {
  stop("Colonnes manquantes dans meta.data : ", paste(missing_needed, collapse = ", "))
}

if (!"umap" %in% names(seuratObj@reductions)) {
  stop("UMAP absent dans seuratObj. Lance RunUMAP() avant ce bloc.")
}

# Ajouter coordonnées UMAP (sans conflit avec embeddings)
if (!all(c("UMAP_1_plot", "UMAP_2_plot") %in% colnames(seuratObj@meta.data))) {
  umap_df <- as.data.frame(Embeddings(seuratObj, "umap"))
  colnames(umap_df) <- c("UMAP_1_plot", "UMAP_2_plot")
  
  seuratObj@meta.data <- seuratObj@meta.data %>%
    rownames_to_column("cell_barcode_tmp") %>%
    left_join(
      umap_df %>% rownames_to_column("cell_barcode_tmp"),
      by = "cell_barcode_tmp"
    ) %>%
    column_to_rownames("cell_barcode_tmp")
}

# Harmoniser types
seuratObj@meta.data$TCR <- as.character(seuratObj@meta.data$TCR)
seuratObj@meta.data$Condition <- as.character(seuratObj@meta.data$Condition)
seuratObj@meta.data$Clonotype <- as.character(seuratObj@meta.data$Clonotype)

# Nombre de couleurs
ncolors <- length(unique(seuratObj@meta.data$TCR))

# Sortie PDF
outfile <- file.path(figuresDir, paste0(lib, "_TCR.pdf"))
if (file.exists(outfile)) {
  try(file.remove(outfile), silent = TRUE)
}

pdf(
  file = outfile,
  width = 9,
  height = 6,
  title = "TCR - no.filt"
)


## 01) UMAP global par TCR
p1 <- DimPlot(
  seuratObj,
  reduction = "umap",
  group.by = "TCR",
  cols = colPaletteGenerator(ncolors),
  pt.size = 0.7,
  alpha = 0.8
) +
  labs(title = "TCR - V gene", colour = "") +
  My_umap_theme()

print(p1)

## 02) UMAP TCR facetté par Condition
legend_rows <- ifelse(ncolors <= 6, 1, ifelse(ncolors <= 12, 2, 3))

p2 <- seuratObj@meta.data %>%
  ggplot(aes(x = UMAP_1_plot, y = UMAP_2_plot)) +
  geom_point(
    aes(colour = TCR),
    size = 0.7,
    alpha = 0.8
  ) +
  facet_grid(~ Condition) +
  My_umap_theme() +
  ggtitle("TCR - V gene") +
  scale_colour_manual(values = colPaletteGenerator(ncolors)) +
  guides(
    colour = guide_legend(
      position = "bottom",
      title = "",
      nrow = legend_rows,
      override.aes = list(size = 3)
    )
  ) +
  theme(
    plot.margin = margin(t = 0.05, b = 0.05, r = 0.05, l = 0.05, unit = "in"),
    plot.title = element_text(margin = margin(b = 0.1, unit = "in"))
  )

print(p2)


## 03) Distribution TCR par Condition (stacked bar)
t <- 150

p3 <- seuratObj@meta.data %>%
  count(Condition, TCR) %>%
  ggplot(aes(x = Condition, y = n, fill = TCR)) +
  geom_bar(
    position = "stack",
    stat = "identity",
    width = 0.8
  ) +
  geom_text(
    aes(label = ifelse(n >= t, n, "")),
    position = position_stack(vjust = 0.5, reverse = FALSE),
    colour = "white",
    size = 4,
    fontface = "bold"
  ) +
  scale_fill_manual(values = colPaletteGenerator(ncolors)) +
  labs(
    title = "TCR Vgenes distribution",
    y = "",
    x = "",
    caption = paste0("Label count min threshold : ", t)
  ) +
  My_bplot_theme()

print(p3)


## 04) Distribution TCR / VDJ_like_class2 (si dispo)
if ("VDJ_like_class2" %in% colnames(seuratObj@meta.data)) {
  p4 <- seuratObj@meta.data %>%
    count(VDJ_like_class2, TCR, Condition) %>%
    ggplot(aes(x = VDJ_like_class2, y = n, fill = TCR)) +
    geom_bar(
      position = "stack",
      stat = "identity",
      width = 0.8
    ) +
    geom_text(
      aes(label = ifelse(n >= t, n, "")),
      position = position_stack(vjust = 0.5, reverse = FALSE),
      colour = "white",
      size = 3,
      fontface = "bold"
    ) +
    scale_fill_manual(values = colPaletteGenerator(ncolors)) +
    labs(
      title = "TCR Vgenes distribution / VDJ-like classification",
      y = "",
      x = "",
      caption = paste0("Label count min threshold : ", t)
    ) +
    My_bplot_theme() +
    RotatedAxis() +
    facet_grid(~ Condition)
  
  print(p4)
} else {
  message("VDJ_like_class2 absent : plot 04 ignoré.")
}

## 05) Distribution TCR / Clonotype
p5 <- seuratObj@meta.data %>%
  count(Clonotype, TCR, Condition) %>%
  ggplot(aes(x = Clonotype, y = n, fill = TCR)) +
  geom_bar(
    position = "stack",
    stat = "identity",
    width = 0.8
  ) +
  geom_text(
    aes(label = ifelse(n >= t, n, "")),
    position = position_stack(vjust = 0.5, reverse = FALSE),
    colour = "white",
    size = 3,
    fontface = "bold"
  ) +
  scale_fill_manual(values = colPaletteGenerator(ncolors)) +
  labs(
    title = "TCR Vgenes distribution / Clonotypes",
    y = "",
    x = "",
    caption = paste0("Label count min threshold : ", t)
  ) +
  My_bplot_theme() +
  RotatedAxis() +
  facet_grid(~ Condition)

print(p5)

dev.off()


## Plot T-ALL status ---------------------------------

# Compatibilité avec ton ancien style
lib <- lib3

# Assurer dossier de sortie
dir.create(figuresDir, recursive = TRUE, showWarnings = FALSE)

# Palette (si pas déjà définie)
if (!exists("colPaletteGenerator")) {
  colPaletteGenerator <- function(n) {
    grDevices::colorRampPalette(c(
      "#1b9e77", "#d95f02", "#7570b3", "#e7298a",
      "#66a61e", "#e6ab02", "#a6761d", "#666666"
    ))(n)
  }
}

# Thème UMAP (si pas déjà défini)
if (!exists("My_umap_theme")) {
  My_umap_theme <- function() {
    theme_bw() +
      theme(
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "bottom",
        strip.background = element_rect(fill = "grey95"),
        strip.text = element_text(face = "bold")
      )
  }
}

# Thème barplot (si pas déjà défini)
if (!exists("My_bplot_theme")) {
  My_bplot_theme <- function() {
    theme_bw() +
      theme(
        panel.grid.major.x = element_blank(),
        legend.position = "bottom",
        plot.title = element_text(face = "bold"),
        axis.text.x = element_text(face = "bold")
      )
  }
}

# Vérifs essentielles
needed_cols <- c("Condition", "Clonal_status", "TCR_chain", "VDJ_like_class2")
missing_needed <- setdiff(needed_cols, colnames(seuratObj@meta.data))
if (length(missing_needed) > 0) {
  stop("Colonnes manquantes dans meta.data : ", paste(missing_needed, collapse = ", "))
}

if (!"umap" %in% names(seuratObj@reductions)) {
  stop("UMAP absent dans seuratObj. Lance RunUMAP() avant ce bloc.")
}

# Ajouter coordonnées UMAP de plot si absentes (sans conflit avec Embeddings)
if (!all(c("UMAP_1_plot", "UMAP_2_plot") %in% colnames(seuratObj@meta.data))) {
  umap_df <- as.data.frame(Embeddings(seuratObj, "umap"))
  colnames(umap_df) <- c("UMAP_1_plot", "UMAP_2_plot")
  
  seuratObj@meta.data <- seuratObj@meta.data %>%
    rownames_to_column("cell_barcode_tmp") %>%
    left_join(
      umap_df %>% rownames_to_column("cell_barcode_tmp"),
      by = "cell_barcode_tmp"
    ) %>%
    column_to_rownames("cell_barcode_tmp")
}

# Harmoniser types (évite bugs de facteurs)
seuratObj@meta.data$Condition       <- as.character(seuratObj@meta.data$Condition)
seuratObj@meta.data$Clonal_status   <- as.character(seuratObj@meta.data$Clonal_status)
seuratObj@meta.data$TCR_chain       <- as.character(seuratObj@meta.data$TCR_chain)
seuratObj@meta.data$VDJ_like_class2 <- as.character(seuratObj@meta.data$VDJ_like_class2)

# Sortie PDF
outfile <- file.path(figuresDir, paste0(lib, "_TALL_status_v1.pdf"))
if (file.exists(outfile)) {
  try(file.remove(outfile), silent = TRUE)
}

pdf(
  file = outfile,
  width = 9,
  height = 6,
  title = "TALL status - no.filt"
)

## 01) UMAP - Clonal_status
ncolors1 <- length(unique(seuratObj@meta.data$Clonal_status))

p1 <- DimPlot(
  seuratObj,
  reduction = "umap",
  group.by = "Clonal_status",
  cols = colPaletteGenerator(ncolors1),
  pt.size = 0.7,
  alpha = 0.8
) +
  labs(title = "Clonal status", colour = "") +
  My_umap_theme()

print(p1)


## 02) UMAP - TCR_chain
ncolors2 <- length(unique(seuratObj@meta.data$TCR_chain))

p2 <- DimPlot(
  seuratObj,
  reduction = "umap",
  group.by = "TCR_chain",
  cols = colPaletteGenerator(ncolors2),
  pt.size = 0.7,
  alpha = 0.8
) +
  labs(title = "TCR chain", colour = "") +
  My_umap_theme()

print(p2)

## 03) UMAP - VDJ_like_class2
ncolors3 <- length(unique(seuratObj@meta.data$VDJ_like_class2))

p3 <- DimPlot(
  seuratObj,
  reduction = "umap",
  group.by = "VDJ_like_class2",
  cols = colPaletteGenerator(ncolors3),
  pt.size = 0.7,
  alpha = 0.8
) +
  labs(title = "TALL VDJ-like classification", colour = "") +
  My_umap_theme()

print(p3)


## 04) UMAP facetté par Condition - VDJ_like_class2
legend_rows <- ifelse(ncolors3 <= 6, 1, ifelse(ncolors3 <= 12, 2, 3))

p4 <- seuratObj@meta.data %>%
  ggplot(aes(x = UMAP_1_plot, y = UMAP_2_plot)) +
  geom_point(
    aes(colour = VDJ_like_class2),
    size = 0.7,
    alpha = 0.8
  ) +
  facet_grid(~ Condition) +
  My_umap_theme() +
  ggtitle("TALL VDJ-like classification") +
  scale_colour_manual(values = colPaletteGenerator(ncolors3)) +
  guides(
    colour = guide_legend(
      position = "bottom",
      title = "",
      nrow = legend_rows,
      override.aes = list(size = 3)
    )
  ) +
  theme(
    plot.margin = margin(t = 0.05, b = 0.05, r = 0.05, l = 0.05, unit = "in"),
    plot.title = element_text(margin = margin(b = 0.1, unit = "in"))
  )

print(p4)

## 05) Distribution VDJ_like_class2 par Condition
t <- 150

p5 <- seuratObj@meta.data %>%
  count(Condition, VDJ_like_class2) %>%
  ggplot(aes(x = Condition, y = n, fill = VDJ_like_class2)) +
  geom_bar(
    position = "stack",
    stat = "identity",
    width = 0.8
  ) +
  geom_text(
    aes(label = ifelse(n >= t, n, "")),
    position = position_stack(vjust = 0.5, reverse = FALSE),
    colour = "white",
    size = 4,
    fontface = "bold"
  ) +
  scale_fill_manual(values = colPaletteGenerator(ncolors3)) +
  labs(
    title = "TALL VDJ-like classification distribution",
    y = "",
    x = "",
    caption = paste0("Label count min threshold : ", t)
  ) +
  My_bplot_theme()

print(p5)

dev.off()

#verifier proportions de late_post_B_checkpoint
tab <- table(seuratObj$Condition, seuratObj$VDJ_like_class2, useNA="ifany")
tab
round(prop.table(tab, margin=1), 3)



## Plot Markers v1 ----------------------------------
lib <- lib3

# Assurer dossier de sortie
dir.create(figuresDir, recursive = TRUE, showWarnings = FALSE)

# Packages utiles
if (!requireNamespace("cowplot", quietly = TRUE)) {
  install.packages("cowplot")
}
library(cowplot)

# Palette (si pas déjà définie)
if (!exists("colPaletteGenerator")) {
  colPaletteGenerator <- function(n) {
    grDevices::colorRampPalette(c(
      "#1b9e77", "#d95f02", "#7570b3", "#e7298a",
      "#66a61e", "#e6ab02", "#a6761d", "#666666"
    ))(n)
  }
}

# Thème UMAP (si pas déjà défini)
if (!exists("My_umap_theme")) {
  My_umap_theme <- function() {
    theme_bw() +
      theme(
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "bottom",
        strip.background = element_rect(fill = "grey95"),
        strip.text = element_text(face = "bold")
      )
  }
}
# 1) Définir le panel de marqueurs
markers_panel <- list(
  "Mature T"         = c("TRAC","CD3D","CD3E","IL7R","CCR7"),
  "T-ALL/immature"   = c("CD7","CD99","HES1","LMO2","MKI67","TOP2A"),
  "Myeloid"          = c("LYZ","S100A8","FCGR3A"),
  "B cells"          = c("MS4A1","CD79A","CD74"),
  "NK/cytotoxic"     = c("NKG7","GNLY","PRF1")
)

# 2) Groupings à tester
markersBy <- c(
  "VDJ_like_class2", "SCT_snn_res.0.1", "SCT_snn_res.0.3", "SCT_snn_res.0.5",
  "SCT_snn_res.1.1"
)

# Garder uniquement les colonnes qui existent vraiment
markersBy <- markersBy[markersBy %in% colnames(seuratObj@meta.data)]

if (length(markersBy) == 0) {
  stop("Aucune colonne de 'markersBy' n'existe dans seuratObj@meta.data.")
}

cat("Groupings utilisés pour les markers:\n")
print(markersBy)


# 3) Préparer les gènes présents dans l'objet
# On prend les features de l'assay actif (SCT si actif, sinon RNA)
assay_used <- DefaultAssay(seuratObj)
cat("Assay actif pour DotPlot :", assay_used, "\n")

all_features <- rownames(seuratObj[[assay_used]])

markers_panel_filtered <- lapply(markers_panel, function(g) intersect(g, all_features))
markers_panel_filtered <- markers_panel_filtered[sapply(markers_panel_filtered, length) > 0]

if (length(markers_panel_filtered) == 0) {
  stop("Aucun gène du panel n'a été trouvé dans l'assay actif ('", assay_used, "').")
}

cat("Marqueurs conservés par panel:\n")
print(markers_panel_filtered)

# Pour fallback DotPlot Seurat (format named list accepté)
markers_vec <- unique(unlist(markers_panel_filtered))

# 4) Fonction fallback si dotPlotCustom n'existe pas
make_marker_dotplot <- function(obj, markersList, groupBy, plotTitle = "Markers average expression") {
  
  if (exists("dotPlotCustom", mode = "function")) {
    # utilise ta fonction custom si elle existe
    p <- dotPlotCustom(
      obj,
      markersList = markersList,
      groupBy = groupBy,
      rotateX = FALSE,
      plotTitle = plotTitle
    )
    return(p)
  }
  
  # Fallback Seurat::DotPlot
  p <- DotPlot(
    object = obj,
    features = markersList,   # liste nommée => facettes par groupe de gènes
    group.by = groupBy,
    assay = DefaultAssay(obj)
  ) +
    RotatedAxis() +
    labs(
      title = plotTitle,
      x = "",
      y = groupBy
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      plot.title = element_text(face = "bold"),
      legend.position = "right"
    )
  
  return(p)
}

# RotatedAxis fallback si absent
if (!exists("RotatedAxis")) {
  RotatedAxis <- function() {
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }
}
# 5) PDF de sortie
outfile <- file.path(figuresDir, paste0(lib, "_Markers_v1.pdf"))
if (file.exists(outfile)) {
  try(file.remove(outfile), silent = TRUE)
}

pdf(
  file = outfile,
  width = 9,
  height = 10,
  title = "Markers-v1 - no.filt"
)
# 6) Boucle sur les groupings (VDJ_like_class2 + résolutions)
for (i in markersBy) {
  
  cat("Plotting markers for grouping:", i, "\n")
  
  # gérer NA pour éviter soucis de plotting
  seuratObj@meta.data[[i]] <- as.character(seuratObj@meta.data[[i]])
  seuratObj@meta.data[[i]][is.na(seuratObj@meta.data[[i]])] <- "missing"
  
  ncolors <- length(unique(seuratObj@meta.data[[i]]))
  
  # UMAP du grouping
  p1 <- DimPlot(
    seuratObj,
    reduction = "umap",
    group.by = i,
    cols = colPaletteGenerator(ncolors),
    pt.size = 0.7,
    alpha = 0.8
  ) +
    labs(title = i, colour = "") +
    My_umap_theme() +
    theme(
      plot.margin = margin(
        l = 0.9,
        r = ifelse(i == "VDJ_like_class2", 0.3, 1),
        unit = "in"
      ),
      plot.title = element_text(size = 14, face = "bold")
    ) +
    guides(
      colour = guide_legend(
        override.aes = list(size = 3),
        ncol = ifelse(ncolors < 11, 1, 2)
      )
    )
  
  # Dotplot des marqueurs
  p2 <- make_marker_dotplot(
    seuratObj,
    markersList = markers_panel_filtered,
    groupBy = i,
    plotTitle = paste("Markers average expression -", i)
  )
  
  # Combiner
  print(
    cowplot::plot_grid(p1, p2, nrow = 2, rel_heights = c(1, 2))
  )
}

dev.off()


#plot markers UMAP bien lisible --------------

outDir <- file.path(figuresDir, "Markers_pretty")
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)

# 1) Palette + thèmes (legend à droite)
colPaletteGenerator <- function(n) {
  grDevices::colorRampPalette(c(
    "#1b9e77", "#d95f02", "#7570b3", "#e7298a",
    "#66a61e", "#e6ab02", "#a6761d", "#666666"
  ))(n)
}

My_umap_theme_right <- function() {
  theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.text  = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(size = 14, face = "bold"),
      legend.position = "right",
      legend.title = element_text(size = 11, face = "bold"),
      legend.text  = element_text(size = 10)
    )
}

My_dotplot_theme_right <- function() {
  theme_bw() +
    theme(
      plot.title = element_text(size = 13, face = "bold"),
      axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 10),
      legend.position = "right",
      legend.title = element_text(size = 10, face = "bold"),
      legend.text  = element_text(size = 9)
    )
}
# 2) Panel de marqueurs
markers_panel <- list(
  "Mature T"         = c("TRAC","CD3D","CD3E","IL7R","CCR7"),
  "T-ALL/immature"   = c("CD7","CD99","HES1","LMO2","MKI67","TOP2A"),
  "Myeloid"          = c("LYZ","S100A8","FCGR3A"),
  "B cells"          = c("MS4A1","CD79A","CD74"),
  "NK/cytotoxic"     = c("NKG7","GNLY","PRF1")
)
# 3) Fonction : créer le plot (UMAP dessus + DotPlot dessous)
plot_umap_dot <- function(obj, groupBy, markers_panel) {
  if (!groupBy %in% colnames(obj@meta.data)) {
    stop("groupBy absent dans meta.data : ", groupBy)
  }
  if (!"umap" %in% names(obj@reductions)) {
    stop("UMAP absent : lance RunUMAP() avant.")
  }
  
  assay_used <- DefaultAssay(obj)
  genes_present <- rownames(obj[[assay_used]])
  
  # Filtrer les gènes absents
  markers_panel2 <- lapply(markers_panel, function(g) intersect(g, genes_present))
  markers_panel2 <- markers_panel2[sapply(markers_panel2, length) > 0]
  if (length(markers_panel2) == 0) stop("Aucun gène du panel trouvé dans l'assay actif.")
  
  ncolors <- length(unique(obj@meta.data[[groupBy]]))
  
  p_umap <- DimPlot(
    obj,
    reduction = "umap",
    group.by = groupBy,
    cols = colPaletteGenerator(ncolors),
    pt.size = 0.6
  ) +
    ggtitle(paste0("UMAP - ", groupBy)) +
    My_umap_theme_right() +
    guides(colour = guide_legend(override.aes = list(size = 3)))
  
  p_dot <- DotPlot(
    obj,
    features = markers_panel2,   # liste => panels
    group.by = groupBy,
    assay = assay_used
  ) +
    ggtitle(paste0("Markers average expression - ", groupBy)) +
    My_dotplot_theme_right()
  
  # UMAP au-dessus, DotPlot dessous
  p_umap / p_dot + plot_layout(heights = c(1, 1.8))
}

# 4) Les groupBy à tracer
markersBy <- c("VDJ_like_class2", "SCT_snn_res.0.1", "SCT_snn_res.0.3", "SCT_snn_res.0.5", "SCT_snn_res.1.1")
markersBy <- markersBy[markersBy %in% colnames(seuratObj@meta.data)]
if (length(markersBy) == 0) stop("Aucun groupBy trouvé dans meta.data.")

# 5) A) Enregistrer 1 PDF multi-pages
pdf_all <- file.path(outDir, paste0(lib3, "_Markers_UMAP_DotPlot_pretty_ALL.pdf"))
pdf(pdf_all, width = 11, height = 13)
library(patchwork)
for (gb in markersBy) {
  p <- plot_umap_dot(seuratObj, gb, markers_panel)
  print(p)
}

dev.off()



# Clonotype analysis en tenant moins compte du cycle cellulaire
# régression S.Score / G2M.Score)----------------------------------

# 0) Copier l'objet pour ne pas écraser l'original
seuratObj_ccreg <- seuratObj

# 1) Vérifier / calculer percent.mt si absent
if (!"percent.mt" %in% colnames(seuratObj_ccreg@meta.data)) {
  seuratObj_ccreg <- PercentageFeatureSet(
    seuratObj_ccreg,
    pattern = "^MT-",
    col.name = "percent.mt",
    assay = "RNA"
  )
}

# 2) Calculer le cycle cellulaire (si Phase absent)
#    -> sur RNA (bonne pratique)
if (!all(c("S.Score", "G2M.Score", "Phase") %in% colnames(seuratObj_ccreg@meta.data))) {
  
  DefaultAssay(seuratObj_ccreg) <- "RNA"
  
  seuratObj_ccreg <- CellCycleScoring(
    seuratObj_ccreg,
    s.features = Seurat::cc.genes.updated.2019$s.genes,
    g2m.features = Seurat::cc.genes.updated.2019$g2m.genes,
    set.ident = FALSE
  )
}

# 3) SCTransform en régressant mt + cycle cellulaire (S.Score et G2M.Score)
# Si un assay SCT existe déjà, on peut le remplacer proprement
if ("SCT" %in% Assays(seuratObj_ccreg)) {
  # Optionnel: supprimer l'ancien SCT si tu veux repartir proprement
  # seuratObj_ccreg[["SCT"]] <- NULL
}

seuratObj_ccreg <- SCTransform(
  seuratObj_ccreg,
  assay = "RNA",
  new.assay.name = "SCT",
  vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"),
  verbose = TRUE,
  ncells = min(5000, ncol(seuratObj_ccreg)),
  variable.features.n = 3000
)

DefaultAssay(seuratObj_ccreg) <- "SCT"

# 4) PCA / Neighbors / Clustering / UMAP
seuratObj_ccreg <- RunPCA(seuratObj_ccreg, npcs = 50, verbose = FALSE)

ndims <- 30

seuratObj_ccreg <- FindNeighbors(
  seuratObj_ccreg,
  dims = 1:ndims,
  verbose = FALSE
)

# Résolution principale = 0.5
seuratObj_ccreg <- FindClusters(
  seuratObj_ccreg,
  resolution = 0.5,
  verbose = FALSE
)

# Multi-résolutions 
seuratObj_ccreg <- FindClusters(
  seuratObj_ccreg,
  resolution = c(seq(0.1, 0.6, by = 0.1), seq(0.7, 1.5, by = 0.2)),
  verbose = FALSE
)

seuratObj_ccreg <- RunUMAP(
  seuratObj_ccreg,
  dims = 1:ndims,
  verbose = FALSE
)

# Choisir la résolution principale = 0.5
Idents(seuratObj_ccreg) <- "SCT_snn_res.0.5"
seuratObj_ccreg$seurat_clusters <- seuratObj_ccreg$SCT_snn_res.0.5

# 5) Ajouter coordonnées UMAP à meta.data
umap_df <- as.data.frame(Embeddings(seuratObj_ccreg, "umap"))
colnames(umap_df) <- c("umap_1", "umap_2")

seuratObj_ccreg@meta.data <- seuratObj_ccreg@meta.data %>%
  tibble::rownames_to_column("cell_barcode_tmp") %>%
  left_join(
    umap_df %>% tibble::rownames_to_column("cell_barcode_tmp"),
    by = "cell_barcode_tmp"
  ) %>%
  tibble::column_to_rownames("cell_barcode_tmp")


# 6) Vérifs
cat("Assays:", paste(Assays(seuratObj_ccreg), collapse = ", "), "\n")
cat("Colonnes cycle cellulaire présentes ? ",
    all(c("S.Score","G2M.Score","Phase") %in% colnames(seuratObj_ccreg@meta.data)),
    "\n")
cat("Colonnes clustering:\n")
print(colnames(seuratObj_ccreg@meta.data)[grep("SCT_snn_res|seurat_clusters", colnames(seuratObj_ccreg@meta.data))])

# 7) Plots clonotype (sur UMAP corrigé du cycle cellulaire)

# Petite palette si tu l'as déjà sinon on la crée
if (!exists("colPaletteGenerator")) {
  colPaletteGenerator <- function(n) {
    grDevices::colorRampPalette(c(
      "#1b9e77", "#d95f02", "#7570b3", "#e7298a",
      "#66a61e", "#e6ab02", "#a6761d", "#666666"
    ))(n)
  }
}

# Thème UMAP si absent
if (!exists("My_umap_theme")) {
  My_umap_theme <- function() {
    theme_bw() +
      theme(
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "bottom",
        strip.background = element_rect(fill = "grey95"),
        strip.text = element_text(face = "bold")
      )
  }
}

# Dossier de sortie
qcDir <- file.path(figuresDir, "QC")
dir.create(qcDir, recursive = TRUE, showWarnings = FALSE)

outfile <- file.path(qcDir, paste0(lib3, "_Clonotype_noCellCycle_UMAP.pdf"))
pdf(outfile, width = 10, height = 8)

# Plot 1: UMAP clonotype 
if ("Clonotype" %in% colnames(seuratObj_ccreg@meta.data)) {
  ncolors_clono <- length(unique(seuratObj_ccreg@meta.data$Clonotype))
  print(
    DimPlot(
      seuratObj_ccreg,
      reduction = "umap",
      group.by = "Clonotype",
      cols = colPaletteGenerator(ncolors_clono),
      pt.size = 0.6,
      alpha = 0.8
    ) +
      ggtitle("Clonotype (UMAP avec régression du cycle cellulaire)") +
      My_umap_theme()
  )
}

# Plot 2: Clonotype facetté par condition
if (all(c("Clonotype","Condition","umap_1","umap_2") %in% colnames(seuratObj_ccreg@meta.data))) {
  ncolors_clono <- length(unique(seuratObj_ccreg@meta.data$Clonotype))
  print(
    seuratObj_ccreg@meta.data %>%
      ggplot(aes(x = umap_1, y = umap_2)) +
      geom_point(aes(colour = Clonotype), size = 0.6, alpha = 0.8) +
      facet_grid(~ Condition) +
      scale_colour_manual(values = colPaletteGenerator(ncolors_clono)) +
      ggtitle("Clonotype par condition (UMAP corrigé cycle)") +
      My_umap_theme() +
      guides(colour = guide_legend(
        title = "",
        nrow = ifelse(ncolors_clono <= 16, 1, 2),
        override.aes = list(size = 3)
      ))
  )
}

# Plot 3: Vérification visuelle du cycle après régression 
if ("Phase" %in% colnames(seuratObj_ccreg@meta.data)) {
  n_phase <- length(unique(seuratObj_ccreg@meta.data$Phase))
  print(
    DimPlot(
      seuratObj_ccreg,
      reduction = "umap",
      group.by = "Phase",
      cols = colPaletteGenerator(n_phase),
      pt.size = 0.6,
      alpha = 0.8
    ) +
      ggtitle("Phase (après régression du cycle)") +
      My_umap_theme()
  )
}

dev.off()

#DEG global ------------------------

DefaultAssay(seuratObj) <- "SCT"
seuratObj$Condition <- dplyr::recode(seuratObj$Condition,
                                     "M187r" = "Relapse",
                                     "M187"  = "Diagnosis")
Idents(seuratObj) <- "Condition"

deg_global <- FindMarkers(
  seuratObj,
  ident.1 = "Relapse",
  ident.2 = "Diagnosis",
  logfc.threshold = 0.25,
  min.pct = 0.1
)

head(deg_global)

deg_global$gene <- rownames(deg_global)
deg_global$neglog10p <- -log10(deg_global$p_val_adj + 1e-300)

ggplot(deg_global, aes(x=avg_log2FC, y=neglog10p)) +
  geom_point(size=0.7) +
  theme_bw() +
  labs(title="Volcano DEG Relapse vs Diagnosis", x="avg_log2FC", y="-log10(adj p)")
top_genes <- rownames(deg_global)[1:30]
DoHeatmap(seuratObj, features = top_genes, group.by="Condition")

degFigDir <- file.path(figuresDir, "DEG_global")
dir.create(degFigDir, recursive = TRUE, showWarnings = FALSE)
p_volcano <- ggplot(deg_global, aes(x=avg_log2FC, y=neglog10p)) +
  geom_point(size=0.7) +
  theme_bw() +
  labs(title="Volcano DEG Relapse vs Diagnosis", x="avg_log2FC", y="-log10(adj p)")

ggsave(
  filename = file.path(degFigDir, paste0(lib3, "_Volcano_DEG_Relapse_vs_Diagnosis.png")),
  plot = p_volcano,
  width = 7, height = 6, units = "in",
  dpi = 300
)
top_genes <- rownames(deg_global)[1:30]

p_heat <- DoHeatmap(seuratObj, features = top_genes, group.by="Condition") +
  ggtitle("Top 30 DEG - Heatmap")

ggsave(
  filename = file.path(degFigDir, paste0(lib3, "_Heatmap_Top30_DEG.png")),
  plot = p_heat,
  width = 10, height = 8, units = "in",
  dpi = 300
)



#DEG par cluster ----------------------------------------

degDir <- file.path(objectsDir, "DEG_by_cluster")
figDir <- file.path(figuresDir, "DEG_by_cluster")
dir.create(degDir, recursive = TRUE, showWarnings = FALSE)
dir.create(figDir, recursive = TRUE, showWarnings = FALSE)

# 1) Réglages
cluster_col <- "SCT_snn_res.0.5"
min_cells_per_condition <- 30
topN_heat <- 20         # heatmap: top 20 gènes/cluster
topN_feat <- 4          # featureplot: top 4 gènes/cluster

DefaultAssay(seuratObj) <- "SCT"
seuratObj$Condition <- factor(seuratObj$Condition, levels = c("Diagnosis", "Relapse"))

stopifnot(cluster_col %in% colnames(seuratObj@meta.data))
stopifnot("Condition" %in% colnames(seuratObj@meta.data))

# Clusters en identités
Idents(seuratObj) <- seuratObj[[cluster_col, drop = TRUE]]
clusters <- levels(Idents(seuratObj))
# 2) Filtre gènes "techniques" : RPL/RPS + MT-

is_tech_gene <- function(g) {
  grepl("^RPL", g) | grepl("^RPS", g) | grepl("^MT-", g)
}

# 3) Volcano helper

make_volcano <- function(deg, title, p_cut=0.05, fc_cut=0.25) {
  df <- deg
  df$gene <- rownames(df)
  df$neglog10p <- -log10(df$p_val_adj + 1e-300)
  df$sig <- (df$p_val_adj < p_cut) & (abs(df$avg_log2FC) >= fc_cut)
  
  ggplot(df, aes(x = avg_log2FC, y = neglog10p)) +
    geom_point(aes(alpha = sig), size = 0.8) +
    geom_vline(xintercept = c(-fc_cut, fc_cut), linetype = "dashed") +
    geom_hline(yintercept = -log10(p_cut), linetype = "dashed") +
    theme_bw() +
    theme(legend.position = "none") +
    labs(title = title,
         x = "avg_log2FC (Relapse vs Diagnosis)",
         y = "-log10(adj p)")
}

# 4) PDFs multi-pages (volcano + heatmaps + featureplots)
pdf_volcano <- file.path(figDir, paste0(lib3, "_DEG_byCluster_volcano.pdf"))
pdf_heat    <- file.path(figDir, paste0(lib3, "_DEG_byCluster_heatmaps.pdf"))
pdf_feat    <- file.path(figDir, paste0(lib3, "_DEG_byCluster_featureplots_top4.pdf"))

pdf(pdf_volcano, width = 7, height = 6)
pdf(pdf_heat,    width = 10, height = 7)
pdf(pdf_feat,    width = 11, height = 7)

all_deg <- list()

for (cl in clusters) {
  
  obj_cl <- subset(seuratObj, idents = cl)
  
  tab <- table(obj_cl$Condition)
  n_diag <- ifelse("Diagnosis" %in% names(tab), tab[["Diagnosis"]], 0)
  n_rel  <- ifelse("Relapse"  %in% names(tab), tab[["Relapse"]], 0)
  
  if (n_diag < min_cells_per_condition || n_rel < min_cells_per_condition) {
    message("Skip cluster ", cl, " (diag=", n_diag, ", relapse=", n_rel, ")")
    next
  }
  
  # DEG dans ce cluster
  Idents(obj_cl) <- "Condition"
  DefaultAssay(obj_cl) <- "SCT"
  
  deg <- FindMarkers(
    obj_cl,
    ident.1 = "Relapse",
    ident.2 = "Diagnosis",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    test.use = "wilcox"
  )
  
  if (nrow(deg) == 0) next
  
  deg$cluster <- cl
  deg$n_diag <- n_diag
  deg$n_relapse <- n_rel
  deg$gene <- rownames(deg)
  
  # Sauvegarde CSV cluster
  write.csv(deg,
            file.path(degDir, paste0("DEG_cluster_", cl, "_Relapse_vs_Diagnosis.csv")),
            row.names = FALSE)
  
  all_deg[[cl]] <- deg
  
  # (A) Volcano (PDF + PNG)
  pV <- make_volcano(deg, paste0("Cluster ", cl, " (diag=", n_diag, ", rel=", n_rel, ")"))
  print(pV)
  ggsave(file.path(figDir, paste0(lib3, "_cluster_", cl, "_volcano.png")),
         plot = pV, width = 7, height = 6, dpi = 300)
  

  # (B) Top genes "propres" (sans RPL/RPS/MT)
  #     -> pour heatmap + featureplot
  deg_clean <- deg %>%
    arrange(p_val_adj) %>%
    filter(!is_tech_gene(gene))
  
  # si trop strict et vide, fallback sur deg complet
  if (nrow(deg_clean) < 5) {
    deg_clean <- deg %>% arrange(p_val_adj)
  }
  
  top_genes_heat <- deg_clean %>%
    slice_head(n = min(topN_heat, nrow(deg_clean))) %>%
    pull(gene)
  
  top_genes_feat <- deg_clean %>%
    slice_head(n = min(topN_feat, nrow(deg_clean))) %>%
    pull(gene)
  
  # (C) Heatmap (PDF + PNG)

  Idents(obj_cl) <- "Condition"
  
  pH <- DoHeatmap(obj_cl, features = top_genes_heat, group.by = "Condition", size = 3) +
    ggtitle(paste0("Cluster ", cl, " - Top ", length(top_genes_heat), " DEG (no RPL/RPS/MT)"))
  print(pH)
  
  png(file.path(figDir, paste0(lib3, "_cluster_", cl, "_heatmap.png")),
      width = 2800, height = 1800, res = 300)
  print(pH)
  dev.off()
  
  # (D) FeaturePlot top4 (UMAP) (PDF + PNG)
  # FeaturePlot fonctionne bien si UMAP existe dans l'objet (hérité de seuratObj)
  pF <- FeaturePlot(obj_cl, features = top_genes_feat, reduction = "umap",
                    ncol = 2, order = TRUE) +
    ggtitle(paste0("Cluster ", cl, " - Top ", length(top_genes_feat), " DEG (FeaturePlot)"))
  print(pF)
  
  png(file.path(figDir, paste0(lib3, "_cluster_", cl, "_featureplot_top", length(top_genes_feat), ".png")),
      width = 2600, height = 1800, res = 300)
  print(pF)
  dev.off()
}

dev.off() # volcano pdf
dev.off() # heatmap pdf
dev.off() # featureplot pdf


# 5) CSV global + résumé
deg_all_df <- bind_rows(all_deg)
write.csv(deg_all_df,
          file.path(degDir, paste0(lib3, "_DEG_ALLclusters_Relapse_vs_Diagnosis.csv")),
          row.names = FALSE)

deg_summary <- deg_all_df %>%
  mutate(sig = (p_val_adj < 0.05) & (abs(avg_log2FC) >= 0.25)) %>%
  group_by(cluster) %>%
  summarise(
    n_genes = n(),
    n_sig = sum(sig, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n_sig))

write.csv(deg_summary,
          file.path(degDir, paste0(lib3, "_DEG_summary_by_cluster.csv")),
          row.names = FALSE)

# 150 UP + 150 DOWN pour CLUE.io (Relapse vs Diagnosis) -------------

# 0) Sécurité : il faut deg_global déjà calculé
stopifnot(exists("deg_global"))
stopifnot(all(c("avg_log2FC","p_val_adj") %in% colnames(deg_global)))

deg_global$gene <- rownames(deg_global)

# 1) Paramètres (tu peux ajuster)
N <- 150
padj_cut <- 0.05
min_fc <- 0.25   # optionnel: cohérent avec ton FindMarkers

# 2) Filtrer significatifs (et éviter les artefacts techniques)
is_tech_gene <- function(g) grepl("^RPL", g) | grepl("^RPS", g) | grepl("^MT-", g)

deg_filt <- deg_global |>
  dplyr::filter(!is.na(p_val_adj), !is.na(avg_log2FC)) |>
  dplyr::filter(p_val_adj < padj_cut, abs(avg_log2FC) >= min_fc) |>
  dplyr::filter(!is_tech_gene(gene))

# 3) Up = logFC positif (Relapse > Diagnosis), Down = logFC négatif
up_tbl <- deg_filt |>
  dplyr::filter(avg_log2FC > 0) |>
  dplyr::arrange(p_val_adj, dplyr::desc(avg_log2FC))

down_tbl <- deg_filt |>
  dplyr::filter(avg_log2FC < 0) |>
  dplyr::arrange(p_val_adj, avg_log2FC)

up_150 <- head(up_tbl$gene, N)
down_150 <- head(down_tbl$gene, N)

cat("UP genes:", length(up_150), "\n")
cat("DOWN genes:", length(down_150), "\n")

# 4) Si tu n’as pas assez de gènes après filtres, fallback moins strict
if (length(up_150) < N || length(down_150) < N) {
  message("Pas assez de gènes avec ces filtres -> fallback (p_val_adj seulement, sans min_fc)")
  
  deg_fallback <- deg_global |>
    dplyr::filter(!is.na(p_val_adj), !is.na(avg_log2FC)) |>
    dplyr::filter(p_val_adj < padj_cut) |>
    dplyr::filter(!is_tech_gene(gene))
  
  up_150 <- head(
    deg_fallback |>
      dplyr::filter(avg_log2FC > 0) |>
      dplyr::arrange(p_val_adj, dplyr::desc(avg_log2FC)) |>
      dplyr::pull(gene),
    N
  )
  down_150 <- head(
    deg_fallback |>
      dplyr::filter(avg_log2FC < 0) |>
      dplyr::arrange(p_val_adj, avg_log2FC) |>
      dplyr::pull(gene),
    N
  )
}

# 5) Sauvegarder en fichiers texte (format le plus simple pour CLUE)
clueDir <- file.path(figuresDir, "CLUE_lists")
dir.create(clueDir, recursive = TRUE, showWarnings = FALSE)

writeLines(up_150,   con = file.path(clueDir, paste0(lib3, "_CLUE_UP_150.txt")))
writeLines(down_150, con = file.path(clueDir, paste0(lib3, "_CLUE_DOWN_150.txt")))

# 6) Optionnel : un CSV récapitulatif
clue_df <- data.frame(
  gene = c(up_150, down_150),
  direction = c(rep("UP", length(up_150)), rep("DOWN", length(down_150)))
)
write.csv(clue_df, file.path(clueDir, paste0(lib3, "_CLUE_UP_DOWN_150.csv")), row.names = FALSE)

cat("Fichiers générés dans :", clueDir, "\n")

  

    # choix du cluster automatique pour clue ------------
  cluster_col <- "SCT_snn_res.0.5"
  
  tab <- table(seuratObj$Condition, seuratObj[[cluster_col, drop=TRUE]])
  prop <- prop.table(tab, margin = 1)
  
  # delta de proportion (Relapse - Diagnosis)
  delta <- prop["Relapse",] - prop["Diagnosis",]
  
  df <- data.frame(
    cluster = colnames(tab),
    n_Diagnosis = as.numeric(tab["Diagnosis",]),
    n_Relapse   = as.numeric(tab["Relapse",]),
    prop_Diagnosis = as.numeric(prop["Diagnosis",]),
    prop_Relapse   = as.numeric(prop["Relapse",]),
    delta = as.numeric(delta)
  )
  
  # garde ceux où on peut faire DEG par cluster (assez de cellules dans les 2)
  df_ok <- df |> dplyr::filter(n_Diagnosis >= 30, n_Relapse >= 30) |> dplyr::arrange(dplyr::desc(abs(delta)))
  
  df_ok
  head(df_ok, 5)   # top 5 clusters à tester dans CLUE
  
  # Générer UP150 + DOWN150 (CLUE) pour plusieurs clusters (0/2/4) ------------

  library(dplyr)
  
  # Clusters à traiter (selon ton df_ok)
  clusters_to_run <- c("0", "2", "4")   # tu peux ajouter "5" etc.
  
  # Paramètres
  N <- 150
  padj_cut <- 0.05
  min_fc <- 0.25
  
  # Filtre gènes techniques
  is_tech_gene <- function(g) grepl("^RPL", g) | grepl("^RPS", g) | grepl("^MT-", g)
  
  # Dossier sortie
  clueDir <- file.path(figuresDir, "CLUE_lists_byCluster")
  dir.create(clueDir, recursive = TRUE, showWarnings = FALSE)
  
  # Résumé
  summary_list <- list()
  
  for (cl in clusters_to_run) {
    
    deg_path <- file.path(objectsDir, "DEG_by_cluster",
                          paste0("DEG_cluster_", cl, "_Relapse_vs_Diagnosis.csv"))
    
    if (!file.exists(deg_path)) {
      warning("Fichier introuvable pour cluster ", cl, " : ", deg_path)
      next
    }
    
    deg_cl <- read.csv(deg_path)
    
    # Vérifs colonnes
    need_cols <- c("gene","avg_log2FC","p_val_adj")
    if (!all(need_cols %in% colnames(deg_cl))) {
      warning("Colonnes manquantes dans ", deg_path, " -> trouvé: ",
              paste(colnames(deg_cl), collapse=", "))
      next
    }
    
    # Filtrage strict
    deg_filt <- deg_cl |>
      filter(!is.na(p_val_adj), !is.na(avg_log2FC)) |>
      filter(p_val_adj < padj_cut, abs(avg_log2FC) >= min_fc) |>
      filter(!is_tech_gene(gene))
    
    # UP/DOWN
    up_tbl <- deg_filt |>
      filter(avg_log2FC > 0) |>
      arrange(p_val_adj, desc(avg_log2FC))
    
    down_tbl <- deg_filt |>
      filter(avg_log2FC < 0) |>
      arrange(p_val_adj, avg_log2FC)
    
    up_150 <- head(up_tbl$gene, N)
    down_150 <- head(down_tbl$gene, N)
    
    # Fallback si pas assez
    if (length(up_150) < N || length(down_150) < N) {
      message("Cluster ", cl, ": pas assez avec filtres stricts -> fallback (p_val_adj seulement)")
      
      deg_fb <- deg_cl |>
        filter(!is.na(p_val_adj), !is.na(avg_log2FC)) |>
        filter(p_val_adj < padj_cut) |>
        filter(!is_tech_gene(gene))
      
      up_150 <- head(
        deg_fb |>
          filter(avg_log2FC > 0) |>
          arrange(p_val_adj, desc(avg_log2FC)) |>
          pull(gene),
        N
      )
      
      down_150 <- head(
        deg_fb |>
          filter(avg_log2FC < 0) |>
          arrange(p_val_adj, avg_log2FC) |>
          pull(gene),
        N
      )
    }
    
    # Sauvegarde txt pour CLUE
    up_file <- file.path(clueDir, paste0(lib3, "_cluster", cl, "_CLUE_UP_", N, ".txt"))
    dn_file <- file.path(clueDir, paste0(lib3, "_cluster", cl, "_CLUE_DOWN_", N, ".txt"))
    
    writeLines(up_150, con = up_file)
    writeLines(down_150, con = dn_file)
    
    # CSV récap cluster
    clue_df <- data.frame(
      gene = c(up_150, down_150),
      direction = c(rep("UP", length(up_150)), rep("DOWN", length(down_150))),
      cluster = cl
    )
    write.csv(
      clue_df,
      file.path(clueDir, paste0(lib3, "_cluster", cl, "_CLUE_UP_DOWN_", N, ".csv")),
      row.names = FALSE
    )
    
    # Résumé
    summary_list[[cl]] <- data.frame(
      cluster = cl,
      up_n = length(up_150),
      down_n = length(down_150),
      up_file = up_file,
      down_file = dn_file
    )
    
    cat("Cluster ", cl, " OK | UP=", length(up_150), " DOWN=", length(down_150), "\n", sep="")
  }
  
  summary_df <- bind_rows(summary_list)
  print(summary_df)
  
  write.csv(
    summary_df,
    file.path(clueDir, paste0(lib3, "_CLUE_lists_summary.csv")),
    row.names = FALSE
  )
  
  cat("\nTous les fichiers CLUE sont dans :", clueDir, "\n")
  
  # generer les genes pour cluster 2 et 4 --------------
  library(dplyr)
    # Helper: filtre gènes techniques
  is_tech_gene <- function(g) grepl("^RPL", g) | grepl("^RPS", g) | grepl("^MT-", g)
  
  # Fonction: génère UP/DOWN pour CLUE depuis un CSV DEG cluster

  make_clue_lists_from_deg <- function(deg_df, N = 150,
                                       try_steps = list(
                                         list(padj = 0.05, min_fc = 0.25),
                                         list(padj = 0.10, min_fc = 0.15),
                                         list(padj = 0.20, min_fc = 0.10),
                                         list(padj = 0.50, min_fc = 0.00)  # ultra large
                                       )) {
    stopifnot(all(c("gene","avg_log2FC","p_val_adj") %in% colnames(deg_df)))
    
    # nettoyer
    deg_df <- deg_df %>%
      filter(!is.na(avg_log2FC), !is.na(p_val_adj)) %>%
      filter(!is_tech_gene(gene))
    
    # base de ranking (utile pour compléter si pas assez)
    up_rank_all <- deg_df %>% filter(avg_log2FC > 0) %>% arrange(p_val_adj, desc(avg_log2FC))
    dn_rank_all <- deg_df %>% filter(avg_log2FC < 0) %>% arrange(p_val_adj, avg_log2FC)
    
    best_up <- character(0)
    best_dn <- character(0)
    used_step <- NA_character_
    
    # essayer plusieurs niveaux de filtres
    for (st in try_steps) {
      padj_cut <- st$padj
      min_fc   <- st$min_fc
      
      df_f <- deg_df %>%
        filter(p_val_adj < padj_cut, abs(avg_log2FC) >= min_fc)
      
      up <- df_f %>% filter(avg_log2FC > 0) %>% arrange(p_val_adj, desc(avg_log2FC)) %>% pull(gene)
      dn <- df_f %>% filter(avg_log2FC < 0) %>% arrange(p_val_adj, avg_log2FC) %>% pull(gene)
      
      if (length(up) >= N && length(dn) >= N) {
        best_up <- head(up, N)
        best_dn <- head(dn, N)
        used_step <- paste0("padj<", padj_cut, " & |logFC|>=", min_fc)
        break
      }
      
      # garde le meilleur trouvé jusqu’ici
      if (length(up) > length(best_up)) best_up <- up
      if (length(dn) > length(best_dn)) best_dn <- dn
      used_step <- paste0("padj<", padj_cut, " & |logFC|>=", min_fc, " (insufficient, will fill)")
    }
    
    # compléter si pas assez (avec ranking global)
    if (length(best_up) < N) {
      fill <- setdiff(up_rank_all$gene, best_up)
      best_up <- c(best_up, head(fill, N - length(best_up)))
    }
    if (length(best_dn) < N) {
      fill <- setdiff(dn_rank_all$gene, best_dn)
      best_dn <- c(best_dn, head(fill, N - length(best_dn)))
    }
    
    # sécurité finale
    best_up <- head(unique(best_up), N)
    best_dn <- head(unique(best_dn), N)
    
    list(up = best_up, down = best_dn, rule = used_step,
         up_available = nrow(up_rank_all), down_available = nrow(dn_rank_all))
  }
  
  # Génération pour cluster 2 et 4
  clusters <- c("2","4")
  N <- 150
  
  clueDir <- file.path(figuresDir, "CLUE_lists_byCluster_moreGenes")
  dir.create(clueDir, recursive = TRUE, showWarnings = FALSE)
  
  for (cl in clusters) {
    deg_path <- file.path(objectsDir, "DEG_by_cluster",
                          paste0("DEG_cluster_", cl, "_Relapse_vs_Diagnosis.csv"))
    
    if (!file.exists(deg_path)) {
      warning("Fichier DEG introuvable : ", deg_path)
      next
    }
    
    deg_cl <- read.csv(deg_path)
    
    # Générer listes CLUE (avec relaxation + complétion)
    res <- make_clue_lists_from_deg(deg_cl, N = N)
    
    # Sauver
    up_file <- file.path(clueDir, paste0(lib3, "_cluster", cl, "_CLUE_UP_", N, ".txt"))
    dn_file <- file.path(clueDir, paste0(lib3, "_cluster", cl, "_CLUE_DOWN_", N, ".txt"))
    
    writeLines(res$up, con = up_file)
    writeLines(res$down, con = dn_file)
    
    # petit résumé
    cat("\nCluster ", cl, ":\n", sep="")
    cat("  Rule used: ", res$rule, "\n", sep="")
    cat("  UP genes written: ", length(res$up), " (available UP in DEG: ", res$up_available, ")\n", sep="")
    cat("  DOWN genes written: ", length(res$down), " (available DOWN in DEG: ", res$down_available, ")\n", sep="")
    cat("  Files:\n   - ", up_file, "\n   - ", dn_file, "\n", sep="")
  }
  
  cat("\nTous les fichiers sont dans : ", clueDir, "\n", sep="")
  
  
# clue pour clonotypes 1  et 2 -------------
  # ── Créer colonne combinée ──────────────────────────────────────
  seuratObj$Condition_Clonotype <- paste0(seuratObj$Condition, "_", seuratObj$Clonotype)
  table(seuratObj$Condition_Clonotype)
  
  # Setter les identités
  Idents(seuratObj) <- "Condition_Clonotype"
  DefaultAssay(seuratObj) <- "SCT"
  
  # ── DEG Clonotype1 : Rechute vs Diagnostic ─────────────────────
  deg_clono1 <- FindMarkers(
    seuratObj,
    ident.1 = "Relapse_clonotype1",
    ident.2 = "Diagnosis_clonotype1",
    logfc.threshold = 0.25,
    min.pct = 0.1
  )
  cat("DEG clonotype1 :", nrow(deg_clono1), "gènes\n")
  
  # ── DEG Clonotype2 : Rechute vs Diagnostic ─────────────────────
  deg_clono2 <- FindMarkers(
    seuratObj,
    ident.1 = "Relapse_clonotype2",
    ident.2 = "Diagnosis_clonotype2",
    logfc.threshold = 0.25,
    min.pct = 0.1
  )
  cat("DEG clonotype2 :", nrow(deg_clono2), "gènes\n")
  
  # ── Extraire gènes UP/DOWN pour CLUE ───────────────────────────
  extract_clue_genes <- function(deg, name) {
    
    genes_up <- deg %>%
      filter(p_val_adj < 0.05 & avg_log2FC > 0.25) %>%
      arrange(desc(avg_log2FC)) %>%
      head(150) %>%
      rownames()
    
    genes_down <- deg %>%
      filter(p_val_adj < 0.05 & avg_log2FC < -0.25) %>%
      arrange(avg_log2FC) %>%
      head(150) %>%
      rownames()
    
    cat("\n========", name, "========\n")
    cat("UP   :", length(genes_up), "gènes\n")
    cat("DOWN :", length(genes_down), "gènes\n")
    
    cat("\n--- GENES UP ---\n")
    cat(genes_up, sep = "\n")
    
    cat("\n--- GENES DOWN ---\n")
    cat(genes_down, sep = "\n")
    
    # Retourner une liste pour usage ultérieur
    list(up = genes_up, down = genes_down)
  }
  
  clue_clono1 <- extract_clue_genes(deg_clono1, "Clonotype1 Relapse vs Diagnosis")
  clue_clono2 <- extract_clue_genes(deg_clono2, "Clonotype2 Relapse vs Diagnosis")

  
  # ── Dossier de sortie ──────────────────────────────────────────
  clueDir <- file.path(figuresDir, "CLUE")
  dir.create(clueDir, recursive = TRUE, showWarnings = FALSE)
  
  # ── Fonction extraction + enregistrement ──────────────────────
  extract_and_save_clue <- function(deg, name, outDir) {
    
    genes_up <- deg %>%
      filter(p_val_adj < 0.05 & avg_log2FC > 0.25) %>%
      arrange(desc(avg_log2FC)) %>%
      head(150) %>%
      rownames()
    
    genes_down <- deg %>%
      filter(p_val_adj < 0.05 & avg_log2FC < -0.25) %>%
      arrange(avg_log2FC) %>%
      head(150) %>%
      rownames()
    
    cat("\n========", name, "========\n")
    cat("UP   :", length(genes_up), "gènes\n")
    cat("DOWN :", length(genes_down), "gènes\n")
    
    # Enregistrer UP
    writeLines(
      genes_up,
      file.path(outDir, paste0(name, "_UP.txt"))
    )
    
    # Enregistrer DOWN
    writeLines(
      genes_down,
      file.path(outDir, paste0(name, "_DOWN.txt"))
    )
    
    # Enregistrer tableau DEG complet
    write.csv(
      deg %>% tibble::rownames_to_column("gene"),
      file.path(outDir, paste0(name, "_DEG_full.csv")),
      row.names = FALSE
    )
    
    cat("Fichiers sauvegardés dans :", outDir, "\n")
    
    list(up = genes_up, down = genes_down)
  }
  
  # ── Lancer pour chaque clonotype ──────────────────────────────
  clue_clono1 <- extract_and_save_clue(
    deg_clono1,
    name   = "Clonotype1_Relapse_vs_Diagnosis",
    outDir = clueDir
  )
  
  clue_clono2 <- extract_and_save_clue(
    deg_clono2,
    name   = "Clonotype2_Relapse_vs_Diagnosis",
    outDir = clueDir
  )
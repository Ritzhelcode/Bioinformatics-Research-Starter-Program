if (!require("BiocManager", quietly = TRUE)) 
  {install.packages("BiocManager")}
BiocManager::install(c("GEOquery", "limma"),ask=FALSE, update=FALSE)
#Install annotation package sesuai platform
#GPL96 = Affymetrix Human Genome U133A
BiocManager::install("hgu133a.db", ask = FALSE, update = FALSE)
install.packages(c("pheatmap", "ggplot2", "dplyr"))
#umap: grafik (plot UMAP)
if (!requireNamespace("umap",quietly=TRUE)){install.packages("umap")}
library(GEOquery)
library(limma)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(hgu133a.db)
library(AnnotationDbi)
library(umap)
gset <- getGEO("GSE10072",GSEMatrix=TRUE, AnnotGPL=TRUE)[[1]]
ex <-exprs(gset)
qx <-as.numeric(quantile(ex, c(0, 0.25, 0.5, 0.75, 0.99, 1), na.rm =TRUE))
LogTransform <-(qx[5] > 100)||(qx[6]-qx[1]>50&&qx[2]>0)
#IF statement:
if (LogTransform) {
  # Nilai <= 0 tidak boleh di-log, maka diubah menjadi NA
  ex[ex <= 0] <- NA
  ex <- log2(ex)}
#pData(): metadata sampel
#source_name_ch1 berisi informasi kondisi biologis sampel
group_info <- pData(gset)[["source_name_ch1"]]  
#make.names(): mengubah teks menjadi format valid untuk R
groups <- make.names(group_info)
gset$group <- factor(groups)
nama_grup <- levels(gset$group)
print(nama_grup)

#model.matrix():
#Membuat matriks desain untuk model linear
#~0 berarti TANPA intercept (best practice limma)
design <- model.matrix(~0 + gset$group)
#colnames(): memberi nama kolom agar mudah dibaca
colnames(design) <- levels(gset$group)
#Menentukan perbandingan biologis
grup_kanker <- nama_grup[1]
grup_normal <- nama_grup[2]
contrast_formula <- paste(grup_kanker, "-", grup_normal)
print(paste("Kontras yang dianalisis:", contrast_formula))

#lmFit():
#Membangun model linear untuk setiap gen
fit <- lmFit(ex, design)
#makeContrasts(): mendefinisikan perbandingan antar grup
contrast_matrix <- makeContrasts(contrasts = contrast_formula, levels
                                 = design)
#contrasts.fit(): menerapkan kontras ke model
fit2 <- contrasts.fit(fit, contrast_matrix)
#eBayes():
#Empirical Bayes untuk menstabilkan estimasi varians
fit2 <- eBayes(fit2)
#topTable():
#Mengambil hasil akhir DEG
#adjust = "fdr" -> koreksi multiple testing
#p.value = 0.01 -> gen sangat signifikan
topTableResults <- topTable(
  fit2,
  adjust = "fdr",
  sort.by = "B",
  number = Inf,
  p.value = 0.01)
head(topTableResults)

probe_ids <- rownames(topTableResults)
#Mapping probe -> gene symbol & gene name
gene_annotation <- AnnotationDbi::select(
  hgu133a.db,
  keys = probe_ids,
  columns = c("SYMBOL", "GENENAME"),
  keytype = "PROBEID"
)
topTableResults$PROBEID <- rownames(topTableResults)
topTableResults <- merge(
  topTableResults,
  gene_annotation,
  by = "PROBEID",
  all.x = TRUE
)
#Hasil Anotasi
head(topTableResults[, c("PROBEID", "SYMBOL", "GENENAME")])

#Boxplot
#Set warna berdasarkan grup
group_colors <- as.numeric(gset$group)
boxplot(
  ex,
  col = group_colors,
  las = 2,
  outline = FALSE,
  main = "Boxplot Distribusi Nilai Ekspresi per Sampel",
  ylab = "Expression Value (log2)"
)
legend(
  "topright",
  legend = levels(gset$group),
  fill = unique(group_colors),
  cex = 0.8
)

#Densityplot
#Gabungkan ekspresi & grup ke data frame
expr_long <- data.frame(
  Expression = as.vector(ex),
  Group = rep(gset$group, each = nrow(ex))
)
if (LogTransform) {
  # Nilai <= 0 tidak boleh di-log, maka diubah menjadi NA
  ex[ex <= 0] <- NA
  ex <- log2(ex)
}
library(ggplot2)
ggplot(expr_long, aes(x = Expression, color = Group)) +
  geom_density(linewidth = 1) +
  theme_minimal() +
  labs(
    title = "Distribusi Nilai Ekspresi Gen",
    x = "Expression Value (log2)",
    y = "Density"
  )
#UMAP
#Transpose matriks ekspresi:
#UMAP bekerja pada OBSERVATION = sampel
umap_input <- t(ex)

#Jalankan UMAP
library(umap)
umap_result<-umap(umap_input)

#Simpan hasil ke data frame
umap_df <- data.frame(
  UMAP1 = umap_result$layout[, 1],
  UMAP2 = umap_result$layout[, 2],
  Group = gset$group
)
library(ggplot2)
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Group)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal() +
  labs(
    title = "UMAP Plot Sampel Berdasarkan Ekspresi Gen",
    x = "UMAP 1",
    y = "UMAP 2"
  )
#VolcanoPlot
#Volcano plot menggabungkan:
#- Log fold change (efek biologis)
#- Signifikansi statistik

volcano_data <- data.frame(
  logFC = topTableResults$logFC,
  adj.P.Val = topTableResults$adj.P.Val,
  Gene = topTableResults$SYMBOL
)
#Klasifikasi status gen
volcano_data$status <- "NO"
volcano_data$status[volcano_data$logFC > 1 & volcano_data$adj.P.Val <
                      0.01] <- "UP"
volcano_data$status[volcano_data$logFC < -1 & volcano_data$adj.P.Val <
                      0.01] <- "DOWN"
#Visualisasi
library(ggplot2)
ggplot(volcano_data, aes(x = logFC, y = -log10(adj.P.Val), color =
                           status)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("DOWN" = "blue", "NO" = "grey", "UP" =
                                  "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed") +
  theme_minimal() +
  ggtitle("Volcano Plot DEG Kanker Paru")

#Heatmap
#Heatmap digunakan untuk melihat pola ekspresi gen
#antar sampel berdasarkan gen-gen paling signifikan

#Pilih 50 gen paling signifikan berdasarkan adj.P.Val
topTableResults <- topTableResults[
  order(topTableResults$adj.P.Val),
]
top50 <- head(topTableResults, 50)
#Ambil matriks ekspresi untuk gen terpilih
mat_heatmap <- ex[top50$PROBEID, ]

#Gunakan Gene Symbol (fallback ke Probe ID)
gene_label <- ifelse(
  is.na(top50$SYMBOL) | top50$SYMBOL == "",
  top50$PROBEID, # jika SYMBOL kosong → probe ID
  top50$SYMBOL # jika ada → gene symbol
)
rownames(mat_heatmap) <- gene_label
#Pembersihan data (WAJIB agar tidak error hclust)
#Hapus baris dengan NA
mat_heatmap <- mat_heatmap[
  rowSums(is.na(mat_heatmap)) == 0,
]
#Hapus gen dengan varians nol
gene_variance <- apply(mat_heatmap, 1, var)
mat_heatmap <- mat_heatmap[gene_variance > 0, ]
#Anotasi kolom (kelompok sampel)
annotation_col <- data.frame(
  Group = gset$group
)
rownames(annotation_col) <- colnames(mat_heatmap)

#Visualisasi heatmap
library(pheatmap)
pheatmap(
  mat_heatmap,
  scale = "row", # Z-score per gen
  annotation_col = annotation_col,
  show_colnames = FALSE, # nama sampel dimatikan
  show_rownames = TRUE,
  fontsize_row = 7,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  main = "Top 50 Differentially Expressed Genes"
)

# write.csv(): menyimpan hasil analisis ke file CSV
write.csv(topTableResults, "Hasil_GSE10072_DEG.csv")
message("Analisis selesai. File hasil telah disimpan.")



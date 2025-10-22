# =============================================================================
# Project 2: TCGA-COAD miRNA — Tumor Primary (TP) vs Normal (NT)
# Pipeline: Query → Prepare → Normalize (TMM) → edgeR DE → Volcano → Exports
# =============================================================================

# ------------------------------- 0) Libraries --------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(SummarizedExperiment)
  library(TCGAbiolinks)
  library(edgeR)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap) # optional heatmap later
  # library(httpgd) # optional browser plotting
})

# ------------------------------- 1) Settings ---------------------------------
project_id   <- "TCGA-COAD"
workflow     <- "BCGSC miRNA Profiling"
fdr_cutoff   <- 0.05          # significance threshold after multiple testing
lfc_cutoff   <- 1.0           # |log2FC| >= 1 (≈2x change)
label_top_n  <- 15            # number of miRNAs to label on the plot
set.seed(42)

# Helper: NULL-coalescing
`%||%` <- function(a, b) if (!is.null(a)) a else b

# ------------------------ 2) Query/Download/Prepare --------------------------
message(">> Querying GDC for TCGA-COAD miRNA (TP & NT) ...")
query <- GDCquery(
  project       = project_id,
  data.category = "Transcriptome Profiling",
  data.type     = "miRNA Expression Quantification",
  workflow.type = workflow,
  sample.type   = c("Primary Tumor","Solid Tissue Normal") # ensure both groups
)
GDCdownload(query)
data <- GDCprepare(query)

# ---------------------- 3) Build Counts Matrix (Robust) ----------------------
make_counts_matrix <- function(obj) {
  if (inherits(obj, "SummarizedExperiment")) {
    a_names <- names(assays(obj))
    counts  <- if ("read_count" %in% a_names) assay(obj, "read_count") else assays(obj)[[1]]
    rownames(counts) <- make.unique(rownames(counts))
    return(counts)
  }
  if (is.data.frame(obj) || is.matrix(obj)) {
    df <- as.data.frame(obj)
    bc_cols <- grep("^TCGA-", colnames(df), value = TRUE)
    if (!length(bc_cols)) bc_cols <- colnames(df)[sapply(df, is.numeric)]
    counts <- as.matrix(df[, bc_cols, drop = FALSE])
    if ("miRNA_ID" %in% colnames(df)) {
      rownames(counts) <- make.unique(df$miRNA_ID)
    } else if (is.null(rownames(counts))) {
      id_col <- setdiff(colnames(df), bc_cols)[1]
      rownames(counts) <- make.unique(df[[id_col]] %||% paste0("miR_", seq_len(nrow(counts))))
    }
    return(counts)
  }
  stop("Unsupported object type from GDCprepare().")
}

counts <- make_counts_matrix(data)

# ---------------------- 4) Define Groups & Sanity Checks ---------------------
all_barcodes    <- colnames(counts)
tumor_barcodes  <- TCGAquery_SampleTypes(barcode = all_barcodes, typesample = "TP")
normal_barcodes <- TCGAquery_SampleTypes(barcode = all_barcodes, typesample = "NT")

grp <- ifelse(all_barcodes %in% normal_barcodes, "NT",
              ifelse(all_barcodes %in% tumor_barcodes,  "TP", NA))
keep_cols <- !is.na(grp)

counts <- counts[, keep_cols, drop = FALSE]
group  <- factor(grp[keep_cols], levels = c("NT", "TP"))
barcodes <- colnames(counts)

message(">> Samples per group:")
print(table(group))
stopifnot(all(table(group) > 0)) # both NT & TP must exist

# --------------------------- 5) edgeR DE Analysis ----------------------------
message(">> edgeR: filtering, normalization, GLM (TP vs NT) ...")
dge <- DGEList(counts = counts, group = group)

# Filter lowly-expressed miRNAs
keep_genes <- filterByExpr(dge, group = group)
dge <- dge[keep_genes, , keep.lib.sizes = FALSE]

# TMM normalization
dge <- calcNormFactors(dge)

# Design: no intercept, contrast TP - NT
design   <- model.matrix(~ 0 + group) # columns "NT","TP"
colnames(design) <- levels(group)
contrast <- makeContrasts(TPvsNT = TP - NT, levels = design)

# Fit GLM with quasi-likelihood
dge  <- estimateDisp(dge, design)
fit  <- glmQLFit(dge, design)
qlf  <- glmQLFTest(fit, contrast = contrast)

# Differential table
de <- topTags(qlf, n = Inf)$table %>%
  as.data.frame()
de$FDR <- p.adjust(de$PValue, "BH")

# ------------------ 6) Volcano Table + Direction Labels ---------------------
message(">> Building volcano data & calling Up/Down ...")
volc <- de %>%
  mutate(
    neglog10FDR = -log10(pmax(FDR, .Machine$double.xmin)),
    status = case_when(
      FDR < fdr_cutoff & logFC >=  lfc_cutoff  ~ "Up (TP)",    # higher in tumor
      FDR < fdr_cutoff & logFC <= -lfc_cutoff  ~ "Down (TP)",  # lower in tumor
      TRUE                                     ~ "NS"          # not significant
    )
  )

up_tp    <- volc %>% filter(status == "Up (TP)")   %>% arrange(FDR, desc(abs(logFC)))
down_tp  <- volc %>% filter(status == "Down (TP)") %>% arrange(FDR, desc(abs(logFC)))
sig_all  <- volc %>% filter(status != "NS")        %>% arrange(FDR)

message(">> Significant counts:")
print(with(volc, table(status)))

# ----------------------------- 7) Volcano Plot -------------------------------
message(">> Plotting volcano ...")
to_label <- volc %>% filter(status != "NS") %>% arrange(FDR) %>% head(label_top_n)

p_volcano <- ggplot(volc, aes(x = logFC, y = neglog10FDR)) +
  geom_point(aes(color = status), alpha = 0.75, size = 1.8) +
  geom_hline(yintercept = -log10(fdr_cutoff), linetype = 2) +
  geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = 2) +
  geom_text_repel(
    data = to_label,
    aes(label = rownames(to_label)),
    size = 3,
    max.overlaps = 30,
    box.padding = 0.3
  ) +
  scale_color_manual(values = c("Up (TP)"="red", "Down (TP)"="green", "NS"="grey")) +
  labs(
    title = "Volcano: Tumor Primary (TP) vs Normal (NT) — TCGA-COAD miRNA (edgeR)",
    x = "log2 Fold Change (TP - NT)",
    y = "-log10(FDR)",
    color = "Direction"
  ) +
  theme_minimal(base_size = 12)

print(p_volcano)
ggsave("Volcano_TP_vs_NT_COAD.png", p_volcano, width = 8, height = 6, dpi = 300)

# ----------------------- 8) Optional: Fold-Change Summary --------------------
# log2CPM group means (helpful for intuition/report tables)
logCPM <- cpm(dge, log = TRUE, prior.count = 1)
mean_tp <- rowMeans(logCPM[, group == "TP", drop = FALSE])
mean_nt <- rowMeans(logCPM[, group == "NT", drop = FALSE])
fc_tbl <- tibble(
  miRNA = rownames(logCPM),
  log2FC_meanCPM = mean_tp - mean_nt
) %>% arrange(desc(abs(log2FC_meanCPM)))
write.csv(fc_tbl, "COAD_miRNA_log2FC_from_logCPM_TP_minus_NT.csv", row.names = FALSE)

# ------------------------------ 9) Save Outputs ------------------------------
write.csv(up_tp,   "COAD_miRNA_upregulated_TP_vs_NT.csv",   row.names = TRUE)
write.csv(down_tp, "COAD_miRNA_downregulated_TP_vs_NT.csv", row.names = TRUE)
write.csv(sig_all, "COAD_miRNA_significant_TP_vs_NT.csv",   row.names = TRUE)

# -------------------- 10) (Optional) Heatmap of Top DE -----------------------
# topN <- 40
# top_features <- rownames(sig_all)[seq_len(min(topN, nrow(sig_all)))]
# Z <- scale(t(logCPM[top_features, ]))  # samples x features
# ann <- data.frame(Group = group)
# rownames(ann) <- colnames(logCPM)
# pheatmap(t(Z), annotation_col = ann, show_rownames = FALSE,
#          main = sprintf("Top %d DE miRNAs — TP vs NT", length(top_features)))

# =============================================================================
# End of pipeline
# =============================================================================

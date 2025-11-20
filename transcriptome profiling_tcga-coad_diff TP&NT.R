# =============================================================================
# Project 2: Transcriptome Profiling TCGA-COAD miRNA — Tumor Primary (TP) vs Normal (NT)
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
})

# ------------------------------- 1) Settings ---------------------------------
fdr_cutoff   <- 0.05          # significance threshold after multiple testing
lfc_cutoff   <- 1.0           # sets the log₂ fold change threshold # |log2FC| >= 1 (≈2x change)
label_top_n  <- 20            # number of miRNAs to label on the plot
set.seed(42)                  # fixes randomization for reproducibility

# Helper: NULL-coalescing
`%||%` <- function(a, b) if (!is.null(a)) a else b # prevent NULL from breaking code

# ------------------------ 2) Query/Download/Prepare --------------------------
message(">> Querying GDC for TCGA-COAD miRNA (TP & NT) ...")
query <- GDCquery(
  project       = "TCGA-COAD",
  data.category = "Transcriptome Profiling",
  data.type     = "miRNA Expression Quantification",
  workflow.type = "BCGSC miRNA Profiling",
  sample.type   = c("Primary Tumor","Solid Tissue Normal") # ensure both groups
)
GDCdownload(query)
data <- GDCprepare(query)

# ---------------------- 3) Build Counts Matrix  ---------------
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

counts  <- counts[, keep_cols, drop = FALSE]
group   <- factor(grp[keep_cols], levels = c("NT", "TP")) # NT = reference
barcodes <- colnames(counts)

message(">> Samples per group:")
print(table(group))
stopifnot(all(table(group) > 0)) # both NT & TP must exist

# --------------------------- 5) edgeR DE Analysis (Simplified) ---------------
message(">> edgeR: filtering, normalization, GLM (TP vs NT) ...")

# Wrap counts into a DGEList object
dge <- DGEList(counts = counts, group = group)

# Filter lowly-expressed miRNAs
keep_genes <- filterByExpr(dge, group = group)
dge <- dge[keep_genes, , keep.lib.sizes = FALSE]

# TMM normalization
dge <- calcNormFactors(dge)

# Design matrix with intercept (NT as baseline, groupTP = TP - NT effect)
design <- model.matrix(~ group)  # columns: "(Intercept)", "groupTP"

# Fit GLM with quasi-likelihood
dge  <- estimateDisp(dge, design)
fit  <- glmQLFit(dge, design)
qlf  <- glmQLFTest(fit, coef = 2)   # coef=2 corresponds to "groupTP"

# Differential expression table
de <- topTags(qlf, n = Inf)$table %>%
  as.data.frame()
de$FDR <- p.adjust(de$PValue, "BH")

# ------------------ 6) Volcano Data + Status Labels --------------------------
message(">> Building volcano data ...")

volc <- de %>%
  mutate(
    neglog10FDR = -log10(pmax(FDR, .Machine$double.xmin)),
    status = case_when(
      FDR < fdr_cutoff & logFC >=  lfc_cutoff  ~ "Up (TP)",    # higher in tumor
      FDR < fdr_cutoff & logFC <= -lfc_cutoff  ~ "Down (NT)",  # lower in tumor
      TRUE                                     ~ "NS"          # not significant
    )
  )

# ----------------------------- 7) Volcano Plot -------------------------------
message(">> Plotting volcano ...")

y_cut <- -log10(fdr_cutoff)

p_volcano <- ggplot(volc, aes(x = logFC, y = neglog10FDR)) +
  # (bottom strip: FDR non-significant area)
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = y_cut,
           fill = "#F0F0F0", alpha = 1) +
  # (middle vertical strip: small |logFC| region)
  annotate("rect", xmin = -lfc_cutoff, xmax =  lfc_cutoff, ymin = 0, ymax = Inf,
           fill = "#F5F5F5", alpha = 1) +
  # points
  geom_point(aes(color = status), size = 1.8, alpha = 0.8) +
  # threshold lines
  geom_hline(yintercept = y_cut, linetype = "dashed", linewidth = 0.4) +
  geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff),
             linetype = "dashed", color = "grey60", linewidth = 0.4) +
  # colors & labels
  scale_color_manual(values = c("Up (TP)"="#FF3030",
                                "Down (NT)"="#54FF9F",
                                "NS"="#A2B5CD")) +
  labs(
    title = "Tumor Primary (TP) vs Normal (NT) — TCGA-COAD miRNA (edgeR)",
    x = "log2 Fold Change (TP - NT)",
    y = "-log10(FDR)",
    color = "Direction Result"
  ) +
  coord_cartesian(
    xlim = c(-10, 10),
    ylim = c(0, max(volc$neglog10FDR, na.rm = TRUE) * 1.05)
  ) +
  theme_light(base_size = 13) +
  theme(
    plot.title  = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right",
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    axis.line    = element_line(linewidth = 0.5),
    axis.ticks   = element_line(linewidth = 0.4)
  )

print(p_volcano)
# ggsave("Volcano_TP_vs_NT_COAD.png", p_volcano, width = 7, height = 6, dpi = 300)

# =============================================================================
# End of pipeline (minimal version for final volcano plot)
# =============================================================================

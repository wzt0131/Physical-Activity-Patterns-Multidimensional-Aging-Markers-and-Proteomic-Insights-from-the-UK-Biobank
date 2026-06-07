# Pathway enrichment

library(dplyr)
library(impute)
library(broom)
results_list <- list()
success_count <- 0
warning_count <- 0
error_count <- 0
target_variable <- "PhenoAgeAccel" 
protein_columns <- colnames(sports_phenoagehigh_pro)[35:2954]  
covariates <- c("smoke","alcohol" ," Sex" ,"Townsend.deprivation.index.at.recruitment "," Body.mass.index..BMI....Instance.0" ," high_bloodpressure" ," ethnic" ,"education","employment") 

for(i in seq_along(protein_columns)) {
  protein <- protein_columns[i]
  cat("Processing protein", i, "/", length(protein_columns), ":", protein, "\n")
  tryCatch({
    if(all(is.na(sports_phenoagehigh_pro[[protein]])) || 
       sd(sports_phenoagehigh_pro[[protein]], na.rm = TRUE) == 0) {
      stop("Protein data is constant or all NA")
    }
    formula_str <- paste0(target_variable, " ~ ", protein, " + ", 
                          paste(covariates, collapse = " + "))
    model <- lm(as.formula(formula_str), 
                data = sports_phenoagehigh_pro)
    model_summary <- tidy(model, conf.int = TRUE) %>% 
      filter(term == protein) %>% 
      mutate(protein = protein,
             formula = formula_str,
             converged = TRUE)  
    results_list[[protein]] <- model_summary
    success_count <- success_count + 1
    
  }, warning = function(w) {
    cat("Warning:", w$message, "\n")
    warning_result <- data.frame(
      protein = protein,
      term = protein,
      estimate = NA,
      std.error = NA,
      statistic = NA,
      p.value = NA,
      conf.low = NA,
      conf.high = NA,
      formula = paste0(target_variable, " ~ ", protein, " + ", 
                       paste(covariates, collapse = " + ")),
      converged = FALSE,
      warning_message = w$message,
      stringsAsFactors = FALSE
    )
    results_list[[protein]] <- warning_result
    warning_count <- warning_count + 1
  }, error = function(e) {
    cat("Error:", e$message, "\n")
    error_result <- data.frame(
      protein = protein,
      term = protein,
      estimate = NA,
      std.error = NA,
      statistic = NA,
      p.value = NA,
      conf.low = NA,
      conf.high = NA,
      formula = paste0(target_variable, " ~ ", protein, " + ", 
                       paste(covariates, collapse = " + ")),
      converged = FALSE,
      error_message = e$message,
      stringsAsFactors = FALSE
    )
    results_list[[protein]] <- error_result
    error_count <- error_count + 1
  })
}
cat("Analysis completed! Success:", success_count, "Warnings:", warning_count, "Errors:", error_count, "\n")
final_results <- bind_rows(results_list)
final_results <- final_results %>% 
  mutate(p.adj = p.adjust(p.value, method = "fdr"))
significant_proteins <- final_results %>% 
  filter(p.adj < 0.05) %>% 
  arrange(p.adj)
head(significant_proteins, 10)
cat("Number of significant proteins after FDR correction:", nrow(significant_proteins), "\n")

library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyverse)

final_results <- read.csv("protein_aging_association_results.csv", stringsAsFactors = FALSE)

volcano_data <- final_results %>%
  filter(!is.na(estimate) & !is.na(p.adj)) %>%
  mutate(
    log_pvalue = -log10(p.adj),  # Using adjusted p-value
    log_odds = estimate,
    significance = case_when(
      p.adj < 0.05 & estimate > 0 ~ "Up",
      p.adj < 0.05 & estimate < 0 ~ "Down",
      TRUE ~ "Not Significant"
    )
  )

sig_proteins <- final_results %>% filter(p.adj < 0.05)
cat("Significantly associated protein names:\n")
print(sig_proteins$term)

top_proteins <- volcano_data %>%
  arrange(p.value) %>%
  head(20)

p <- ggplot(volcano_data, aes(x = log_odds, y = log_pvalue, color = significance)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(
    values = c(
      "Up" = "red",
      "Down" = "blue", 
      "Not Significant" = "grey"
    ),
    name = "Significance"
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "darkgrey") +
  geom_hline(yintercept = -log10(0.05/nrow(volcano_data)), linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "darkgrey") +
  geom_label_repel(
    data = top_proteins,
    aes(label = term),
    size = 2,
    box.padding = 0.3,
    point.padding = 0.1,
    max.overlaps = 50,
    show.legend = FALSE
  )+
  labs(
    title = "Protein Association with PhenoAgeAccel",
    x = "Estimate (Effect Size)",
    y = "-log10(P-value)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title = element_text(size = 12),
    legend.position = "right"
  )
print(p)

library(BiocManager)
required_packages <- c("clusterProfiler", "org.Hs.eg.db", "enrichplot", "ggplot2", "dplyr", "stringr")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg)
  }
  library(pkg, character.only = TRUE)
}

if (!dir.exists("enrichment_results")) {
  dir.create("enrichment_results")
}
if (!dir.exists("enrichment_plots")) {
  dir.create("enrichment_plots")
}

cat("Data frame column names:", colnames(significant_proteins), "\n")
cat("Number of significant proteins:", nrow(significant_proteins), "\n")

protein_list <- significant_proteins$term
protein_list <- trimws(protein_list)
protein_list <- protein_list[protein_list != "" & protein_list != "NA" & !is.na(protein_list)]
cat("Number of valid proteins:", length(protein_list), "\n")

if (length(protein_list) < 5) {
  stop("Insufficient number of valid proteins (n =", length(protein_list), "), cannot perform enrichment analysis")
}

install.packages("shadowtext", type = "binary")
BiocManager::install("clusterProfiler", force = TRUE)
library(clusterProfiler)

if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  BiocManager::install("org.Hs.eg.db")
}
library(org.Hs.eg.db)

gene_df <- tryCatch({
  bitr(protein_list, fromType = "SYMBOL", 
       toType = "ENTREZID", 
       OrgDb = org.Hs.eg.db)
}, error = function(e) {
  cat("Gene ID conversion error:", e$message, "\n")
  return(data.frame())
})

if (nrow(gene_df) == 0) {
  stop("Unable to convert gene symbols, please check protein name format")
}
cat("Number of successfully converted genes:", nrow(gene_df), "\n")

# GO
go_bp <- tryCatch({
  enrichGO(gene = gene_df$ENTREZID,
           OrgDb = org.Hs.eg.db,
           keyType = "ENTREZID",
           ont = "BP",
           pAdjustMethod = "BH",
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.2,
           readable = TRUE)
}, error = function(e) {
  cat("GO enrichment analysis error:", e$message, "\n")
  return(NULL)
})

if (!is.null(go_bp) && nrow(go_bp) > 0) {
  write.csv(go_bp, file = "enrichment_results/GO_BP_enrichment.csv", row.names = FALSE)
  cat("GO results saved, number of significant terms:", nrow(go_bp), "\n")
  tryCatch({
    p1 <- barplot(go_bp, showCategory = 15) + 
      ggtitle("GO Biological Process Enrichment")
    ggsave("enrichment_plots/GO_BP_barplot.png", plot = p1, width = 12, height = 8)
    p2 <- dotplot(go_bp, showCategory = 15) + 
      ggtitle("GO Biological Process Enrichment")
    ggsave("enrichment_plots/GO_BP_dotplot.png", plot = p2, width = 12, height = 8)
    p3 <- cnetplot(go_bp, showCategory = 10)
    ggsave("enrichment_plots/GO_BP_network.png", plot = p3, width = 12, height = 8)
    
  }, error = function(e) {
    cat("Error generating GO plots:", e$message, "\n")
  })
} else {
  cat("No significant GO terms found\n")
}

# KEGG
kegg <- tryCatch({
  enrichKEGG(gene = gene_df$ENTREZID,
             organism = "hsa",
             pAdjustMethod = "BH",
             pvalueCutoff = 0.05,
             qvalueCutoff = 0.2)
}, error = function(e) {
  cat("KEGG enrichment analysis error:", e$message, "\n")
  return(NULL)
})

if (!is.null(kegg) && nrow(kegg) > 0) {
  kegg <- setReadable(kegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  write.csv(kegg, file = "enrichment_results/KEGG_enrichment.csv", row.names = FALSE)
  cat("KEGG results saved, number of significant pathways:", nrow(kegg), "\n")
  tryCatch({
    p1 <- barplot(kegg, showCategory = 15) + 
      ggtitle("KEGG Pathway Enrichment")
    ggsave("enrichment_plots/KEGG_barplot.png", plot = p1, width = 12, height = 8)
    p2 <- dotplot(kegg, showCategory = 15) + 
      ggtitle("KEGG Pathway Enrichment")
    ggsave("enrichment_plots/KEGG_dotplot.png", plot = p2, width = 12, height = 8)
    
  }, error = function(e) {
    cat("Error generating KEGG plots:", e$message, "\n")
  })
} else {
  cat("No significant KEGG pathways found\n")
}

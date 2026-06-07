# Druggability genome tier enrichment

library(tidyverse)
library(readxl)
library(clusterProfiler)
library(enrichplot)
library(openxlsx)
library(ggplot2)

setwd("D:/OneDrive/Desktop/drug_target")

mr_data <- read.csv("mr_results_summary.csv", stringsAsFactors = FALSE, 
                    fileEncoding = "UTF-8")

druggable_data <- read_excel("NIHMS80906-supplement-Table_S1.xlsx")

head(mr_data)
head(druggable_data)

mr_data <- mr_data %>%
  mutate(gene_symbol = str_extract(exposure, "^[^_]+"))

table(mr_data$gene_symbol)

background_genes <- unique(mr_data$gene_symbol)

significant_genes <- mr_data %>%
  filter(n_significant > 0) %>%
  pull(gene_symbol) %>%
  unique()

cat("Background gene count:", length(background_genes), "\n")
cat("Significant gene count:", length(significant_genes), "\n")

term2gene <- druggable_data %>%
  select(term = druggability_tier, gene = hgnc_names) %>%
  na.omit()

table(term2gene$term)

enrich_result <- enricher(
  gene = significant_genes,
  universe = background_genes,
  TERM2GENE = term2gene,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  qvalueCutoff = 1,
  minGSSize = 1,
  maxGSSize = 1000
)

enrich_df <- as.data.frame(enrich_result)
print(enrich_df)

if(nrow(enrich_df) > 0) {
  p1 <- barplot(enrich_result, showCategory = 15, 
                font.size = 10, 
                title = "Druggability Genome Tier Enrichment Analysis - Barplot")
  print(p1)
}

if(nrow(enrich_df) > 0) {
  p2 <- dotplot(enrich_result, showCategory = 15, 
                font.size = 10,
                title = "Druggability Genome Tier Enrichment Analysis - Dotplot")
  print(p2)
}

if(nrow(enrich_df) > 0) {
  enrich_score_plot <- enrich_df %>%
    mutate(Enrichment_Score = -log10(pvalue)) %>%
    arrange(desc(Enrichment_Score)) %>%
    head(15) %>%
    ggplot(aes(x = reorder(Description, Enrichment_Score), y = Enrichment_Score)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    coord_flip() +
    labs(x = "Druggability Tier", y = "-log10(p-value)", 
         title = "Druggability Tier Enrichment Score") +
    theme_minimal()
  print(enrich_score_plot)
}

write.csv(enrich_df, "druggability_tier_enrichment_results.csv", 
          row.names = FALSE, fileEncoding = "UTF-8")

tier_stats <- term2gene %>%
  group_by(term) %>%
  summarise(
    total_genes_in_tier = n(),
    significant_genes_in_tier = sum(gene %in% significant_genes),
    background_genes_in_tier = sum(gene %in% background_genes)
  ) %>%
  mutate(
    enrichment_ratio = (significant_genes_in_tier / length(significant_genes)) / 
      (background_genes_in_tier / length(background_genes))
  )

print(tier_stats)

write.csv(tier_stats, "druggability_tier_detailed_stats.csv", 
          row.names = FALSE, fileEncoding = "UTF-8")

gene_tier_mapping <- significant_genes %>%
  as.data.frame() %>%
  setNames("gene_symbol") %>%
  left_join(druggable_data %>% select(hgnc_names, druggability_tier), 
            by = c("gene_symbol" = "hgnc_names")) %>%
  left_join(mr_data %>% select(gene_symbol, n_significant, signif_rate, mean_or), 
            by = "gene_symbol")

write.csv(gene_tier_mapping, "significant_genes_tier_mapping.csv", 
          row.names = FALSE, fileEncoding = "UTF-8")

if(nrow(enrich_df) > 0) {
  volcano_plot <- enrich_df %>%
    mutate(log10pval = -log10(pvalue),
           log2OR = log2(as.numeric(gsub("/.*", "", GeneRatio)) / 
                           as.numeric(gsub(".*/", "", GeneRatio)) /
                           (as.numeric(gsub("/.*", "", BgRatio)) / 
                              as.numeric(gsub(".*/", "", BgRatio))))) %>%
    ggplot(aes(x = log2OR, y = log10pval, size = Count, color = Description)) +
    geom_point(alpha = 0.7) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
    labs(x = "log2(Odds Ratio)", y = "-log10(p-value)",
         title = "Druggability Tier Enrichment Analysis - Volcano Plot") +
    theme_minimal() +
    theme(legend.position = "bottom")
  print(volcano_plot)
}

cat("Analysis completed! Results saved to current working directory.\n")
cat("Main output files:\n")
cat("1. druggability_tier_enrichment_results.csv - Main enrichment results\n")
cat("2. druggability_tier_detailed_stats.csv - Detailed statistics\n")
cat("3. significant_genes_tier_mapping.csv - Gene-Tier mapping\n")

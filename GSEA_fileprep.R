setwd("E:/KY-INBRE PostDoc/Yvonne/Matt/RNAseq_Run1/GSEA/DESEQ2")

#list of folders for the comparisons
folder_list <- c("CS2_8WK__VS__CS1_8WK", "CS2_17WK__VS__CS1_17WK", "CS2_28WK__VS__CS1_28WK", "CS3_8WK__VS__CS1_8WK", "CS3_17WK__VS__CS1_17WK", "CS3_28WK__VS__CS1_28WK", "CS4_17WK__VS__CS1_17WK", "CS4_28WK__VS__CS1_28WK", "CS7_8WK__VS__CS1_8WK", "CS7_17WK__VS__CS1_17WK", "CS7_28WK__VS__CS1_28WK",
                 "CS12_8WK__VS__CS1_8WK", "CS12_17WK__VS__CS1_17WK", "CS12_28WK__VS__CS1_28WK", "CS14_8WK__VS__CS13_8WK", "CS14_17WK__VS__CS13_17WK", "CS14_28WK__VS__CS13_28WK", "CS13_8WK__VS__CS1_8WK", "CS13_17WK__VS__CS1_17WK", "CS13_28WK__VS__CS1_28WK", "CS14_8WK__VS__CS2_8WK", 
                 "CS14_17WK__VS__CS2_17WK", "CS14_28WK__VS__CS2_17WK")

#Divide the DESEQ2 results in to upregulated and downregulated gene lists
for (folder in folder_list){
  readfile <- paste0("./", folder, "/", folder, "_DESeq2_gene_exp.diff")
  df <- read.table(readfile, header = TRUE, sep = "\t")
  
  positive_logfc <- df[df$log2FoldChange > 0, ]  
  negative_logfc <- df[df$log2FoldChange < 0, ] 
  
  write.table(positive_logfc, paste0("./", folder, "/", folder, "_positive_logfc_UP.txt"), row.names = FALSE, quote = FALSE, sep = "\t")
  write.table(negative_logfc, paste0("./", folder, "/", folder, "_negative_logfc_DN.txt"), row.names = FALSE, quote = FALSE, sep = "\t")
 
}

#get the rank for UP
for (folder in folder_list){
  readfile <- paste0("./", folder, "/", folder, "_positive_logfc_UP.txt")
  df <- read.table(readfile, header = TRUE, sep = "\t")
  
  df$Rank <- df$log2FoldChange  * -log10(df$pvalue)
  df_sorted <- df[order(df$Rank, decreasing = TRUE), ]
  df_selected <- df_sorted[, c("ENSEMBL", "Rank")]
  df_selected <- df_selected[!is.na(df_selected$Rank), ]
  
  write.table(df_selected, paste0("./", folder, "/", folder, "_positive_logfc_UP_Ranked.rnk"), row.names = FALSE, quote = FALSE, sep = "\t")
}

#get the rank for DN
for (folder in folder_list){
  readfile <- paste0("./", folder, "/", folder, "_negative_logfc_DN.txt")
  df <- read.table(readfile, header = TRUE, sep = "\t")
  
  df$Rank <- df$log2FoldChange  * -log10(df$pvalue)
  #df$Rank <- format(df$Rank, scientific = FALSE, digits = 6)
  df_sorted <- df[order(df$Rank, decreasing = TRUE), ]
  df_selected <- df_sorted[, c("ENSEMBL", "Rank")]
  df_selected <- df_selected[!is.na(df_selected$Rank), ]
  
  write.table(df_selected, paste0("./", folder, "/", folder, "_negative_logfc_DN_Ranked.rnk"), row.names = FALSE, quote = FALSE, sep = "\t")
}

#convert the ensembl ids into gene symbol
biomart <- read.table("E:/KY-INBRE PostDoc/Yvonne/Matt/RNAseq_Run1/GSEA/hg38.p12_ENSEMBL93_GeneToSymbolAndDescription.txt", header = TRUE, sep = "\t", fill = TRUE, quote = "")


for (folder in folder_list){
  readfile <- paste0("./", folder, "/", folder, "_positive_logfc_UP_Ranked.rnk")
  df <- read.table(readfile, header = TRUE, sep = "\t")
  
  df_merged <- merge(df, biomart,  by.x = "ENSEMBL", by.y = "Gene.stable.ID", all.x = TRUE)
  df_selected <- df_merged[, c("HGNC.symbol", "Rank")]
  df_selected <- df_selected[!is.na(df_selected$Rank), ]
  
  write.table(df_selected, paste0("./", folder, "/", folder, "_positive_logfc_UP_Ranked.rnk"), row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
}

for (folder in folder_list){
  readfile <- paste0("./", folder, "/", folder, "_negative_logfc_DN_Ranked.rnk")
  df <- read.table(readfile, header = TRUE, sep = "\t")
  
  df_merged <- merge(df, biomart, by.x = "ENSEMBL", by.y = "Gene.stable.ID", all.x = TRUE)
  df_selected <- df_merged[, c("HGNC.symbol", "Rank")]
  df_selected <- df_selected[!is.na(df_selected$Rank), ]
  
  write.table(df_selected, paste0("./", folder, "/", folder, "_negative_logfc_DN_Ranked.rnk"), row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
}


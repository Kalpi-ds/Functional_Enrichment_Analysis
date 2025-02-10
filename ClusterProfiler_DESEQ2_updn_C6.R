library(GO.db)
library(graph)
library(clusterProfiler)
library(enrichplot)
library(AnnotationDbi)
library(GSEABase)
library(dplyr)
library("org.Hs.eg.db")
library(ggplot2)
library(msigdbr)


setwd("E:/KY-INBRE PostDoc/Yvonne/Matt/ClusterProfiler_1updn_C6")




compareList <- c("CS2_8WK__VS__CS1_8WK", "CS2_17WK__VS__CS1_17WK", "CS2_28WK__VS__CS1_28WK", "CS3_8WK__VS__CS1_8WK", "CS3_17WK__VS__CS1_17WK", "CS3_28WK__VS__CS1_28WK", "CS4_17WK__VS__CS1_17WK", "CS4_28WK__VS__CS1_28WK", "CS7_8WK__VS__CS1_8WK", "CS7_17WK__VS__CS1_17WK", "CS7_28WK__VS__CS1_28WK", 
 "CS12_8WK__VS__CS1_8WK", "CS12_17WK__VS__CS1_17WK", "CS12_28WK__VS__CS1_28WK", "CS14_8WK__VS__CS13_8WK", "CS14_17WK__VS__CS13_17WK", "CS14_28WK__VS__CS13_28WK", "CS13_8WK__VS__CS1_8WK", "CS13_17WK__VS__CS1_17WK", "CS13_28WK__VS__CS1_28WK", "CS14_8WK__VS__CS2_8WK", "CS14_17WK__VS__CS2_17WK",
 "CS14_28WK__VS__CS2_17WK")

#compareList <- c("CS2_17WK__VS__CS1_17WK", "CS2_28WK__VS__CS1_28WK", "CS3_17WK__VS__CS1_17WK", "CS3_28WK__VS__CS1_28WK", "CS4_17WK__VS__CS1_17WK", "CS4_28WK__VS__CS1_28WK", "CS12_17WK__VS__CS1_17WK",
#                 "CS12_28WK__VS__CS1_28WK", "CS13_17WK__VS__CS1_17WK", "CS14_28WK__VS__CS2_17WK", "CS14_28WK__VS__CS13_28WK")
######################################
## READ IN DIFFERENTIALLY EXPRESSED ##
## ENTREZ IDENTIFIERS AND CREATE    ##
## GENE LISTS FOR CATEGORICAL       ##
## ENRICHMENT                       ##
######################################
numCompare <- 23
for(i in 1:numCompare) {
  currInFN <- paste("./", compareList[[i]], "_DEG_UP_P",
                    0.05, "_Q0.05", "_FC0", "_FPKM1", "_FPKMAVG1",
                    "_FPKMSAM3", "_MINCOUNT10", "_ENTREZ.txt", sep="")
  C6_TXT_FN <- paste("./", compareList[[i]], "_DEG_UP_P",
                       0.05, "_Q0.05", "_FC0", "_FPKM1", "_FPKMAVG1",
                       "_FPKMSAM3", "_MINCOUNT10", "_TOP200C6.txt", sep="")
  C6_PDF_FN <- paste("./", compareList[[i]], "_DEG_UP_P",
                       0.05, "_Q0.05", "_FC0", "_FPKM1", "_FPKMAVG1",
                       "_FPKMSAM3", "_MINCOUNT10", "_TOP20C6.pdf", sep="")
  
  
  tmp <- readLines(currInFN)
  numLines <- length(tmp) - 1
  if(numLines > 0) { ## MAKE SURE THERE IS DATA TO READ ##
    data <- as.matrix(read.table(currInFN, header=FALSE, sep="\t"))
    data <- as.character(data)
    
    ## Convert Entrez IDs to Gene Symbols ##
    gene_symbols <- mapIds(org.Hs.eg.db, keys=data, column="SYMBOL", keytype="ENTREZID", multiVals="first")
    gene_symbols <- na.omit(gene_symbols)  # Remove NAs (unmapped IDs)
    
    if (length(gene_symbols) == 0) {
      message("No genes could be mapped to Gene Symbols. Skipping enrichment analysis.")
      next
    }
    
    # Load C6 gene sets
    c6_gene_sets <- msigdbr(species = "Homo sapiens", category = "C6")
    # Prepare gene set list for enricher
    c6_list <- split(c6_gene_sets$gene_symbol, c6_gene_sets$gs_name)
    
    # Perform enrichment analysis
    enrichedListC6 <- enricher(gene = gene_symbols, TERM2GENE = c6_gene_sets[, c("gs_name", "gene_symbol")], pvalueCutoff = 1.0, qvalueCutoff = 1.0)
    if (!is.null(enrichedListC6) && nrow(enrichedListC6@result) > 0) { 
      topC6 <- head(enrichedListC6@result, 200)
      names(topC6)[1] <- "C6 Pathway"
      topC6$GeneRatio <- gsub("/", " out of ", topC6$GeneRatio)
      topC6$BgRatio <- gsub("/", " out of ", topC6$BgRatio)
      topC6$geneID <- gsub("/", "///", topC6$geneID)
      
      write.table(topC6, file=C6_TXT_FN, sep="\t", row.names=FALSE)
      
      p <- barplot(enrichedListC6, showCategory=20)
      pdf(C6_PDF_FN, width=9, height=7)
      print(p)
      dev.off()
    }
  }
}



for(i in 1:numCompare) {
  currInFN <- paste("./", compareList[[i]], "_DEG_DN_P",
                    0.05, "_Q0.05", "_FC0", "_FPKM1", "_FPKMAVG1",
                    "_FPKMSAM3", "_MINCOUNT10", "_ENTREZ.txt", sep="")
  C6_TXT_FN <- paste("./", compareList[[i]], "_DEG_DN_P",
                     0.05, "_Q0.05", "_FC0", "_FPKM1", "_FPKMAVG1",
                     "_FPKMSAM3", "_MINCOUNT10", "_TOP20C6.txt", sep="")
  C6_PDF_FN <- paste("./", compareList[[i]], "_DEG_DN_P",
                     0.05, "_Q0.05", "_FC0", "_FPKM1", "_FPKMAVG1",
                     "_FPKMSAM3", "_MINCOUNT10", "_TOP20C6.pdf", sep="")
  
  tmp <- readLines(currInFN)
  numLines <- length(tmp) - 1
  if(numLines > 0) { ## MAKE SURE THERE IS DATA TO READ ##
    data <- as.matrix(read.table(currInFN, header=FALSE, sep="\t"))
    data <- as.character(data)
    
    ## Convert Entrez IDs to Gene Symbols ##
    gene_symbols <- mapIds(org.Hs.eg.db, keys=data, column="SYMBOL", keytype="ENTREZID", multiVals="first")
    gene_symbols <- na.omit(gene_symbols)  # Remove NAs (unmapped IDs)
    
    if (length(gene_symbols) == 0) {
      message("No genes could be mapped to Gene Symbols. Skipping enrichment analysis.")
      next
    }
    
    # Load C6 gene sets
    c6_gene_sets <- msigdbr(species = "Homo sapiens", category = "C6")
    # Prepare gene set list for enricher
    c6_list <- split(c6_gene_sets$gene_symbol, c6_gene_sets$gs_name)
    
    # Perform enrichment analysis
    enrichedListC6 <- enricher(gene = gene_symbols, TERM2GENE = c6_gene_sets[, c("gs_name", "gene_symbol")], pvalueCutoff = 1.0, qvalueCutoff = 1.0)
    if (!is.null(enrichedListC6) && nrow(enrichedListC6@result) > 0) { 
      topC6 <- head(enrichedListC6@result, 200)
      names(topC6)[1] <- "C6 Pathway"
      topC6$GeneRatio <- gsub("/", " out of ", topC6$GeneRatio)
      topC6$BgRatio <- gsub("/", " out of ", topC6$BgRatio)
      topC6$geneID <- gsub("/", "///", topC6$geneID)
      
      write.table(topC6, file=C6_TXT_FN, sep="\t", row.names=FALSE)
      
      p <- barplot(enrichedListC6, showCategory=20)
      pdf(C6_PDF_FN, width=9, height=7)
      print(p)
      dev.off()
    }
  }
}


sessionInfo()

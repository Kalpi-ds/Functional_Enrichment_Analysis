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


setwd("E:/KY-INBRE PostDoc/Yvonne/Matt/ClusterProfiler_1updnfix")




#compareList <- c("CS2_8WK__VS__CS1_8WK", "CS2_17WK__VS__CS1_17WK", "CS2_28WK__VS__CS1_28WK", "CS3_8WK__VS__CS1_8WK", "CS3_17WK__VS__CS1_17WK", "CS3_28WK__VS__CS1_28WK", "CS4_17WK__VS__CS1_17WK", "CS4_28WK__VS__CS1_28WK", "CS7_8WK__VS__CS1_8WK", "CS7_17WK__VS__CS1_17WK", "CS7_28WK__VS__CS1_28WK", 
                # "CS12_8WK__VS__CS1_8WK", "CS12_17WK__VS__CS1_17WK", "CS12_28WK__VS__CS1_28WK", "CS14_8WK__VS__CS13_8WK", "CS14_17WK__VS__CS13_17WK", "CS14_28WK__VS__CS13_28WK", "CS13_8WK__VS__CS1_8WK", "CS13_17WK__VS__CS1_17WK", "CS13_28WK__VS__CS1_28WK", "CS14_8WK__VS__CS2_8WK", "CS14_17WK__VS__CS2_17WK",
                # "CS14_28WK__VS__CS2_17WK")

compareList <- c("CS2_17WK__VS__CS1_17WK", "CS2_28WK__VS__CS1_28WK", "CS3_17WK__VS__CS1_17WK", "CS3_28WK__VS__CS1_28WK", "CS4_17WK__VS__CS1_17WK", "CS4_28WK__VS__CS1_28WK", "CS12_17WK__VS__CS1_17WK",
                 "CS12_28WK__VS__CS1_28WK", "CS13_17WK__VS__CS1_17WK", "CS14_28WK__VS__CS2_17WK", "CS14_28WK__VS__CS13_28WK")
######################################
## READ IN DIFFERENTIALLY EXPRESSED ##
## ENTREZ IDENTIFIERS AND CREATE    ##
## GENE LISTS FOR CATEGORICAL       ##
## ENRICHMENT                       ##
######################################
numCompare <- 11
for(i in 1:numCompare) {
  currInFN <- paste("./", compareList[[i]], "_DEG_UP_P",
                    0.05, "_Q0.05", "_FC0", "_FPKM1", "_FPKMAVG1",
                    "_FPKMSAM3", "_MINCOUNT10", "_ENTREZ.txt", sep="")
  KEGG_TXT_FN <- paste("./", compareList[[i]], "_DEG_UP_P",
                       0.05, "_Q0.05", "_FC0", "_FPKM1", "_FPKMAVG1",
                       "_FPKMSAM3", "_MINCOUNT10", "_TOP200KEGG.txt", sep="")
  KEGG_PDF_FN <- paste("./", compareList[[i]], "_DEG_UP_P",
                       0.05, "_Q0.05", "_FC0", "_FPKM1", "_FPKMAVG1",
                       "_FPKMSAM3", "_MINCOUNT10", "_TOP20KEGG.pdf", sep="")
  
  
  tmp <- readLines(currInFN)
  numLines <- length(tmp) - 1
  if(numLines > 0) { ## MAKE SURE THERE IS DATA TO READ ##
    data <- as.matrix(read.table(currInFN, header=FALSE, sep="\t"))
    data <- as.character(data)
    
    ## FIRST DO KEGG ANALYSIS ##
    enrichedListKEGG <- enrichMKEGG(gene=data, organism='hsa', pvalueCutoff=1.0, qvalueCutoff=1.0)
    if(length(enrichedListKEGG) > 0) { 
      topKEGG <- head(enrichedListKEGG, 200)
      names(topKEGG)[1] <- "KEGG ID"
      topKEGG$GeneRatio <- gsub("/", " out of ", topKEGG$GeneRatio)
      topKEGG$BgRatio <- gsub("/", " out of ", topKEGG$BgRatio)
      topKEGG$geneID <- gsub("/", "///", topKEGG$geneID)
      
      write.table(topKEGG, file=KEGG_TXT_FN, sep="\t", row.names=FALSE)
      
      p <- barplot(enrichedListKEGG, showCategory=20)
      pdf(KEGG_PDF_FN, width=9, height=7)
      print(p)
      dev.off()
    }
  }
}

compareList <- c("CS7_17WK__VS__CS1_17WK")

numCompare <- 1
for(i in 1:numCompare) {
  currInFN <- paste("./", compareList[[i]], "_DEG_UP_P",
                    0.05, "_Q0.05", "_FC0", "_FPKM1", "_FPKMAVG1",
                    "_FPKMSAM3", "_MINCOUNT10", "_ENTREZ.txt", sep="")
    GO_TXT_FN <- paste("./", compareList[[i]], "_DEG_UP_P",
                     0.05, "_Q0.05", "_FC0", "_FPKM1", "_FPKMAVG1",
                     "_FPKMSAM3", "_MINCOUNT10", "_TOP200GO.txt", sep="")
  GO_PDF_FN <- paste("./", compareList[[i]], "_DEG_UP_P",
                     0.05, "_Q0.05", "_FC0", "_FPKM1", "_FPKMAVG1",
                     "_FPKMSAM3", "_MINCOUNT10", "_TOP20GO.pdf", sep="")
  
  tmp <- readLines(currInFN)
  numLines <- length(tmp) - 1
  if(numLines > 0) { ## MAKE SURE THERE IS DATA TO READ ##
    data <- as.matrix(read.table(currInFN, header=FALSE, sep="\t"))
    data <- as.character(data)
    
    
    ## DO GO ANALYSIS ##
    enrichedListGO <- enrichGO(gene=data, OrgDb='org.Hs.eg.db', ont="BP", pvalueCutoff=1.0, qvalueCutoff=1.0)
    if(length(enrichedListGO) > 0) { 
      topGO <- head(enrichedListGO, 200)
      names(topGO)[1] <- "GO ID"
      topGO$GeneRatio <- gsub("/", " out of ", topGO$GeneRatio)
      topGO$BgRatio <- gsub("/", " out of ", topGO$BgRatio)
      topGO$geneID <- gsub("/", "///", topGO$geneID)
      
      write.table(topGO, file=GO_TXT_FN, sep="\t", row.names=FALSE)
      p <- barplot(enrichedListGO, showCategory=20)
      pdf(GO_PDF_FN, width=9, height=7)
      print(p)
      dev.off()
    }
  }
}


sessionInfo()

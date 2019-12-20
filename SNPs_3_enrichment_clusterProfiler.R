setwd("data/")

lal_snps <- unlist(read.csv("lookalike_genes.txt", stringsAsFactors = F, header = F))
bck_snps <- unlist(read.csv("all_infinium_genes.txt", stringsAsFactors = F, header = F))

library(clusterProfiler)
library(org.Hs.eg.db)


egobp <- enrichGO(gene = lal_snps, universe = bck_snps,  keyType = "ENSEMBL", OrgDb = org.Hs.eg.db,
                  minGSSize = 1, maxGSSize = 22000,
                  ont = "BP",
                  pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.5, readable = T)

egomf <- enrichGO(gene = lal_snps, universe = bck_snps,  keyType = "ENSEMBL", OrgDb = org.Hs.eg.db,
                  minGSSize = 1, maxGSSize = 22000,
                  ont = "MF",
                  pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.5, readable = T)

egocc <- enrichGO(gene = lal_snps, universe = bck_snps,  keyType = "ENSEMBL", OrgDb = org.Hs.eg.db,
                  minGSSize = 1, maxGSSize = 22000,
                  ont = "CC",
                  pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.5, readable = T)


for(go in c("bp", "mf", "cc")) {
  ego <- get(paste0("ego", go))
  
  write.table(ego, file = paste0("clusterProfiler_results/", toupper(go), "cluster_profiler_enrichGO_table.tsv"), 
              row.names = F,
              quote = F, sep = "\t")
}

#### PLOTS ####
# Biological Process
  pdf("clusterProfiler_results/emapplot_all_BP_no_names.pdf")
  emapplot(egobp)
  dev.off()
  
  # TOP 10
  pdf("clusterProfiler_results/emapplot_top10_BP_no_names.pdf")
  ep <- emapplot(egobp, showCategory = 10)
  ep$data[, "name"] <- ""
  plot(ep)
  dev.off()

# Molecular function
  pdf("clusterProfiler_results/emapplot_all_MF_no_names.pdf")
  emapplot(egomf)
  dev.off()
  
  # TOP 10
  pdf("clusterProfiler_results/emapplot_top10_MF_no_names.pdf")
  ep <- emapplot(egomf, showCategory = 10)
  ep$data[, "name"] <- ""
  plot(ep)
  dev.off()
  
  
# Cellular compartment
  pdf(paste0("clusterProfiler_results/emapplot_all_CC_no_names.pdf"))
  emapplot(egocc)
  dev.off()
  
  # TOP 10
  pdf(paste0("clusterProfiler_results/emapplot_top10_CC_no_names.pdf"))
  ep <- emapplot(egocc, showCategory = 10 )
  ep$data[, "name"] <- ""
  plot(ep)
  dev.off()



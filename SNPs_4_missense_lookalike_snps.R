## SNP analysis
setwd("data")

# Load table with Look-alike snps (the ones with shared genotype in all pairs)
load("merged_snp9_infinium_table.RData", verbose = T)

## rs snps
rs_snps <- unique(all_snp_information$Name[grepl("rs",all_snp_information$Name)])


# Some lines are duplicated except for the name where the snp is annotated (kgp or rs)
# I remove these lines
all_snp_information <- all_snp_information[!duplicated(all_snp_information[, colnames(all_snp_information) != "Name"]),]

# Types of snp
vartypes <- unlist(all_snp_information$Mutation.s.)
vartypes <- unique(vartypes)
vartypes <- unlist(strsplit(vartypes,  split = ","))
vartypes <- unique(vartypes)
# Silent, Missense, Synonymous

# Missense 
missense <- all_snp_information[grepl(pattern = "Missense", x = all_snp_information$Mutation.s.), ]

# write.table(missense, file = "missense_SNPs.tsv", row.names = F, quote = F, sep = "\t")


genes_with_missense <- c()
for(i in 1:nrow(missense)) {
  mut_type <- strsplit(missense[i, "Mutation.s."], split = ",")[[1]]
  genes    <- strsplit(missense[i, "Gene.s."], split = ",")[[1]]
  genes_with_missense <- c(genes_with_missense, unique(genes[grepl("Missense", mut_type)]))
}
genes_with_missense <- unique(genes_with_missense)



## Enrichment of genes with missense 

# Defining background (all SNPs in SNP array)
InfiniumOmni5 <- read.table("InfiniumOmni5-4v1-2_A1.annotated.txt", 
                            stringsAsFactors = F,
                            header = T, 
                            fill = T, # I need to do this
                            sep = "\t")

all_genes_infinium <- unique(InfiniumOmni5$Gene.s.)
all_genes_infinium <- unlist(strsplit(all_genes_infinium, split = ","))
all_genes_infinium <- unique(all_genes_infinium)


infinium_ranges <- makeGRangesFromDataFrame(InfiniumOmni5, seqnames.field = "Chr", start.field = "MapInfo",
                                            end.field = "MapInfo")

background_SNPs <- infinium_ranges[!duplicated(ranges(infinium_ranges))]
background_SNPs <- background_SNPs[!seqnames(background_SNPs) %in% c("M", "0")]

seqlevelsStyle(background_SNPs) <- "UCSC"


output_folder <- "missense_SNPs_enrichment/"

if (!file.exists(output_folder)) {
  dir.create(output_folder)
}
 
# cat(paste0(sort(genes_with_missense)), sep = "\n", 
#     file = paste0(output_folder, "genes_with_lookalike_missense_SNPs.txt"))
# cat(paste0(sort(all_genes_infinium)), sep = "\n",  
#     file = paste0(output_folder,  "background_for_genes_with_lookalike_missense_SNPs_all_infinium_genes.txt"))


ens_to_hugo <- read.table("mart_export_75_ENS_to_HUGO.csv", header = T, sep = ",",
                          row.names = 1,
                          stringsAsFactors = F)

ENS_genes_with_missense <- rownames(ens_to_hugo)[ens_to_hugo$HGNC.symbol %in% genes_with_missense]
ENS_all_genes_infinium <- rownames(ens_to_hugo)[ens_to_hugo$HGNC.symbol %in% all_genes_infinium]

# cat(paste0(sort(ENS_genes_with_missense)), sep = "\n", 
#     file = paste0(output_folder, "ENS_genes_with_lookalike_missense_SNPs.txt"))
# cat(paste0(sort(ENS_all_genes_infinium)), sep = "\n",  
#     file = paste0(output_folder,  "ENS_background_for_genes_with_lookalike_missense_SNPs_all_infinium_genes.txt"))


#### MISSENSE ENRICHMENT WITH CLUSTERPROFILER ####
library(clusterProfiler)
library(org.Hs.eg.db)


egobp <- enrichGO(gene = ENS_genes_with_missense, universe = ENS_all_genes_infinium,  keyType = "ENSEMBL", OrgDb = org.Hs.eg.db,
                  minGSSize = 1, maxGSSize = 22000,
                  ont = "BP",
                  pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1, readable = T)

egomf <- enrichGO(gene = ENS_genes_with_missense, universe = ENS_all_genes_infinium,  keyType = "ENSEMBL", OrgDb = org.Hs.eg.db,
                  minGSSize = 1, maxGSSize = 22000,
                  ont = "MF",
                  pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1, readable = T)

egocc <- enrichGO(gene = ENS_genes_with_missense, universe = ENS_all_genes_infinium,  keyType = "ENSEMBL", OrgDb = org.Hs.eg.db,
                  minGSSize = 1, maxGSSize = 22000,
                  ont = "CC",
                  pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1, readable = T)



# TOP 10 PROCESSES
for(go in c("bp", "mf", "cc")) {
  ego <- get(paste0("ego", go))
  
  write.table(head(ego, 10), file = paste0(output_folder, "/missense_", toupper(go), "cluster_profiler_enrichGO_table.tsv"), 
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
# ep$data[, "name"] <- ""
plot(ep)
dev.off()

# Molecular function
pdf("clusterProfiler_results/emapplot_all_MF_no_names.pdf")
emapplot(egomf)
dev.off()

# TOP 10
pdf("clusterProfiler_results/emapplot_top10_MF_no_names.pdf")
ep <- emapplot(egomf, showCategory = 10)
# ep$data[, "name"] <- ""
plot(ep)
dev.off()


# Cellular compartment
pdf(paste0("clusterProfiler_results/emapplot_all_CC_no_names.pdf"))
emapplot(egocc)
dev.off()

# TOP 10
pdf(paste0("clusterProfiler_results/emapplot_top10_CC_no_names.pdf"))
ep <- emapplot(egocc, showCategory = 10 )
# ep$data[, "name"] <- ""
plot(ep)
dev.off()










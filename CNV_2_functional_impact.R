# FUNCTIONAL IMPACT OF CNVs
library(GenomicRanges)
library(gplots)

setwd("data/")


samples <- substr(list.files("tmp"), 1,10)
samples <- unique(samples[grepl("LAL", samples)])

pair_ids <- sort(unique(substr(samples,1, 8)))

list.files("tmp/")


# Biomart not possible (libcurl not available for R verion 3.6.1)
# library("biomaRt")
# ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)
# # now we select a Dataset
# genes <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol", "chromosome_name","start_position","end_position", "description"), filters = "ensembl_gene_id", values = names(asdf), mart = ensembl)
# rownames(genes) <- genes$hgnc_symbol

# I download gene coordinates manually:
# Ensembl v75, build GRCh37
# http://feb2014.archive.ensembl.org/biomart/martview/47ddba8032a07f5fb3f01ffea22e9a49
all_genes <- read.table("mart_export_genes_grch37_v75.txt", sep = "\t", header = T)


### ANALYSIS OF PROTEIN-CODING GENES
pc_genes <- makeGRangesFromDataFrame(all_genes[all_genes$Gene.Biotype == "protein_coding", ], seqnames.field = "Chromosome.Name",
                                     start.field = "Gene.Start..bp.", end.field = "Gene.End..bp.")

names(pc_genes) <- all_genes[all_genes$Gene.Biotype == "protein_coding", "Ensembl.Gene.ID"]



# Remove regions in which CN is shared in all individuals
# Keep regions that coincide by pair
all_genes_CN <- data.frame(Gene = names(pc_genes), stringsAsFactors = F)
rownames(all_genes_CN) <- all_genes_CN$Gene


for(pair in pair_ids) {
  for(AB in c("A", "B")) {
    load(paste0("tmp/", pair, ".", AB, "breakpoint_table.RData"), verbose = F) 
    # I change chromosome XY (a pseudoautosomal region) for X,  (Only coordinates from PAR1, 
    # which are interchangeable with Y, as it is located at the begining of the chromosome)
    table_per_sample$chr[table_per_sample$chr == "XY"] <- "X"
    granges <- makeGRangesFromDataFrame(table_per_sample, keep.extra.columns = T)
    # granges <- granges[granges$CN != 2]
    fo <- findOverlaps(pc_genes, granges, type = "within")
    cn <- granges$CN[subjectHits(fo)]
    names(cn) <- names(pc_genes)[queryHits(fo)]
    all_genes_CN[names(cn), paste(pair, AB, sep = ".")] <- cn
  }
}

# Remove if all NAs
all_genes_CN <- all_genes_CN[apply(all_genes_CN[, 2:ncol(all_genes_CN)], 1, function(x) !all(is.na(x))),]

# Remove genes with same CN in all samples
all_genes_CN <-all_genes_CN[apply(all_genes_CN[,2:ncol(all_genes_CN)], 1, function(x) length(unique(x))) != 1, ]


shared_in_pair <- c()
for(gene in rownames(all_genes_CN)) {
  x <- all_genes_CN[gene, 2:ncol(all_genes_CN)]
  x[is.na(x)] <- 999
  line_shared_in_pair <- c()
  for(pair in pair_ids) {
    line_shared_in_pair <- c(line_shared_in_pair,
                             x[, paste0(pair, ".A")] == x[, paste0(pair, ".B")])

  }
  shared_in_pair <- rbind(shared_in_pair,
                          line_shared_in_pair)
}
rownames(shared_in_pair) <- rownames(all_genes_CN)

# Select genes shared in all pairs
shared_in_pair <- shared_in_pair[apply(shared_in_pair, 1, all),] 

# remove genes in X or Y chromosomes 
shared_in_pair <- shared_in_pair[ !rownames(shared_in_pair) %in% names(pc_genes)[as.character(seqnames(pc_genes)) %in% c("X", "Y")],]


genes_with_shared_CN_in_all_pairs <- all_genes_CN[rownames(shared_in_pair), ]
genes_with_shared_CN_in_all_pairs <- genes_with_shared_CN_in_all_pairs[complete.cases(genes_with_shared_CN_in_all_pairs), ]

### lists for Gene set enrichment analysis
# write.table(matrix(names(pc_genes)[as.character(seqnames(pc_genes)) %in% 1:22], ncol = 1), 
#             quote = F, row.names = F, col.names = F,
#             file = "autosomal_pc_genes_grch37.txt")
# write.table(matrix(rownames(genes_with_shared_CN_in_all_pairs), ncol = 1), 
#             quote = F, row.names = F, col.names = F,
#             file = "genes_with_shared_CN_in_all_pairs.txt")

cn_lal_genes <- read.table("genes_with_shared_CN_in_all_pairs.txt", stringsAsFactors = F)[, 1]
bck          <- read.csv("../snp_analysis/all_infinium_genes.txt", stringsAsFactors = F, header = F)[,1]

## OR 
# cn_lal_genes <- rownames(genes_with_shared_CN_in_all_pairs)
# bck <- names(pc_genes)[as.character(seqnames(pc_genes)) %in% 1:22]

## ENRICHMENT ANALYSIS WITH clusterProfiler
library(clusterProfiler)
library(org.Hs.eg.db)
egobp <- enrichGO(gene = cn_lal_genes, universe = bck,  keyType = "ENSEMBL", OrgDb = org.Hs.eg.db,
                  minGSSize = 1, maxGSSize = 22000,
                  ont = "BP",
                  pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.5, readable = T)

egomf <- enrichGO(gene = cn_lal_genes, universe = bck,  keyType = "ENSEMBL", OrgDb = org.Hs.eg.db,
                  minGSSize = 1, maxGSSize = 22000,
                  ont = "MF",
                  pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.5, readable = T)

egocc <- enrichGO(gene = cn_lal_genes, universe = bck,  keyType = "ENSEMBL", OrgDb = org.Hs.eg.db,
                  minGSSize = 1, maxGSSize = 22000,
                  ont = "CC",
                  pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.5, readable = T)

emapplot(egobp)
emapplot(egomf)
emapplot(egocc)



# PLOT (HEATMAP)
m <- as.matrix(genes_with_shared_CN_in_all_pairs[, 2:ncol(genes_with_shared_CN_in_all_pairs)])

m <- m[order(as.integer(seqnames(pc_genes[rownames(genes_with_shared_CN_in_all_pairs)]))),]

ens_to_hugo <- read.table("~/_projects/A_FINAL_VERSION_CNVs/BSC/mart_export_75_ENS_to_HUGO.csv", stringsAsFactors = F,
                          sep = ",", row.names = 1,
                          header = T)

new_names <- ens_to_hugo[rownames(m), "HGNC.symbol"]
new_names[is.na(new_names)] <- rownames(m)[is.na(new_names)]
rownames(m) <- new_names

# pdf(file = "heatmap_shared_CN.pdf", height = 10, width = 11)
heatmap.2(m,
          Rowv = F,
          cexCol = 1.05,
          cexRow = 1.05,
          # Colv = F,
          dendrogram = "none",
          colsep = 1:32, rowsep = 1:24,
          col = colorRampPalette(c("white", "royalblue3"))(5)[c(2,5)],
          sepwidth = c(0.01, 0.01),
          margins = c(10,10),
          sepcolor = "white", trace = "none")
# dev.off()





##### Analysis for non-coding genes ####
# RNAs 
rnas <- as.character(unique(all_genes$Gene.Biotype)[grepl("rna", unique(all_genes$Gene.Biotype), ignore.case = T)])

NC_genes <- makeGRangesFromDataFrame(all_genes[all_genes$Gene.Biotype %in% rnas, ], seqnames.field = "Chromosome.Name",
                                     start.field = "Gene.Start..bp.", end.field = "Gene.End..bp.")

names(NC_genes) <- all_genes[all_genes$Gene.Biotype %in% rnas, "Ensembl.Gene.ID"]



# Remove regions in which CN is shared in all individuals
# Keep regions that coincide by pair
all_genes_CN <- data.frame(Gene = names(NC_genes), stringsAsFactors = F)
rownames(all_genes_CN) <- all_genes_CN$Gene


for(pair in pair_ids) {
  for(AB in c("A", "B")) {
    load(paste0("tmp/", pair, ".", AB, "breakpoint_table.RData"), verbose = F) 
    # I change chromosome XY (a pseudoautosomal region) for X,  (Only coordinates from PAR1, 
    # which are interchangeable with Y, as it is located at the begining of the chromosome)
    table_per_sample$chr[table_per_sample$chr == "XY"] <- "X"
    granges <- makeGRangesFromDataFrame(table_per_sample, keep.extra.columns = T)
    # granges <- granges[granges$CN != 2]
    fo <- findOverlaps(NC_genes, granges, type = "within")
    cn <- granges$CN[subjectHits(fo)]
    names(cn) <- names(NC_genes)[queryHits(fo)]
    all_genes_CN[names(cn), paste(pair, AB, sep = ".")] <- cn
  }
}

# remove lines with all NAs
all_genes_CN <- all_genes_CN[apply(all_genes_CN[, 2:ncol(all_genes_CN)], 1, function(x) !all(is.na(x))),]

# Remove genes with same CN in all samples
all_genes_CN <-all_genes_CN[apply(all_genes_CN[,2:ncol(all_genes_CN)], 1, function(x) length(unique(x))) != 1, ]


shared_in_pair <- c()
for(gene in rownames(all_genes_CN)) {
  x <- all_genes_CN[gene, 2:ncol(all_genes_CN)]
  x[is.na(x)] <- 999
  line_shared_in_pair <- c()
  for(pair in pair_ids) {
    line_shared_in_pair <- c(line_shared_in_pair,
                             x[, paste0(pair, ".A")] == x[, paste0(pair, ".B")])
    
  }
  shared_in_pair <- rbind(shared_in_pair,
                          line_shared_in_pair)
}
rownames(shared_in_pair) <- rownames(all_genes_CN)

# Select genes shared in all pairs
shared_in_pair <- shared_in_pair[apply(shared_in_pair, 1, all),] 

# remove genes in X or Y chromosomes 
shared_in_pair <- shared_in_pair[ !rownames(shared_in_pair) %in% names(NC_genes)[as.character(seqnames(NC_genes)) %in% c("X", "Y")],]


genes_with_shared_CN_in_all_pairs <- all_genes_CN[rownames(shared_in_pair), ]
genes_with_shared_CN_in_all_pairs <- genes_with_shared_CN_in_all_pairs[complete.cases(genes_with_shared_CN_in_all_pairs), ]


# write.table(matrix(names(NC_genes)[as.character(seqnames(NC_genes)) %in% 1:22], ncol = 1), 
#             quote = F, row.names = F, col.names = F,
#             file = "autosomal_RNA_genes_grch37.txt")
# write.table(matrix(rownames(genes_with_shared_CN_in_all_pairs), ncol = 1), 
#             quote = F, row.names = F, col.names = F,
#             file = "RNA_genes_with_shared_CN_in_all_pairs.txt")




## PLOT 
m <- as.matrix(genes_with_shared_CN_in_all_pairs[, 2:ncol(genes_with_shared_CN_in_all_pairs)])
m <- m[order(as.integer(seqnames(NC_genes[rownames(genes_with_shared_CN_in_all_pairs)]))),]
# pdf(file = "heatmap_shared_CN_RNA_genes.pdf", height = 10, width = 11)
heatmap.2(m, 
          Rowv = F, 
          cexCol = 1.05,
          cexRow = 1.05,
          # Colv = F, 
          dendrogram = "none", 
          colsep = 1:ncol(m), rowsep = 1:nrow(m), 
          col = colorRampPalette(c("white", "royalblue3"))(5)[c(2,5)],
          sepwidth = c(0.01, 0.01),
          margins = c(10,10),
          sepcolor = "white", trace = "none")
# dev.off()

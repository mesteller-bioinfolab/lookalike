## SNP analysis with gwas catalog
setwd("data/")

library(GenomicRanges)

load("merged_snp9_infinium_table.RData", verbose = T)

# Some lines are duplicated except for the name where the snp is annotated (kgp or rs)
# I remove these lines
all_snp_information <- all_snp_information[!duplicated(all_snp_information[, colnames(all_snp_information) != "Name"]),]

lookalike_snps <- makeGRangesFromDataFrame(df = all_snp_information, 
                                           seqnames.field = "V1", 
                                           start.field = "V2", 
                                           end.field = "V3", 
                                           keep.extra.columns = T)
seqlevelsStyle(lookalike_snps) <- "UCSC"

load("gwas_catalog_filtered_and_lifted_to_GRCh37.RData")

gwas <- makeGRangesFromDataFrame(df = gwas_catalog_filtered[, c("gwas_id_maria",
                                                                "chr_build37",
                                                                "pos_build37",
                                                                "DISEASE.TRAIT",
                                                                "MAPPED_TRAIT",
                                                                "CONTEXT",
                                                                "MAPPED_GENE",
                                                                "UPSTREAM_GENE_ID",
                                                                "DOWNSTREAM_GENE_ID",
                                                                "SNP_GENE_IDS",
                                                                "LINK",
                                                                "CONTEXT")],
                                 seqnames.field = "chr_build37", start.field = "pos_build37", end.field = "pos_build37", 
                                 keep.extra.columns = T)

# Find overlaps between look-alike snps and GWAS catalog
fo <- findOverlaps(lookalike_snps, gwas)

snps    <- as.data.frame(lookalike_snps[queryHits(fo)], row.names = NULL)
gwas_df <- as.data.frame(gwas[subjectHits(fo)])

snps    <- snps[, c("Name", "Chr", "start", "end", "Alleles")]
gwas_df <- gwas_df[, c("DISEASE.TRAIT", "MAPPED_TRAIT", "MAPPED_GENE", "LINK")]
table_to_save <- cbind(snps, gwas_df)
table_to_save <- table_to_save[!duplicated(table_to_save),]
table_to_save$LINK <- gsub(table_to_save$LINK, pattern = "www.ncbi.nlm.nih.gov/pubmed/", replacement = "")

# Some rows are repeated, with only PMID changing
table_to_save$repeated <- as.integer(as.factor(apply(table_to_save[, 1:8], 1, function(x) paste0(x, collapse = "_"))))


pmids <- unlist(lapply(by(table_to_save$LINK, table_to_save$repeated, as.character), function(x) paste(x, collapse = ", ")))
table_to_save$PMID <- pmids[as.character(table_to_save$repeated)]

table_to_save <- table_to_save[, !colnames(table_to_save) %in% c("end", "LINK", "repeated", "DISEASE.TRAIT")]
table_to_save <- table_to_save[!duplicated(table_to_save), ]

colnames(table_to_save)[colnames(table_to_save) == "start"] <- "Position"
colnames(table_to_save)[colnames(table_to_save) == "MAPPED_TRAIT"] <- "Mapped.trait"
colnames(table_to_save)[colnames(table_to_save) == "MAPPED_GENE"] <- "Mapped.gene"

write.table(table_to_save, row.names = F, quote = F, sep = "\t", file = "SNPS_in_GWAS_catalog.tsv" )


#### custom ENRICHMENT analysis
InfiniumOmni5 <- read.table("InfiniumOmni5-4v1-2_A1.annotated.txt", 
                            stringsAsFactors = F,
                            header = T, 
                            fill = T, # I need to do this
                            sep = "\t")

all_genes_infinium <- unique(InfiniumOmni5$Gene.s.)
all_genes_infinium <- unlist(strsplit(all_genes_infinium, split = ","))
all_genes_infinium <- unique(all_genes_infinium)

InfiniumOmni5 <- InfiniumOmni5[!duplicated(InfiniumOmni5[, 2:3]), ]


infinium_ranges <- makeGRangesFromDataFrame(InfiniumOmni5, seqnames.field = "Chr", start.field = "MapInfo",
                                            end.field = "MapInfo")
seqlevelsStyle(infinium_ranges) <- "UCSC"

rm(InfiniumOmni5)


#### ENRICHMENT
traits <- unique(gwas$MAPPED_TRAIT)

### Loop to make contingency tables ### 

# table_by_trait <- list()
# length(traits)
# for(trait in traits) {
#   snp_traits <- gwas[gwas$MAPPED_TRAIT == trait]
#   df <- data.frame(lookalike_snp = overlapsAny(infinium_ranges, lookalike_snps),
#                    trait_snp     = overlapsAny(infinium_ranges, snp_traits))
#   table_by_trait[[trait]] <- table(df)
#   message(paste(which(traits == trait), "of", length(traits)))
# }
# 
# save(table_by_trait, file = "contingency_tables_by_trait_background_all_snps.RData")


load("contingency_tables_by_trait_background_all_snps.RData")

for(i in 1:length(table_by_trait)) {
  tt <- table_by_trait[[i]]
  if(!"TRUE" %in% colnames(tt)) {
    tt <- cbind(tt, "TRUE" = c(0,0))
  }
  table_by_trait[[i]] <- tt
}

pvals <- unlist(lapply(table_by_trait, function(x) fisher.test(x)$p.value))

# The text (very long names) is giving me problems
pvals_unnamed <- unname(pvals)
adjusted_pvals <- p.adjust(pvals_unnamed, method = "BH")

df <- data.frame(trait = names(pvals),
                 pval = pvals, 
                 adj_pval = adjusted_pvals)[pvals < 0.05, ]
df <- df[order(df$pval, decreasing = F), ]

# write.table(df, quote = F, sep = "\t", row.names = F,
#             file = "significant_enrichment_all_snps_background.tsv")






### Manual annotation of what traits can be seen at first sight

gwas_common <- gwas[overlapsAny(gwas, infinium_ranges)]

gwas_common$likely_to_be_seen <- as.logical(NA)
# 2 will mean very likely to be visually obvious trait, 0 not observable, 1, doubt

all_traits <- rep(as.logical(NA), length(unique(gwas_common$MAPPED_TRAIT)))
names(all_traits) <- unique(gwas_common$MAPPED_TRAIT)

# Before really putting a 0, I check them manually
all_traits[grepl("Tumor", names(all_traits))] <- F
all_traits[grepl("tumor", names(all_traits))] <- F
all_traits[grepl("cancer", names(all_traits))] <- F
all_traits[grepl("carcinoma", names(all_traits))] <- F
all_traits[grepl("benign", names(all_traits))] <- F
all_traits[grepl("count", names(all_traits))] <- F
all_traits[grepl("disease", names(all_traits))] <- F
all_traits[grepl("impairment", names(all_traits))] <- F
all_traits[grepl("measurement", names(all_traits))] <- F
all_traits[grepl("syndrome", names(all_traits))] <- F
all_traits[grepl("Affymetrix", names(all_traits))] <- F
all_traits[grepl("response", names(all_traits))] <- F

# Face
all_traits[grepl("morphology", names(all_traits))] <- T
all_traits[grepl("height", names(all_traits))] <- T
all_traits[grepl("waist-hip ratio", names(all_traits))] <- T
all_traits[grepl("facial", names(all_traits))] <- T
all_traits[grepl("anthropometric", names(all_traits))] <- T
all_traits[grepl("periodontal measurement", names(all_traits))] <- T
all_traits[grepl("hair shape measurement", names(all_traits))] <- T
all_traits[grepl("balding measurement", names(all_traits))] <- T
all_traits[grepl("cartilage thickness measurement", names(all_traits))] <- T
all_traits[grepl("optic disc size measurement", names(all_traits))] <- T
all_traits[grepl("eye measurement", names(all_traits))] <- T
all_traits[grepl("cleft palate", names(all_traits))] <- T
all_traits[grepl("Alopecia", names(all_traits))] <- T
all_traits[grepl("alopecia", names(all_traits))] <- T
all_traits[grepl("eye color", names(all_traits))] <- T
all_traits[grepl("tooth agenesis", names(all_traits))] <- T
all_traits[grepl("suntan", names(all_traits))] <- T
all_traits[grepl("sunburn", names(all_traits))] <- T
all_traits[grepl("cleft lip", names(all_traits))] <- T
all_traits[grepl("skin sensitivity to sun", names(all_traits))] <- T
all_traits[grepl("skin pigmentation", names(all_traits))] <- T
all_traits[grepl("lobe attachment", names(all_traits))] <- T
all_traits[grepl("braces", names(all_traits))] <- T
all_traits[grepl("earlobe", names(all_traits))] <- T
all_traits[grepl("hair", names(all_traits))] <- T

# No chemotherapy-induced alopecia 
all_traits[grepl("chemotherapy", names(all_traits))] <- F


all_traits[is.na(all_traits)] <- F


gwas_common$likely_to_be_seen <- unname(all_traits[gwas_common$MAPPED_TRAIT])


df <- data.frame(lookalike_snp    = overlapsAny(gwas_common, lookalike_snps),
                 visible_gwas     = gwas_common$likely_to_be_seen)
table(df)
fisher.test(table(df))



### TRAITS with "morphology" 
all_traits_morphology <- all_traits
all_traits_morphology[] <- F 
all_traits_morphology[grepl("morphology", names(all_traits_morphology))] <- T

gwas_common$morphology <- unname(all_traits_morphology[gwas_common$MAPPED_TRAIT])


df <- data.frame(lookalike_snp    = overlapsAny(gwas_common, lookalike_snps),
                 morphology_gwas  = gwas_common$morphology)
table(df)
fisher.test(table(df))

gwas_common[overlapsAny(gwas_common, lookalike_snps) & gwas_common$morphology]


#######################
 
 # FaceBase dataset
 # https://www.facebase.org/chaise/recordset/#1/isa:dataset/*::facets::N4IghgdgJiBcDaoDOB7ArgJwMYFM6JAEsIAjdafIpMEAGhCjABcwkcmB9FDAc0kKQBbDoxZtOhKBwBmAaxwBPEAF0AvrVDomZNBQRUa9Ua3Zde-IWb4QBwuYpXqiMZfSwALFIVxJKAOQBhACEASQAVAEEADQB5P1gATgA2AAYkxzUgA@sort(accession::desc::)@after(FB00000988)
 

 facebase <- read.csv("FaceBase_Dataset.csv", stringsAsFactors = F)
 
 candidate_genes <- strsplit(facebase$title, fixed =  T, split = "Candidate Gene")
 candidate_genes <- unlist(lapply(candidate_genes[unlist(lapply(candidate_genes, length)) == 2], function(x) x[[2]]))
 candidate_genes <- gsub("s: ", "", candidate_genes)
 candidate_genes <- gsub(": ", "", candidate_genes)
 candidate_genes <- unique(unlist(strsplit(candidate_genes, fixed =  T, split = " & ")))
 
 
   
  ens_to_hugo <- read.table("mart_export_75_ENS_to_HUGO.csv", header = T, sep = ",",
                            row.names = 1,
                            stringsAsFactors = F)
 
 all_text_facebase <- unlist(c(facebase[, c("title", "summary", "X_keywords", "description")]))
 
 # Function to capitalize words
 simpleCap <- function(x) {
   s <- strsplit(x, " ")[[1]]
   paste(toupper(substring(s, 1, 1)), substring(s, 2),
         sep = "", collapse = " ")
 }
 
 
 gene_in_facebase <- c()
 for(x in 1:nrow(ens_to_hugo)) {
   hugo <- ens_to_hugo$HGNC.symbol[x]
   if((sum(grepl(hugo, all_text_facebase)) + sum(grepl(simpleCap(tolower(hugo)), all_text_facebase))) > 0) {
     gene_in_facebase <- c(gene_in_facebase, hugo)
   }
   
   if(x%%100 == 0 ) print(x)
 }
 
 gene_in_facebase <- gene_in_facebase[!gene_in_facebase == ""]
 
 candidate_genes[candidate_genes %in% gene_in_facebase]
 # IRF6, MAFB, ARHGAP29, 8q24, PAX7, VAX1, NTN1, NOG, FOXE1, MSX1, BMP4, FGFR2, PTCH1.
 
 # Manually remove what are not genes
 gene_in_facebase <- gene_in_facebase[!gene_in_facebase %in% c("ACR", "AGL", "ANG", "APP", "ARC",
                                                               "C2", "C3", "C5", "CAT", "CHIA", "CIT", "CP", "CS",
                                                               "DES", "DR1", "F11", "GDF1", 
                                                               "HAL", "HR", "IDE", "INA", "KL", "KY", "LAT", "LXN", "MAF", "MAL",
                                                               "MAX", "MB", "MET", "NIN", "NOV", "OLR1", "PC", 
                                                               "RAN", "REL", "RELA", "REM2", "REN", "RET", 
                                                               "SAG", "SET", "SHE", "SI", "SMO", "T", "TBX1", "TH", "TRIO")]
 

 df <- data.frame(lookalike = all_genes_infinium %in% lookalike_genes,
                  facial   = all_genes_infinium %in% unique(c(candidate_genes, gene_in_facebase)))
 
 fisher.test(table(df))
 
 
 # Exploratory genotype–phenotype correlations of facial form and asymmetry in 
 # unaffected relatives of children with non-syndromic cleft lip and/or palate
 # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4025896/
 
 pmc4025896 <- c("PAX7", "ABCA4-ARHGAP292", "IRF6", "MSX1", "rs987525", "FOXE1", "TGFB3", "MAFB")
 
 # Three novel genes tied to mandibular prognathism in eastern Mediterranean families.
 # https://www.ncbi.nlm.nih.gov/pubmed/31256822
 pubmed_31256822 <- c("C1orf167", "NBPF8", "NBPF9")
 
 
 # Genetic Screening in Patients with Craniofacial Malformations
 # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5123894/#JR1616-1
 pmc5123894 <- c("FOXE1", "IRF6", "TGFA", "TBX22", "SATB2", "FGFR3", "FGFR2", "TCOF1", "POLR1C", "POLR1D", "SF3B4", "PRX1", "PRX2")
 
 literature_genes <- unique(c(gene_in_facebase, 
                              pmc4025896, 
                              pubmed_31256822, 
                              pmc5123894))
 literature_genes <- literature_genes[literature_genes %in% ens_to_hugo$HGNC.symbol]
 
 
 df <- data.frame(lookalike = all_genes_infinium %in% lookalike_genes,
                  facial   = all_genes_infinium %in% literature_genes)
 
 fisher.test(table(df))
 
 # Claes et al. 2018
 # TABLE 1
 strong_evidence <- c("TBX15", "ASPM", "PKDCC", 
                      "HOXD1", "HOXD3", "HOXD4", "HOXD8", "HOXD9", "HOXD10", "HOXD11", "HOXD12", "HOXD13", 
                      "PAX3", "RAB7A", "ACAD9", "EPHB3", "DVL3", "DCHS2", "SUPT3H", "RPS12", 
                      "EYA4", "DLX6", "DYNC1L1", "BC039327", "SOX9", "KCTD15")
 
 weak_evidence <- c("CDH18", "NHP2", "ZNF354A", "PAX1")
 
 
 strong_evidence <- strong_evidence[strong_evidence%in% ens_to_hugo$HGNC.symbol]
 weak_evidence   <- weak_evidence[weak_evidence%in% ens_to_hugo$HGNC.symbol]
 
 df <- data.frame(lookalike = all_genes_infinium %in% lookalike_genes,
                  facial   = all_genes_infinium %in% c(strong_evidence, weak_evidence))
 
 fisher.test(table(df))
 
 ######################
 #### GWAS CENTRAL ####
 ######################
 
 # Two files, one for the small datasets, one for the huge one
 gwas_central <- read.table("gwas_central_martquery_0712105315_100.txt.gz", 
                              sep = "\t", stringsAsFactors = F, header = T)
 
 gwas_central_ranges <- makeGRangesFromDataFrame(gwas_central, keep.extra.columns = T, seqnames.field = "Chromosome", start.field = "Marker.Start", end.field = "Marker.Stop")
 
 seqlevelsStyle(gwas_central_ranges) <- "UCSC"
 
   df <- data.frame(lookalike_snp = overlapsAny(infinium_ranges, lookalike_snps),
                    trait_snp     = overlapsAny(infinium_ranges, gwas_central_ranges))
 
   fisher.test(table(df))

   
   ## For each study separately
   central <- c()
   for(study in unique(gwas_central$GC.Study.Identifier)) {
     snps_in_study <- gwas_central_ranges[gwas_central_ranges$GC.Study.Identifier == study]
     
     df <- data.frame(lookalike_snp = overlapsAny(infinium_ranges, lookalike_snps),
                      trait_snp     = overlapsAny(infinium_ranges, snps_in_study))
     
     central <- rbind(central, c("Study ID" = study,
                       "Study Name" = unique(gwas_central_ranges$Study.Name[gwas_central_ranges$GC.Study.Identifier == study]),
                       "Total SNPs" = length(snps_in_study),
                       "OR" = fisher.test(table(df))$estimate,
                       "P-value" = fisher.test(table(df))$p.value))
   }
   
 # Print
  for(i in 1:nrow(central)) {
    if(i == 1)    cat(paste0(colnames(central), sep = "\t\t"), "\n")
    cat(paste0(central[i, ], sep = "\t\t"), "\n")
  }
   
   
   ## For each phenotype separately
   central <- c()
   for(phenotype in unique(gwas_central$Phenotype.Annotation.Identifier)) {
     snps_in_phenotype <- gwas_central_ranges[gwas_central_ranges$Phenotype.Annotation.Identifier  == phenotype]
    
     df <- data.frame(lookalike_snp = overlapsAny(infinium_ranges, lookalike_snps),
                      trait_snp     = overlapsAny(infinium_ranges, snps_in_phenotype))
    
     central <- rbind(central, c("phenotype ID" = phenotype,
                                 "phenotype Name" = unique(gwas_central_ranges$Annotation.Name[gwas_central_ranges$Phenotype.Annotation.Identifier == phenotype]),
                                 "Total SNPs" = length(snps_in_phenotype),
                                 "OR" = fisher.test(table(df))$estimate,
                                 "P-value" = fisher.test(table(df))$p.value))
   }
   
   
   # Print
   central <- cbind(central, "Adj.pval" = p.adjust(as.numeric(central[, "P-value"])))
   central <- central[order(central[, "P-value"]), ]
   central[, "P-value"] <- round(as.numeric(central[, "P-value"]), 6)
   
   central[, "OR.odds ratio"] <- round(as.numeric(central[, "OR.odds ratio"]), 4)
   cat(paste0(colnames(central), sep = "\t\t"), "\n")
   for(i in 1:nrow(central)) {
     if(i == 1) cat(paste0(colnames(central), sep = "\t\t"), "\n")
     cat(paste0(central[i, ], sep = "\t\t"), "\n")
   }
   
   
      
  #######################
  # Gene lists for further enrichment analysis ####

 
 lookalike_genes    <- unique(rownames(ens_to_hugo)[ens_to_hugo$HGNC.symbol %in% lookalike_genes])
 all_genes_infinium <- unique(rownames(ens_to_hugo)[ens_to_hugo$HGNC.symbol %in% all_genes_infinium])
 
 
 write.table(matrix(lookalike_genes, ncol = 1), col.names = F, row.names = F, quote = F, file = "lookalike_genes.txt")
 write.table(matrix(all_genes_infinium, ncol = 1), col.names = F, row.names = F, quote = F, file = "all_infinium_genes.txt")

 # Run DAVID enrichment with these lists (all_genes_infinium is background)
 # OMIM
 # UP_keywords
 # GOTERM_BP_direct
 # GOTERM_CC_direct
 # GOTERM_MF_direct
 # PUBMED_ID
 # KEGG_PATHWAY
 # REACTOME_PATHWAY
 # UP_TISSUE
 
 ## GOrilla enrichment
 
 

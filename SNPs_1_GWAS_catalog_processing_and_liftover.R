## Filtering and liftover of the GWAS catalog
setwd("data/")

load("merged_snp9_infinium_table.RData", verbose = T)

## rs snps
rs_snps <- unique(all_snp_information$Name[grepl("rs",all_snp_information$Name)])

# Some lines are duplicated except for the name where the snp is annotated (kgp or rs)
# I remove these lines
all_snp_information <- all_snp_information[!duplicated(all_snp_information[, colnames(all_snp_information) != "Name"]),]


gwas_catalog <- read.table("gwas_catalog_v1.0.2-associations_e95_r2019-03-01.tsv", sep = "\t",
                           header = T, stringsAsFactors = F, fill = T)


# Some chromosome and position are elswere in the matrix, I will create new chr and position columns and fix this
rownames(gwas_catalog) <- paste0("gwas_catalog_", 1:nrow(gwas_catalog))
gwas_catalog$chr <- paste0("chr", gwas_catalog$CHR_ID)
gwas_catalog$pos <- as.integer(gwas_catalog$CHR_POS)

# Fixable if they have position in column SNPS
fixable <- which(is.na(gwas_catalog$pos) & grepl("chr", gwas_catalog$SNPS, ignore.case = T) )
gwas_catalog$SNPS_corrected <- gsub(pattern = "_", ":", gwas_catalog$SNPS)
chr <- unlist(lapply(strsplit(gwas_catalog[fixable, "SNPS_corrected"], split = ":"), function(x) x[1]))
pos <- as.integer(unlist(lapply(strsplit(gwas_catalog[fixable, "SNPS_corrected"], split = ":"), function(x) x[2])))

gwas_catalog$chr[fixable] <- chr
gwas_catalog$pos[fixable] <- pos

### 
maybe_fixable <-  which(is.na(gwas_catalog$pos) & !grepl("chr", gwas_catalog$SNPS) & grepl(":", gwas_catalog$SNPS) & 
                          !grepl("\\*", gwas_catalog$SNPS) &  !grepl("European", gwas_catalog$SNPS))

chr <- paste0("chr", unlist(lapply(strsplit(gwas_catalog[maybe_fixable, "SNPS_corrected"], split = ":"), function(x) x[1])))
pos <- as.integer(unlist(lapply(strsplit(gwas_catalog[maybe_fixable, "SNPS_corrected"], split = ":"), function(x) x[2])))

gwas_catalog$chr[maybe_fixable] <- chr
gwas_catalog$pos[maybe_fixable] <- pos

# Still 2816 not fixed
not_fixed <- gwas_catalog[!gwas_catalog$chr %in% paste0("chr", c(1:22, "X", "Y")), ]

gwas_catalog_filtered <- gwas_catalog[gwas_catalog$chr %in% paste0("chr", c(1:22, "X", "Y")),  ]
# Also filter out NAs (from not reported)
gwas_catalog_filtered <- gwas_catalog_filtered[!is.na(gwas_catalog_filtered$pos), ]

gwas_ranges <- GRanges(seqnames = gwas_catalog_filtered$chr, 
                       ranges   = IRanges(start = gwas_catalog_filtered$pos,
                                          end   = gwas_catalog_filtered$pos,
                                          names = rownames(gwas_catalog_filtered)))

library(rtracklayer)
ch = import.chain("hg38ToHg19.over.chain")

gwas_ranges_37 = unlist(liftOver(gwas_ranges, ch))

gwas_catalog_filtered$chr_build37 <- as.character(NA)
gwas_catalog_filtered$pos_build37 <- as.integer(NA)

gwas_catalog_filtered[names(gwas_ranges_37), "chr_build37"] <- as.character(seqnames(gwas_ranges_37))
gwas_catalog_filtered[names(gwas_ranges_37), "pos_build37"] <- start(gwas_ranges_37)

# Remove the ones that were not lifted
gwas_catalog_filtered <- gwas_catalog_filtered[!is.na(gwas_catalog_filtered$chr_build37),]
gwas_catalog_filtered <- cbind("gwas_id_maria" = as.character(rownames(gwas_catalog_filtered)), gwas_catalog_filtered)
gwas_catalog_filtered$gwas_id_maria <- as.character(gwas_catalog_filtered$gwas_id_maria)

gwas_catalog_filtered[is.na(gwas_catalog_filtered$pos_build37), "REGION" ]


### Save
save(gwas_catalog_filtered, file = "gwas_catalog_filtered_and_lifted_to_GRCh37.RData")


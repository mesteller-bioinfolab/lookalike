setwd("data/")
library(gplots)
library(pvclust)


ultra_lals <- c("03", "04", "07", "08", "13", "15", "16", "20", "25", "29", "33", "45",  "46", "47", "49", "58")
ultra_lals <- paste0("LAL.00", ultra_lals)
total_ultra_lals <- 32


#### SNPs ####
genome_studio_output <- read.table("LAL_Omni5_CNV_16p.csv", sep = "\t", stringsAsFactors = F, header = T)
# last line is empty 
genome_studio_output <- genome_studio_output[1:(nrow(genome_studio_output)-1), ]

snps <- genome_studio_output[, grepl(pattern = "GType", colnames(genome_studio_output))]
colnames(snps) <- gsub(pattern = ".GType", replacement = "", colnames(snps))
snps <- t(snps)
snps[snps == "NC"] <- NA
snps[snps == "AA"] <- 0
snps[snps == "AB"] <- 1
snps[snps == "BB"] <- 2

snps_numeric <- c()
for(i in 1:nrow(snps)) {
  snps_numeric <- rbind(snps_numeric, as.integer(snps[i, ]))
  message(i)
}
rownames(snps_numeric) <- rownames(snps)
 
dist_snps <- dist(snps_numeric, method = "euclidean")
save(dist_snps, file = "dist_snps.RData")


#### CNVs ####
cnvs <- genome_studio_output[, grepl("CNV.Value", colnames(genome_studio_output))]
colnames(cnvs) <- gsub(pattern = ".CNV.Value", replacement = "", colnames(cnvs))
cnvs <- t(cnvs)

dist_cnvs <- dist(cnvs, method = "euclidean")
save(dist_cnvs, file = "dist_cnvs.RData")


#### Methylation 
bvalues <- read.csv("b_values.csv", sep = ",", stringsAsFactors = F, row.names = 1)
bvalues <- as.matrix(bvalues)

dist_bvals_all <- dist(t(bvalues), 
                      method = "euclidean")
save(dist_bvals_all, file = "dist_bvals_all.RData")

#### Microbiome ####
# Qualitative
qualit_meta <- read.table("table_diff_alpha_metag.csv", sep = ",", stringsAsFactors = F, header = T, row.names = 1)
qualit_meta <- as.matrix(t(qualit_meta))
dist_qualit <- dist(qualit_meta, method = "euclidean")
save(dist_qualit, file = "dist_qualit.RData")

# Quantitative
quant_meta  <- read.table("table_metag16p_proportions.csv", sep = ",", stringsAsFactors = F, header = T, row.names = 1)
quant_meta <- as.matrix(quant_meta)
dist_quanti <- dist(t(quant_meta), method = "euclidean")
save(dist_quanti, file = "dist_quanti.RData")

# Questionnaire (normalized to 1 numeric data)
# quest <- read.table("~/_projects/Look_alike_multiomics/SAMPLE_INFORMATION_LOOKALIKE/norm_to_1_curated_questionaire_all_numeric_with_Age_Height_Weight_ULTRALALs_only.csv",
# sep = ",", stringsAsFactors = F, header = T)
# 
# quest <- read.table("OUTPUT_lookalike_form_numeric.csv", 
#                     sep = ",", stringsAsFactors = F, header = T)

load("lookalike_form_normalized_to_1.RData")

norm1_df <- norm1_df[order(norm1_df$ID_PEBC), ]
rownames(norm1_df) <- gsub(norm1_df$ID_PEBC, pattern = "-", replacement = ".")
norm1_df$pair_id   <- gsub(norm1_df$pair_id, pattern = "-", replacement = ".")
norm1_df <- norm1_df[norm1_df$pair_id %in% ultra_lals, ]

norm1_df <- norm1_df[, 3:ncol(norm1_df)]
norm1_df <- as.matrix(norm1_df)

dist_quest <- dist(norm1_df, method = "euclidean")
save(dist_quest, file = "dist_quest.RData")





pdf("boxplots_similarity_6_data_types.pdf", width = 10)
par(mfrow = c(2,3), mar = c(4,4,3,1))
for(y in c("qualit", "quanti",  "snps", "cnvs", "bvals_all", "quest")) {
  
  load(paste0("dist_", y ,".RData"))
  
  distances <- get(paste0("dist_", y))
  
  plot_title <- c("qualit" = "Metagenome (qualitative)",
                  "quanti" = "Metagenome (quantitative)",
                  "snps" = "SNPs",
                  "cnvs" = "CNVs",
                  "bvals_all" = "Methylome",
                  "quest" = "Questionnaire")[y]
  
  
  valores_por_columna <- list()
  
  x <- 0
  for(col_num in 1:(total_ultra_lals -1)) {
    first_value <- x+1
    last_value  <- x+total_ultra_lals-col_num
    
    valores_por_columna[[col_num]]  <- c(rep(NA, (total_ultra_lals -1)-length(first_value:last_value)), 
                                         distances[first_value:last_value])
    
    x <- x+total_ultra_lals-col_num
  }
  
  dist_matrix <- do.call("cbind", valores_por_columna)
  colnames(dist_matrix) <- labels(distances)[1:(total_ultra_lals-1)]
  rownames(dist_matrix) <- labels(distances)[2:total_ultra_lals]
  
  
  ###  Are samples more similar to their than with the others?
  matrix_extra <- dist_matrix
  for(lal in ultra_lals) {
    lal_A <- paste0(lal, ".A")
    lal_B <- paste0(lal, ".B")
    if(lal_A %in% rownames(matrix_extra) & lal_A %in% colnames(matrix_extra)) matrix_extra[lal_A, lal_A] <- NA
    if(lal_B %in% rownames(matrix_extra) & lal_B %in% colnames(matrix_extra)) matrix_extra[lal_B, lal_B] <- NA
    if(lal_A %in% rownames(matrix_extra) & lal_B %in% colnames(matrix_extra)) matrix_extra[lal_A, lal_B] <- NA
    if(lal_B %in% rownames(matrix_extra) & lal_A %in% colnames(matrix_extra)) matrix_extra[lal_B, lal_A] <- NA
  }
  
  
  # intra and extra scores
  intra_pair <- c()
  for(lal in ultra_lals) {
    a <- NA
    b <- NA
    lal_A <- paste0(lal, ".A")
    lal_B <- paste0(lal, ".B")
    if(lal_B %in% rownames(matrix_extra) & lal_A %in% colnames(matrix_extra)) a <- dist_matrix[lal_B, lal_A]
    if(lal_A %in% rownames(matrix_extra) & lal_B %in% colnames(matrix_extra)) {
      b <- dist_matrix[lal_A, lal_B]
    } else { b <- NA }
    # print(c(a,b))
    intra_pair <- c(intra_pair, c(a,b)[!is.na(c(a,b))])
  }
  
  extra_pair <- c(matrix_extra[!is.na(matrix_extra)])
  
  # pdf("intra_extra_pair_distances_boxplot_prueba.pdf")
  lista <- list("Distance\nwithin pairs" = intra_pair, "Distance with\nother individuals" = extra_pair)
  ymax  <- max(boxplot(lista, plot = F)$stats)
  ymin  <- min(boxplot(lista, plot = F)$stats)
  
  
  bp <- boxplot(lista,
                main = plot_title,
                ylab = "Euclidean distance",
                xaxt = "n",
                outline = F,
                ylim = c(ymin-0.05*(ymax-ymin), ymax + 0.1*(ymax-ymin)), 
                lwd = 1.25, col = "gray70")
  axis(side = 1, at = c(1,2), labels = c("Distance\nwithin pairs\n(intra-pair)", "Distance with\nother individuals\n(extra-pair)"), 
       tick = F, line = 1.5)
  pval <- t.test(intra_pair, extra_pair)$p.value 
  
  if(pval < 0.05) {
    arrows(1, ymax + 0.05*(ymax-ymin),2, ymax + 0.05*(ymax-ymin), code = 3, angle = 90, length = 0.07, lwd = 1.5)
    label <- ifelse(pval < 0.05, yes = paste0("P = ", round(pval, 5)), no = "n.s.")
    text(x = 1.5, y = ymax + 0.1*(ymax-ymin), labels = label, cex = 1.1)
  }
  # dev.off()
  
  
}
dev.off()




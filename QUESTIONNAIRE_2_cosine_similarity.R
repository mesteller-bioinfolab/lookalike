setwd("data/")

ultra_lals <- c("03", "04", "07", "08", "13", "15", "16", "20", "25", "29", "33", "45",  "46", "47", "49", "58")
ultra_lals <- paste0("LAL-00", ultra_lals)


# Load numeric table
df <- read.table("OUTPUT_lookalike_form_numeric.csv",sep = ",", stringsAsFactors = F, header = T)

df <- df[df$pair_id %in% ultra_lals,]

# I remove physical traits and of place of origin 
categoric_physical_or_origin_traits <- c("Country",
                                         "Gender",
                                         "Language",
                                         "Race",
                                         "Current.Recidence.Category",
                                         "Hair.Color",
                                         "Hair.Shape",
                                         "Eyes.Color",
                                         "Country_born",
                                         "City_village_born",
                                         "Country_residence",
                                         "City_village_residence",
                                         "Lived_elsewhere")

numeric_physical_or_origin_traits   <- c("Age_turning_in_2016",
                                         "Weight",
                                         "Height")

# Remove categoric variables (all related to origin or physical traits)
df <- df[, !colnames(df) %in% categoric_physical_or_origin_traits]

# Remove columns with all 0
all_zeros_cols <- colnames(df)[3:ncol(df)][apply(df[,3:ncol(df)], 2, sum) == 0]
df             <- df[, !colnames(df) %in% all_zeros_cols]
num_m <- df[, 3:ncol(df)]
rownames(num_m) <- df$ID_PEBC


# Normalize columns to 1
for(col_x in 3:ncol(df)) {
  df[, col_x] <- df[, col_x]/max(df[, col_x], na.rm = T)
}


# order matrix
df <- df[order(df$ID_PEBC), ]

norm1_df <- df
rm(df)


cosine_sim <- function(a, b) crossprod(a,b)/sqrt(crossprod(a)*crossprod(b))

cosine_matrix <- matrix(as.numeric(NA), ncol = nrow(norm1_df), nrow = nrow(norm1_df))
for(x in 1:nrow(norm1_df)) {
  for(y in 1:nrow(norm1_df)) {
    a <- as.numeric(norm1_df[x, 3:ncol(norm1_df)])
    b <- as.numeric(norm1_df[y, 3:ncol(norm1_df)])
    
    cosine_matrix[x, y] <- cosine_sim(a,b)
  }
}
colnames(cosine_matrix) <- norm1_df$ID_PEBC
rownames(cosine_matrix) <- norm1_df$ID_PEBC
library(gplots)
heatmap.2(cosine_matrix,
          Rowv = F, 
          Colv = F,
          dendrogram = "none", 
          colsep = 1:ncol(cosine_matrix), rowsep = 1:nrow(cosine_matrix), 
          density.info = "none",
          cexCol = 0.8,
          cexRow = 0.8,
          col = colorRampPalette(c("white", "lightblue", "darkblue"))(80)[5:80],
          sepwidth = c(0.01, 0.01),
          margins = c(10,10),
          sepcolor = "white", trace = "none")

# Ranking of scores
ranking_scores <- c()
for(lal in unique(norm1_df$pair)) {
  lal_A <- paste0(lal, "-A")
  lal_B <- paste0(lal, "-B")
  ranking_scores <- c(ranking_scores, cosine_matrix[lal_A, lal_B])
  names(ranking_scores)[length(ranking_scores)] <- lal
}

ranking_scores <- sort(ranking_scores, decreasing = T)

ranking_scores <- matrix(round(ranking_scores, 4), ncol = 1, dimnames = list(names(ranking_scores), "cosine_score"))


ultraLAL <- paste0("LAL-00", c("07", "13", "15", "16", "25", "33", "45", "47", "58"))

distance_score   <- cbind(ranking_scores,   order = order(ranking_scores[, "cosine_score"], decreasing = T))

distance_score   <- cbind(distance_score,   ultraLAL = rownames(distance_score) %in% ultraLAL)
distance_score[, "ultraLAL"] <- ifelse(distance_score[, "ultraLAL"] == 1, yes = "*", no = "")
distance_score  <- cbind(rownames(distance_score), distance_score)

# write.table(distance_score, file = "cosine_ranking_with_age_height_weight_ULTRALALs_only.csv", quote = F, sep = ",", row.names = F)
# write.table(cosine_matrix , file = "cosine_matrix_with_age_height_weight_ULTRALALs_only.csv", quote = F, sep = ",", row.names = T)


# Filtering out origin height, weight and age

norm1_df <- norm1_df[, !colnames(norm1_df) %in% numeric_physical_or_origin_traits]

save(norm1_df, file = "lookalike_form_normalized_to_1.RData")


cosine_matrix <- matrix(as.numeric(NA), ncol = nrow(norm1_df), nrow = nrow(norm1_df))
for(x in 1:nrow(norm1_df)) {
  for(y in 1:nrow(norm1_df)) {
    a <- as.numeric(norm1_df[x, 3:ncol(norm1_df)])
    b <- as.numeric(norm1_df[y, 3:ncol(norm1_df)])
    
    cosine_matrix[x, y] <- cosine_sim(a,b)
  }
}
colnames(cosine_matrix) <- norm1_df$ID_PEBC
rownames(cosine_matrix) <- norm1_df$ID_PEBC





library(gplots)

# pdf("cosine_matrix_filtering_out_age_height_weight_ULTRALALs_only.pdf", height = 12, width = 12)
heatmap.2(cosine_matrix,
          Rowv = F, 
          Colv = F,
          dendrogram = "none", 
          colsep = 1:ncol(cosine_matrix), rowsep = 1:nrow(cosine_matrix), 
          density.info = "none",
          cexCol = 0.7,
          cexRow = 0.7,
          col = colorRampPalette(c("white", "lightblue", "darkblue"))(80)[5:80],
          sepwidth = c(0.01, 0.01),
          margins = c(10,10),
          sepcolor = "white", trace = "none")
# dev.off()

# Ranking of scores
ranking_scores <- c()
for(lal in unique(norm1_df$pair)) {
  lal_A <- paste0(lal, "-A")
  lal_B <- paste0(lal, "-B")
  ranking_scores <- c(ranking_scores, cosine_matrix[lal_A, lal_B])
  names(ranking_scores)[length(ranking_scores)] <- lal
}

ranking_scores <- sort(ranking_scores, decreasing = T)

ranking_scores <- matrix(round(ranking_scores, 4), ncol = 1, dimnames = list(names(ranking_scores), "cosine_score"))


ultraLAL <- paste0("LAL-00", c("07", "13", "15", "16", "25", "33", "45", "47", "58"))

distance_score   <- cbind(ranking_scores,   order = order(ranking_scores[, "cosine_score"], decreasing = T))

distance_score   <- cbind(distance_score,   ultraLAL = rownames(distance_score) %in% ultraLAL)
distance_score[, "ultraLAL"] <- ifelse(distance_score[, "ultraLAL"] == 1, yes = "*", no = "")
distance_score  <- cbind(rownames(distance_score), distance_score)

# write.table(distance_score, file = "cosine_ranking_no_age_height_weight_ULTRALALs_only.csv", quote = F, sep = ",", row.names = F)
# write.table(cosine_matrix, file = "cosine_matrix_no_age_height_weight_ULTRALALs_only.csv", quote = F, sep = ",", row.names = T)

###  Are samples more similar to their than with the others?

cosine_matrix_extra <- cosine_matrix
for(lal in unique(norm1_df$pair)) {
  lal_A <- paste0(lal, "-A")
  lal_B <- paste0(lal, "-B")
  cosine_matrix_extra[lal_A, lal_A] <- NA
  cosine_matrix_extra[lal_B, lal_B] <- NA
  cosine_matrix_extra[lal_A, lal_B] <- NA
  cosine_matrix_extra[lal_B, lal_A] <- NA
}

for(x in 3:nrow(cosine_matrix_extra)) {
  for(y in 1:x) {
    cosine_matrix_extra[x, y] <- NA
  }
}

# intra and extra scores
intra_pair <- c(ranking_scores)
extra_pair <- c(cosine_matrix_extra[!is.na(cosine_matrix_extra)])
# pdf("intra_extra_pair_cosine_similarity_boxplot.pdf", width = 5, height = 5)
bp <- boxplot(list(intra_pair, extra_pair),
              main = "Cosine similarity matrix",
              ylab = "Cosine distance",
              ylim = c(0.2,1), lwd = 1.25, col = "gray70")

arrows(1,0.95,2, 0.95, code = 3, angle = 90, length = 0.07, lwd = 1.5)
# label <- ifelse(t.test(intra_pair, extra_pair)$p.value < 0.05, yes = "*", no = "n.s.")
# text(x = 1.5, y = 0.98, labels = label, cex = 1.5)

pval <- t.test(intra_pair, extra_pair)$p.value
text(x = 1.5, y = 0.98, labels = paste0("P = ", round(pval, 5)), cex = 1)
# dev.off()




# Generate boxplots of distance for each variable idependently
setwd("data/")


# Load 
load("lookalike_form_normalized_to_1.RData")

# Append Height and Weight 
df <- read.table("INPUT_curated_questionaire.csv",sep = ",", stringsAsFactors = F, header = T)
# First normalize the variables to 1
for(col_x in c("Height", "Weight")) {
  df[, col_x] <- df[, col_x]/max(df[, col_x], na.rm = T)
}

# Put same rownames in both tables
rownames(df) <- df$ID_PEBC
rownames(norm1_df) <- norm1_df$ID_PEBC
# Append
norm1_df$Height <- df[rownames(norm1_df), "Height"]
norm1_df$Weight <- df[rownames(norm1_df), "Weight"]

# Intra-pair and Extra-pair combinations of individuals  
all_combinations <- t(combn(norm1_df$ID_PEBC, m = 2))
intra_combinations <- all_combinations[substr(all_combinations[, 1], start = 1,8) == substr(all_combinations[, 2], start = 1,8),]
extra_combinations <- all_combinations[substr(all_combinations[, 1], start = 1,8) != substr(all_combinations[, 2], start = 1,8),]

pdf("boxplots_distance_independent_variables.pdf", width = 10)
par(mfrow = c(2,3), mar = c(4,4,3,1))
layout(mat = rbind(c(1,2,3), c(5,6,4)))
for(i in 3:ncol(norm1_df)) {
  intra <- c()
  for(x in 1:nrow(intra_combinations)) {
    intra <- c(intra, abs(diff(norm1_df[intra_combinations[x, ],i])))
    
  }
  extra <- c()
  for(x in 1:nrow(extra_combinations)) {
    extra <- c(extra, abs(diff(norm1_df[extra_combinations[x, ],i])))
    
  }
  
  if(t.test(x = intra, y = extra)$p.value < 0.05) {
    message(colnames(norm1_df)[i])
    message(t.test(x = intra, y = extra)$p.value)
    
    lista <- list("Distance\nwithin pairs" = intra, "Distance with\nother individuals" = extra)
    ymax  <- max(boxplot(lista, plot = F)$stats)
    ymin  <- min(boxplot(lista, plot = F)$stats)
    
    bp <- boxplot(lista,
                  main = colnames(norm1_df)[i],
                  ylab = "Distance",
                  xaxt = "n",
                  outline = T,
                  ylim = c(ymin-0.05*(ymax-ymin), ymax + 0.1*(ymax-ymin)),
                  lwd = 1.25, col = "gray70")
    axis(side = 1, at = c(1,2), labels = c("Distance\nwithin pairs\n(intra-pair)", "Distance with\nother individuals\n(extra-pair)"), 
         tick = F, line = 1.5)
    pval <- wilcox.test(intra, extra)$p.value 

    # Add p.value
    if(pval < 0.05) {
      if(pval < 0.01) {
        pval_char <- formatC(pval, digits = 2, format = "E") 
      } else {
        pval_char <- round(pval, 2)
      }
      arrows(1, ymax + 0.05*(ymax-ymin),2, ymax + 0.05*(ymax-ymin), code = 3, angle = 90, length = 0.07, lwd = 1.5)
      text(x = 1.5, y = ymax + 0.1*(ymax-ymin), 
           labels = paste0("P = ",pval_char ), 
           cex = 1.1)
    }

  }
}
dev.off()

# Script to determine the breakpoints of the CNVs in each sample
library(GenomicRanges)

setwd("data/")

genome_studio_output <- read.table("LAL_Omni5_CNV_16p.csv", sep = "\t", stringsAsFactors = F, header = T)

# Sample names
samples <- substr(colnames(genome_studio_output), 1,10)
samples <- unique(samples[grepl("LAL", samples)])


# There is a chromosome 0, with no CN, I remove it
genome_studio_output <- genome_studio_output[genome_studio_output$Chr != "0",]

chromosomes <- c(1:22, "X", "XY", "MT")

#### Determine breakpoints ####

# Folder for output
if(!file.exists("tmp")) dir.create("tmp")

for(indiv in samples) {
  print(st <- Sys.time())
  # print(indiv)
  table_per_sample <- c()
  for(chr in chromosomes) {
    # print(chr)
    chr_table <- genome_studio_output[genome_studio_output$Chr == chr, c("Chr", "Position", paste0(indiv, ".CNV.Value"))]
    chr_table <- chr_table[order(chr_table$Position),] # Order all positions in a chromosome
    chr_table[is.na(chr_table[, paste0(indiv, ".CNV.Value")]), paste0(indiv, ".CNV.Value")] <- 999 
    
    # First position is always a breakpoint 
    breakpoints <- T
    
    # Central positions
    for(i in 2:(nrow(chr_table)-1)) {
      if( length(unique(chr_table[c(i-1, i, i+1), paste0(indiv, ".CNV.Value")])) == 1) {
        breakpoints <- c(breakpoints, F)
      } else {
        breakpoints <- c(breakpoints, T)
      }
    }
    # Final position is always a breakpoint 
    breakpoints <- c(breakpoints, T)
    
    chr_table <- chr_table[breakpoints, ]
    
    
    if(nrow(chr_table) == 2) {
      starts <- 1
      ends   <- 2
    } else {
      # Starts and ends
      starts <- 1
      ends  <- c()
      
      # If starting position is CNV (independent of i, then it is an end too)
      if(chr_table[1, paste0(indiv, ".CNV.Value")] != chr_table[2, paste0(indiv, ".CNV.Value")]) ends <- c(ends, 1)
      
      # Central positions
      for(i in 2:(nrow(chr_table)-1)) {
        if(chr_table[i, paste0(indiv, ".CNV.Value")] != chr_table[i-1, paste0(indiv, ".CNV.Value")]) starts <- c(starts, i)
        if(chr_table[i, paste0(indiv, ".CNV.Value")] != chr_table[i+1, paste0(indiv, ".CNV.Value")]) ends   <- c(ends, i)
      }
      
      # Final position
      # If last position is CNV (independent of prelast, then it is a start too)
      if(chr_table[nrow(chr_table), paste0(indiv, ".CNV.Value")] != chr_table[nrow(chr_table)-1, paste0(indiv, ".CNV.Value")]) starts <- c(starts, nrow(chr_table))
      # Last position is always an end
      ends <- c(ends, nrow(chr_table))
    }
    
    
    table_per_sample <- rbind(table_per_sample, 
                              cbind(chr   = chr_table[starts, "Chr"], 
                                    start = chr_table[starts, "Position"],
                                    end   = chr_table[ends, "Position"],
                                    CN    = Rle(chr_table[, paste0(indiv, ".CNV.Value")])@values))
    
    
  }
  
  table_per_sample        <- as.data.frame(table_per_sample, stringsAsFactors = F)
  table_per_sample$start  <- as.integer(table_per_sample$start)
  table_per_sample$end    <- as.integer(table_per_sample$end)
  table_per_sample$CN     <- as.integer(table_per_sample$CN)
  
  
  save(table_per_sample, file = paste0("tmp/", paste0(indiv, "breakpoint_table.RData" )))
  print(Sys.time() - st)
}




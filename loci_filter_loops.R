library(vcfR)

# read in
vcf <- read.vcfR("pruned_pinelake.vcf")

# what filtering thresholds you want to test
locus_thresh <- c(.1,.15,.2,.25,.3,.35)
sample_thresh <- c(.1,.15,.2,.25,.3,.35,.4,.45,.5)

# filtering function
filter_vcf <- function(vcf, locus_thresh, sample_thresh) {
  # Get all sample columns (keep FORMAT)
  gt_full <- vcf@gt
  
  # Extract just the GT part from each sample
  gt_matrix <- apply(gt_full[, -1, drop = FALSE], c(1,2), function(x) {
    if (is.na(x)) return(NA)
    strsplit(x, ":")[[1]][1]
  })
  gt_matrix[gt_matrix == "./." | gt_matrix == "." | gt_matrix == ""] <- NA
  
  # Calculate missingness
  loci_missing <- apply(gt_matrix, 1, function(x) mean(is.na(x)))
  samples_missing <- apply(gt_matrix, 2, function(x) mean(is.na(x)))
  
  # Subset VCF by loci
  vcf <- vcf[loci_missing <= locus_thresh, ]
  
  # Recompute GT matrix with updated VCF
  gt_full <- vcf@gt
  gt_matrix <- apply(gt_full[, -1, drop = FALSE], c(1,2), function(x) {
    if (is.na(x)) return(NA)
    strsplit(x, ":")[[1]][1]
  })
  gt_matrix[gt_matrix == "./." | gt_matrix == "." | gt_matrix == ""] <- NA
  
  # Subset VCF by samples (preserving FORMAT)
  keep_samples <- names(which(apply(gt_matrix, 2, function(x) mean(is.na(x))) <= sample_thresh))
  vcf <- vcf[, c(TRUE, colnames(vcf@gt)[-1] %in% keep_samples)]  # Keep FORMAT + desired samples
  
  return(vcf)
}

# store results
results <- list()

# filtering loops
for (i in locus_thresh) {
  for (j in sample_thresh) {
    name <- paste0("locus threshold", i, "and sample threshold", j)
    results[[name]] <- filter_vcf(vcf, i, j)
    cat("Finished", name, "\n")
  }
}

# look at the results object, this will just print each summary one after other, 
# so what i did was just copy and paste it into a text editor and then look through it like that
# (and annotate)
results

# take the filtered file you want
# either grab from the results object
filtered <- results[[54]]

# or if you didn't run the whole loop/just want to be sure
filtered <- filter_vcf(vcf, 0.35, 0.5)

# export
# note: exports as .gz
write.vcf(filtered, file = "filtered.vcf.gz")

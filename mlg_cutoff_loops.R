library(vcfR)
library(poppr)
library(tidyr)
library(tibble)

# file read in and set up to get from vcf -> SNPclone
vcf <- read.vcfR("path/to/file")
meta <- read.table("path/to/population/meta/file", sep = "\t", 
                   header = T)
gl <- vcfR2genlight(vcf)
sc <- as.snpclone(gl)
# if you have done any sample filtering that results in your SNPclone object 
# having less individuals than your original meta file
meta <- meta[meta$Sample %in% sc$ind.names, ]  
strata(sc) <- meta

# calculate distances however you like
# I am doing bitwise non-Euclidean
dist <- bitwise.dist(sc, euclidean = F)

# view the distributions of distances throughout the matrix
fil <- filter_stats(sc, distance = dist, plot=T)
cutoff_predictor(fil$farthest$THRESHOLDS)
cutoff_predictor(fil$nearest$THRESHOLDS)
cutoff_predictor(fil$average$THRESHOLDS)

# create a list of your duplicates' names
# this is how I assess the accurateness of clone correction. 
# I want all my DNA extraction duplicates to collapse together (or as many as possible)
# you may learn from this process some of your duplicates just refuse to collapse 
# and you may want/need to remove those from the dataset because of sequence messiness/error.
dupes <- list(c("sample1_dup", "sample2_dup", "sample3_dup", 
                "sample4_dup", "sample5_dup", "sample6_dup"), 
              c("sample1", "sample2", "sample3", 
                "sample4", "sample5", "sample6")
              )

# concatenated list of collapsing cutoffs you want to use, informed by lines 20-23
distance_cutoffs <- c(1, 2, 3, 4, 5, 6)
mlg_df <- list()
k <- 1

# test farthest neighbour algorithm, at each distance cutoff listed in line 34
for (i in distance_cutoffs) {
  mlg.filter(sc, algorithm = "farthest_neighbor",
             distance = dist) <- i
  # grab sample names and MLG assignments for each individual
  mlg_df[[k]] <- mlg.id(sc) %>%
    lapply(as.character) %>%
    enframe(name = "MLG", value = "Sample") %>%
    unnest(cols = "Sample")
  # reassign a new lcoation in the list for each loop
  k <- k + 1
}

# are duplicates getting collapsed or not
far_results <- replicate(
  k-1,
  logical(length(dupes[[1]])),
  simplify = F
  )
for (n in 1:(k-1)) {
  for (i in seq_along(dupes[[1]])) {
    d1 <- dupes[[1]][i]
    d2 <- dupes[[2]][i]
    
    mlg1 <- mlg_df[[n]]$MLG[mlg_df[[n]]$Sample == d1]
    mlg2 <- mlg_df[[n]]$MLG[mlg_df[[n]]$Sample == d2]
    
    far_results[[n]][i] <- mlg1 == mlg2 
  }
}

far_results


# FYI : all intermediate objects named in each algorithm section are named the same, 
# but the result output has different names

# concatenated list of collapsing cutoffs you want to use, informed by lines 20-23
distance_cutoffs <- c(1, 2, 3, 4, 5, 6)
mlg_df <- list()
k <- 1

# test average neighbour algorithm (UPGMA), at each distance cutoff listed in line 34
for (i in distance_cutoffs) {
  mlg.filter(sc, algorithm = "average_neighbor",
             distance = dist) <- i
  # grab sample names and MLG assignments for each individual
  mlg_df[[k]] <- mlg.id(sc) %>%
    lapply(as.character) %>%
    enframe(name = "MLG", value = "Sample") %>%
    unnest(cols = "Sample")
  # reassign a new lcoation in the list for each loop
  k <- k + 1
}

# are duplicates getting collapsed or nawr
ave_results <- replicate(
  k-1,
  logical(length(dupes[[1]])),
  simplify = F
)
for (n in 1:(k-1)) {
  for (i in seq_along(dupes[[1]])) {
    d1 <- dupes[[1]][i]
    d2 <- dupes[[2]][i]
    
    mlg1 <- mlg_df[[n]]$MLG[mlg_df[[n]]$Sample == d1]
    mlg2 <- mlg_df[[n]]$MLG[mlg_df[[n]]$Sample == d2]
    
    ave_results[[n]][i] <- mlg1 == mlg2 
  }
}

ave_results


# concatenated list of collapsing cutoffs you want to use, informed by lines 20-23
distance_cutoffs <- c(1, 2, 3, 4, 5, 6)
mlg_df <- list()
k <- 1

# test nearest neighbour algorithm, at each distance cutoff listed in line 34
for (i in distance_cutoffs) {
  mlg.filter(sc, algorithm = "nearest_neighbor",
             distance = dist) <- i
  # grab sample names and MLG assignments for each individual
  mlg_df[[k]] <- mlg.id(sc) %>%
    lapply(as.character) %>%
    enframe(name = "MLG", value = "Sample") %>%
    unnest(cols = "Sample")
  # reassign a new lcoation in the list for each loop
  k <- k + 1
}

# are duplicates getting collapsed or nawr
nar_results <- replicate(
  k-1,
  logical(length(dupes[[1]])),
  simplify = F
)
for (n in 1:(k-1)) {
  for (i in seq_along(dupes[[1]])) {
    d1 <- dupes[[1]][i]
    d2 <- dupes[[2]][i]
    
    mlg1 <- mlg_df[[n]]$MLG[mlg_df[[n]]$Sample == d1]
    mlg2 <- mlg_df[[n]]$MLG[mlg_df[[n]]$Sample == d2]
    
    nar_results[[n]][i] <- mlg1 == mlg2 
  }
}

nar_results

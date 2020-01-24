#!/usr/bin/env Rscript

library(ape)
tr <- read.tree(commandArgs(T)[1])

cat(Ntip(tr), "tips:", paste(tr$tip.label, collapse=','), '\n')

cat("is.rooted:", is.rooted(tr), "\n")
cat("is.binary:", is.binary(tr), "\n")

min_thr <- 0.0000001
max_thr <- 10
thr <- min_thr
is.ultr <- FALSE
#thresholds <- c()
results <- c()
while (thr <= max_thr && !is.ultr) {
  is.ultr <- is.ultrametric(tr, thr)
  results <- c(results, is.ultr)
  #thresholds <- c(threshold, thr)
  thr <- thr * 10
}

first_true <- which(results)
if ( !length(first_true) ) {
  cat(paste0("is.ultrametric(", max_thr, "): ", FALSE, "\n"))
} else {
  cat(paste0("is.ultrametric(", min_thr * 10^(first_true-2), "): ", results[first_true - 1], "\n"))
  cat(paste0("is.ultrametric(", min_thr * 10^(first_true-1), "): ", results[first_true], "\n"))
}

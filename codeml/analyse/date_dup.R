#!/usr/bin/env Rscript

library(ape)


load_calibration  <- function(agefile) {
    agefile <- read.delim(agefile)
}

date_line <- function(line){
	tree <- read.tree(text=line)
}


date_all <- function(treefile) {
	f <- file(treefile, "r")
	while ( TRUE ) {
		line <- readLines(f, n=1)
		if ( length(line) == 0 ) {
			break
		}
		date_line(line)
	}
	close(f)
}

#!/usr/bin/env Rscript

library(ape)

readopts <- function(...) {
  progname <- basename(sub('--file=', '',
                           grep('--file=', commandArgs(), value=T, fixed=T),
                           fixed=T)
                       )
  expectopts <- list(...)
  opts <- lapply(expectopts, `[[`, 'default')  # fill with default values
  # Apply default FALSE to all logical options
  opts[sapply(expectopts, `[[`, "type") == "logical"] <- FALSE
  help <- paste(names(expectopts),
                ifelse(as.vector(lapply(opts, is.null)),
                       '', paste0('[', as.vector(opts), ']')),
                lapply(expectopts, `[[`, 'help'),
                collapse='\n')
  given_args <- commandArgs(T)
  if( !length(given_args) ) {
    cat(progname, '\n\n', help, '\n', sep='')
    quit(status=2)
  }
  pos_args <- c()
  i <- 1
  while (i <= length(given_args)) {
    arg <- given_args[i]
    if( arg %in% names(expectopts) ) {
      type <- expectopts[[arg]]$type
      if( type == "logical" ) {
        opts[[arg]] <- TRUE
      } else {
        opts[[arg]] <- do.call(paste0("as.", type), list(given_args[i+1]))
        i <- i+1
      }
    } else {
      pos_args <- c(pos_args, arg)
    }
    i <- i+1
  }
  return(list(args=pos_args, opts=opts))
}

argopts <- readopts("-b"=list(type="logical",
                              help="Return 1 if tree not binary."),
                    "-r"=list(type="logical",
                              help="return 1 if tree not rooted"),
                    "-B"=list(type="logical",
                              help="1 if tree is binary"),
                    "-R"=list(type="logical",
                              help="1 if tree is rooted."),
                    "-u"=list(type="double", default=Inf,
                              help="1 if tree not ultrametric above given precision"),
                    "-q"=list(type="logical",
                              help="quiet")
                    )
tr <- read.tree(argopts$args[1])
opts <- argopts$opts



if( !opts[['-q']] )
  cat(Ntip(tr), "tips:", paste(tr$tip.label, collapse=','), '\n',
      "is.rooted:", is.rooted(tr), "\n",
      "is.binary:", is.binary(tr), "\n", sep='')

fail <- (opts[['-r']] && !is.rooted(tr)) || ( opts[['-R']] && is.rooted(tr))
fail <- fail || (opts[['-b']] && !is.binary(tr)) || (opts[['-B']] && is.binary(tr))

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
if ( !opts[['-q']] ) { 
  if ( !length(first_true) ) {
    cat(paste0("is.ultrametric(", max_thr, "): ", FALSE, "\n"))
    fail <- fail || (opts[['-u']] < Inf)
  } else {
    cat(paste0("is.ultrametric(", min_thr * 10^(first_true-2), "): ", results[first_true - 1], "\n"))
    cat(paste0("is.ultrametric(", min_thr * 10^(first_true-1), "): ", results[first_true], "\n"))
    fail <- fail || (opts[['-u']] < min_thr * 10^(first_true-1))
  }
} else {
  fail <- fail || !is.ultrametric(tr, opts[['-u']])
}

if(!interactive())
  quit(status=fail)

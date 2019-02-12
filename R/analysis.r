# analysis.r -- functions for empirical analysis

annotate_density <- function(annot, tiles, reduce=TRUE) {
  if (class(annot) != 'GRanges' || class(tiles) != 'GRanges')
    stop("both 'anot' and 'tiles' must be GRanges objects.")
  strand(annot) <- '*'
  if (reduce)
    annot <- GenomicRanges::reduce(annot)
  hits <- GenomicRanges::findOverlaps(annot, tiles)

  overlaps <- GenomicRanges::pintersect(annot[queryHits(hits)], 
                                        tiles[subjectHits(hits)])
  # summarize by width
  widths <- tibble(width=width(overlaps), tile=subjectHits(hits)) 

  out <- tibble(seqnames=as.factor(GenomicRanges::seqnames(tiles)), 
                start=GenomicRanges::start(tiles), 
                end=GenomicRanges::end(tiles), 
                tile=seq_along(tiles)) %>% 
                left_join(widths) %>% 
                mutate(width = ifelse(is.na(width), 0, width)) %>% 
                group_by(seqnames, start, end, tile) %>% 
                summarize(codingbases=sum(width))
   bind_cols(out, as_tibble(mcols(tiles)))
}

#' 
#'
#' 
#' 
weighedOverlapsMean <- function(x, tiles, col) {
  hits <- GenomicRanges::findOverlaps(tiles, x)
  overlaps <- GenomicRanges::pintersect(tiles[queryHits(hits)], 
                                        x[subjectHits(hits)])
     
}


#' Create a GenomicRanges object from a tibble
#'
#' @param x a tibble with columns (chrom | chr) and (start,end | pos). 
#'   Optionally, a strand column. All other columns will be included
#'   as metadata.
#'
#' @export
to_granges <- function(x) {
 if (class(x) == 'GRanges') return(x)
 if (!('chrom' %in% colnames(x))) {
    if ('chr' %in% colnames(x)) {
      # rename column name
      colnames(x)[which(colnames(x) == 'chr')] <- 'chrom'
    } else {
      stop("'x' must contain either a 'chrom' or 'chr' column")
    }
  }
  if (!('start' %in% colnames(x) && 'end' %in% colnames(x))) {
    if ('pos' %in% colnames(x)) {
      x$start <- x$pos
      x$end <- x$pos
      x <- x[,-which(colnames(x) == 'pos')]
    } else {
      msg <- "'x' must have either 'start' and 'end', or 'pos' in column names"
      stop(msg)
    }
  }
  out <- GenomicRanges::GRanges(seqnames=as.character(x$chrom),
                                        IRanges::IRanges(x$start, x$end))
  if ('strand' %in% colnames(out))
    strand(out) <- x$strand
 
  range_cols <- c('chrom', 'start', 'end', 'strand')
  metacols <- x[, !(colnames(x) %in% range_cols)]
  GenomicRanges::mcols(out) <- S4Vectors::DataFrame(metacols)
  out
}


expand_region <- function(x, amount) {
  strand(x) <- '*'
  suppressWarnings(GenomicRanges::start(x) <- GenomicRanges::start(x) - amount)
  suppressWarnings(GenomicRanges::end(x) <- GenomicRanges::end(x) + amount)
  trim(x)
}

#' Reshape long dataframe frequencies into a wide matrix
#'
#' @param x a dataframe with columns \code{key}, \code{freq}, \code{time}
#'
#' @export
reshape_freqs <- function(x) {
  out <- select(x, key, freq, time) %>% 
            distinct(key, time, .keep_all=TRUE) %>% spread(time, freq)
  mat <- as.matrix(out[, -1])
  rownames(mat) <- out$key
  mat
}

#' Sort all entries into low/high quantiles
#'
#' @param x numeric vector
#' @param alpha cutoff level
#' @param suffix a character vector suffix
#'
#' @export
quantile_bin <- function(x, alpha=0.05, suffix='') {
  qs <- quantile(x, c(alpha, 1-alpha))
  ifelse(x < qs[1], paste0('low', suffix),
         ifelse(x > qs[2], paste0('high', suffix), NA))
}




rbinom_prob <- Vectorize(function(size, prob) rbinom(n=1, size, prob), c('size', 'prob'))

#' Bootstrap from level of frequency
bootstrap_freqs <- function(x, fun, ..., B=100) {
  n <- nrow(x)
  vals <- replicate(B, { fun(x[sample(seq_len(n), replace=TRUE), ], ...) }, simplify=FALSE)
  tibble(boot_i = seq_len(B), boot_val = vals)
}

#' Bootstrap a sample
#'
#' @param x dataset
#'
bootstrap_temp_cov <- function(x, fun, B=100, as_df=FALSE, swap=TRUE, N=NULL) {
  n <- nrow(x)
  vals <- replicate(B, {
      temp_cov(x[sample(seq_len(n), replace=TRUE), ], as_df=as_df, swap=swap, N=N)
  }, simplify=FALSE)
  tibble(boot_i = seq_len(B), boot_val = vals)
}



get_rate <- Vectorize(function(chrom, pos, map_funs) {
  chrom <- as.character(chrom)
  if (!(chrom %in% names(map_funs)))
    stop(sprintf("Chromosome '%s' not found in linkage map predict functions.",
                 chrom))
  map_funs[[chrom]](pos)
}, c('chrom', 'pos'))


get_genlen <- Vectorize(function(chrom, start, end, map_funs) {
  chrom <- as.character(chrom)
  if (!(chrom %in% names(map_funs)))
    stop(sprintf("Chromosome '%s' not found in linkage map predict functions.",
                 chrom))
  map_funs[[chrom]](end) - map_funs[[chrom]](start)
}, c('chrom', 'start', 'end'))


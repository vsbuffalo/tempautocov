# syncfile.r -- load in syncfiles

SYNC_ALLES <- c('A', 'T', 'C', 'G', 'N', 'del')
TDF_COLS <- c('chrom', 'pos', 'ref', 'alt', 'sample', 
              'time', 'rep', 'treatment', 'ref_count', 'depth')

read_sync <- function(file) {
  out <-  read_tsv(file, 
                   col_types=cols(Chrom='c', Pos='i', Ref='c', .default='c'))
  colnames(out) <- tolower(colnames(out))
  out %>% gather(sample, counts, -(chrom:ref)) %>% 
          separate(counts, into=SYNC_ALLES, sep=':', convert=TRUE)
}

separate_sample <- function(.data, .remove=FALSE) {
  .data %>% extract(sample, into=c('time', 'rep', 'treatment'), 
                    'f(\\d+)_r(\\d+)_?(.*)?', 
                    convert=TRUE, remove=.remove) %>% 
  mutate(treatment=ifelse(treatment == 'h', 'hot', 'control'))
}


# we need to add 'chr' to chromosome names to work with genetic maps
add_chr <- function(x) paste0('chr', x)

filter_loci <- function(.data, biallelic=TRUE, dels=FALSE, 
                        add_alt=TRUE)  {
  ALLELES <- c('A', 'T', 'C', 'G')
  out <- .data
  if (biallelic) {
    out$biallelic <- apply(out[, ALLELES], 1, function(x) sum(x > 0) == 2)
    message(sprintf('removing %d/%d non-biallelic loci', sum(!out$biallelic), length(out$biallelic)))
    out <- out %>% filter(biallelic) %>% select(-biallelic)
  } 
  if (add_alt) {
    if (!biallelic)
      stop("alt alleles can only be added if non-biallelic loci are filtered out.")
    alleles <- t(apply(out[, ALLELES], 1, function(x) ALLELES[which(x > 0)]))
    alt <- Map(function(a, r) setdiff(a, r),
                      split(alleles, 1:nrow(alleles)), out$ref)
    # remove variants that don't agree with the ref allele
    keep <- sapply(alt, length) == 1
    if (any(!keep)) {
      out <- out[keep, ]
      warning(sprintf('removing %d/%d variants that disagree with ref allele',
                      sum(!keep), length(keep)))
    }
    out$alt <- unlist(alt[keep])
  }
  if (!dels) {
    out <- out %>% filter(del == 0)
  }
  out
}

syncfile_calc_freqs <- function(.data) {
  ALLELES <- c('A', 'T', 'C', 'G')
  x <- as.matrix(.data[, ALLELES])
  .data$ref_count <- as.integer(x[cbind(1:nrow(x), match(.data$ref, ALLELES))])
  .data %>% mutate(tot=(A + T + C + G), freq=ref_count / tot)
}


#' Read a VCF File
#' 
#' @param file Path to input VCF file.
read_vcffile <- function(file) {
  vcf_cols <- cols(.default=col_character())
  res <- read_tsv(file, col_types=vcf_cols, comment="##")
  colnames(res) <- tolower(sub('^#', '', colnames(res)))
  mutate(res, pos=as.integer(pos))
}



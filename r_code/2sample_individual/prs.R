## Set of functions for calculating the prs score

get_ivs <- function(gwas_fname, p_threshold, clump_threshold){

  # Read in the gwas data
  gwas = read.csv(gwas_fname, stringsAsFactors=FALSE)

  # choose ivs greater than 5*10^-8
  iv = gwas[gwas$pval<p_threshold,]

  # Clump IVs - remove correlated IVs
  iv.clump = ieugwasr::ld_clump(iv, clump_r2 = clump_threshold)

  return(iv.clump)
}

preprocessing <- function(geno, snps, minor_allele) {

  # re-order the genotype data
  geno <- geno[, snps]

  # convert to a numeric matrix while preserving the eids
  eids <- rownames(geno)
  geno <- as.matrix(sapply(geno, as.numeric))
  rownames(geno) <- eids

  # convert all 0s to NAs
  geno[geno == 0] <- NA

  # convert to minor homozygote to 2
  # the heterozygote to 1
  # and the major homozygote to 0
  recode_alleles <- function(x) {abs(x - 3)}

  # apply to the geno dataset
  geno <- apply(geno, MARGIN = 2, FUN = recode_alleles)

  return(geno)
}


sum_prs <- function(geno, log_or) {
  
  # multiply each SNP by the log odds ratio
  rs <- t(t(geno) * log_or)
  
  # sum
  prs <- rowSums(rs, na.rm = TRUE)
  
  return(prs)
}

mean_prs <- function(geno, log_or) {

  # multiply each SNP by the log odds ratio
  rs <- t(t(geno) * log_or)

  # mean
  prs <- rowMeans(rs, na.rm = TRUE)

  return(prs)
}

get_prs <- function(genotype_fname, gwas_fname, p_threshold, clump_threshold) {

  # read in the data
  geno <- readRDS(genotype_fname)
  ivs <- get_ivs(gwas_fname, p_threshold, clump_threshold)

  # preprocess
  geno <- preprocessing(geno, ivs$rsid)

  # get the prs score
  prs <- sum_prs(geno, log(ivs$OR))

  return(prs)
}

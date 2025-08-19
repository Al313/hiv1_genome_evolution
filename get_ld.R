#=============================
# Load packages
#=============================
library(Rsamtools)
library(dplyr)
library(GenomicRanges)
library(parallel)
library(pbapply)
pboptions(type="txt")

if (file.exists("/home/amovas/")){
  print("Remote HPC Connection!")
  wd <- "/home/amovas/data/genome-evo-proj/"
} else {
  print("Local PC Connection!")
  wd <- "/Users/alimos313/Documents/studies/phd/hpc-research/genome-evo-proj/"
}


#=============================
# Inputs
#=============================
bam_file <- "/home/amovas/shared/genome-evo-proj/data/processed-data/mappings/pipeline-outputs/iii/13/13MT2EXPIIIVP100seq07062018_S4_L001_sorted.bam"
# "/Users/alimos313/Desktop/scrap/bam-portal/13/13MT2EXPIIIVP100seq07062018_S4_L001_sorted.bam"
# load variant data
source(paste0(wd, "src/analysis/1st-man/readin_data.R"))

variants <- variants_expiii %>%
  filter(exp_line == "MT-2_1" & passage == 100 & allele_freq >= 0.01 & genomic_pos <= 1000) %>%
  select(genomic_pos, ref_allele, alt_allele)

max_dist <- 400   # maximum distance between SNPs to consider

#=============================
# Safe null-coalesce
#=============================
`%||%` <- function(a, b) if (!is.null(a)) a else b

#=============================
# Helper: extract haplotypes from a read pair
#=============================
get_haplotype <- function(read_start, read_seq, mate_start, mate_seq, pos1, pos2) {
  cov_read <- read_start:(read_start + nchar(read_seq) - 1)
  cov_mate <- mate_start:(mate_start + nchar(mate_seq) - 1)

  seq_read <- setNames(strsplit(read_seq, "")[[1]], cov_read)
  seq_mate <- setNames(strsplit(mate_seq, "")[[1]], cov_mate)

  seq_pair <- c(seq_read, seq_mate)

  if (all(c(pos1, pos2) %in% names(seq_pair))) {
    allele1 <- seq_pair[as.character(pos1)]
    allele2 <- seq_pair[as.character(pos2)]
    return(paste0(allele1, allele2))
  } else {
    return(NA)
  }
}

#=============================
# Helper: safe count from table
#=============================
get_count <- function(counts, key) {
  if (!(key %in% names(counts))) return(0L)
  val <- counts[[key]]
  if (is.null(val) || length(val) == 0 || is.na(val)) return(0L)
  return(as.integer(val))
}

#=============================
# Function: compute LD (r^2) between two variants
#=============================
compute_ld <- function(var1, var2, bam_file) {
  region <- GRanges("NL43_ann_wk0virusPassRef_plasmid",
                    IRanges(start=min(var1$genomic_pos, var2$genomic_pos) - 10,
                            end=max(var1$genomic_pos, var2$genomic_pos) + 10))

  param <- ScanBamParam(which=region, what=c("qname","pos","seq","flag"))
  aln <- scanBam(bam_file, param=param)[[1]]

  if (length(aln$qname) == 0) return(NA)

  pairs <- split(seq_along(aln$qname), aln$qname)

  haplotypes <- lapply(pairs, function(ix) {
    if (length(ix) < 2) return(NA)

    read_start <- aln$pos[ix[1]]
    read_seq   <- as.character(aln$seq[ix[1]])
    mate_start <- aln$pos[ix[2]]
    mate_seq   <- as.character(aln$seq[ix[2]])

    if (is.na(read_start) || is.na(mate_start) ||
        nchar(read_seq) == 0 || nchar(mate_seq) == 0) return(NA)

    get_haplotype(read_start, read_seq, mate_start, mate_seq,
                  var1$genomic_pos, var2$genomic_pos)
  })

  haplotypes <- unlist(haplotypes)
  haplotypes <- haplotypes[!is.na(haplotypes)]
  if (length(haplotypes) == 0) return(NA)

  counts <- table(haplotypes)

  nAB <- get_count(counts, paste0(var1$alt_allele, var2$alt_allele))
  nAb <- get_count(counts, paste0(var1$alt_allele, var2$ref_allele))
  naB <- get_count(counts, paste0(var1$ref_allele, var2$alt_allele))
  nab <- get_count(counts, paste0(var1$ref_allele, var2$ref_allele))

  N <- nAB + nAb + naB + nab
  if (N == 0) return(NA)

  fAB <- nAB/N; fAb <- nAb/N; faB <- naB/N; fab <- nab/N
  pA <- fAB + fAb
  pB <- fAB + faB
  D  <- fAB - pA * pB
  r2 <- (D^2) / (pA*(1-pA)*pB*(1-pB))
  return(r2)
}

#=============================
# Prepare variant pairs within max_dist
#=============================
variant_pairs <- combn(nrow(variants), 2, simplify = FALSE)
variant_pairs <- Filter(function(pair) {
  abs(variants$genomic_pos[pair[1]] - variants$genomic_pos[pair[2]]) <= max_dist
}, variant_pairs)

#=============================
# Parallel computation with progress
#=============================
n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")) - 1  # detectCores() - 1
print(paste0("Cores available are ", n_cores))
cl <- makeCluster(n_cores)
clusterExport(cl, c("variants", "bam_file", "compute_ld", "get_haplotype", "get_count"))
clusterEvalQ(cl, {library(Rsamtools); library(GenomicRanges)})

print("starting ...")

# Use pbapply::pblapply to wrap parLapply with progress
ld_results <- pbapply::pblapply(variant_pairs, cl=cl, FUN=function(pair) {
  var1 <- variants[pair[1], ]
  var2 <- variants[pair[2], ]
  r2 <- compute_ld(var1, var2, bam_file)
  data.frame(pos1 = var1$genomic_pos,
             pos2 = var2$genomic_pos,
             r2   = r2)
})

stopCluster(cl)

# Combine results
ld_df <- do.call(rbind, ld_results)


#=============================
# Save output
#=============================

write.table(ld_df, paste0(wd, "results/tables/ld/100_1000.tsv"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)



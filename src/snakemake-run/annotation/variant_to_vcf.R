

# load the libraries

library(Biostrings)
library(stringr)

# read in fasta file to retrieve ref information for adjacent positions
fastaFile <- readDNAStringSet("/home/amovas/data/genome-evo-proj/data/reference/plasmid/hiv_plasmid_ref_genome.fasta")
seq_name = names(fastaFile)

sequence = paste(fastaFile)


# read in the variant file
variants <- read.table(file = "/home/amovas/data/genome-evo-proj/results/tables/2-p/all_variants.csv.gz", sep = ",", stringsAsFactors = F, header = T)

# remove the large deletion variations
variants<-variants[variants$mut_type!= "LD",]
# remove variants with AF below 0.01
variants <- variants[variants$fraction > 0.01,]

# 1st field of VCF in CHROM
variants$CHROM <- "AF324493.2"

# 2nd filed of VCF is POS and 3rd and 4th are REF and ALT fileds (remember that they are all in 1-based index coordinate system)
## for single base substitutions the conversion is pretty simple and straightforward.
variants$POS[variants$mut_type=="M"] <- variants$end[variants$mut_type=="M"] #end and start positions for SBS are the same and can be used interchangeably
variants$REF[variants$mut_type=="M"] <- variants$ref[variants$mut_type=="M"]
variants$ALT[variants$mut_type=="M"] <- variants$alt[variants$mut_type=="M"]

## for small insertions the conversion is done as follow:
### first the position of the insertion is the 1-based position of the base before the actual inseriton
variants$POS[variants$mut_type=="I"] <- variants$start[variants$mut_type=="I"]
### second the REF field in a VCF has to have the sequence of one base prior to the insetions (as opposed to being empty)
variants$REF[variants$mut_type=="I"] <- stringr::str_sub(sequence, start = variants$start[variants$mut_type=="I"], end = variants$start[variants$mut_type=="I"])
### lastly the ALT field of a VCF contains the one base prior to insetion and full sequence of insertion
variants$ALT[variants$mut_type=="I"] <- paste0(variants$REF[variants$mut_type=="I"],variants$alt[variants$mut_type=="I"])

## for small deletions the conversion is done as follows:
### first the position of deletion is the 1-based index of the base right before the deletion occurs
variants$POS[variants$mut_type=="D"] <- variants$start[variants$mut_type=="D"]
### second the REF field in a VCF has to have the sequence of one base prior to deletion and the full sequence of the bases that are deleted 
variants$REF[variants$mut_type=="D"] <- paste0(stringr::str_sub(sequence, start = variants$start[variants$mut_type=="D"], end = variants$start[variants$mut_type=="D"]),variants$ref[variants$mut_type=="D"])
### lastly the ALT field of a VCF for small deletions only contains the sequence of the one base prior to the deletion
variants$ALT[variants$mut_type=="D"] <- stringr::str_sub(sequence, start = variants$start[variants$mut_type=="D"], end = variants$start[variants$mut_type=="D"])

# sort the dataframe based on the position of variants and reassign the row numbers to be used as ids
variants <- variants[order(variants$POS),]
row.names(variants) <- 1:nrow(variants)

# 5th field of a VCF is the variant ID. In case the ID is not known it can be assigned to "."
variants$ID <- row.names(variants)
# 6th field of a VCF is the variant call quality. In this case since it is not known and irrelevant I just set to an arbitrary number.
variants$QUAL <- "."
# 7th field of a VCF is the filtering citeria. It's "PASSED" in our case
variants$FILTER <- "."
# 8th field of a VCF is the INFO field. It can contain any information that are not stored by other fields. It has to be in the correct format. Later on when I do the annotation
# with snpEff the annotation will be added to this field.
#variants$INFO <- "."
variants$INFO <- paste0("AF=",variants$fraction, ";MT=", variants$mut_type, ";PASSAGE=", variants$passage, ";LINE=", variants$line)

# Last and 9th field of a VCF is the FORMAT which can contain some technical information. In our case it is empty.
variants$FORMAT <- "."

# remove the redundant fields to get a df in VCF format
variants <- variants[,-c(1:10)]

vcf_variants <- variants
vcf_variants <- vcf_variants[,c(1,2,5,3,4,6,7,8,9)]

# Let's save the dataframe in VCF format. I'd try to   do so using vcfR package.

vcf_header <- '##fileformat=VCFv4.3
##contig=<ID=AF324493.2,length=14825>
##INFO=<ID=AF,Number=1,Type=Float,Description="Allele Freq.">
##INFO=<ID=MT,Number=1,Type=String,Description="Allele Freq.">
##INFO=<ID=PASSAGE,Number=1,Type=Integer,Description="Allele Freq.">
##INFO=<ID=LINE,Number=1,Type=Integer,Description="Allele Freq.">'


colnames(vcf_variants)[1] <- paste0("#",colnames(vcf_variants)[1])

#vcf_variants <- vcf_variants[1:25,]
#vcf_variants <- vcf_variants[1:11660,]
#209850
#sed '11673q;d' all.variants.ann.vcf


cat(vcf_header, file = '/home/amovas/data/genome-evo-proj/results/tables/2-p/all_variants.vcf', sep = '\n')

write.table(vcf_variants, file = '/home/amovas/data/genome-evo-proj/results/tables/2-p/all_variants.vcf', sep = "\t", row.names = F, quote = F, append = T)  


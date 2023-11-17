For beulding the snpEff database I followed the following steps:

First and foremost install the snpEff program.
Then download the genebank file of NL4-3 genome from NCBI.
Find the location of snpEff program installation running following command:
$ find ~ -name snpEff.config
Make a copy of the config file at the local directory.
Edit the config file as indicated in the documentation.

add the following to the config file after # databses and genomes section
 
# my added genes
hiv_plasmid_ref_genome.genome : HIV_1
        hiv_plasmid_ref_genome.chromosomes : AF324493.2


Run the below command to build the database:
$ snpEff build -c snpEff.config -genbank -v hiv_plasmid_ref_genome > snpEff.stderr 2> snpEff.stdout
Run the below command to take a peek into the database:
$ snpEff dump hiv_plasmid_ref_genome | less

To do the annotation after the reference database is prepared:
$snpEff -v hiv_plasmid_ref_genome ../../../results/tables/3-p/all_variants.vcf > ../../../results/tables/3-p/all.variants.ann.vcf
$java -Xmx8g -jar /Users/alimos313/snpEff/snpEff.jar -v hiv_plasmid_ref_genome ../../../results/tables/3-p/all_variants.vcf > ../../../results/tables/3-p/all_variants.ann.vcf

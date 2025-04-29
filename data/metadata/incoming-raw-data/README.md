Instructions:


Check for new data on the shared folder. Transfer them to local pc through SwitchDrive along with an NGS sample list excel version that contains those newly sequenced samples.

Create a new intermediary directory in incoming-raw-data (name corresponds to date of transfer). Carefully transfer them to their relevant folder in freezed data and run the snakemake pipeline (only until the production of bam/bam.bai files).

Check the coverage of repeated samples and decide whether to combine them or not.

If the need to be combined combine them in a separate directory (zcat a.R1.fastq.gz b.R1.fastq.gz | gzip > acombined.R1.fastq.gz ... then rename to new sequencing sample name including "combined before seq in the name). Do this for both reads of samples. Then transfer these samples to their respective directories.


For repeated samples whether combined or replaced make sure to move the bam/bam.bai of previous versions to archived.


Last but not least update the metadata following the script in the src directory.



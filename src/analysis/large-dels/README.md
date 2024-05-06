I tried to use different SV callers (GRIDSS, manta) but they did not call any SVs in our dataset.

I used LUMPY and it did output some SV calls, but not all of them. So, for example for sample 13_540 it correctly called the large deletions 8789-9039 which was supported by a lot of split read pieces, but the other large deletions 9102-9246 was not not called.

The following commands were used to invoke lumpy on those samples:


java -jar build/libs/picard.jar AddOrReplaceReadGroups I=/home/amovas/shared/genome-evo-proj/data/processed-data/mappings/pipeline-outputs/iii/13/13MT2EXPIIIVP550seq13072023_S4_L001_sorted.bam O=/home/amovas/data/scrap/out_with_rg.bam SORT_ORDER=coordinate RGID=13 RGLB=bar RGPL=illumina RGPU=10 RGSM=1 CREATE_INDEX=True

samtools view -b -F 1294 /home/amovas/data/scrap/out_with_rg.bam > ./13MT2EXPIIIVP550seq13072023_S4_L001.discordants.unsorted.bam

samtools view -h /home/amovas/data/scrap/out_with_rg.bam     | ../scripts/extractSplitReads_BwaMem -i stdin     | samtools view -Sb -     > ./13MT2EXPIIIVP550seq13072023_S4_L001.splitters.unsorted.bam

samtools sort 13MT2EXPIIIVP550seq13072023_S4_L001.discordants.unsorted.bam > 13MT2EXPIIIVP550seq13072023_S4_L001.discordants.bam
samtools sort 13MT2EXPIIIVP550seq13072023_S4_L001.splitters.unsorted.bam > 13MT2EXPIIIVP550seq13072023_S4_L001.splitters.bam 

./bin/lumpyexpress -B /home/amovas/data/scrap/out_with_rg.bam -S ./own_data/13MT2EXPIIIVP550seq13072023_S4_L001.splitters.bam -D ./own_data/13MT2EXPIIIVP550seq13072023_S4_L001.discordants.bam -o sample.vcf


Also when I tried to use the svtyper tool which is in the downstream of the lumpy pipeline, I could not make it work:

the main reason was that the source code of this package is in python 2 format and the package has long been left unsupported by its developers.


svtyper -i sample.vcf -B /home/amovas/data/scrap/out_with_rg.bam -S ./own_data/13MT2EXPIIIVP550seq13072023_S4_L001.splitters.bam sample.gt.vcf


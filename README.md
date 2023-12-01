# TRYOCC metagenomics
Repository containing bioinformatic pipelines and supplementary material for the project **"Facultatively intra-bacterial localization of a planthopper symbiont as an apparent adaptation to its vertical transmission."**

_Here we demonstrate an unusual transmission strategy adapted by one of the endosymbionts of the planthopper Trypetimorpha occidentalis (Hemiptera: Tropiduchidae) from Bulgaria. In this species, an Acetobacteraceae endosymbiont is transmitted transovarially within deep invaginations of cellular membranes of an ancient endosymbiont Sulcia - strikingly resembling recently described way of plant virus transmission. However, Acetobacteraceae co-colonizes the same bacteriocytes as Sulcia but does not enter its cells in males.  Then, the unusual endobacterial localization of Acetobacteraceae observed in females appears to be a unique adaptation to maternal transmission._

## Software requirement
* [Trim_galore](https://github.com/FelixKrueger/TrimGalore)
* [NCBI blast v2.11.0+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
* [Anvi'o] (https://anvio.org/)
* [Megahit](https://github.com/voutcn/megahit)

##Removing adapters using trim_galore
trim_galore --gzip --paired TRYOCC_R1.fq.gz TRYOCC_R2.fq.gz -o trimmed_seqs/ -j$NUM_THREADS --trim1 --phred33

megahit -1 TRYOCC-QUALITY_PASSED_R1.fastq -2 TRYOCC-QUALITY_PASSED_R2.fastq --min-contig-len 1000 -o ../TRYOCC_assembly/ --out-prefix TRYOCC -t $NUM_THREADS --k-max 255
anvi-script-reformat-fasta TRYOCC.contigs.fa -o TRYOCC_contigs.fasta --simplify-names --prefix TRYOCC --report-file TRYOCC_contigs_info.txt
bowtie2-build TRYOCC_contigs.fasta contigs
bowtie2 --threads $NUM_THREADS -x contigs -1 ../Quality_control/TRYOCC-QUALITY_PASSED_R1.fastq -2 ../Quality_control/TRYOCC-QUALITY_PASSED_R2.fastq --no-unal -S TRYOCC.sam
samtools view -F 4 -bS -@ $NUM_THREADS TRYOCC.sam > TRYOCC-RAW.bam
anvi-init-bam TRYOCC-RAW.bam -o TRYOCC.bam
anvi-gen-contigs-database -f TRYOCC_contigs.fasta -o TRYOCC.db -n 'Trypetimorpha occidentalis'
anvi-run-hmms -c TRYOCC.db -T $NUM_THREADS
anvi-run-ncbi-cogs -c TRYOCC.db --num-threads $NUM_THREADS
anvi-get-sequences-for-gene-calls -c TRYOCC.db -o TRYOCC-gene-calls.fa

anvi-profile -i TRYOCC.bam -c TRYOCC.db --output-dir TRYOCC_profile --sample-name TRYOCC -T $NUM_THREADS

pseudofinder.py annotate --genome prokka_SD7/SD7.gbk --outprefix SD7_pseudofinder --database ~/db/nr --threads $NUM_THREADS --skip_makedb


## Contigs taxonomic classification against whole NCBI nucleotide database

blastn -query TRYOCC_contigs.fasta -db nt -outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle' -evalue 1e-10 -task blastn -num_threads $NUM_THREADS -max_target_seqs 1 -max_hsps 1 > blastn_TRYOCC_contigs.txt

## Getting Avg. fold, GC percentage and contig size from bam file
pileup.sh in=TRYOCC.bam out=cov_TRYOCC.txt

# Merging the taxonomy and coverage outputs with [Nanotax](https://github.com/diecasfranco/Nanotax/blob/main/NanoTax_v2.2.py)

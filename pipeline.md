## Software requirement
* [Trim_galore](https://github.com/FelixKrueger/TrimGalore)
* [Anvi'o](https://anvio.org/)
* [Megahit](https://github.com/voutcn/megahit)
* [Bowtie2](https://github.com/BenLangmead/bowtie2)
* [Samtools](https://www.htslib.org/)
* [hmmer (version 3.2.1)](http://www.hmmer.org/)
* [NCBI blast v2.11.0+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
* BBMap [pileup.sh](https://github.com/BioInfoTools/BBMap/blob/master/sh/pileup.sh)
* [Bakta](https://github.com/oschwengers/bakta)
* [Pseudofinder](https://github.com/filip-husnik/pseudofinder)

### Removing adapters using trim_galore
```shell
trim_galore --gzip --paired TRYOCC_R1.fq.gz TRYOCC_R2.fq.gz -o trimmed_seqs/ -j $NUM_THREADS --trim1 --phred33
```
### Assembly with Megahit using max_kmer_size=255
```shell
megahit -1 TRYOCC-QUALITY_PASSED_R1.fastq -2 TRYOCC-QUALITY_PASSED_R2.fastq --min-contig-len 1000 -o ../TRYOCC_assembly/ --out-prefix TRYOCC -t $NUM_THREADS --k-max 255
```
### Reformating contig headers
```shell
anvi-script-reformat-fasta TRYOCC.contigs.fa -o TRYOCC_contigs.fasta --simplify-names --prefix TRYOCC --report-file TRYOCC_contigs_info.txt
```
### Creating references for alignment
```shell
bowtie2-build TRYOCC_contigs.fasta contigs
```
### Aligning raw reads against contigs
```shell
bowtie2 --threads $NUM_THREADS -x contigs -1 ../Quality_control/TRYOCC-QUALITY_PASSED_R1.fastq -2 ../Quality_control/TRYOCC-QUALITY_PASSED_R2.fastq --no-unal -S TRYOCC.sam
```
### Sorting sam file and converting to bam
```shell
samtools view -F 4 -bS -@ $NUM_THREADS TRYOCC.sam > TRYOCC-RAW.bam
anvi-init-bam TRYOCC-RAW.bam -o TRYOCC.bam
```
### Creating Anvi'o contigs database
```shell
anvi-gen-contigs-database -f TRYOCC_contigs.fasta -o TRYOCC.db -n 'Trypetimorpha occidentalis'
```
### Running hmm profiles
```shell
anvi-run-hmms -c TRYOCC.db -T $NUM_THREADS
```

## Contigs taxonomic classification against whole NCBI nucleotide database
```shell
blastn -query TRYOCC_contigs.fasta -db nt -outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle' -evalue 1e-10 -task blastn -num_threads $NUM_THREADS -max_target_seqs 1 -max_hsps 1 > blastn_TRYOCC_contigs.txt
```
## Getting Avg. fold, GC percentage and contig size information from bam file
```shell
pileup.sh in=TRYOCC.bam out=cov_TRYOCC.txt
```

### Merging the taxonomy and coverage outputs with [Nanotax](https://github.com/diecasfranco/Nanotax/blob/main/NanoTax_v2.2.py)



### Running bakta (in this case we used bakta on A. bogorensis W19, Ca. Kirkpatrickella and Acetobacteraceae_TRYOCC)
```shell
for file in *fasta; do SampleName=`basename $file .fasta`; bakta "$SampleName".fasta --complete --translation-table 4 --keep-contig-headers --compliant --threads $NUM_THREADS --output a_bakta/bakta_"$SampleName"; done
```

### Running pseudofinder (example for Acetobacteraceae_TRYOCC)
```shell
pseudofinder.py annotate --genome bakta_AcetoTRYOCC/AcetoTRYOCC.gbk --outprefix AcetoTRYOCC_pseudofinder --database ~/db/nr --threads $NUM_THREADS --skip_makedb
```

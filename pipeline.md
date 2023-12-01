## Software requirement
* [Trim_galore](https://github.com/FelixKrueger/TrimGalore)
* [Anvi'o](https://anvio.org/)
* [Megahit](https://github.com/voutcn/megahit)
* [Bowtie2](https://github.com/BenLangmead/bowtie2)
* [Samtools](https://www.htslib.org/)
* [hmmer (version 3.2.1)](http://www.hmmer.org/)
* [NCBI blast v2.11.0+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
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
#anvi-run-ncbi-cogs -c TRYOCC.db --num-threads $NUM_THREADS

anvi-get-sequences-for-gene-calls -c TRYOCC.db -o TRYOCC-gene-calls.fa

anvi-profile -i TRYOCC.bam -c TRYOCC.db --output-dir TRYOCC_profile --sample-name TRYOCC -T $NUM_THREADS

pseudofinder.py annotate --genome prokka_SD7/SD7.gbk --outprefix SD7_pseudofinder --database ~/db/nr --threads $NUM_THREADS --skip_makedb


## Contigs taxonomic classification against whole NCBI nucleotide database
```shell
blastn -query TRYOCC_contigs.fasta -db nt -outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle' -evalue 1e-10 -task blastn -num_threads $NUM_THREADS -max_target_seqs 1 -max_hsps 1 > blastn_TRYOCC_contigs.txt
```
## Getting Avg. fold, GC percentage and contig size information from bam file
```shell
pileup.sh in=TRYOCC.bam out=cov_TRYOCC.txt
```

### Merging the taxonomy and coverage outputs with [Nanotax](https://github.com/diecasfranco/Nanotax/blob/main/NanoTax_v2.2.py)


```shell
# You may need to install pip3 before. #
curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py
python3 get-pip.py

# You may need to install or upgrade setuptools and wheel using pip3 before. #
pip3 install --upgrade setuptools wheel

# Download and install MetaDecoder version 1.0.18 #
pip3 install -U https://github.com/liu-congcong/MetaDecoder/releases/download/v1.0.18/metadecoder-1.0.18-py3-none-any.whl
```

# SNP calling
## Quality control
Data comes back from AgResearch in 2 lanes + 1 key file: 'lane1_fastq.txt.gz', 'lane2_fastq.txt.gz', 'key.txt'

Subsample to look at quality:

```
mkdir mohua_gbs
zcat lane1_fastq.txt.gz | head -n 1000000 > lane1.fq
zcat lane2_fastq.txt.gz | head -n 1000000 > lane2.fq
```

load and run FastQC:

```
module load FastQC
fastqc *fq
```

Adapter removal (adapter is Illumina Universal Adapter - AGATCGGAAGAG):

```
module load cutadapt
```

```
cutadapt -a AGATCGGAAGAG -m 30 -o S lane1_cleaned.fastq lane1_fastq.txt.gz
```

```
cutadapt -a AGATCGGAAGAG -m 30 -o lane2_cleaned.fastq lane2_fastq.txt.gz
```

Run FastQC again on the trimmed/cleaned files:

```
fastqc lane1_cleaned.fastq
fastqc lane2_cleaned.fastq
```

## Demultiplexing

First create a barcode file using key.txt and ID file (send by AgResearch). Match unique barcodes to samples ID to create barcodes.txt file

```
mkdir raw samples source_files
```

```
mv lane*cleaned* ./source_files
```

Going into the raw folder and create the link to the raw datafiles:

```
cd ./raw
```

```
ln -s ../source_files/lane*cleaned* . 
```

Go into /mohua_gbs directory and load Stacks: 

```
module load Stacks
```

then run process_radtags command:

```
process_radtags -p ./raw -b barcodes.txt -e pstI -o ./samples -c -q -r --inline_null
```

## Alignment and variant calling


```
mkdir refmap_output samples_mapped
```

Upload mohua reference genome (NCBI accession number: GCA 013398855.1) & index it:

```
head -n 10 mohua_ref_genome.fna
```
```
bwa index mohua_ref_genome.fna
```

Load needed modules:

```
module load BWA
module load SAMtools
```

Map reads:

```
for filename in samples/*fq.gz
do base=$(basename ${filename} .fq.gz)
echo $base
bwa mem -t 8 mohua_ref_genome.fna samples/${base}.fq.gz | samtools view -b | samtools sort --threads 4 > samples_mapped/${base}.bam
done
```

Upload popmap.txt, then run ref_map.pl:

```
ref_map.pl --samples samples_mapped/ -o refmap_output/ --popmap popmap.txt -T 8
```

Outputs a .vcf file which will be used for SNP filtering

## combining dataset from first and second GBS batch

Repeat demultiplexing steps from first batch, do the same for second batch. 
Then copy data for 1st and 2nd plate into one directory.
Go into folder were demultiplexed samples are (i.e., mohua_gbs_plat1/samples and mohua_gbs_plate2/samples), then copy all into a new folder containing all data:

```
cp *fq.gz ../../mohua_gbs/samples_all/
```


To combine 2 samples (e.g., "1R" and "2R"):


```
zcat ../mohua_gbs_plate1/samples/1R.fq.gz ../mohua_gbs_plate2/samples/1R.fq.gz | gzip -c > ../mohua_gbs/samples_all/1R.fq.gz
```

```
zcat ../mohua_gbs_plate1/samples/2R.fq.gz ../mohua_gbs_plate2/samples/2R.fq.gz | gzip -c > ../mohua_gbs/samples_all/2R.fq.gz
```

From here on follow alignment and variant calling script as done for first batch (but this time with data from both batches)

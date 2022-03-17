## Trimming and Quality Control

### Rename sequences for convenience

Keep individual ID and strand direction, discard sequence ID and direction ID

```
id_loc="/Users/hetzler/Amphipod/MiSeq/fastQ/"

cd $id_loc

a="$(ls *.fastq)"

for i in $a
do
	b="$(echo $i | perl -pe 's/(\S)_S[0-9]+(\S)/$1$2/')"
	
	mv $i $b
done
```



### Trim adapters with cutadapt

From the MiSeq sequences trim COI and 18S primers from one fastq file.

Run before BASH scirpt: conda activate cutadaptenv

```
#run first: conda activate cutadaptenv

cutadapt \
-g CHACWAAYCATAAAGATATYGG \
-g ATCCCTGGTTGATCCTGCCAGT \
-G AWACTTCVGGRTGVCCAAARAATCA \
-G CCTGCGCTCGATACTGACAT \
-o trmFastQ/A1_R1.fastq.gz -p trmFastQ/A1_R2.fastq.gz \
--discard-untrimmed A1_S218_L001_R1_001.fastq.gz A1_S218_L001_R2_001.fastq.gz
```

From the MiSeq sequences trim COI and 18S primers from all fastq files.

Before running BASH script conda environment has to be loaded for cutadat: conda activate cutadaptenv


```
#run first: conda activate cutadaptenv


inp_loc="/Users/hetzler/Amphipod/MiSeq/fastq"
out_loc="/Users/hetzler/Amphipod/MiSeq/trimmed"
cd $inp_loc

#trim adaptor for all fastq files
for i in $(eval echo {A..D}); do
  for x in {1..5}; do
      FWD1=CHACWAAYCATAAAGATATYGG  #forward primer COI
      FWD2=ATCCCTGGTTGATCCTGCCAGT  #forward primer 18S
      REV1=AWACTTCVGGRTGVCCAAARAATCA #reverse primer COI
      REV2=CCTGCGCTCGATACTGACAT  #forward primer 18S
      cutadapt -g $FWD1 -g $FWD2 -G $REV1 -G $REV2 \
      -o $out_loc/${i}${x}\_R1.fastq -p $out_loc/${i}${x}\_R2.fastq \
      ${i}${x}\_L001_R2_001.fastq \
      ${i}${x}\_L001_R1_001.fastq
    done
done
```

### create fastQC reports for multiQC

```
set_loc="/Users/hetzler/Amphipod/MiSeq/trmFastQ"

cd $set_loc

find . -name "*.fastq" | xargs fastqc 
```

Run multiqc . for multiQC report

## Sequence mapping

### Create reference sequence index
Create index reference files for bowtie2 alignment.

reference.fasta contains COI and 18S reference sequence. 

```
bowtie2-build "/Users/hetzler/Amphipod/referenceSeq/reference.fasta" "/Users/hetzler/Amphipod/referenceSeq/reference"
```

### Map reads to reference sequence

Using bowtie2 to map reads to reference sequence. 

How to manipulate sam/bam files:
https://medium.com/@shilparaopradeep/samtools-guide-learning-how-to-filter-and-manipulate-with-sam-bam-files-2c28b25d29e8

#### Script map.sh
```
inp_loc="/Users/hetzler/Amphipod/MiSeq/trimmed"
out_loc="/Users/hetzler/Amphipod/mappedReads"

#map to reference
for i in $(eval echo {A..D}); do
  for x in {1..5}
    do
      bowtie2 --end-to-end -N 0 -L 20 --dpad 15 --gbar 4 --seed 0 --threads 1 -q \
      --rg-id ${i}${x} \
      -x "/Users/hetzler/Amphipod/referenceSeq_test/reference" \
      -U $inp_loc/${i}${x}\_R1.fastq,$inp_loc/${i}${x}\_R2.fastq \
      -S $out_loc/${i}${x}.sam
    done
done
```

### Sort and convert from .sam to .bam

```
set_loc="/Users/hetzler/Amphipod/mappedReads"

cd $set_loc

# samtools:  sort .sam file and convert to .bam file
for i in $(eval echo {A..D}); do
	for x in {1..5}; do
    samtools view -bS ${i}${x}.sam | samtools sort - -o sort_${i}${x}.bam
  done
done
```

#### variant calling

Call variants from mapped and aligned sequences using bcftools

[mpileup] maximum number of reads per input file set to -d 250

```
ref_file="/Users/hetzler/Amphipod/referenceSeq/reference.fasta"
id_loc="/Users/hetzler/Amphipod/mappedReads"
out_loc="/Users/hetzler/Amphipod/variantCalling"
cd $id_loc

#variant calling
for i in $(eval echo {A..D}); do
	for x in {1..5}; do
    bcftools mpileup -f $ref_file sort_${i}${x}.bam | bcftools call -mv -Ob -o $out_loc/calls_${i}${x}.vcf.gz
  done
done
```

Create variant calling index file. 

```
set_loc="/Users/hetzler/Amphipod/variantCalling"

cd $set_loc

#index vcf
for i in $(eval echo {A..D}); do
  for x in {1..5}; do
   bcftools index calls_${i}${x}.vcf.gz
    done
done
```


```
set_loc="/Users/hetzler/Amphipod/variantCalling"

cd $set_loc

ls calls_*.vcf.gz > merge.txt

bcftools merge -l merge.txt -0 -Oz -o pop.vcf.gz

bcftools merge -m indels -l merge.txt -0 -Oz -o indel.vcf.gz

bcftools annotate -x INFO,^FORMAT/GT pop.vcf.gz -Oz -o popAno.vcf.gz
```


Parse to GENO file

```
set_loc="/Users/hetzler/Amphipod/variantCalling"

cd $set_loc

parseVCF.py -i pop_merged.vcf.gz -o pop.geno
```


Create a file containing population information using grep

```
set_loc="/Users/hetzler/Amphipod/variantCalling"
inp_loc="/Users/hetzler/Amphipod/mappedReads"

cd $set_loc

ls $inp_loc/sort* > samples

#CREATE POPULATION IDENTIFYER FILE 
grep  'sort_A' samples > SAL
grep  'sort_B' samples >> SAL
grep 'sort_C' samples > SKJ
grep 'sort_D' samples >> SKJ

awk '{print $1"\tSAL"}' SAL > pop_file
awk '{print $1"\tSKJ"}' SKJ >> pop_file
```

population genomics analysis
```
gzip pop.geno

popgenWindows.py -g pop.geno.gz -o div_stat.csv -f phased -w 1500 -m 5 -s 25000 -p SAL -p SKJ --popsFile pop_file --writeFailedWindow
```

#### Multifasta from consensus files.

https://samtools.github.io/bcftools/howtos/consensus-sequence.html

__NB: try to normalized indels and filter calls__

Create consensus files from reference fasta and indexed VCF file

```
#!/bin/bash
set_loc="/home/jhetzler/MiSeqOSL/VC"

cd $set_loc

for i in $(eval echo {A..D}); do
  for x in {1..5}; do
cat /home/jhetzler/MiSeqOSL/ref/reference.fasta | bcftools consensus -s ${i}${x} ${i}${x}.vcf.gz > multi/${i}${x}.fa
  done
done
```

Add prefix of fasta headers
````
#!/bin/bash
set_loc="/home/jhetzler/MiSeqOSL/VC"
cd $set_loc

for i in $(eval echo {A..D}); do
  for x in {1..50}; do
        perl -pi -e "s/^>/>${i}${x}_/g" multi/${i}${x}.fa
  done
done
```

Add suffix of fasta headers
````
#!/bin/bash
set_loc="/home/jhetzler/MiSeqOSL/VC"
cd $set_loc

for i in $(eval echo {A..D}); do
  for x in {1..50}; do
        perl -pi -e 's/^(>.*)$/$1-'$i''$x'/g' multi/${i}${x}.fa
  done
done
```

create multi.fasta file from fasta files
```
cat *.fa > multi.fas
```

separate multifasta files

```
#!/bin/bash

while read line ; do
  if [ ${line:0:1} == ">" ] ; then
    filename=$(echo "$line" | cut -d "_" -f1 | tr -d ">")
    touch ./"$filename".fasta
    echo "$line" >> ./"${filename}".fasta
  else
    echo "$line" >> ./"${filename}".fasta
  fi
done < $1
```


Masked fasta???

````
1. Call variants.

2. Identify positions with nocalls (say, by emitting all sites) and convert to a BED:

grep "\./\." [YOUR_VCF] | awk '{OFS="\t"; if ($0 !~ /\#/); print $1, $2-1, $2}' > [YOUR_FILTERED_BED]
3. Use BEDtools' maskfasta.
```


## Fasta Alternate Reference Maker

Create reference sequence dictionary

```
java -jar /home/jhetzler/tools/picard/picard.jar CreateSequenceDictionary R=../ref/reference.fasta O=../ref/reference.dict
```

Create index files for vcf

```
gatk IndexFeatureFile -I A1.vcf.gz
```

Create an alternative reference by combining a fasta with a vcf.

https://gatk.broadinstitute.org/hc/en-us/articles/360037594571-FastaAlternateReferenceMaker

```
gatk FastaAlternateReferenceMaker -R /home/jhetzler/MiSeqOSL/ref/reference.fasta -O A1_FARM.fasta -V A1.vcf.gz
```


## Haplotype Network

PopART
DNAsp
Arlequin

# TO DO:

Check Alignment score
Haplotype network



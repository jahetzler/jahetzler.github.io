## Basic Linux commands

Count files in folder using list (ls) and word count (wc)

```
ls /path | wc -l
```

Check folder size (du)

```
du -sh /path
```

Copy files or folder from server using scp (-r copy folder)

```
# scp -r [FROM:PATH] [TO:PATH]

scp -r username@192.168.0.1:/path/ /path

scp username@192.168.0.1:/path/file /path
```

Read length of fastq file

echo $(zcat *.fastq.gz|wc -l)/4|bc

## BASH 

#### Nested loop

Set location and run nested loop example

```
#!/bin/bash

inp_loc="/Users/path"
cd $inp_loc

#trim adaptor for all fastq files
for i in $(eval echo {A..D}); do
  for x in {1..5}; do
      print${i}${x}
    done
done
```

## BCFTools

#### Manual:
https://samtools.github.io/bcftools/bcftools.html#view

#### View

View VCF file

```
bcftools view file.vcf.gz
```

Extract specific chromosome data

```
bcftools view file.vcf.gz --regions [CHROM NAME] > output.vcf

bgzip output.vcf
```

Filter on minor allele frequency MAF 
```
bcftools view -q 0.05:minor pop.vcf.gz > popMAF.vcf
```

#### Merge:

Merge VCF files

```
ls *.vcf.gz > merge.txt

bcftools merge -l merge.txt -0 -Oz -o pop.vcf.gz
```


#### Index:

View per contig stats (SNPs)

```
#multi loci file
bcftools index -s input.vcf.gz

#single loci file
bcftools index -n input.vcf.gz
```

Create VCF index files

```
bcftools index file.vcf.gz
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

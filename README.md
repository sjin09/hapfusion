## <a name="started"></a>Getting started

```sh
## clone github repository and install hapfusion
git clone https://github.com/sjin09/hapfusion
cd hapfusion
/bin/bash install.sh

## use miniamp2 and samtools to align, sort (and merge) PacBio CCS read alignments
minimap2 "@RG\tSM:sample" -ax map-hifi --cs ref.fa pacbio.ccs.fastq.gz | samtools sort -o aln.sorted.bam # SM tag must be provided to retrieve sample ID
samtools merge *.sorted.bam | samtools sort -o aln.mergeSorted.bam - ## if there are multiple BAM files, merge and sort the BAM files 
samtools view -bh -F 0x900 aln.sorted.bam > aln.primary_alignments.sorted.bam # select primary alignments
samtools index aln.sorted.bam
samtools index aln.primary_alignments.sorted.bam 

## use deepvariant to call germline mutations
deepvariant.simg /opt/deepvariant/bin/run_deepvariant --model_type=PACBIO --ref ref.fa --reads=aln.primary_alignments.sorted.bam --output_vcf=germline.vcf

## use hapfusion to phase germline mutations
hapfusion phase --bam aln.primary_alignments.sorted.bam --vcf germline.vcf.bgz -o germline.phased.vcf 
bgzip -c germline.phased.vcf > germline.phased.vcf.bgz
tabix -p vcf germline.phased.vcf.bgz

## call crossovers and gene conversions
hapfusion call -i aln.primary_alignments.sorted.bam --vcf germline.phased.vcf.bgz -o recombinations.txt

## plot crossovers and gene conversions
hapfusion plot -i aln.primary_alignments.sorted.bam --vcf germline.phased.vcf.bgz --fusion recombinations.txt -o pdf
```

## <a name="uguide"></a>Users' Guide

hapfusion leverages CCS base accuracy and read length to crossovers (CO), gene conversions (NCO) and complex gene conversions (CNCO).

#### <a name="alignment"></a>minimap2 CCS read alignment

Currently, hapfusion only accepts minimap2 CCS read alignment SAM/BAM file with a `--cs` tag, like the `cigar` string provides information about the matches and mismatches between the read and the reference genome. pbmm2 and ngmlr does not provide the `--cs` tag and as a result, pbmm2 and ngmlr SAM/BAM files are incompatible with himut. The `--cs` tag is favored over the `cigar` string because it offers a more elegant representation of the differences between the read and the reference genome (source: https://lh3.github.io/2018/03/27/the-history-the-cigar-x-operator-and-the-md-tag).

#### <a name="general"></a>General Usage

hapfusion accepts as input a BAM file with primary CCS read alignments and a VCF file with phased germline mutations and returns a file with crossovers and gene conversions. 

```sh
## phase germline mutations
hapfusion phase --bam aln.primary_alignments.sorted.bam --vcf germline.vcf.bgz -o germline.phased.vcf

## call crossovers and gene conversions
hapfusion call -i aln.primary_alignments.sorted.bam --vcf germline.phased.vcf.bgz -o recombinations.txt

## plot crossovers and gene conversions
hapfusion plot -i aln.primary_alignments.sorted.bam --vcf germline.phased.vcf.bgz --fusion recombinations.txt -o pdf
```

#### <a name="limits"></a>Limitations

The number of crossovers and gene conversions detcted from the sample is dependent on the number of phased hetSNPs, length of haplotype blocks and CCS read length.

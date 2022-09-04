## <a name="started"></a>Getting started

```sh
## clone github repository and install hapmash
git clone https://github.com/sjin09/hapmash
cd hapmash
/bin/bash install.sh

## use miniamp2 and samtools to align, sort (and merge) PacBio CCS read alignments
minimap2 "@RG\tSM:sample" -ax map-hifi --cs ref.fa pacbio.ccs.fastq.gz | samtools sort -o aln.sorted.bam # SM tag must be provided to retrieve sample ID
samtools view -bh -F 0x900 aln.sorted.bam > aln.primary_alignments.sorted.bam # select primary alignments
samtools index aln.primary_alignments.sorted.bam
samtools merge *.primary_alignments.sorted.bam | samtools sort -o aln.primary_alignments.mergeSorted.bam -
samtools index aln.primary_alignments.mergeSorted.bam 

## use deepvariant to call germline mutations
deepvariant.simg /opt/deepvariant/bin/run_deepvariant --model_type=PACBIO --ref ref.fa --reads=aln.primary_alignments.sorted.bam --output_vcf=germline.vcf

## use hapmash to phase germline mutations
hapmash phase --bam aln.primary_alignments.sorted.bam --vcf germline.vcf.bgz -o germline.phased.vcf 

## call crossovers and gene conversions
hapmash call -i aln.primary_alignments.sorted.bam --vcf germline.phased.vcf -o recombinants.vcf 
```

## <a name="uguide"></a>Users' Guide

hapmash (hapmash) is leverages the base accuracy and read length of Pacific Biosciences (PacBio) CCS reads to call crossovers and gene conversions. 

### <a name="general"></a>General Usage

hapmash accepts as input BAM file with CCS read alignments and VCF file with phased germline mutations and returns a VCF file with crossovers and gene conversions. hapmash, currently, only accepts minimap2 generated BAM files as `--cs` tag is required for somatic mutation calling. The read alignments must be sorted, compressed and indexed using SAMtools. In addition, secondary and supplementary alignment needs to be removed from the BAM file as deepvariant can only process primary alignments.

```sh
minimap2 "@RG\tSM:sample" -ax map-hifi --cs ref.fa pacbio.ccs.fastq.gz | samtools sort -o aln.sorted.bam # SM tag must be provided to retrieve sample ID
samtools view -bh -F 0x900 aln.sorted.bam > aln.primary_alignments.sorted.bam # select primary alignments
samtools merge *.primary_alignments.sorted.bam | samtools sort -o aln.primary_alignments.mergeSorted.bam -
samtools index aln.primary_alignments.sorted.bam
samtools index aln.primary_alignments.mergeSorted.bam 
```

#### Haplotype phase heterozygous single nucleotide polymorphisms 

CCS reads have an average read length of ~10-20kb and can span a number of heterozygous single nucleotide polymorphisms (hetSNPs). To haplotype phase the hetSNPs, we treat each hetSNP as a node in a graph and as CCS reads are able to span multiple hetSNPs, we are able to count the number of edges between two nodes and determine whether the two nodes belong to the same haplotype. Two or more haplotype consistent nodes are connected to construct a haplotype block across the chromosome and we are able to assign CCS reads to a haplotype block based on their hetSNP composition. If the CCS read belongs to multiple haplotype blocks or if the CCS read does not have a hetSNP, CCS read is determined to be not phased. We would like to highlight that the length of haplotype blocks and the number of phased hetSNPs will be dependent on SNP density, heterozygosity and CCS read length.

```sh
hapmash phase --bam aln.primary_alignments.sorted.bam --vcf germline.vcf.bgz -o germline.phased.vcf 
```

hetSNPs belonging to the same phase set (PS) in the VCF file belongs to the same haplotype block.

#### call crossovers and gene conversions 

```sh
hapmash call -i aln.primary_alignments.sorted.bam --vcf germline.vcf -o recombinant.vcf 
```

### Limitations

As mentioned above, the number of called crossovers and gene conversions will be dependent on the number of phased hetSNPs and length of haplotype blocks, which will be dependent on SNP density, heterozygosity and CCS read length.


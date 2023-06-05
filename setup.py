# -*- coding: utf-8 -*-
from setuptools import setup

package_dir = \
{'': 'src'}

packages = \
['hapfusion']

package_data = \
{'': ['*']}

install_requires = \
['argparse>=1.4.0,<2.0.0',
 'biopython>=1.78,<2.0',
 'cyvcf2>=0.30.18,<0.31.0',
 'natsort>=8.0.0,<9.0.0',
 'numpy>=1.20.2,<2.0.0',
 'plotnine>=0.9.0,<0.10.0',
 'psutil>=5.8.0,<6.0.0',
 'pyfastx>=0.8.4,<0.9.0',
 'pysam>=0.18.0,<0.19.0',
 'pytabix>=0.1,<0.2']

entry_points = \
{'console_scripts': ['hapfusion = hapfusion.__main__:main']}

setup_kwargs = {
    'name': 'hapfusion',
    'version': '0.0.1',
    'description': 'hapfusion: Gene conversion and crossover detection using PacBio CCS reads',
    'long_description': '## <a name="started"></a>Getting started\n\n```sh\n## clone github repository and install hapfusion\ngit clone https://github.com/sjin09/hapfusion\ncd hapfusion\n/bin/bash install.sh\n\n## use miniamp2 and samtools to align, sort (and merge) PacBio CCS read alignments\nminimap2 "@RG\\tSM:sample" -ax map-hifi --cs ref.fa pacbio.ccs.fastq.gz | samtools sort -o aln.sorted.bam # SM tag must be provided to retrieve sample ID\nsamtools view -bh -F 0x900 aln.sorted.bam > aln.primary_alignments.sorted.bam # select primary alignments\nsamtools index aln.primary_alignments.sorted.bam\nsamtools merge *.primary_alignments.sorted.bam | samtools sort -o aln.primary_alignments.mergeSorted.bam -\nsamtools index aln.primary_alignments.mergeSorted.bam \n\n## use deepvariant to call germline mutations\ndeepvariant.simg /opt/deepvariant/bin/run_deepvariant --model_type=PACBIO --ref ref.fa --reads=aln.primary_alignments.sorted.bam --output_vcf=germline.vcf\n\n## use hapfusion to phase germline mutations\nhapfusion phase --bam aln.primary_alignments.sorted.bam --vcf germline.vcf.bgz -o germline.phased.vcf \n\n## call crossovers and gene conversions\nhapfusion call -i aln.primary_alignments.sorted.bam --vcf germline.phased.vcf -o recombinants.vcf \n```\n\n## <a name="uguide"></a>Users\' Guide\n\nhapfusion (hapfusion) is leverages the base accuracy and read length of Pacific Biosciences (PacBio) CCS reads to call crossovers and gene conversions. \n\n### <a name="general"></a>General Usage\n\nhapfusion accepts as input BAM file with CCS read alignments and VCF file with phased germline mutations and returns a VCF file with crossovers and gene conversions. hapfusion, currently, only accepts minimap2 generated BAM files as `--cs` tag is required for somatic mutation calling. The read alignments must be sorted, compressed and indexed using SAMtools. In addition, secondary and supplementary alignment needs to be removed from the BAM file as deepvariant can only process primary alignments.\n\n```sh\nminimap2 "@RG\\tSM:sample" -ax map-hifi --cs ref.fa pacbio.ccs.fastq.gz | samtools sort -o aln.sorted.bam # SM tag must be provided to retrieve sample ID\nsamtools view -bh -F 0x900 aln.sorted.bam > aln.primary_alignments.sorted.bam # select primary alignments\nsamtools merge *.primary_alignments.sorted.bam | samtools sort -o aln.primary_alignments.mergeSorted.bam -\nsamtools index aln.primary_alignments.sorted.bam\nsamtools index aln.primary_alignments.mergeSorted.bam \n```\n\n#### Haplotype phase heterozygous single nucleotide polymorphisms \n\nCCS reads have an average read length of ~10-20kb and can span a number of heterozygous single nucleotide polymorphisms (hetSNPs). To haplotype phase the hetSNPs, we treat each hetSNP as a node in a graph and as CCS reads are able to span multiple hetSNPs, we are able to count the number of edges between two nodes and determine whether the two nodes belong to the same haplotype. Two or more haplotype consistent nodes are connected to construct a haplotype block across the chromosome and we are able to assign CCS reads to a haplotype block based on their hetSNP composition. If the CCS read belongs to multiple haplotype blocks or if the CCS read does not have a hetSNP, CCS read is determined to be not phased. We would like to highlight that the length of haplotype blocks and the number of phased hetSNPs will be dependent on SNP density, heterozygosity and CCS read length.\n\n```sh\nhapfusion phase --bam aln.primary_alignments.sorted.bam --vcf germline.vcf.bgz -o germline.phased.vcf \n```\n\nhetSNPs belonging to the same phase set (PS) in the VCF file belongs to the same haplotype block.\n\n#### call crossovers and gene conversions \n\n```sh\nhapfusion call -i aln.primary_alignments.sorted.bam --vcf germline.vcf -o recombinantions.txt\n```\n\n### Limitations\n\nAs mentioned above, the number of called crossovers and gene conversions will be dependent on the number of phased hetSNPs and length of haplotype blocks, which will be dependent on SNP density, heterozygosity and CCS read length.\n\n',
    'author': 'Sangjin Lee',
    'author_email': 'sl17@sanger.ac.uk',
    'maintainer': 'None',
    'maintainer_email': 'None',
    'url': 'https://github.com/sjin09/hapfusion',
    'package_dir': package_dir,
    'packages': packages,
    'package_data': package_data,
    'install_requires': install_requires,
    'entry_points': entry_points,
    'python_requires': '>=3.8,<3.10',
}


setup(**setup_kwargs)


# PAVFinder_transcriptome


Author:	Readman	Chiu (rchiu@bcgsc.ca)

PAVFinder_transcript (**PVT**) is a Python package written to identify structural variants in transcriptome assemblies.
In a nutshell, the algorithm infers variants from non-contiguous (split or gapped) alignments of assembled contig sequences to the reference genome.
With the aid of gene-model annotation(s), diversified classes of variants such as gene fusions, read-throughs, internal and partial tandem duplications, indels and novel splice variants are classified.

The program is usually preceded by de novo assembly of RNAseq sequences followed by alignment to the reference genome.
As such, a pipeline that bundles the 3 analysis	steps called **TAP** (**T**ransabyss-**A**lignment-**P**AVFinder) is provided as a standalone application.  TAP can also be run in targeted mode on selected genes.  This requires a Bloom Filter of target gene sequences to be created beforehand.  Whereas the full assembly of a RNAseq library with over 100 million read pairs requires more than 24 hours to complete, a target assembly and analysis of a gene list (e.g. COSMIC) of several hundred can be completed within half an hour.

## Requirements
1. External software

 * [BWA-mem](http://bio-bwa.sourceforge.net/)
 * [GMAP](http://research-pub.gene.com/gmap/)
 * [samtools](http://samtools.sourceforge.net/)
 * [BBT](http://www.bcgsc.ca/platform/bioinfo/software/biobloomtools) (if TAP is run in targeted mode)

2. Reference files

 * single reference genome FASTA indexed by samtools faidx and GMAP
 * gene model(s) in VCF format with chromosome names matching reference genome

See INSTALL for more details

## Usage
1. Run PVT (for structural variants)

 ```python
 find_sv.py --gbam <contigs_to_genome_bam> --tbam <contigs_to_transcripts_bam> --transcripts_fasta <indexed_transcripts_fasta> --genome_index <GMAP index genome directory and name> --r2c <reads_to_contigs_bam> <contigs_fasta> <gtf> <genome_fasta> <outdir>
 ```

2. Run TAP

  ```tap.py <sample> <outdir> --bf <target_genes.bf> --fq_list <file_listing_FASTQ_pairs> --k <space-delimited k values> --readlen <read_length>  --nprocs <number_of_processes> --params <parameters_file>```

See USAGE for more details
* Run PVT to detect structural variants (fusions, read-throughs, ITDs, PTDs, InDels)

  ```
  find_sv.py --gbam <contigs_to_genome_bam> --tbam <contigs_to_transcripts_bam> --transcripts_fasta <indexed_transcripts_fasta> --genome_index <GMAP index genome directory and name> --r2c <reads_to_contigs_bam> <contigs_fasta> <gtf> <genome_fasta> <outdir>
  ```

* Run PVT to detect novel splice variants (exon_skipping, novel_exon, novel_intron, novel_donor, novel_acceptor, retained_intron)

  ```
  map_splice.py <contigs_to_genome_bam> <contigs_fasta> <gtf> <genome_fasta> <outdir> --r2c <reads_to_contigs_bam> --suppl_annot <supplmental_annotations>
  ```

* Run full (assembly + analysis) TAP in targeted model

  ```
  tap.py <sample> <outdir> --bf <target_genes.bf> --fq_list <file_listing_FASTQ_pairs> --k <space-delimited k values> --readlen <read_length>  --nprocs <number_of_processes> --params <parameters_file>
  ```

* Run full (assembly + analysis) TAP for entire transcriptome

  ```
  tap.py <sample> <outdir> --fq_list <file_listing_FASTQ_pairs> --k <space-delimited k values> --readlen <read_length> --nprocs <number_of_processes> --params <parameters_file>
  ```

* Run TAP for just de novo assembly

  ```
  tap.py <sample> <outdir> --fq_list <file_listing_FASTQ_pairs> --k <space-delimited k values> --readlen <read_length> --nprocs <number_of_processes> --params <parameters_file> --only_assembly
  ```
  
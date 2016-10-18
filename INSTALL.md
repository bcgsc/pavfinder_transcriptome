# Installation

## External software

* [BWA-mem](http://bio-bwa.sourceforge.net/) v0.7.12
* [GMAP](http://research-pub.gene.com/gmap/) 2014-12-28
* [samtools](http://samtools.sourceforge.net/) v0.1.19
* [BBT](http://www.bcgsc.ca/platform/bioinfo/software/biobloomtools) v3.0.0b (required) (if TAP is run in targeted mode)

*tested versions indicated, may not be the most recent version.  For BBT, v3.0.0b is required.


## External files
To run PAVFinder, the following reference sequence and annotation files are required:

* genome FASTA and index files
 * For example, to prepare hg19 files from UCSC:
 
   ```bash
   wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr*.fa.gz
   zcat chr*.fa.gz > hg19.fa && rm chr*.fa.gz
   samtools faidx hg19.fa
   bwa index hg19.fa
   gmap_build -D . -d hg19 hg19.fa
   ```
* annotation VCF and index files
 * For example, to prepare Refseq genes from UCSC:

   ```bash
   wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/refGene.txt.gz
   zcat refGene.txt.gz | cut -f2- | ~rchiu/bin/genePredToGtf file stdin refGene.gtf
   sort -k1,1 -k4,4n refGene.gtf > refGene.sorted.gtf
   bgzip refGene.sorted.gtf
   tabix -p gff refGene.sorted.gtf.gz
   ```

*  transcriptome FASTA and index files
 * To create transcriptome FASTA file PAVFinder provides a script for doing it (after Pip virtualenv install it should be present under '/bin/')

   ```extract_transcript_sequence.py <gtf> <fasta_output> <genome_fasta> --index```

   A transcript FASTA and a corresponding BWA index will be genearated

If running TAP, after installation of PVT(below) the full paths of the above reference and annotation files can be specified in /path/to/virtualenv/config/tap.cfg which can be passed to TAP via --params

  ```bash
  genome_index = /real/path/to/gmapdb hg19
  transcripts_fasta = /real/path/to/transcripts.fa
  genome_fasta = /real/path/to/genome.fa
  gtf = /real/path/to/gtf
  suppl_annot = /real/path/to/supplementary_gtf (optional)
  ```

## Installing the Python PVT package

1. ```pip install virtualenv```
2. ```virtualenv <DIR>```
3. ```source <DIR>/bin/activate```
4. ```pip install -U cython```
5. ```pip install git+https://github.com/bcgsc/pavfinder_transcriptome.git#egg=pavfinder_transcriptome```

After successful installation, the following will be made avaiable in different sub-directories under the virtualenv directory:

1. Python scripts for detecting structural (find\_sv.py) and splice (map\_splice.py) variants, and the TAP pipeline script (tap.py) will be in "bin"
2. A template config file (tap.cfg) for running TAP will be in "config"
3. Sample data will be in "test"
# Installation

## External software

* [BWA-mem](http://bio-bwa.sourceforge.net/) v0.7.12
* [GMAP](http://research-pub.gene.com/gmap/) 2014-12-28
* [samtools](http://samtools.sourceforge.net/) v0.1.19
* [BBT](http://www.bcgsc.ca/platform/bioinfo/software/biobloomtools) v3.0.0 (required) (if TAP is run in targeted mode)

*tested version is provided, may not be the most recent version.  For BBT, v3.0.0 is required.


## External files

## PAVFinder_transcriptome

1. ```pip install virtualenv```
2. ```virtualenv <DIR>```
3. ```source <DIR>/bin/activate```
4. ```pip install -U cython```
5. ```pip install git+https://github.com/bcgsc/pavfinder_transcriptome.git#egg=pavfinder_transcriptome```

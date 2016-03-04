Installation
------------

Install the following softwares and put there executables in ``PATH``:

-  `BWA 0.7.4
   <http://sourceforge.net/projects/bio-bwa/files/>`_
-  `Samtools
   <http://sourceforge.net/projects/samtools/files/samtools/>`_
-  `GMAP 2014-12-28
   <http://research-pub.gene.com/gmap/src/gmap-gsnap-2014-01-21.tar.gz>`_

Install pavfinder-transcriptome:

::

   $ pip install virtualenv
   $ virtualenv <DIR>
   $ source <DIR>/bin/activate
   $ pip install -U cython
   $ pip install git+https://github.com/bcgsc/pavfinder_transcriptome.git#egg=pavfinder_transcriptome

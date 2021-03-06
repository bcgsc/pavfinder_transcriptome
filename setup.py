import os
from setuptools import setup, find_packages
from pavfinder_transcriptome import __version__

setup(
    name='pavfinder_transcriptome',
    version=__version__,
    description='Post Assembly Variant Finder - transcriptome',
    long_description='Identifies transcriptomic structural variants from sequence assembly',
    url='https://github.com/bcgsc/pavfinder_transcriptome.git',
    author='Readman Chiu',
    author_email='rchiu@bcgsc.ca',
    license='BCCA',
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Programming Language :: Python :: 2.7',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        ],
    packages=find_packages(),
    install_requires = [
        'pysam==0.8.2.1',
        'pybedtools==0.6.2',
        'intspan>=0.701',
        'numpy==1.9.2',
        'pandas',
        'biopython',
        ],
    package_data = {'pavfinder_transcriptome': ["test/*", "scripts/*.py"]},
    data_files = [('config', ['pavfinder_transcriptome/cfg/tap.cfg']), ('test', ['pavfinder_transcriptome/test/sample_1.fastq.gz', 'pavfinder_transcriptome/test/sample_2.fastq.gz', 'pavfinder_transcriptome/test/sample_tap_output.tar.gz', 'pavfinder_transcriptome/test/cancer_census.tar.gz', 'pavfinder_transcriptome/test/README'])],
    scripts = ['pavfinder_transcriptome/scripts/find_sv.py',
               'pavfinder_transcriptome/scripts/map_splice.py',
               'pavfinder_transcriptome/scripts/tap.py',
               'pavfinder_transcriptome/scripts/extract_transcript_sequence.py',
               'pavfinder_transcriptome/scripts/rescue_fusion.py',
               ],
)

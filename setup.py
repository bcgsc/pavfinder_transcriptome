import os
from setuptools import setup, find_packages
execfile(os.path.dirname(os.path.realpath(__file__)) + "/pavfinder_transcriptome/version.py")
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
        'pysam>=0.8.1,<0.9.0',
        'pybedtools>=0.6.2',
        'intspan>=0.701',
        'numpy',
        'pandas',
        'biopython',
        ],
    scripts = ['pavfinder_transcriptome/scripts/find_events.py',
               ],
)

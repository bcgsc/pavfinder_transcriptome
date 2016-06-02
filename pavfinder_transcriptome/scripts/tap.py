#!/usr/bin/env python
from ruffus import *
import ruffus.cmdline as cmdline
import subprocess
from collections import defaultdict
import gzip
import os
import sys
import glob
import re
import fileinput
import datetime
import itertools

def run_cmd(cmd, force=False):
    process = subprocess.Popen(cmd,
                               stdout = subprocess.PIPE,
                               stderr = subprocess.PIPE,
                               shell = True)

    stdout_str, stderr_str = process.communicate()

    if process.returncode != 0 and not force:
        raise Exception("Failed to run '%s'\n%s%sNon-zero exit status %s" %
                        (cmd, stdout_str, stderr_str, process.returncode))

    return stdout_str, stderr_str
    
def format_read_pairs(fqs=None, list_file=None):
    fastqs = []
    if list_file is not None:
        with open(list_file, 'r') as ff:
            for line in ff:
                fastqs.append(line.rstrip('\n'))
                
    elif fastqs is not None:
        fastqs = fqs
            
    fastqs1 = sorted([fq for fq in fastqs if '1.' in fq])
    fastqs2 = sorted([fq for fq in fastqs if '2.' in fq])
    
    if len(fastqs1) != len(fastqs2):
        raise Exception('input fastqs not paired')
    
    return fastqs1, fastqs2


parser = cmdline.get_argparse(description='TAP pipeline')
parser.add_argument('sample', type=str, help='sample name')
parser.add_argument('outdir', type=str, help='output directory')
parser.add_argument('--bf', type=str, help='path to bloomfilter')
parser.add_argument('--fq', type=str, nargs='+', help='input gzipped fastqs')
parser.add_argument('--fq_list', type=str, help='text file of input fastq paths')
parser.add_argument('--nprocs', type=int, default=64, help='number of threads/processes')
parser.add_argument('--k', type=int, nargs='+', help='k sizes for assembly')
parser.add_argument('--readlen', type=int, help='read length')
parser.add_argument('--gmap_index', type=str, nargs=2, help='gmap index')
parser.add_argument('--bwa_index', type=str, help='bwa index of transcript sequences')
parser.add_argument('--gtf', type=str, help='gtf')
parser.add_argument('--sort_mem', type=str, help='samtools sort memory. Default:5G', default='5G')
parser.add_argument('--genome_fasta', type=str, help='genome fasta')

args = parser.parse_args()

logs_dir = args.outdir + '/logs'
if not os.path.exists(logs_dir):
    os.makedirs(logs_dir)

log_file = '%s/log.%s.txt' % (logs_dir, datetime.datetime.now().strftime("%Y-%m-%d_%H:%M:%S"))
logger, logging_mutex = cmdline.setup_logging(__name__,
                                              log_file,
                                              args.verbose)
print 'log_file:', log_file

cmdline.run(args)

read_pairs = []
if args.fq:
    read_pairs = format_read_pairs(fqs=args.fq)
elif args.fq_list:
    read_pairs = format_read_pairs(list_file=args.fq_list)

bbt_outdir = '%s/bbt_v3.0.0b' % args.outdir
assembly_outdir = '%s/transabyss_v1.5.4' % args.outdir
pvt_outdir = '%s/pvt_v0.3.0' % args.outdir

bbt_prefix = bbt_outdir + '/' + args.sample

assembly_input = []

@follows(mkdir(bbt_outdir))
@transform([read_pairs],
           formatter('_1.fastq.gz$'),
           bbt_outdir + '/' + args.sample + '_reads.fq',
           bbt_prefix,
           args.bf,
           args.nprocs)
def classify(paired_fqs, bbt_output, prefix, bf, nthreads):  
    if len(paired_fqs[0]) == 1:
        input_fq1 = paired_fqs[0][0]
    else:
        input_fq1 = '<(zcat %s)' % ' '.join(paired_fqs[0])
    if len(paired_fqs[1]) == 1:
        input_fq2 = paired_fqs[1][0]
    else:
        input_fq2 = '<(zcat %s)' % ' '.join(paired_fqs[1])
    
    cmd = 'biobloomcategorizer --fq -i -p %s -a 2 -t %d -e -f %s %s %s' % (prefix,
                                                                           nthreads,
                                                                           bf,
                                                                           input_fq1,
                                                                           input_fq2)
        
    run_cmd('/bin/bash -c "%s"' % cmd)
    
@split(classify,
       bbt_outdir + '/*.fastq')
def split_input(bbt_fastq, split_fastqs, genes=None):
    seqs1 = defaultdict(list)
    seqs2 = defaultdict(list)

    with open(bbt_fastq[0], 'r') as fq:
        for lines in itertools.izip_longest(*[fq]*8):
            target = lines[0].rstrip().split(' ')[1].split('.')[0]
            seqs1[target].extend(lines[:4])
            seqs2[target].extend(lines[-4:])

    split_fastqs = output_split_pairs(seqs1, seqs2, bbt_outdir)

def output_split_pairs(seqs1, seqs2, outdir):
    fqs = []
    for target in sorted(seqs1.keys()):
        fq1 = '%s/%s_1.fastq' % (outdir, target)
        fq2 = '%s/%s_2.fastq' % (outdir, target)
        with open(fq1, 'w') as out1:
            out1.writelines(seqs1[target])
                
        with open(fq2, 'w') as out2:
            out2.writelines(seqs2[target])

        fqs.append(fq1)
        fqs.append(fq2)
        
    return fqs

@collate(split_input,
         formatter(".+/(.+)_[12].fastq"),
         ["{path[0]}/{1[0]}_1.fastq",
          "{path[0]}/{1[0]}_2.fastq"])
def pair_reads(fastqs, pair):
    pass

@follows(mkdir(assembly_outdir))
@subdivide(pair_reads,
           formatter(".+/(.+)_1.fastq"),
           assembly_outdir + "/{1[0]}/k*",
           args.k)
def symlink_assembly_input(read_pairs, k_dirs, ks):  
    gene = os.path.basename(read_pairs[0]).split('.')[0].split('_')[0]
    linked_reads = []
    
    for k in ks:
        k_dir = '%s/%s/k%d' % (assembly_outdir, gene, k)
        k_dirs.append(k_dir)
        if not os.path.exists(k_dir):
            os.makedirs(k_dir)
        
        for i in range(2):
            source = os.path.relpath(read_pairs[i], k_dir)
            target = '%s/%s_%d.fastq' % (k_dir, gene, i + 1)
            linked_reads.append(target)
        
            if not os.path.exists(target):
                os.symlink(source, target)
                
@transform(symlink_assembly_input,
           formatter(".+/(.+)/(k\d+)"),
           assembly_outdir + "/{1[0]}/{2[0]}/{1[0]}-final.fa",
           logger, logging_mutex)
def assemble_single_gene(k_dir, contigs_file, logger, logging_mutex):
    gene, k = filter(None, k_dir.split(os.sep))[-2:]
    
    cmd = 'transabyss --kmer %s --pe %s %s --outdir %s --name %s --cleanup 3' % (k.lstrip('k'),
                                                                                 '%s/%s_1.fastq' % (k_dir, gene),
                                                                                 '%s/%s_2.fastq' % (k_dir, gene),
                                                                                 k_dir,
                                                                                 gene)
    stdout_str, stderr_str = run_cmd(cmd, force=True)

    if stderr_str and\
       re.search('error', stderr_str, re.IGNORECASE):
        run_cmd('touch %s/%s-final.fa %s/%s-FINAL.COMPLETE' % (k_dir, gene, k_dir, gene))
        with logging_mutex:
            logger.info(stderr_str)
    
@collate(assemble_single_gene,
         formatter(".+/k\d+/(.+)-final.fa$"),
         assembly_outdir + "/{1[0]}/{1[0]}-merged.fa",
         args.readlen)
def merge_assemblies(k_assemblies, merged_fasta, readlen):
    ks = []
    for k_assembly in k_assemblies:
        gene, k = filter(None, k_assembly.split(os.sep))[-3:-1]
        ks.append(int(k.lstrip('k')))
    
    prefixes = ' '.join(['%s.k%d.' % (gene, k) for k in ks])
    
    merged_fasta = '/'.join(filter(None, k_assemblies[0].split(os.sep))[:-2]) + '/%s-merged.fa' % gene
    
    # check if we have all empty assemblies
    num_empty_assemblies = 0
    for assembly in k_assemblies:
        if os.path.getsize(assembly) == 0:
            num_empty_assemblies += 1
    
    if num_empty_assemblies < len(k_assemblies):
        cmd = 'transabyss-merge --mink %d --maxk %d --prefixes %s --length %d %s --out %s --force' % (min(ks),
                                                                                                      max(ks),
                                                                                                      prefixes,
                                                                                                      readlen,
                                                                                                      ' '.join(k_assemblies),
                                                                                                      merged_fasta)
    else:
        cmd = 'touch %s' % merged_fasta

    run_cmd(cmd)
    
@transform(merge_assemblies,
           formatter(".+-merged.fa"),
           "{0[0]}.bwt")
def r2c_bwa_index(merged_fasta, index):
    if os.path.getsize(merged_fasta) > 0:
        cmd = 'bwa index %s' % merged_fasta

    else:
        cmd = 'touch %s' % index
        
    run_cmd(cmd)
    
@transform(r2c_bwa_index,
           formatter(),
           "{path[0]}/r2c.bam",
           args.nprocs,
           args.sort_mem)
def r2c(index, r2c_bam, nthreads, sort_mem):
    gene = filter(None, index.split(os.sep))[-2]
    reads1 = '%s/%s_1.fastq' % (bbt_outdir, gene)
    reads2 = '%s/%s_2.fastq' % (bbt_outdir, gene)
    
    if os.path.getsize(index) > 0:
        cmd = 'bwa mem -t %d %s %s %s | samtools view -uhS - | samtools sort -m %s - %s' % (nthreads,
                                                                                            os.path.splitext(index)[0],
                                                                                            reads1,
                                                                                            reads2,
                                                                                            sort_mem,
                                                                                            os.path.splitext(r2c_bam)[0])
        run_cmd('/bin/bash -c "%s"' % cmd)
        
    else:
        cmd = 'touch %s' % r2c_bam
        
    run_cmd('/bin/bash -c "%s"' % cmd)
    
@merge(merge_assemblies,
       '%s/%s.fa' % (assembly_outdir, args.sample))
def concat_fasta(gene_fastas, single_merged_fasta):
    fin = fileinput.input(gene_fastas)
    with open(single_merged_fasta, 'w') as fout:
        for line in fin:
            fout.write(line)
    fin.close()
    
@merge(r2c,
       '%s/r2c_cat.bam' % assembly_outdir,
       args.sort_mem)
def r2c_concat(r2c_bams, r2c_cat_bam, sort_mem):
    header_file = '%s/r2c_cat.header' % assembly_outdir
    sam_file = '%s/r2c_cat.sam' % assembly_outdir
    with open(header_file, 'w') as header, open(sam_file, 'w') as sam:
        for i in range(len(r2c_bams)):
            if os.path.getsize(r2c_bams[i]) == 0:
                    continue
            cmd = subprocess.Popen('samtools view -H %s' % r2c_bams[i], shell=True, stdout=subprocess.PIPE)
            for line in cmd.stdout:
                if i == 0 and line[:3] == '@HD':
                    header.write(line)
                elif line[:3] == '@SQ':
                    header.write(line)

            cmd = subprocess.Popen('samtools view %s' % r2c_bams[i], shell=True, stdout=subprocess.PIPE)
            sam.writelines(cmd.stdout)

    cmd = 'cat %s %s | samtools view -Su - | samtools sort -m %s - %s' % (header_file,
                                                                          sam_file,
                                                                          sort_mem,
                                                                          os.path.splitext(r2c_cat_bam)[0])
    
    run_cmd('/bin/bash -c "%s"' % cmd)

def r2c_cleanup():
    temp_files = ['%s/r2c_cat.header' % assembly_outdir, '%s/r2c_cat.sam' % assembly_outdir]
    for ff in glob.glob('%s/*/*-merged.fa.*' % assembly_outdir):
        temp_files.append(ff)
    for ff in glob.glob('%s/*/r2c.bam' % assembly_outdir):
        temp_files.append(ff)

    for temp_file in temp_files:
        if os.path.exists(temp_file):
            os.remove(temp_file)

@posttask(r2c_cleanup)
@transform(r2c_concat,
           suffix('.bam'),
           '.bam.bai')
def r2c_index_concat(r2c_cat_sorted_bam, r2c_cat_sorted_bam_index):
    cmd = 'samtools index %s' % r2c_cat_sorted_bam
    run_cmd(cmd)
           
@transform(concat_fasta,
           formatter(".+.fa$"),
           "{path[0]}/c2g.bam",
           args.gmap_index,
           args.nprocs)
def c2g(contigs_fasta, c2g_bam, gmap_index, nthreads):
    if os.path.getsize(contigs_fasta) > 0:
        cmd = 'gmap -D %s -d %s %s -t %d -f samse -n 0 -x 10 | samtools view -bhS - -o %s' % (gmap_index[0],
                                                                                              gmap_index[1],
                                                                                              contigs_fasta,
                                                                                              nthreads,
                                                                                              c2g_bam)
        run_cmd('/bin/bash -c "%s"' % cmd)
    else:
        run_cmd('touch %s' % c2g_bam)
    
@transform(concat_fasta,
           formatter(".+.fa$"),
           "{path[0]}/c2t.bam",
           args.bwa_index,
           args.nprocs)
def c2t(contigs_fasta, c2t_bam, bwa_index, nthreads):
    if os.path.getsize(contigs_fasta) > 0:
        cmd = 'bwa mem -t %d %s %s | samtools view -bhS - -o %s' % (nthreads,
                                                                    bwa_index,
                                                                    contigs_fasta,
                                                                    c2t_bam)
        run_cmd('/bin/bash -c "%s"' % cmd)
    else:
        run_cmd('touch %s' % c2t_bam)
    
@follows(mkdir(pvt_outdir))
@merge([concat_fasta, c2g, c2t, r2c_index_concat],
       pvt_outdir + '/events.bedpe',
       args.bwa_index,
       args.gmap_index,
       args.gtf,
       args.genome_fasta)
def find_events(inputs, events_output, transcripts_index, gmap_index, gtf, genome_fasta):
    merged_fasta, c2g_bam, c2t_bam, r2c_index = inputs
    if os.path.getsize(merged_fasta) > 0:
        cmd = 'find_events.py --gbam %s --tbam %s --transcripts_fasta %s --genome_index %s --r2c %s %s %s %s %s' % (c2g_bam,
                                                                                                                    c2t_bam,
                                                                                                                    transcripts_index,
                                                                                                                    ' '.join(gmap_index),
                                                                                                                    os.path.splitext(r2c_index)[0],
                                                                                                                    merged_fasta,
                                                                                                                    gtf,
                                                                                                                    genome_fasta,
                                                                                                                    os.path.dirname(events_output)
                                                                                                                    )
        run_cmd(cmd)
    else:
        run_cmd('touch %s' % events_output)

pipeline_printout(sys.stdout, verbose=3)
pipeline_run(verbose=3, multiprocess=args.nprocs, logger=logger)

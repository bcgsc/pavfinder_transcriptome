#!/usr/bin/env python

from ruffus import *
import ruffus.cmdline as cmdline
import subprocess
from collections import defaultdict
import gzip
import os
#from ruffus.combinatorics import *
import sys

def run_cmd(cmd):
    process = subprocess.Popen(cmd,
                               stdout = subprocess.PIPE,
                               stderr = subprocess.PIPE,
                               shell = True)
    stdout_str, stderr_str = process.communicate()
    if process.returncode != 0:
        raise Exception("Failed to run '%s'\n%s%sNon-zero exit status %s" %
                        (cmd, stdout_str, stderr_str, process.returncode))
    
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
parser.add_argument('--fq_list', type=str, nargs='+', help='text file of input fastq paths')
parser.add_argument('--nprocs', type=int, default=64, help='number of threads/processes')
parser.add_argument('--k', type=int, nargs='+', help='k sizes for assembly')
parser.add_argument('--readlen', type=int, help='read length')

args = parser.parse_args()
logger, logger_mutex = cmdline.setup_logging (__name__,
                                              args.log_file,
                                              args.verbose)

cmdline.run(args)

read_pairs = []
if args.fq:
    read_pairs = format_read_pairs(fqs=args.fq)
elif args.fq_list:
    read_pairs = format_read_pairs(list_file=args.fq_list)
#print 'www', read_pairs

bbt_outdir = '%s/bbt_v3.0.0b' % args.outdir
assembly_outdir = '%s/transabyss_v1.5.4' % args.outdir
pvt_outdir = '%s/pvt_v0.3.0' % args.outdir
if not os.path.exists(bbt_outdir):
    os.mkdir(bbt_outdir)

bbt_prefix = bbt_outdir + '/' + args.sample

assembly_input = []

#@follows('check_dir', mkdir(args.outdir))
@transform([read_pairs],
           formatter('_1.fastq.gz$'),
           bbt_outdir + '/' + args.sample + '_reads.fq',
           bbt_prefix,
           args.bf,
           args.nprocs)
def classify(paired_fqs, bbt_output, prefix, bf, nthreads):  
    #bbt_output = ['%s_reads.fq' % prefix, '%s_summary.tsv' % prefix]
    
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
    #with open(bbt_output[0], 'w') as out:
        #out.write('aa\n')
        
    print cmd
    run_cmd('/bin/bash -c "%s"' % cmd)
    
    print 'eee', bbt_output
    
@split(classify,
       bbt_outdir + '/*.fastq.gz')
def split_input(bbt_fastq, split_fastqs, genes=None):
    print 'zzzzz', bbt_fastq
    seqs1 = defaultdict(list)
    seqs2 = defaultdict(list)
    
    count = 0
    with open(bbt_fastq[0], 'r') as ff:
        for line in ff:
            if line[0] == '@' and count % 4 == 0:
                if len(line.split(' ')) > 1:
                    target = line.split(' ')[1].split('.')[0]
                    if genes is not None and not target in genes:
                        target = None
                    if count == 0 or count == 8:
                        seqs = seqs1
                        if count == 8:
                            count = 0
                    elif count == 4:
                        seqs = seqs2

            if target is not None:
                seqs[target].append(line)

            count += 1

    split_fastqs = output_split_pairs(seqs1, seqs2, bbt_outdir)
    
def output_split_pairs(seqs1, seqs2, outdir):
    fqs = []
    for target in seqs1.keys():
        fq1 = '%s/%s_1.fastq.gz' % (outdir, target)
        fq2 = '%s/%s_2.fastq.gz' % (outdir, target)
        with gzip.open(fq1, 'wb') as out1:
            for line in seqs1[target]:
                out1.write(line)
                
        with gzip.open(fq2, 'wb') as out2:
            for line in seqs2[target]:
                out2.write(line)
               
        fqs.append(fq1)
        fqs.append(fq2)
        
    return fqs

@collate(split_input,
         formatter(".+/(.+)_[12].fastq.gz$"),
         ["{path[0]}/{1[0]}_1.fastq.gz",
          "{path[0]}/{1[0]}_2.fastq.gz"])
def pair_reads(fastqs, pair):
    print 'ggg', fastqs
    pass

@subdivide(pair_reads,
           formatter(".+/(.+)_1.fastq.gz"),
           assembly_outdir + "/{1[0]}/k*",
           #regex(r".+/(.+)_1.fastq.gz"),
           #assembly_outdir + "/*/k*",
           args.k)
           #[assembly_outdir + "/" + r"\1" + "/k" + str(args.k[0]),
            #assembly_outdir + "/" + r"\1" + "/k" + str(args.k[1])])
def symlink_assembly_input(read_pairs, k_dirs, ks):  
    #print 'yyy', read_pairs
    gene = os.path.basename(read_pairs[0]).split('.')[0].split('_')[0]
    linked_reads = []
    #print 'ccc', read_pairs, ks, gene
    
    for k in ks:
        k_dir = '%s/%s/k%d' % (assembly_outdir, gene, k)
        k_dirs.append(k_dir)
        #print '777', k_dir
        if not os.path.exists(k_dir):
            os.makedirs(k_dir)
        
        for i in range(2):
            source = os.path.relpath(read_pairs[i], k_dir)
            target = '%s/%s_%d.fastq.gz' % (k_dir, gene, i + 1)
            linked_reads.append(target)
            #print source, target
        
            if not os.path.exists(target):
                os.symlink(source, target)
                
    for d in k_dirs:
        print 'ggg', d
                
@transform(symlink_assembly_input,
           formatter(".+/(.+)/k\d+"),
           "{0[0]}/{1[0]}-final.fa")
def assemble_single_gene(k_dir, contigs_file):
    gene, k = filter(None, k_dir.split(os.sep))[-2:]
    print 'ccc', gene, k
    
    cmd = 'transabyss --kmer %s --pe %s %s --outdir %s --name %s --cleanup 3' % (k.lstrip('k'),
                                                                                 '%s/%s_1.fastq.gz' % (k_dir, gene),
                                                                                 '%s/%s_2.fastq.gz' % (k_dir, gene),
                                                                                 k_dir,
                                                                                 gene)
    print cmd
    run_cmd(cmd)    

pipeline_printout(sys.stdout, verbose=3)
pipeline_run(verbose=3, multiprocess = 15)
#pipeline_run([assemble], verbose=2, multiprocess = 5)
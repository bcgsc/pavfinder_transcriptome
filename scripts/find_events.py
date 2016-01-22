import argparse
import os
import pysam
from sets import Set
from shared.transcript import Transcript
from event_finder import EventFinder
from exon_mapper import ExonMapper
from shared.adjacency import Adjacency
from read_support import find_support
from shared.translate import check_frame

def combine_events(events, mappings):
    """Combine events via genome and transcripts alignment on contig level"""
    def same_event(event_t, event_g, window=100):
        if event_t.rearrangement == event_g.rearrangement or\
           event_t.event == event_g.event:
            genome_breaks_t = sorted(event_t.genome_breaks)
            genome_breaks_g = sorted(event_g.genome_breaks)
            seq_breaks_t = sorted(event_t.seq_breaks)
            seq_breaks_g = sorted(event_g.seq_breaks)
            
            if abs(seq_breaks_t[0] - seq_breaks_g[0]) <= window and\
               abs(seq_breaks_t[1] - seq_breaks_g[1]) <= window:
                if event_g.exon_bounds and\
                   event_g.exon_bounds[0] and\
                   event_g.exon_bounds[1]:
                    return event_g
                elif event_t.exon_bounds and\
                     event_t.exon_bounds[0] and\
                     event_t.exon_bounds[1]:
                    return event_t
                else:
                    return event_g
        return False
            
    def same_mapping(query):
        passed = False
        if mappings['via_transcripts'] and mappings['via_genome']:
            if mappings['via_transcripts'].has_key(query) and\
               mappings['via_genome'].has_key(query):
                if not mappings['via_transcripts'][query] or not mappings['via_transcripts'][query][0]:
                    print '%s: mapping disagreed transcripts:None genome:%s' % (query,
                                                                                mappings['via_genome'][query])
                elif not mappings['via_genome'][query] or not mappings['via_genome'][query][0]:
                    print '%s: mapping disagreed transcripts:%s genome:None' % (query,
                                                                                mappings['via_transcripts'][query])
                elif mappings['via_transcripts'][query][0] & mappings['via_genome'][query][0]:
                    passed = True
                else:
                    print '%s: mapping disagreed transcripts:%s genome:%s' % (query,
                                                                              mappings['via_transcripts'][query],
                                                                              mappings['via_genome'][query])
            elif not mappings['via_transcripts'].has_key(query):
                print '%s: mapping disagreed transcripts:None genome:%s' % (query,
                                                                            mappings['via_genome'][query])
            else:
                print '%s: mapping disagreed transcripts:%s genome:None' % (query,
                                                                            mappings['via_transcripts'][query])
                
        return passed

    combined_events = []
    for query in Set(events['via_genome'].keys()) | Set(events['via_transcripts'].keys()):
        if events['via_genome'] and events['via_transcripts'] and not same_mapping(query):
            continue
        
        if not events['via_transcripts'].has_key(query):
            print 'onlygenome', query
            for event in events['via_genome'][query]:
                #if event.rearrangement == 'fusion' or event.rearrangement == 'read_through':
                combined_events.append(event)
            #else:
                #combined_events.extend(events['via_genome'][query])
        elif not events['via_genome'].has_key(query):
            print 'onlytranscripts', query
            for event in events['via_transcripts'][query]:
                combined_events.append(event)
        else:
            events_t = list(events['via_transcripts'][query])
            events_g = list(events['via_genome'][query])
            used_t = Set()
            used_g = Set()
            for i in range(len(events_t)):
                if i in used_t:
                    continue
                for j in range(len(events_g)):
                    if j in used_g:
                        continue
                    event = same_event(events_t[i], events_g[j])
                    if event:
                        combined_events.append(event)
                    #if same_event(events_t[i], events_g[j]):
                        #combined_events.append(events_t[i])
                        used_t.add(i)
                        used_g.add(j)
                        
            for i in range(len(events_t)):
                if not i in used_t:
                    if not events_t[i].event in ('ins', 'del'):
                        combined_events.append(events_t[i])
            for i in range(len(events_g)):
                if not i in used_g:
                    if not events_g[i].event in ('ins', 'del'):
                        combined_events.append(events_g[i])
                        
    return combined_events
     
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("query_fasta", type=str, help="path to query sequences")
    parser.add_argument("gtf", type=str, help="gtf file")
    parser.add_argument("genome_fasta", type=str, help="genome_fasta")
    parser.add_argument("outdir", type=str, help="path to output directory")
    parser.add_argument("--debug", dest="debug", action="store_true", help="debug mode")
    parser.add_argument("--tbam", type=str, help="query to transcript bam")
    parser.add_argument("--gbam", type=str, help="query to genome bam")
    parser.add_argument("--transcripts_fasta", type=str, help="path to target sequences")
    parser.add_argument("--r2c", type=str, help="reads to genome bam")
    parser.add_argument("--nproc", type=int, help="number of processes. Default:4", default=4)
    parser.add_argument("--genome_index", type=str, help="genome index path and name", nargs=2)
    filtering = parser.add_argument_group('filtering')
    filtering.add_argument("--min_support", type=int, help="minimum read support. Default:4", default=4)
    filtering.add_argument("--min_indel_size", type=int, help="minimum indel size. Default:3", default=3)
    filtering.add_argument("--no_utr", action="store_true", help="don't report events in UTR")
    filtering.add_argument("--include_nonsense_fusion", action="store_true", help="include non-sense fusions")
    filtering.add_argument("--include_non_exon_bound_fusion", action="store_true", help="include non-exon-bound fusions")
    filtering.add_argument("--max_homol_len", type=int, help="maximum homology sequence length. Default:5", default=5)
    filtering.add_argument("--max_novel_len", type=int, help="maximum novel sequence length. Default:5", default=5)
    filtering.add_argument("--subseq_len", type=int, help="subsequence length for filtering. Default:50", default=50)
    args = parser.parse_args()
    return args

def create_pysam_bam(path):
    if path is not None and os.path.exists(path):
        return pysam.AlignmentFile(path)
    return None

def create_pysam_fasta(path):
    if path is not None and os.path.exists(path):
        return pysam.FastaFile(path)
    return None

def create_pysam_tabix(path):
    if path is not None and os.path.exists(path):
        return pysam.Tabixfile(path, parser=pysam.asGTF())
    return None

def main():
    args = parse_args()
        
    gbam = create_pysam_bam(args.gbam)
    tbam = create_pysam_bam(args.tbam)
    query_fasta = create_pysam_fasta(args.query_fasta)
    genome_fasta = create_pysam_fasta(args.genome_fasta)
    transcripts_fasta = create_pysam_fasta(args.transcripts_fasta)
    transcripts_dict = Transcript.extract_transcripts(args.gtf)
    annot_tabix = create_pysam_tabix(args.gtf)
            
    ef = EventFinder(genome_fasta, annot_tabix, transcripts_dict, args.outdir, debug=args.debug)
    events = {'via_genome': {}, 'via_transcripts': {}}
    mappings = {'via_genome': {}, 'via_transcripts': {}}
    gene_hits = None
    if gbam and annot_tabix:
        events['via_genome'], mappings['via_genome'] = ef.find_events(gbam,
                                                                      query_fasta,
                                                                      genome_fasta,
                                                                      'genome',
                                                                      min_indel_size=args.min_indel_size,
                                                                      no_utr=args.no_utr,
                                                                      max_homol_len=args.max_novel_len,
                                                                      max_novel_len=args.max_homol_len,
                                                                      only_sense_fusion=not args.include_nonsense_fusion,
                                                                      only_exon_bound_fusion=not args.include_non_exon_bound_fusion
                                                                      )
        
    if tbam:
        events['via_transcripts'], mappings['via_transcripts'] = ef.find_events(tbam,
                                                                                query_fasta,
                                                                                transcripts_fasta,
                                                                                'transcripts',
                                                                                external_mappings=mappings['via_genome'],
                                                                                min_indel_size=args.min_indel_size,
                                                                                no_utr=args.no_utr,
                                                                                no_indels=True,
                                                                                max_homol_len=args.max_novel_len,
                                                                                max_novel_len=args.max_homol_len,
                                                                                only_sense_fusion=not args.include_nonsense_fusion,
                                                                                only_exon_bound_fusion=not args.include_non_exon_bound_fusion
                                                                                )        

    events_combined = combine_events(events, mappings)
    events_merged = Adjacency.merge(events_combined)
    
    if events_merged and args.genome_index and len(args.genome_index) == 2:
        ef.filter_probes(events_merged, args.genome_index[0], args.genome_index[1], args.outdir, debug=args.debug)
        ef.filter_subseqs(events_merged, query_fasta, args.genome_index[0], args.genome_index[1], args.outdir,
                          subseq_len=args.subseq_len, debug=args.debug)

    if args.r2c:
        find_support(events_merged, args.r2c, args.query_fasta, num_procs=args.nproc, debug=args.debug)
        events_filtered = [event for event in events_merged if event.support >= args.min_support]
    else:
        events_filtered = events_merged
        
    ef.set_frame(events_filtered, query_fasta, genome_fasta)

    Adjacency.report_events(events_filtered, '%s/events.tsv' % args.outdir)

main()
    

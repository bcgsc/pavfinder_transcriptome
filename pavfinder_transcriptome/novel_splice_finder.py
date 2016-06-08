from sets import Set
from alignment import reverse_complement
from adjacency import Adjacency
from transcript import Transcript
from collections import OrderedDict
from copy import deepcopy

report_items = OrderedDict(
        [('ID', None),
         ('event', 'event'),
         ('chrom', 'chroms,0'),
         ('start', 'genome_breaks,0'),
         ('end', 'genome_breaks,1'),
         ('gene', None),
         ('transcript', None),
         ('size', 'size'),
         ('exon1', 'exons,0',),
         ('exon2', 'exons,1',),
         ('seq_id', 'seq_id'),
         ('seq_breaks', 'seq_breaks'),
         ('in_frame', 'in_frame'),
         ('probe', 'probe'),
         ('support_reads', 'support'),
         ]
    )

def filter_events(events, min_support):
    failed = Set()
    event_to_index = dict((events[i], i) for i in range(len(events)))
    for i in range(len(events)):
	if events[i].spanning < min_support:
	    failed.add(i)
	    if events[i].link:
		print 'failed link', events[i].seq_id, events[i].event
		if event_to_index.has_key(events[i].link):
		    failed.add(event_to_index[events[i].link])
		else:
		    print 'cannot find link', events[i].seq_id, events[i].event, events[i].link
    for i in sorted(list(failed), reverse=True):
	del events[i]
    
def extract_features(gtfs, feature_types=('exon', 'junction')):
    annotated_features = {}
    for feature_type in feature_types:
	annotated_features[feature_type] = Set()
    for gtf in gtfs:
	features = Transcript.extract_features(gtf)
	for feature_type in feature_types:
	    annotated_features[feature_type] = annotated_features[feature_type].union(features[feature_type])
    return annotated_features

def find_novel_junctions(matches, align, transcript, query_seq, ref_fasta, accessory_known_features=None):
    """Find novel junctions within a single gene/transcript
    
    Args:
	block_matches: (list) dictionaries where 
					  key=transcript name, 
					  value=[match1, match2, ...] where
						match1 = matches of each alignment block
							 i.e.
							 [(exon_id, '=='), (exon_id, '==')] 
	align: (Alignment) alignment object
	transcripts_dict: (dictionary) key=transcript_id value=Transcript object
	ref_fasta: (Pysam.Fastafile) Pysam handle to access reference sequence, for checking splice motif
    Returns:
	List of event (dictionary storing type of event, exon indices, genomic coordinates and size)
    """
    def event_to_adj(event):
	adj = Adjacency(align.query,
	                (align.target, align.target),
	                (event['seq_breaks'][0], event['seq_breaks'][1]),
	                (event['pos'][0], event['pos'][1]),
	                transcripts = (transcript, transcript),
	                chroms = (align.target, align.target),
	                genome_breaks = (event['pos'][0], event['pos'][1]),
	                event = event['event'],
	                size = event['size']
	                )
	if event.has_key('exons') and event['exons']:
	    if len(event['exons']) == 1:
		adj.exons = (transcript.exon_num(event['exons'][0]),
		             transcript.exon_num(event['exons'][0]))
	    else:
		adj.exons = (transcript.exon_num(event['exons'][0]),
		             transcript.exon_num(event['exons'][1]))
	if adj.event == 'retained_intron':
	    adj.set_probe(query_seq)
	else:
	    adj.set_probe(query_seq)
	return adj
    
    def split_event(event):
	event1 = deepcopy(event)
	event2 = deepcopy(event)
	if event['event'] == 'retained_intron':
	    event1['pos'] = (event['pos'][0], event['pos'][0] + 1)
	    event1['seq_breaks'] = (event['seq_breaks'][0], event['seq_breaks'][0] + 1)
	    event2['pos'] = (event['pos'][1] - 1, event['pos'][1])
	    event2['seq_breaks'] = (event['seq_breaks'][1] - 1, event['seq_breaks'][1])
	else:
	    event1['pos'] = (event['pos'][0] - 1, event['pos'][0])
	    event1['seq_breaks'] = (event['seq_breaks'][0] - 1, event['seq_breaks'][0])
	    event2['pos'] = (event['pos'][1], event['pos'][1] + 1)
	    event2['seq_breaks'] = (event['seq_breaks'][1], event['seq_breaks'][1] + 1)
	print 'qqqq', event1
	print 'qqqq', event2
	adj1 = event_to_adj(event1)
	adj2 = event_to_adj(event2)
	adj1.link = adj2
	adj2.link = adj1
	return (adj1, adj2)
    
    #print 'aaa', align.blocks
    #print 'aaa', align.query_blocks
    #print 'aaa', align.strand

    # find annotated junctions
    annotated = Set()
    # sort multiple exon matches for single exon by exon num
    [m.sort(key=lambda mm: int(mm[0])) for m in matches if m is not None]	
    for i in range(len(matches) - 1):
	j = i + 1
	    
	if matches[i] == None or matches[j] == None:
	    continue
	
	if (i, j) in annotated:
	    continue
	
	if is_junction_annotated(matches[i][-1], matches[j][0]):
	    annotated.add((i, j))
	    continue
	    
    known_exons = known_juncs = None
    if accessory_known_features is not None:
	known_exons = accessory_known_features['exon']
	known_juncs = accessory_known_features['junction']

    all_events = []
    for i in range(len(matches) - 1):
	j = i + 1
		
	if matches[i] is None and matches[j] is not None:
	    # special case where the first 2 blocks is the utr and there's an insertion separating the 2 blocks
	    if i == 0:
		events = classify_novel_junction(matches[i], 
	                                         matches[j][0], 
	                                         align.target, 
	                                         align.blocks[i:j+1], 
	                                         transcript,
	                                         ref_fasta,
	                                         known_juncs = known_juncs,
	                                         known_exons = known_exons,
	                                         )
		for e in events:
		    e['blocks'] = (i, j)
		    e['transcript'] = transcript.id
		    e['seq_breaks'] = [align.query_blocks[i][1], align.query_blocks[j][0]]
		#all_events.extend([event_to_adj(event) for event in events])
		all_events.extend(events)
	    continue
	
	# for retained intron, not a 'junction'
	if matches[i] is not None and len(matches[i]) > 1:
	    events = classify_novel_junction(matches[i], 
	                                     None,
	                                     align.target, 
	                                     align.blocks[i], 
	                                     transcript,
	                                     ref_fasta,
	                                     known_juncs = known_juncs,
	                                     known_exons = known_exons,
	                                     )
	    if events:
		for e in events:
		    e['blocks'] = (i, j)
		    e['transcript'] = transcript.id
		    e['seq_breaks'] = [align.query_blocks[i][0], align.query_blocks[i][1]]
		    print 'zzz4', e
		#all_events.extend([event_to_adj(event) for event in events])
		all_events.extend(events)
	# skip if junction is annotated
	if (i, j) in annotated:
	    continue
	
	# for novel exon
	if matches[j] == None and matches[i] is not None:
	    for k in range(j + 1, len(matches)):
		if matches[k] is not None and matches[k][0][0] - matches[i][-1][0] == 1:
		    j = k
		    break
	
	if matches[i] is not None and matches[j] is not None:
	    # matches (i and j) are the flanking matches, blocks are the middle "novel" blocks
	    events = classify_novel_junction(matches[i][-1], 
	                                     matches[j][0], 
	                                     align.target, 
	                                     align.blocks[i:j+1], 
	                                     transcript,
	                                     ref_fasta,
	                                     known_juncs = known_juncs,
	                                     known_exons = known_exons,
	                                     )
	    if events:
		for e in events:
		    e['blocks'] = range(i + 1, j)
		    e['seq_breaks'] = [align.query_blocks[i][1], align.query_blocks[j][0]]
		    e['transcript'] = transcript.id
		#all_events.extend([event_to_adj(event) for event in events])
		all_events.extend(events)	
    #for event in all_events:
	#print 'ggg', event.details()
	
    adjs = []
    for event in all_events:
	if event['event'] in ('novel_exon', 'retained_intron'):
	    if event['event'] == 'retained_intron':
		if align.strand == '+':
		    event['seq_breaks'][0] += event['seq_break_offsets'][0]
		    event['seq_breaks'][1] -= event['seq_break_offsets'][1]
		else:
		    event['seq_breaks'][0] -= event['seq_break_offsets'][0]
		    event['seq_breaks'][1] += event['seq_break_offsets'][1]
	    print 'qqq', event
	    adjs.extend(split_event(event))
	else:
	    adjs.append(event_to_adj(event))
	
    return adjs
	    	
def report(event, event_id=None):
    data = []
    for item, label in report_items.iteritems():
	value = 'na'
	if item == "ID":
	    if event_id is not None:
		value = event_id
		
	elif label is not None and ',' in label:
	    attr, index = label.split(',')
	    if hasattr(event, attr):
		values = getattr(event, attr)
		if (type(values) is tuple or type(values) is list) and\
	           len(values) > int(index):
		    value = values[int(index)]
			
	elif item == 'gene':
	    value = event.transcripts[0].gene
	    
	elif item == 'transcript':
	    value = event.transcripts[0].id
		    
	elif hasattr(event, label):
	    val = getattr(event, label)
	    if val is not None:
		value = val
		
	data.append(str(value))

    return '\t'.join(data)
    
def classify_novel_junction(match1, match2, chrom, blocks, transcript, ref_fasta, min_intron_size=20, 
                            known_juncs=None, known_exons=None, no_indel=True):
    """Classify given junction into different splicing events or indel
    
    Args:
	match1: (tuple or list) single tuple: exon_index, 2-char match result e.g. '==', '>=', etc
				list of 2 tuples for retained_intron (special case)
	match2: (tuple) exon_index, 2-char match result e.g. '==', '>=', etc
	chrom: (str) chromosome, for checking splice motif
	blocks: (list) list of list (block coordinates)
	transcript: (Transcript) object of transcript, for getting exon coordinates
	ref_fasta: (Pysam.Fastafile) Pysam handle to access reference sequence, for checking splice motif
    Returns:
	List of event (dictionary storing type of event, exon indices, genomic coordinates and size)
    """
    events = []
    # set default values for 'pos'
    if type(blocks[0]) is int:
	pos = (blocks[0], blocks[1])
	if known_exons and (chrom, pos[0], pos[1]) in known_exons:
	    return events
    else:
	pos = (blocks[0][1], blocks[1][0])
	if known_juncs and (chrom, pos[0], pos[1]) in known_juncs:
	    return events

    if match2 is None:
	if len(match1) == 2:
	    exons = [m[0] for m in match1]
	    print 'zzz', blocks, exons
	    if match1[0][1] == '=>' and\
	       match1[-1][1] == '<=' and\
	       len([(a, b) for a, b in zip(exons, exons[1:]) if b == a + 1]) == len(match1) - 1:
		size = transcript.exons[exons[1]][0] - transcript.exons[exons[0]][1] - 1
		if not known_exons or not (chrom, pos[0], pos[1]) in known_exons:
		    events.append({'event': 'retained_intron',
		                   'exons': exons,
		                   'pos': (transcript.exons[exons[0]][1], transcript.exons[exons[1]][0]),
		                   'seq_break_offsets': (transcript.exons[exons[0]][1] - pos[0],
		                                         pos[1] - transcript.exons[exons[1]][0]),
		                   'size':size})
	      
    # genomic deletion or novel_intron
    elif match1 is None and match2 is not None:
	if match2[0] == 0 and match2[1] == '<=':
	    if transcript.strand == '+':
		donor_start = blocks[0][1] + 1
		acceptor_start = blocks[1][0] - 2
	    else:
		donor_start = blocks[1][0] - 2
		acceptor_start = blocks[0][1] + 1
		
	    gap_size = blocks[1][0] - blocks[0][1] - 1
	    pos = (blocks[0][1], blocks[1][0])
	    event = None
	    if gap_size > 0:
		if gap_size > min_intron_size and\
	           check_splice_motif(chrom, donor_start, acceptor_start, transcript.strand, ref_fasta):
		    event = 'novel_intron'
		else:
		    event = 'del'
	    if event is not None:
		events.append({'event': event, 'exons': [match2[0]], 'pos':pos, 'size':gap_size})
	    
    else:
	if match2[0] > match1[0] + 1 and\
           '=' in match1[1] and\
           '=' in match2[1]:
	    size = 0
	    for e in range(match1[0] + 1, match2[0]):
		exon = transcript.exons[e]
		size += exon[1] - exon[0] + 1
	    events.append({'event': 'skipped_exon', 'exons': range(match1[0] + 1, match2[0]), 'pos':pos, 'size':size})
	   
	if match1[0] == match2[0] and\
           match1[1][1] == '<' and match2[1][0] == '>':		
	    if transcript.strand == '+':
		donor_start = blocks[0][1] + 1
		acceptor_start = blocks[1][0] - 2
	    else:
		donor_start = blocks[1][0] - 2
		acceptor_start = blocks[0][1] + 1
	    
	    gap_size = blocks[1][0] - blocks[0][1] - 1
	    pos = (blocks[0][1], blocks[1][0])
	    event = None
	    if gap_size > 0:
		if gap_size > min_intron_size and\
	           check_splice_motif(chrom, donor_start, acceptor_start, transcript.strand, ref_fasta):
		    event = 'novel_intron'
		else:
		    event = 'del'
	    elif gap_size == 0:
		event = 'ins'
	    
	    if event is not None:
		if event != 'ins':
		    events.append({'event': event, 'exons': [match1[0]], 'pos':pos, 'size':gap_size})
		else:
		    events.append({'event': event, 'exons': [match1[0]], 'pos':pos})
	    
	# novel donor and acceptor
	if match1[1][1] == '=' and match2[1][0] != '=':
	    size = abs(blocks[1][0] - transcript.exons[match2[0]][0])
	    if transcript.strand == '+':
		event = 'novel_acceptor'
		donor_start = blocks[0][1] + 1
		acceptor_start = blocks[1][0] - 2
	    else:
		event = 'novel_donor'
		donor_start = blocks[1][0] - 2
		acceptor_start = blocks[0][1] + 1
	    # check splice motif
	    if check_splice_motif(chrom, donor_start, acceptor_start, transcript.strand, ref_fasta):
		events.append({'event': event, 'exons': [match1[0], match2[0]], 'pos':pos, 'size':size})
		
	if match1[1][1] != '=' and match2[1][0] == '=':
	    size = abs(blocks[0][1] - transcript.exons[match1[0]][1])
	    if transcript.strand == '+':
		event = 'novel_donor'
		donor_start = blocks[0][1] + 1
		acceptor_start = blocks[1][0] - 2
	    else:
		event = 'novel_acceptor'
		donor_start = blocks[1][0] - 2
		acceptor_start = blocks[0][1] + 1
	    # check splice motif
	    if check_splice_motif(chrom, donor_start, acceptor_start, transcript.strand, ref_fasta):
		events.append({'event': event, 'exons': [match1[0], match2[0]], 'pos':pos, 'size':size})
		    
	if match2[0] == match1[0] + 1 and\
           match1[1][1] == '=' and\
           match2[1][0] == '=': 
	    pos = (blocks[1][0], blocks[-2][1])
	    if not known_exons or not (chrom, pos[0], pos[1]) in known_exons:
		size = blocks[-2][1] - blocks[1][0] + 1
		if transcript.strand == '+':
		    donor_start = pos[1] + 1
		    acceptor_start = pos[0] - 2
		else:
		    donor_start = pos[0] - 2
		    acceptor_start = pos[1] + 1
		# check splice motif
		if check_splice_motif(chrom, donor_start, acceptor_start, transcript.strand, ref_fasta):
		    events.append({'event': 'novel_exon', 'exons': [], 'pos':pos, 'size':size})
	    
    # set size to None for event that doesn't have size i.e. 'ins'
    for event in events:
	if not event.has_key('size'):
	    event['size'] = None
    
    if no_indel:
	return [e for e in events if e['event'] not in ('ins', 'del')]
    else:
	return events

def check_splice_motif(chrom, donor_start, acceptor_start, strand, ref_fasta):
    """Check if the 4-base splice motif of a novel junction is canonical (gtag)
	
    Right now only considers 'gtag' as canonical
    
    Args:
	chrom: (str) chromosome
	donor_start: (int) genomic position of first(smallest) base of donor site (1-based)
	acceptor_start: (int) genomic position of first(smallest) base of acceptor site (1-based)
	strand: (str) transcript strand '+' or '-'
	ref_fasta: (Pysam.Fastafile) Pysam handle to access reference sequence
    Returns:
	True if it's canonical, False otherwise
    """
    canonical_motifs = Set()
    canonical_motifs.add('gtag')
    
    donor_seq = ref_fasta.fetch(chrom, donor_start - 1, donor_start - 1 + 2)
    acceptor_seq = ref_fasta.fetch(chrom, acceptor_start - 1, acceptor_start - 1 + 2)
    
    if strand == '+':
	motif = donor_seq + acceptor_seq
    else:
	motif = reverse_complement(acceptor_seq + donor_seq)
    
    if motif.lower() in canonical_motifs:
	return True
    else:
	return False
    
def is_junction_annotated(match1, match2):
    """Checks if junction is in gene model
    
    Args:
	match1: (tuple) exon_index, 2-char match result e.g. '==', '>=', etc
	match2: (tuple) exon_index, 2-char match result e.g. '==', '>=', etc
    Returns:
	True or False
    """
    if match2[0] == match1[0] + 1 and\
       match1[1][1] == '=' and\
       match2[1][0] == '=':
	return True
    
    return False

#def annotate_ref_junctions(events, junction_depths, transcripts_dict):
    #"""Annotate 5' and 3' gene reference junction depth/coverage

    #Arguments:
        #events: (list) of events (Event)
        #juncton_depths: (dict) {chrom[(start, end)] = depth}
        #transcripts_dict: (dict) transcript ID to object mapping (for extracting exon coordinates)
    #"""
    #for event in events:
	#txt = transcripts_dict[event.transcripts[0]]
	#jn5 = None
	#jn3 = None
	#if event.rna_event == 'novel_donor':
	    #jn5 = txt.exon(event.exons[0])[1], txt.exon(event.exons[1])[0]

	#elif event.rna_event == 'novel_acceptor':
	    #jn5 = txt.exon(event.exons[0])[1], txt.exon(event.exons[1])[0]

	#elif event.rna_event == 'skipped_exon':
	    #exon5 = min(event.exons) - 1
	    #exon3 = max(event.exons) + 1
	    #exon5_coord = txt.exon(exon5)
	    #exon3_coord = txt.exon(exon3)
	    #if exon5_coord is not None:
		#exon_coord = txt.exon(min(event.exons))
		#if txt.strand == '+':
		    #jn5 = exon5_coord[1], exon_coord[0]
		#else:
		    #jn5 = exon_coord[1], exon5_coord[0]

	    #if exon3_coord is not None:
		#exon_coord = txt.exon(max(event.exons))
		#if txt.strand == '+':
		    #jn3 = exon_coord[1], exon3_coord[0]
		#else:
		    #jn3 = exon3_coord[1], exon_coord[0]

	#elif event.rna_event == 'retained_intron':
	    #jn5 = txt.exon(event.exons[0])[1], txt.exon(event.exons[1])[0]

	#elif event.rna_event == 'novel_exon':
	    #for i in range(len(txt.exons) - 1):
		#intron = txt.exons[i][1] + 1, txt.exons[i + 1][0] - 1
		#if event.breaks[0] > intron[0] and event.breaks[1] < intron[1]:
		    #jn5 = txt.exons[i][1], txt.exons[i + 1][0]
		    #break

	#elif event.rna_event == 'novel_intron':
	    #exon = event.exons[0]
	    #exon_coord = txt.exon(exon)
	    #exon5 = exon - 1
	    #exon3 = exon + 1
	    #exon5_coord = txt.exon(exon5)
	    #exon3_coord = txt.exon(exon3)
	    #if exon5_coord is not None:
		#if txt.strand == '+':
		    #jn5 = exon5_coord[1], exon_coord[0]
		#else:
		    #jn5 = exon_coord[1], exon5_coord[0]

	    #if exon3_coord is not None:
		#if txt.strand == '+':
		    #jn3 = exon_coord[1], exon3_coord[0]
		#else:
		    #jn3 = exon3_coord[1], exon_coord[0]

	#if jn5 is not None or jn3 is not None:
	    #if junction_depths.has_key(event.chroms[0]):
		#if jn5 is not None and junction_depths[event.chroms[0]].has_key(jn5):
		    #event.ref5_depth = junction_depths[event.chroms[0]][jn5]
		    #event.ref5_coord = '%s:%s-%s' % (event.chroms[0], jn5[0], jn5[1])
		#if jn3 is not None and junction_depths[event.chroms[0]].has_key(jn3):
		    #event.ref3_depth = junction_depths[event.chroms[0]][jn3]
		    #event.ref3_coord = '%s:%s-%s' % (event.chroms[0], jn3[0], jn3[1])

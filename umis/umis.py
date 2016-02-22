#!/usr/bin/env python

import itertools
import collections
import regex as re
import json
import gzip
import sys
import logging
import time
import multiprocessing
from functools import partial
import toolz as tz

import click

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def stream_fastq(file_handler):
    ''' Generator which gives all four lines if a fastq read as one string
    '''
    next_element = ''
    for i, line in enumerate(file_handler):
        next_element += line
        if i % 4 == 3:
            yield next_element
            next_element = ''

@click.command()
@click.argument('transform', required=True)
@click.argument('fastq1', required=True)
@click.argument('fastq2', default=None, required=False)
@click.option('--separate_cb', is_flag=True, help="Keep dual index barcodes separate.")
@click.option('--demuxed_cb', default=None)
@click.option('--dual_index', is_flag=True)
@click.option('--cores', default=1)
# @profile
def fastqtransform(transform, fastq1, fastq2, separate_cb, demuxed_cb, dual_index, cores):
    ''' Transform input reads to the tagcounts compatible read layout using
    regular expressions as defined in a transform file. Outputs new format to
    stdout.
    '''
    if dual_index and separate_cb:
        read_template = '{name}:CELL_{CB1}-{CB2}:UMI_{MB}\n{seq}\n+\n{qual}\n'
    else:
        read_template = '{name}:CELL_{CB}:UMI_{MB}\n{seq}\n+\n{qual}\n'

    transform = json.load(open(transform))
    read1_regex = re.compile(transform['read1'])
    read2_regex = re.compile(transform['read2']) if fastq2 else None

    fastq1_fh = open(fastq1)
    if fastq1.endswith('gz'):
        fastq1_fh = gzip.GzipFile(fileobj=fastq1_fh)

    fastq_file1 = stream_fastq(fastq1_fh)

    if fastq2:
        fastq2_fh = open(fastq2)
        if fastq2.endswith('gz'):
            fastq2_fh = gzip.GzipFile(fileobj=fastq2_fh)

        fastq_file2 = stream_fastq(fastq2_fh)

    else:
        fastq_file2 = itertools.cycle((None,))

    transform = partial(transformer, read1_regex=read1_regex,
                          read2_regex=read2_regex, paired=fastq2)
    p = multiprocessing.Pool(cores)

    chunks = tz.partition_all(10000, itertools.izip(fastq_file1, fastq_file2))
    bigchunks = tz.partition_all(cores, chunks)
    for bigchunk in bigchunks:
        for chunk in p.map(transform, list(bigchunk)):
            for read1_dict in chunk:
                if dual_index:
                    if not separate_cb:
                        read1_dict['CB'] = read1_dict['CB1'] + read1_dict['CB2']

                if demuxed_cb:
                    read1_dict['CB'] = demuxed_cb

                # Deal with spaces in read names
                read1_dict['name'] = read1_dict['name'].partition(' ')[0]
                sys.stdout.write(read_template.format(**read1_dict))

def transformer(chunk, read1_regex, read2_regex, paired):
    # Parse the reads with the regexes
    reads = []
    for read1, read2 in chunk:
        read1_match = read1_regex.search(read1)
        if not read1_match:
            continue

        read1_dict = read1_match.groupdict()

        if paired:
            read2_match = read2_regex.search(read2)
            if not read2_match:
                continue

            read2_dict = read2_match.groupdict()

        else:
            read2_dict = dict()

        read1_dict.update(read2_dict)

        # Output the restrutured read
        reads.append(read1_dict)
    return reads

@click.command()
@click.argument('sam')
@click.argument('out')
@click.option('--genemap', required=False, default=None)
@click.option('--output_evidence_table', default=None)
@click.option('--positional', default=False)
@click.option('--minevidence', required=False, default=1.0, type=float)
# @profile
def tagcount(sam, out, genemap, output_evidence_table, positional, minevidence):
    ''' Count up evidence for tagged molecules
    '''
    from pysam import AlignmentFile
    from cStringIO import StringIO
    import pandas as pd

    from utils import weigh_evidence

    logger.info('Reading optional files')

    gene_map = None
    if genemap:
        with open(genemap) as fh:
            gene_map = dict(p.strip().split() for p in fh)

    if positional:
        tuple_template = '{0},{1},{2},{3}'
    else:
        tuple_template = '{0},{1},{3}'

    parser_re = re.compile('.*:CELL_(?P<CB>.*):UMI_(?P<MB>.*)')

    logger.info('Tallying evidence')
    start_tally = time.time()

    evidence = collections.defaultdict(int)

    sam_file = AlignmentFile(sam, mode='r')
    track = sam_file.fetch(until_eof=True)
    for i, aln in enumerate(track):
        if aln.is_unmapped:
            continue

        match = parser_re.match(aln.qname)
        CB = match.group('CB')
        MB = match.group('MB')

        txid = sam_file.getrname(aln.reference_id)
        if gene_map:
            target_name = gene_map[txid]

        else:
            target_name = txid

        e_tuple = tuple_template.format(CB, target_name, aln.pos, MB)

        # Scale evidence by number of hits
        evidence[e_tuple] += weigh_evidence(aln.tags)

    tally_time = time.time() - start_tally
    logger.info('Tally done - {:.3}s, {:,} alns/min'.format(tally_time, int(60. * i / tally_time)))
    logger.info('Collapsing evidence')

    buf = StringIO()
    for key in evidence:
        line = '{},{}\n'.format(key, evidence[key])
        buf.write(line)

    buf.seek(0)
    evidence_table = pd.read_csv(buf)
    evidence_query = 'evidence >= %f' % minevidence
    if positional:
        evidence_table.columns=['cell', 'gene', 'umi', 'pos', 'evidence']
        collapsed = evidence_table.query(evidence_query).groupby(['cell', 'gene'])['umi', 'pos'].size()

    else:
        evidence_table.columns=['cell', 'gene', 'umi', 'evidence']
        collapsed = evidence_table.query(evidence_query).groupby(['cell', 'gene'])['umi'].size()

    expanded = collapsed.unstack().T

    if gene_map:
        # This Series is just for sorting the index
        genes = pd.Series(index=set(gene_map.values()))
        genes = genes.sort_index()
        # Now genes is assigned to a DataFrame
        genes = expanded.ix[genes.index]

    else:
        genes = expanded

    genes.replace(pd.np.nan, 0, inplace=True)

    logger.info('Output results')

    if output_evidence_table:
        import shutil
        buf.seek(0)
        with open(output_evidence_table, 'w') as etab_fh:
            shutil.copyfileobj(buf, etab_fh)

    genes.to_csv(out)


@click.command()
@click.argument('fastq', type=click.File('r'))
def cb_histogram(fastq):
    ''' Counts the number of reads for each cellular barcode

    Expects formatted fastq files.
    '''
    parser_re = re.compile('(.*):CELL_(?P<CB>.*):UMI_(.*)\\n(.*)\\n\\+\\n(.*)\\n')

    counter = collections.Counter()
    for read in stream_fastq(fastq):
        match = parser_re.search(read).groupdict()
        counter[match['CB']] += 1

    for bc, count in counter.most_common():
        sys.stdout.write('{}\t{}\n'.format(bc, count))

def cb_filterer(chunk, bc1, bc2):
    parser_re = re.compile('(.*):CELL_(?P<CB>.*):UMI_(.*)\\n(.*)\\n\\+\\n(.*)\\n')
    kept = []
    for read in chunk:
        match = parser_re.search(read).groupdict()
        cb1 = match['CB']
        if bc2:
            cb1, cb2 = cb1.split("-")
        if cb1 not in bc1:
            continue
        if bc2 and cb2 not in bc2:
            continue
        kept.append(read)
    return kept

@click.command()
@click.argument('fastq', type=click.File('r'))
@click.option('--bc1', type=click.File('r'))
@click.option('--bc2', type=click.File('r'), required=False)
@click.option('--cores', default=1)
def cb_filter(fastq, bc1, bc2, cores):
    ''' Filters reads with non-matching barcodes
    Expects formatted fastq files.
    '''

    bc1 = set(cb.strip() for cb in bc1)
    if bc2:
        bc2 = set(cb.strip() for cb in bc2)

    filter_cb = partial(cb_filterer, bc1=bc1, bc2=bc2)
    p = multiprocessing.Pool(cores)

    chunks = tz.partition_all(10000, stream_fastq(fastq))
    bigchunks = tz.partition_all(cores, chunks)
    for bigchunk in bigchunks:
        for chunk in p.map(filter_cb, list(bigchunk)):
            for read in chunk:
                sys.stdout.write(read)

@click.group()
def umis():
    pass

umis.add_command(fastqtransform)
umis.add_command(tagcount)
umis.add_command(cb_histogram)
umis.add_command(cb_filter)

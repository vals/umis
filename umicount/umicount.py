#!/usr/bin/env python

import sys
import itertools
import collections
import copy
import re
from re import findall
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import json
import gzip


def stream_fastq(file_handler):
    ''' Generator which gives all four lines if a fastq read as one string
    '''
    next_element = ''
    for i, line in enumerate(file_handler):
        next_element += line
        if i % 4 == 3:
            yield next_element
            next_element =''


def fastq_transform(args):
    ''' Transform input reads to the umicount compatible read layout using regular expressions
    as defined in a transform file. [To be described]
    '''
    read_template = '{name}:CELL_{CB}:UMI_{MB}\n{seq}\n+\n{qual}\n'

    transform = json.load(open(args.transform))
    read1_regex = re.compile(transform['read1'])
    read2_regex = re.compile(transform['read2']) if args.fastq2 else None

    fastq1_fh = open(args.fastq1)
    if args.fastq1.endswith('gz'):
        fastq1_fh = gzip.GzipFile(fileobj=fastq1_fh)

    fastq_file1 = stream_fastq(fastq1_fh)

    if args.fastq2:
        fastq2_fh = open(args.fastq2)
        if args.fastq2.endswith('gz'):
            fastq2_fh = gzip.GzipFile(fileobj=fastq2_fh)

        fastq_file2 = stream_fastq(fastq2_fh)

    else:
        fastq_file2 = itertools.cycle((None,))

    fastq_out = open(args.outfastq, "w")
    for read1, read2 in itertools.izip(fastq_file1, fastq_file2):
        # Parse the reads with the regexes
        read1_match = read1_regex.search(read1)
        if not read1_match:
            continue

        read1_dict = read1_match.groupdict()

        if args.fastq2:
            read2_match = read2_regex.search(read2)
            if not read2_match:
                continue

            read2_dict = read2_match.groupdict()

        else:
            read2_dict = dict()

        read1_dict.update(read2_dict)

        if args.demuxed_cb:
            read1_dict['CB'] = args.demuxed_cb

        # Output the restrutured read
        fastq_out.write(read_template.format(**read1_dict))

    fastq_out.close()

"""
Tools for making a UMI count table (genes by cells) in the single-merged SAM file
"""

def tag_count(args):
    ''' Count up evidence for tagged molecules
    '''
    from simplesam import Reader
    from cStringIO import StringIO
    import pandas as pd

    sam_file = Reader(open(args.sam))

    gene_map = None
    if args.geneMap:
        with open(args.geneMap) as fh:
            gene_map = dict(p.strip().split() for p in fh)

    parser_re = re.compile('(.*):CELL_(?P<CB>.*):UMI_(?P<MB>.*)')

    evidence = collections.defaultdict(int)

    for i, aln in enumerate(sam_file):
        if aln.mapped:
            match = parser_re.search(aln.qname).groupdict()
            CB = match['CB']
            MB = match['MB']

            if gene_map:
                target_name = gene_map[aln.rname]
            else:
                target_name = aln.rname

            if args.positional:
                e_tuple = (CB, target_name, aln.pos, MB)
            else:
                e_tuple = (CB, target_name, MB)
            
            # TODO: Parsing NH should be more robust.
            nh = float(aln._tags[-1].split('NH:i:')[-1])  # Number of hits per read
            evidence[e_tuple] += 1. / nh

    buf = StringIO()
    for key in evidence:
        line = ','.join(map(str, key)) + ',' + str(evidence[key]) + '\n'
        buf.write(line)

    buf.seek(0)
    evidence_table = pd.read_csv(buf)
    evidence_table.columns=['cell', 'gene', 'umi', 'evidence']

    # TODO: Make required amount of evidence for a tagged molecule a parameter
    # TODO: How to use positional information?
    collapsed = evidence_table.query('evidence > 1').groupby(['cell', 'gene'])['umi'].size()
    expanded = collapsed.unstack().T

    genes = pd.Series(index=set(gene_map.values()))
    genes = genes.sort_index()
    genes = expanded.ix[genes.index]
    genes.replace(pd.np.nan, 0, inplace=True)

    if args.evidence_table:
        import shutil
        buf.seek(0)
        with open(args.evidence_table, 'w') as etab_fh:
            shutil.copyfileobj(buf, etab_fh)

    genes.to_csv(args.out)


def sam_spike_count(sam_file, cell_barcodes, gene_cell_umi_sets, gene_umi_sets, minaqual, umilen):
    for aln in sam_file:
        if aln.aligned and aln.aQual >= minaqual:
            umi = extract_umi(aln.read.name)
            if umilen is not None:
                umi = umi[:umilen]
            cell = extract_cellbarcode(aln.read.name)
            if cell in cell_barcodes:
                if cell not in gene_cell_umi_sets:
                    gene_cell_umi_sets[cell] = copy.deepcopy(gene_umi_sets)
                gene_cell_umi_sets[cell][aln.iv.chrom].add(umi)

    gene_cell_counts = collections.defaultdict(dict)
    for cell in gene_cell_umi_sets:
        for gene in gene_cell_umi_sets[cell]:
            gene_cell_counts[cell][gene] = len(gene_cell_umi_sets[cell][gene])
    
    return gene_cell_counts


def main():
    parser = ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(help="subcommad help")

    subparser_fastqtransform = subparsers.add_parser("fastqtransform", description="Reformat fastq reads to umicount compatible format",
                                                                  formatter_class=ArgumentDefaultsHelpFormatter,
                                                                  help="trim cell and molecular barcodes and incorporate them into read name")
    subparser_fastqtransform.add_argument("--fastq1", metavar="FASTQ1", help="input FASTQ file 1", required=True)
    subparser_fastqtransform.add_argument("--fastq2", metavar="FASTQ1", help="input FASTQ file 2 for paired-end reads", required=False)
    subparser_fastqtransform.add_argument("--transform", metavar="TRANSFORM", help="FASTQ Transform JSON file", required=True)
    subparser_fastqtransform.add_argument("--outfastq", metavar="FASTQOUT", help="output FASTQ file for FASTQ1", required=True)
    subparser_fastqtransform.add_argument("--demuxed_cb", metavar="DEMUXED_CB", help="Set CB value to this in the transformed read name. Use this if your files have already been demultiplexed (e.g. STRT-Seq).", required=False)
    subparser_fastqtransform.set_defaults(func=fastq_transform)

    subparser_tagcount = subparsers.add_parser("tagcount", description="Count tag evidence from the SAM file",
                                                        formatter_class=ArgumentDefaultsHelpFormatter,
                                                        help="count reads from the SAM file")
    subparser_tagcount.add_argument("--sam", metavar="SAM", help="SAM file", required=True)
    subparser_tagcount.add_argument("--geneMap", "-g", metavar="GENEMAP",
                                                       help="Mapping of transcripts to genes", required=False)
    subparser_tagcount.add_argument("--positional", help="Consider position in transcript as molecular evidence",
                                                    required=False, action='store_true')
    subparser_tagcount.add_argument("--out", metavar="OUT", help="Output file", required=True)
    subparser_tagcount.add_argument("--evidence_table", metavar="ETAB", help="Save evidence table", required=False)
    subparser_tagcount.set_defaults(func=tag_count)

    
    args = parser.parse_args()
    args.func(args)

if __name__ == "__main__":
    main()

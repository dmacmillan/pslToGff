import argparse
import os
import sys
import logging
from itertools import groupby
import time

def parse_psl(psl, feature, gff, gff_cols, exons={}):
    stime = time.time()
    # Keep track of Qnames so as to avoid looking at
    # secondary alignments
    to_analyze = {}
    with open(psl, 'r') as f:
        logging.debug('Reading psl: {}'.format(psl))
        # Read 5 lines to skip the header
        for i in range(5):
            f.readline()
        for line in f:
            #logging.debug('parsing line: {}'.format(line))
            cols = line.strip().split('\t')
            match = int(cols[0])
            cols[0] = match
            mrna_name = cols[9]
            blocksizes = [int(x) for x in cols[18].strip(',').split(',')]
            nblocks = len(blocksizes)
            if mrna_name in to_analyze:
                if match > to_analyze[mrna_name][0]:
                    to_analyze[mrna_name] = cols
            else:
                to_analyze[mrna_name] = cols
        for mrna_name in to_analyze:
            cols = to_analyze[mrna_name]
            # psl uses a 0-based coordinate system
            # while gff uses a 1-based
            tstarts = [int(x)+1 for x in cols[20].strip(',').split(',')]
            blocksizes = [int(x) for x in cols[18].strip(',').split(',')]
            logging.debug('mrna_name: {}'.format(mrna_name))
            logging.debug('cols: {}'.format(cols))
            logging.debug('tstarts: {}'.format(tstarts))
            logging.debug('blocksizes: {}'.format(blocksizes))
            strand = cols[8]
            if exons or exons == {}:
                # Get the gene name
                #gene_name = cols[9][:cols[9].rindex('-')]
                gene_name = mrna_name.rsplit('-',2)[0]
                if gene_name not in exons:
                    exons[gene_name] = {
                        mrna_name: {
                            'start': tstarts[0],
                            'end': tstarts[-1] + blocksizes[-1],
                            'strand': strand
                        }
                    }
                elif mrna_name not in exons[gene_name]:
                    exons[gene_name][mrna_name] = {
                        'start': tstarts[0],
                        'end': tstarts[-1] + blocksizes[-1],
                        'strand': strand
                    }
                else:
                    exons[gene_name][mrna_name]['start'] = min(
                        exons[gene_name][mrna_name]['start'],
                        tstarts[0]
                    )
                    exons[gene_name][mrna_name]['end'] = max(
                        exons[gene_name][mrna_name]['end'],
                        tstarts[-1] + blocksizes[-1]
                    )
            # Check number of blocks is equal
            if mrna_name not in gff:
                logging.warning(
                    '"{}" not found in reference gff! ' \
                    'Cannot check block sizes!'.format(mrna_name)
                )
            else:
                len_blocksizes = len(blocksizes)
                len_gff = len(gff[mrna_name][feature])
                gff_sizes = set([x[1]-x[0] for x in gff[mrna_name][feature]])
                set_blocksizes = set([x-1 for x in blocksizes])
    #            elif len_blocksizes != len_gff:
    #                gff_sizes = set([x[1]-x[0] for x in gff[mrna_name][feature]])
    #                set_blocksizes = set(blocksizes)
    #                logging.warning(
    #                    'Number of blocks not equal "{}"'.format(mrna_name)
    #                )
                if len_blocksizes < len_gff:
                    logging.warning(
                        '"{}" GFF contains blocksizes {} not in PSL'.format(
                            mrna_name,
                            (', ').join([str(x) for x in gff_sizes - set_blocksizes])
                        )
                    )
                elif len_gff < len_blocksizes:
                    logging.warning(
                        '"{}" GFF is missing blocksizes {} in PSL'.format(
                            mrna_name,
                            (', ').join([str(x) for x in set_blocksizes - gff_sizes])
                        )
                    )
#                    'Number of blocks not equal "{}". ' \
#                    'Ref gff = "{}", psl blocks = "{}"\n' \
#                    '{}'.format(
#                        mrna_name,
#                        gff[mrna_name][feature],
#                        len(blocksizes),
#                        cols
#                    )
            for i,start in enumerate(tstarts):
                out = {key:None for key in gff_cols}
                out['source'] = 'pslToGff'
                out['feature'] = feature
                out['start'] = str(tstarts[i])
                logging.debug('blocksizes: {}'.format(blocksizes))
                logging.debug('tstarts: {}'.format(tstarts))
                logging.debug('i: {}'.format(i))
                out['end'] = str(tstarts[i] + blocksizes[i])
                out['score'] = '0'
                # These boundaries can be fixed in this case
                out['seqname'] = mrna_name[6:26]
                out['strand'] = strand
                out['frame'] = '0'
                attributes = {
                    'ID': '{}:{}:{}'.format(mrna_name, feature, i),
                    'Parent': mrna_name
                }
                out['attribute'] = ('; ').join(
                    [
                        '{}={}'.format(k,v) for k,v in attributes.items()
                    ]
                )
                print(('\t').join(
                    [out[key] for key in gff_cols]
                ))
    logging.info('Parsed "{}" in {}s'.format(psl, time.time()-stime))
    return exons

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Given psl files, return a gff file',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument('gff', help='reference gff file with the "ID" attribute set to the Qname of the psl')
    parser.add_argument('-ep', '--exon_psls', nargs='+', help='exon psl file input')
    parser.add_argument('-cp', '--cds_psls', nargs='+', help='CDS psl file input')
    #parser.add_argument('-f', '--feature', default='exon', help='A feature name for the gff output, must match the feature names from the gff input')
    parser.add_argument("-l", "--log", dest="logLevel", default='warning', choices=['debug', 'info', 'warning', 'error', 'critical'], help='set the logging level. default = "warning"')
    #parser.add_argument('-n', '--name', help='Name for file output')
    #parser.add_argument('-o', '--outdir', default=os.getcwd(), help='Path to output to. Will be created if it does not exist.')
    
    args = parser.parse_args()
    
    # Create output directory if it does not exist
#    if not os.path.isdir(args.outdir):
#        try:
#            os.makedirs(args.outdir)
#        except OSError:
#            logging.warning('Attempted to create {}, but directory may already exist'.format(args.outdir))
    
    # Logging
    logging.basicConfig(
        filename = os.path.join(os.getcwd(), '{}.log'.format(os.path.basename(__file__))),
        level = getattr(logging, args.logLevel.upper()),
        filemode = 'w'
    )

    # Keep track of gene sizes
    genes = {}

    # GFF format columns
    gff_cols = (
        'seqname',
        'source',
        'feature',
        'start',
        'end',
        'score',
        'strand',
        'frame',
        'attribute'
    )

    # Load gff into memory
    gff = {}

    with open(args.gff, 'r') as f:
        for line in f:
            if line[0] == '#':
                continue
            cols = line.strip().split('\t')
            attributes = cols[8]
            start = int(cols[3])
            end = int(cols[4])
            feature = cols[2]
            # Since this script is very specific I can 
            # hardcode the CDS and exon
            if feature not in ('CDS', 'exon'):
                continue
            attributes = dict([x.split('=') for x in attributes.split(';')])
            names = attributes['Parent'].split(',')
#            try:
#                name = attributes['ID'][:attributes['ID'].index(':')]
#            except ValueError:
#                name = attributes['ID']
            for name in names:
                if name not in gff:
                    gff[name] = {}
                if feature not in gff[name]:
                    gff[name][feature] = []
                gff[name][feature].append((start, end))

    # Keep track of exon boundaries
    exons = {}
    if args.exon_psls:
        for psl in args.exon_psls:
            exons = parse_psl(psl, 'exon', gff, gff_cols, exons=exons)
    if args.cds_psls:
        for psl in args.cds_psls:
            parse_psl(psl, 'CDS', gff, gff_cols)

    for gene in exons:
        out = {key:None for key in gff_cols}
        out['source'] = 'pslToGff'
        out['feature'] = 'gene'
        out['start'] = str(min([exons[gene][x]['start'] for x in exons[gene]]))
        out['end'] = str(max([exons[gene][x]['end'] for x in exons[gene]]))
        out['score'] = '0'
        # These boundaries can be fixed in this case
        out['seqname'] = gene[6:26]
        out['strand'] = exons[gene][list(exons[gene])[0]]['strand']
        out['frame'] = '0'
        out['attribute'] = 'ID={}'.format(gene)
        print(('\t').join(
            [out[key] for key in gff_cols]
        ))
        for mrna in exons[gene]:
            out['feature'] = 'mRNA'
            out['start'] = str(exons[gene][mrna]['start'])
            out['end'] = str(exons[gene][mrna]['end'])
            # These boundaries can be fixed in this case
            attributes = {
                'ID': mrna,
                'Parent': gene
            }
            out['attribute'] = ('; ').join(
                [
                    '{}={}'.format(k,v) for k,v in attributes.items()
                ]
            )
            print(('\t').join(
                [out[key] for key in gff_cols]
            ))

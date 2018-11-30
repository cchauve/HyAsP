#!/usr/bin/env python

# Filters BLAST output file (outfmt 6) by quality (match identity) and length (proportion of query covered).
# Optionally, fragmented matches (i.e. several short high-identity hits that together satisfy length threshold)
# can be detected.


import numpy as np
import os.path
import pandas as pd
import sys

from Bio import SeqIO


# default values / constants
DEF_IDENTITY_THRESHOLD = 0.95
DEF_LENGTH_THRESHOLD = 0.95
DEF_FIND_FRAGMENTED = False
DEF_VERBOSE = False


# read FASTA file and create mapping from id to sequence
def read_queries(queries_file):
    queries = dict()
    with open(queries_file, 'r') as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            queries[record.id] = record.seq

    return queries


# helper method for extracting rows from Pandas data frame associated to given query
def hits_of(df, query):
    return df.loc[df.qseqid == query]


# filters Pandas data frame according to length and identity threshold
# returns (1) reduced data frame, (2) data frame of low-identity hits, and
# (3) data frame of high-identity hits that cover a too small fraction of the query
def filter_data(data, queries, length_threshold, identity_threshold):
    identity_threshold *= 100

    filtered_out_identity = data.loc[data.pident < identity_threshold]
    data = data.loc[data.pident >= identity_threshold]

    align_lengths = np.array(data.length.tolist())
    query_lengths = np.array([len(queries[q]) for q in data.qseqid])
    too_short = align_lengths / query_lengths < length_threshold
    filtered_out_length = data.loc[too_short]
    data = data.loc[~too_short]

    return data, filtered_out_identity, filtered_out_length


# --- collection of methods for finding fragmented hits ---

# search for collection of hits (each with high identity) that together satisfy length threshold
# do not consider containment and overlaps
def find_fragmented_match_simple(data, query_id, query_len, length_threshold, identity_threshold):
    hits = hits_of(data, query_id)
    hits = hits.loc[hits.pident >= 100 * identity_threshold]

    # check total length of remaining fragments
    fragments = []
    if sum(hits.length) >= query_len * length_threshold:
        fragments = list(hits.index)

    return fragments


# search for collection of hits (each with high identity) that - after removing containments - together satisfy length threshold
# do not consider overlaps
def find_fragmented_match_uncontained(data, query_id, query_len, length_threshold, identity_threshold):
    hits = hits_of(data, query_id)
    hits = hits.loc[hits.pident >= 100 * identity_threshold]

    # remove containments
    contained_hits = []
    for index, row in hits.iterrows():
        if len(hits.loc[(hits.index != index) & (hits.qstart <= row.qstart) & (row.qend <= hits.qend)]) > 0:
            contained_hits.append(index)
    hits = hits.loc[[(c not in contained_hits) for c in hits.index]]

    # check total length of remaining fragments
    fragments = []
    if sum(hits.length) >= query_len * length_threshold:
        fragments = list(hits.index)

    return fragments


# search for collection of hits (each with high identity) with limited total overlap (and no containment) that together satisfy length threshold
def find_fragmented_match_limited_overlaps(data, query_id, query_len, length_threshold, identity_threshold, overlap_threshold = 0):
    hits = hits_of(data, query_id)
    hits = hits.loc[hits.pident >= 100 * identity_threshold]

    # remove containments
    contained_hits = []
    for index, row in hits.iterrows():
        if len(hits.loc[(hits.index != index) & (hits.qstart <= row.qstart) & (row.qend <= hits.qend)]) > 0:
            contained_hits.append(index)
    hits = hits.loc[[(c not in contained_hits) for c in hits.index]]

    ## check overlaps
    intervals = list(zip(hits.qstart, hits.qend))  # portion of query matched
    intervals.sort(key=lambda x: x[1])
    intervals.sort(key=lambda x: x[0]) # intervals is now sorted by first component (start) with second component (end) as tie-breaker, both ascending

    # first, check whether enough of the query is covered (positions covered by multiple hits contribute only once)
    align_len_no_overlaps = 0
    last_pos_covered = 0 # last (right-most) position of query covered so far
    for start, end in intervals:
        if end <= last_pos_covered:
            pass
        else:
            align_len_no_overlaps += end - max(last_pos_covered + 1, start) + 1
            last_pos_covered = end
    if align_len_no_overlaps < query_len * length_threshold:
        return []

    # second, check that the total length of pairwise overlaps of the hits is not too big
    total_length_pairwise_overlaps = 0
    for i in range(0, len(intervals)):
        a_start = intervals[i][0]
        a_end = intervals[i][1]
        for j in range(i + 1, len(intervals)):
            b_start = intervals[j][0]
            b_end = intervals[j][1]

            if b_end <= a_end:
                total_length_pairwise_overlaps += b_end - b_start + 1
            elif a_end >= b_start:
                total_length_pairwise_overlaps += a_end - b_start + 1
            else:
                pass
    if total_length_pairwise_overlaps > overlap_threshold:
        return []

    return list(hits.index)


# filters the BLAST output ('mapping') by match quality and length, prints some statistics to stdout and writes results to file(s)
def filter_mapping(queries, blast_file, filtered_blast_file, contigs_output = None, genes_output = None, length_threshold = DEF_LENGTH_THRESHOLD,
                   identity_threshold = DEF_IDENTITY_THRESHOLD, find_fragmented = DEF_FIND_FRAGMENTED, verbose = DEF_VERBOSE):
    col_names = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'] # outfmt 6
    data = pd.read_csv(blast_file, sep = '\t', names = col_names)

    # check consistency of the provided genes and mapping
    inconsistent_genes = set(data.qseqid.unique()).difference(queries)
    if len(inconsistent_genes) > 0:
        raise ValueError('The used gene-contig mapping and gene database are not consistent '
                         + '(%i genes are in the mapping but not in the database). ' % len(inconsistent_genes)
                         + 'Please check whether mapping %s was created using the database stated in the configuration.'
                         % blast_file)

    filtered_data, too_different, too_short = filter_data(data, queries, length_threshold, identity_threshold)

    located_genes = data.qseqid.unique()
    unlocated_genes = list(set(queries.keys()).difference(located_genes))
    genes_located_once = []
    genes_located_multi = []
    groups = data.groupby(['qseqid'])
    for gene, hits in groups:
        genes_located_once.append(gene) if len(hits) == 1 else genes_located_multi.append(gene)
    target_contigs = data.sseqid.unique()

    filtered_located_genes = filtered_data.qseqid.unique()
    filtered_unlocated_genes = list(set(queries.keys()).difference(filtered_located_genes))
    filtered_genes_located_once = []
    filtered_genes_located_multi = []
    groups = filtered_data.groupby(['qseqid'])
    for gene, hits in groups:
        filtered_genes_located_once.append(gene) if len(hits) == 1 else filtered_genes_located_multi.append(gene)
    filtered_target_contigs = filtered_data.sseqid.unique()

    # statistics
    if verbose:
        print('Search results before fragment search (and before filtering):')
        print('Number of queries:\t%i' % len(queries.keys()))
        print('Number of hits:\t%i (%i)' % (len(filtered_data.qseqid), len(data.qseqid)))
        print('Number of genes found:\t%i (%i)' % (len(filtered_located_genes), len(located_genes)))
        print('Number of genes w/ 0 hits:\t%i (%i)' % (len(filtered_unlocated_genes), len(unlocated_genes)))
        print('Number of genes w/ 1 hit:\t%i (%i)' % (len(filtered_genes_located_once), len(genes_located_once)))
        print('Number of genes w/ >1 hits:\t%i (%i)' % (len(filtered_genes_located_multi), len(genes_located_multi)))
        print('Number of gene-containing contigs:\t%i (%i)' % (len(filtered_target_contigs), len(target_contigs)))

        print('\nListing of multiple hits:')
        for gene in filtered_genes_located_multi:
            print(hits_of(filtered_data, gene).to_string())
            print()
        print('\n - End of listing - ')

    # optional search for fragmented hits
    if find_fragmented:
        genes_from_fragments = set() # genes (qseqid) of additional matches
        hits_from_fragments = set() # row numbers (index) of additional matches
        groups = data.groupby(['qseqid'])
        for gene in set(filtered_unlocated_genes).intersection(data.qseqid.unique()): # search only unlocated genes that have at least one hit
            fragments = find_fragmented_match_limited_overlaps(groups.get_group(gene), gene, len(queries[gene]), length_threshold,identity_threshold)
            if (fragments != []):
                genes_from_fragments.add(gene)
                hits_from_fragments.update(fragments)

        hits_set = set(filtered_data.index)
        prev_num_hits = len(hits_set)
        hits_set.update(hits_from_fragments)
        if verbose:
            print()
            print('Number of hits found by fragment search:\t%i (%i new)' % (len(hits_from_fragments), len(hits_set) - prev_num_hits))

        hits_from_fragments = data.loc[hits_from_fragments] # turn indices into full data frame rows

        genes_set = set(filtered_located_genes)
        prev_num_genes = len(genes_set)
        genes_set.update(genes_from_fragments)
        if verbose:
            print('Number of genes found by fragment search:\t%i (%i new)' % (len(genes_from_fragments), len(genes_set) - prev_num_genes))

        contigs_set = set(filtered_target_contigs)
        prev_num_contigs = len(contigs_set)
        contigs_set.update(hits_from_fragments.sseqid)
        if verbose:
            print('Number of contigs found by fragment search:\t%i (%i new)' % (len(hits_from_fragments.sseqid.unique()), len(contigs_set) - prev_num_contigs))

        filtered_data = data.loc[hits_set]
        filtered_located_genes = list(genes_set)
        filtered_target_contigs = list(contigs_set)

    # output generation
    filtered_data.to_csv(filtered_blast_file, sep = '\t', header = False, index = False)

    if contigs_output is not None:
        with open(contigs_output, 'w') as out:
            out.write('\n'.join(list(map(str, filtered_target_contigs))) + '\n')

    if genes_output is not None:
        with open(genes_output, 'w') as out:
            out.write('\n'.join(filtered_located_genes) + '\n')


if __name__ == '__main__':
    import argparse
    argparser = argparse.ArgumentParser()
    argparser.add_argument('queries_file', help = 'FASTA file of BLAST queries (genes)')
    argparser.add_argument('blast_file', help = 'BLAST output (default outfmt 6)')
    argparser.add_argument('filtered_blast_file', help = 'output file for remaining BLAST hits')
    argparser.add_argument('-c', '--contigs_output', default = None, help = 'list of matched contigs in filtered BLAST hits')
    argparser.add_argument('-g', '--genes_output', default = None, help = 'list of query genes in filtered BLAST hits')
    argparser.add_argument('-i', '--identity_threshold', type = float, default = DEF_IDENTITY_THRESHOLD, help = 'minimum identity in a hit to keep it')
    argparser.add_argument('-l', '--length_threshold', type = float, default = DEF_LENGTH_THRESHOLD, help = 'minimum fraction of query that has to be matched to keep a hit')
    argparser.add_argument('-f', '--find_fragmented', action = 'store_true', help = 'search for fragmented hits, i.e. several short high-identity hits that together satisfy the length threshold')
    argparser.add_argument('-v', '--verbose', action = 'store_true', help = 'print more information')
    args = argparser.parse_args()

    if not os.path.isfile(args.queries_file):
        print('ERROR: FASTA file of queries %s does not exist' % args.queries_file, file = sys.stderr)
    elif not os.path.isfile(args.blast_file):
        print('ERROR: BLAST output file %s does not exist.' % args.blast_file, file = sys.stderr)
    else:
        queries = read_queries(args.queries_file)
        try:
            filter_mapping(queries, args.blast_file, args.filtered_blast_file, contigs_output = args.contigs_output,
                           genes_output = args.genes_output, length_threshold = args.length_threshold,
                           identity_threshold = args.identity_threshold, find_fragmented = args.find_fragmented)
        except ValueError as err:
            print('ERROR: %s' % err)
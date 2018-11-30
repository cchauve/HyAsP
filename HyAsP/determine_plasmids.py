#!/usr/bin/env python

# Determines putative plasmids in assembly graph by
# (1) searching given plasmid genes in the contigs of the assembly,
# (2) filtering the resulting gene-contig matches by quality and length,
# (3) greedily contructing plasmids using the assembly graph and the gene-contig matches.
#
# Requirements:
#  - standard UNIX tools (mkdir)
#  # recursively
#  - makeblastdb (makeblastdb / --makeblastdb; tested with BLAST+ v2.6.0)
#  - blastn (blastn / --blastn; tested with BLAST+ v2.6.0)
#  - standard UNIX tools (rm)
#
# makeblastdb and blastn have to be in $PATH or explicitly specified using the respective path options.


import os

from subprocess import call

from HyAsP import map_genes as mg, find_plasmids as fp


# default values / constants
DEF_IDENTITY_THRESHOLD = mg.DEF_IDENTITY_THRESHOLD
DEF_LENGTH_THRESHOLD = mg.DEF_LENGTH_THRESHOLD
DEF_FIND_FRAGMENTED = mg.DEF_FIND_FRAGMENTED
DEF_MAKEBLASTDB_PATH = mg.DEF_MAKEBLASTDB_PATH
DEF_BLASTN_PATH = mg.DEF_BLASTN_PATH


# map genes to contigs, filter the mapping and execute greedy algorithm
def map_and_find(out_dir, genes_file, assembly_graph, assembly_fasta, identity_threshold = mg.DEF_IDENTITY_THRESHOLD,
                 length_threshold = mg.DEF_LENGTH_THRESHOLD, find_fragmented = mg.DEF_FIND_FRAGMENTED, min_gene_density = fp.DEF_MIN_GENE_DENSITY,
                 min_seed_gene_density = fp.DEF_MIN_SEED_GENE_DENSITY, min_length = fp.DEF_MIN_LENGTH, max_length = fp.DEF_MAX_LENGTH,
                 min_read_depth = fp.DEF_MIN_READ_DEPTH, min_plasmid_read_depth = fp.DEF_MIN_PLASMID_READ_DEPTH,
                 max_gc_diff = fp.DEF_MAX_GC_DIFF, max_intermediate_contigs = fp.DEF_MAX_INTERMEDIATE_CONTIGS,
                 max_intermediate_nt = fp.DEF_MAX_INTERMEDIATE_NT, max_score = fp.DEF_MAX_SCORE, score_weights = fp.DEF_SCORE_WEIGHTS,
                 keep_subplasmids = fp.DEF_KEEP_SUBPLASMIDS, overlap_ends = fp.DEF_OVERLAP_ENDS, binning = fp.DEF_BINNING,
                 fanout = fp.DEF_FANOUT, probabilistic = fp.DEF_PROBABILISTIC, use_median = fp.DEF_USE_MEDIAN,
                 use_node_based = fp.DEF_USE_NODE_BASED, verbose = fp.DEF_VERBOSE, makeblastdb = DEF_MAKEBLASTDB_PATH, blastn = DEF_BLASTN_PATH):

    mapping_subdir = os.path.join(out_dir, 'mapping') # subdirectory for files generating while mapping plasmid genes to contigs
    greedy_subdir = os.path.join(out_dir, 'greedy') # subdirectory for outputs of greedy algorithm

    blast_output = os.path.join(mapping_subdir, 'genes_to_contigs.csv') # initial gene-contig mapping
    filtered_blast_output = os.path.join(mapping_subdir, 'filtered_genes_to_contigs.csv') # quality-filtered gene-contig mapping

    call('mkdir -p %s %s' % (mapping_subdir, greedy_subdir), shell = True)

    # search contigs in assembly containing plasmid genes
    if verbose:
        print('\nMapping plasmid genes to assembly contigs...')
    mg.map(blast_output, genes_file, assembly_fasta, verbose = verbose, makeblastdb = makeblastdb, blastn = blastn)

    # filter gene-contig matches by quality and length of match
    if verbose:
        print('\nFiltering gene-contig matches by quality and length to obtain potential seed contigs...')
    mg.filter_blast(genes_file, blast_output, filtered_blast_output, identity_threshold = identity_threshold,
                    length_threshold = length_threshold, find_fragmented = find_fragmented, verbose = verbose)

    # perform greedy plasmid search
    if verbose:
        print('\nSearching greedily for putative plasmids...')
    fp.greedy(assembly_graph, genes_file, filtered_blast_output, greedy_subdir,
              min_gene_density = min_gene_density, min_seed_gene_density = min_seed_gene_density,
              min_length = min_length, max_length = max_length, min_read_depth = min_read_depth,
              min_plasmid_read_depth = min_plasmid_read_depth, max_gc_diff = max_gc_diff,
              max_intermediate_contigs = max_intermediate_contigs, max_intermediate_nt = max_intermediate_nt,
              max_score = max_score, score_weights = score_weights, keep_subplasmids = keep_subplasmids,
              overlap_ends = overlap_ends, binning = binning, fanout = fanout, probabilistic = probabilistic,
              use_median = use_median, use_node_based = use_node_based, verbose = verbose)


if __name__ == '__main__':
    import argparse
    argparser = argparse.ArgumentParser()
    argparser.add_argument('out_dir', help = 'output directory')
    argparser.add_argument('genes_file', help = 'FASTA file of plasmid genes')
    argparser.add_argument('assembly_graph', help = 'assembly graph (GFA format)')
    argparser.add_argument('assembly_fasta', help = 'assembly contigs (FASTA format)')
    argparser.add_argument('-i', '--identity_threshold', type = float, default = mg.DEF_IDENTITY_THRESHOLD, help = 'minimum identity in a hit to keep it')
    argparser.add_argument('-l', '--length_threshold', type = float, default = mg.DEF_LENGTH_THRESHOLD, help = 'minimum fraction of query that has to be matched to keep a hit')
    argparser.add_argument('-f', '--find_fragmented', action = 'store_true', help = 'search for fragmented hits, i.e. several short high-identity hits that together satisfy the length threshold')
    # options of greedy algorithm
    argparser.add_argument('--min_gene_density', type = float, default = fp.DEF_MIN_GENE_DENSITY, help = 'greedy algorithm: minimum gene density of putative plasmids (plasmids with lower gene density are marked as questionable)')
    argparser.add_argument('--min_seed_gene_density', type = float, default = fp.DEF_MIN_SEED_GENE_DENSITY, help = 'greedy algorithm: minimum gene density of a contig to be considered as a seed')
    argparser.add_argument('--min_length', type = int, default = fp.DEF_MIN_LENGTH, help = 'greedy algorithm: minimum length of putative plasmids (shorter plasmids are marked as questionable)')
    argparser.add_argument('--max_length', type = int, default = fp.DEF_MAX_LENGTH, help = 'greedy algorithm: maximum length of putative plasmids (longer plasmids are marked as questionable)')
    argparser.add_argument('--min_read_depth', type = float, default = fp.DEF_MIN_READ_DEPTH, help = 'greedy algorithm: minimum read depth of contigs to be used in plasmids')
    argparser.add_argument('--min_plasmid_read_depth', type = float, default = fp.DEF_MIN_PLASMID_READ_DEPTH, help = 'greedy algorithm: minimum read depth of putative plasmids (plasmids with lower read depth are marked as questionable)')
    argparser.add_argument('--max_gc_diff', type = float, default = fp.DEF_MAX_GC_DIFF, help = 'greedy algorithm: maximum difference in GC content between contig to be added and current plasmid')
    argparser.add_argument('--max_intermediate_contigs', type = int, default = fp.DEF_MAX_INTERMEDIATE_CONTIGS, help = 'greedy algorithm: maximum number of contigs between gene-containing contigs in a plasmid')
    argparser.add_argument('--max_intermediate_nt', type = int, default = fp.DEF_MAX_INTERMEDIATE_NT, help = 'greedy algorithm: maximum sum of lengths of contigs between gene-containing contigs in a plasmid')
    argparser.add_argument('--max_score', type = float, default = fp.DEF_MAX_SCORE, help = 'greedy algorithm: maximum score for potential extensions')
    argparser.add_argument('--score_weights', type = str, default = fp.DEF_SCORE_WEIGHTS, help = 'weights of the score components (comma-separated list of elements <component>=<weight>)')
    argparser.add_argument('--keep_subplasmids', action = 'store_true', help = 'greedy algorithm: do not mark plasmids contained by others as questionable')
    argparser.add_argument('--overlap_ends', type = int, default = fp.DEF_OVERLAP_ENDS, help = 'greedy algorithm: minimum overlap between ends of plasmid sequence for circularisation')
    argparser.add_argument('--binning', type = float, default = fp.DEF_BINNING, help = 'greedy algorithm: number standard deviations the read depth and GC content of plasmids are allowed to differ from the centre of their bin (binning is activated by setting the parameter to any value unequal NaN)')
    argparser.add_argument('--fanout', type = int, default = fp.DEF_FANOUT, help = 'greedy algorithm: maximum number of predecessors / successors of any contig in a plasmid (resp. contig collection)')
    argparser.add_argument('--probabilistic', action = 'store_true', help = 'greedy algorithm: activates probabilistic extensions')
    argparser.add_argument('--use_median', action = 'store_true', help = 'greedy algorithm: determine average read depth of a plasmid using the median instead of the mean')
    argparser.add_argument('--use_node_based', action = 'store_true', help = 'greedy algorithm: activates node-based extensions')
    argparser.add_argument('--verbose', action = 'store_true', help = 'print more information')
    argparser.add_argument('--makeblastdb', default = DEF_MAKEBLASTDB_PATH, help = 'path to makeblastdb executable')
    argparser.add_argument('--blastn', default = DEF_BLASTN_PATH, help = 'path to blastn executable')
    args = argparser.parse_args()

    try:
        map_and_find(args.out_dir, args.genes_file, args.assembly_graph, args.assembly_fasta,
                     identity_threshold = args.identity_threshold, length_threshold = args.length_threshold, find_fragmented = args.find_fragmented,
                     min_gene_density = args.min_gene_density, min_seed_gene_density = args.min_seed_gene_density,
                     min_length = args.min_length, max_length = args.max_length, min_read_depth = args.min_read_depth,
                     min_plasmid_read_depth = args.min_plasmid_read_depth, max_gc_diff = args.max_gc_diff,
                     max_intermediate_contigs = args.max_intermediate_contigs, max_intermediate_nt = args.max_intermediate_nt,
                     max_score = args.max_score, score_weights = args.score_weights, keep_subplasmids = args.keep_subplasmids,
                     overlap_ends = args.overlap_ends, binning = args.binning, fanout = args.fanout, probabilistic = args.probabilistic,
                     use_median = args.use_median, use_node_based = args.use_node_based, verbose = args.verbose,
                     makeblastdb = args.makeblastdb, blastn = args.blastn)
    except ValueError as err:
        print('ERROR: %s' % err)
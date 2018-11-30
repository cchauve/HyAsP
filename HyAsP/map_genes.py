#!/usr/bin/env python

# Implements commands 'map' and 'filter' of HyAsP.
#
# Creates gene-contig mapping by searching the genes in the contigs using BLAST.
# The contigs can be provided directly as a FASTA file or indirectly through as assembly graph (GFA format).
#
# Filtering discards hits whose identity is too low or which do not cover enough of the gene.
#
# Requirements:
#  - makeblastdb (makeblastdb / --makeblastdb; tested with BLAST+ v2.6.0)
#  - blastn (blastn / --blastn; tested with BLAST+ v2.6.0)
#  - standard UNIX tools (rm)
#
# makeblastdb and blastn have to be in $PATH or explicitly specified using the respective path options.


import os.path

from subprocess import call

from HyAsP import filter_blast_mapping as fbm


# default values / constants
DEF_FROM_FASTA = ''
DEF_FROM_GFA = ''
DEF_CONTIGS_FILE = 'assembly_contigs.fasta'
DEF_CONTIGS_DB = 'contigs_blast_db'
DEF_IDENTITY_THRESHOLD = fbm.DEF_IDENTITY_THRESHOLD
DEF_LENGTH_THRESHOLD = fbm.DEF_LENGTH_THRESHOLD
DEF_FIND_FRAGMENTED = fbm.DEF_FIND_FRAGMENTED
DEF_CLEAN = False
DEF_VERBOSE = False
DEF_MAKEBLASTDB_PATH = 'makeblastdb'
DEF_BLASTN_PATH = 'blastn'


# create BLAST database from the given contigs and map the genes to them using blastn
def map_genes_fasta(genes_file, contigs_file, mapping_file, clean = DEF_CLEAN, verbose = DEF_VERBOSE,
                    makeblastdb = DEF_MAKEBLASTDB_PATH, blastn = DEF_BLASTN_PATH):

    blast_db = os.path.join(os.path.dirname(os.path.abspath(mapping_file)), DEF_CONTIGS_DB)

    # search contigs in assembly containing plasmid genes
    if verbose:
        print('\nMapping plasmid genes to assembly contigs...')
    call('%s -in %s -dbtype nucl -out %s' % (makeblastdb, contigs_file, blast_db), shell = True)
    call('%s -task megablast -query %s -db %s -out %s -outfmt 6' % (blastn, genes_file, blast_db, mapping_file), shell = True)

    if clean:
        call('rm contigs_blast_db.nhr contigs_blast_db.nin contigs_blast_db.nsq')


# extract contigs (FASTA) from assembly graph (GFA), create BLAST database from them
# and map the genes to them using blastn
def map_genes_gfa(genes_file, assembly_file, mapping_file, clean = DEF_CLEAN, verbose = DEF_VERBOSE,
                  makeblastdb = DEF_MAKEBLASTDB_PATH, blastn = DEF_BLASTN_PATH):

    contigs_file = os.path.join(os.path.dirname(os.path.abspath(mapping_file)), DEF_CONTIGS_FILE)
    blast_db = os.path.join(os.path.dirname(os.path.abspath(mapping_file)), DEF_CONTIGS_DB)

    with open(assembly_file, 'r') as in_file, open(contigs_file, 'w') as out_file:
        for line in in_file:
            if line.startswith('S'):
                tokens = line.split('\t')

                name = tokens[1]
                seq = str(tokens[2])

                out_file.write('>%s\n%s\n' % (name, seq))

    # search contigs in assembly containing plasmid genes
    if verbose:
        print('\nMapping plasmid genes to assembly contigs...')
    call('%s -in %s -dbtype nucl -out %s' % (makeblastdb, contigs_file, blast_db), shell = True)
    call('%s -task megablast -query %s -db %s -out %s -outfmt 6' % (blastn, genes_file, blast_db, mapping_file), shell = True)

    if clean:
        call('rm %s' % ' '.join([contigs_file, blast_db + '.nhr', blast_db + '.nin', blast_db + '.nsq']), shell = True)


# print configuration of mapping genes to contigs
def show_config_map(genes_file, mapping_file, used_src, from_fasta, from_gfa, clean, verbose, makeblastdb, blastn):
    print('################################')
    print('### Configuration of mapping ###\n')

    print('>>> Input / output')
    print('Genes file: %s' % genes_file)
    print('Contigs source: %s' % used_src)
    print('Mapping file: %s' % mapping_file)

    print('\n>>> Other options')
    print('Clean: %i' % clean)
    print('Verbose: %i' % verbose)
    print('makeblastdb executable: %s' % makeblastdb)
    print('blastn executable: %s' % blastn)

    print('################################\n')


# map the genes to the contigs
def map(mapping_file, genes_file, from_fasta = DEF_FROM_FASTA, from_gfa = DEF_FROM_GFA, clean = DEF_CLEAN,
        verbose = DEF_VERBOSE, makeblastdb = DEF_MAKEBLASTDB_PATH, blastn = DEF_BLASTN_PATH):

    used_src = [src for src in [from_fasta, from_gfa] if src != '']  # exactly one will be != '' (form required but mutually exclusive group)
    if len(used_src) != 1:
        print('ERROR: Only one source (from_fasta, from_gfa) can be used at once.')
        return
    used_src = used_src[0]

    show_config_map(genes_file, mapping_file, used_src, from_fasta, from_gfa, clean, verbose, makeblastdb, blastn)

    if from_fasta != DEF_FROM_FASTA:
        map_genes_fasta(genes_file, from_fasta, mapping_file, clean = clean, verbose = verbose, makeblastdb = makeblastdb,
                        blastn = blastn)
    else:
        map_genes_gfa(genes_file, from_gfa, mapping_file, clean = clean, verbose = verbose, makeblastdb = makeblastdb,
                      blastn = blastn)


# print configuration of filtering of mapping file
def show_config_filter(genes_file, blast_output, filtered_blast_output, identity_threshold, length_threshold, find_fragmented, verbose):
    print('##################################')
    print('### Configuration of filtering ###\n')

    print('>>> Input / output')
    print('Genes file: %s' % genes_file)
    print('Unfiltered mapping: %s' % blast_output)
    print('Filtered mapping: %s' % filtered_blast_output)

    print('\n>>> Filtering options')
    print('Identity threshold: %f' % identity_threshold)
    print('Length threshold: %f' % length_threshold)
    print('Find fragmented hits: %i' % find_fragmented)

    print('\n>>> Other options')
    print('Verbose: %i' % verbose)

    print('##################################\n')


# filter gene-contig matches by quality and length of match
def filter_blast(genes_file, blast_output, filtered_blast_output, identity_threshold = DEF_IDENTITY_THRESHOLD,
                 length_threshold = DEF_LENGTH_THRESHOLD, find_fragmented = DEF_FIND_FRAGMENTED, verbose = DEF_VERBOSE):

    show_config_filter(genes_file, blast_output, filtered_blast_output, identity_threshold, length_threshold, find_fragmented, verbose)

    if verbose:
        print('\nFiltering gene-contig matches by quality and length to obtain potential seed contigs...')
    queries = fbm.read_queries(genes_file)
    fbm.filter_mapping(queries, blast_output, filtered_blast_output, identity_threshold = identity_threshold,
                       length_threshold = length_threshold, find_fragmented = find_fragmented, verbose = verbose)
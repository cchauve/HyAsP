#!/usr/bin/env python

# Implements command 'create' of HyAsP.
#
# Takes a list of plasmids with annotated genes and extracts all genes into a FASTA file.
# The plasmids can be provided as GenBank files or as accession numbers (in this case,
# the GenBank files are temporarily downloaded).
#
# In addition, a plasmid table obtained from NCBI can be used as the source.
# Such a table can be downloaded from https://www.ncbi.nlm.nih.gov/genome/browse#!/plasmids/
# and all possible columns should be selected.
# Additional filtering options are available for this source:
#  - choose only plasmids released before a particular date
#  - blacklisting of plasmids not to be included in the database
#  - choosing the type of sequence database (RefSeq, GenBank)
#
# There are further options allowing to dereplicate the genes (based on their sequence),
# to consider only plasmids in a certain length range, and to collect the plasmid sequences in a FASTA file.
#
# Requirements:
# standard UNIX tools (curl, rm)


import math
import os
import pandas as pd

from subprocess import call

from Bio import SeqIO


# default values / constants
DEF_FROM_ACCESSION = ''
DEF_FROM_GENBANK = ''
DEF_FROM_PLASMID_TABLE = ''
DEF_KEEP_PLASMIDS = ''
DEF_DEREPLICATE = False
DEF_FROM_COMMAND_LINE = False
DEF_EXTEND = False
DEF_RELEASED_BEFORE = ''
DEF_TYPE = 'both'
DEF_BLACKLIST = ''
DEF_MIN_LENGTH = 0
DEF_MAX_LENGTH = math.inf
DEF_MIN_GENE_LENGTH = 0
DEF_NUM_ATTEMPTS = 25
DEF_VERBOSE = False


# reads GenBank file and extracts all genes from it
def extract_all_genes(gb_file):
    gene_collection = []
    for gb_record in SeqIO.parse(open(gb_file, 'r'), 'genbank'):
        # find all genes
        genes = []
        for feature in gb_record.features:
            if feature.type == 'gene':
                genes.append(feature)

        # create FASTA entry per gene
        i = 0
        for g in genes:
            if 'locus_tag' in g.qualifiers:
                gene_collection.append((g.qualifiers['locus_tag'][0], str(g.extract(gb_record.seq))))
            elif 'gene' in g.qualifiers:
                gene_collection.append(('gene%i_%s' % (i, g.qualifiers['gene'][0]), str(g.extract(gb_record.seq))))
            else:
                print('WARNING: No identifier for a gene in %s:' % gb_file)
                print(g)
            i += 1

    return gene_collection


# reads GenBank file and extracts the sequence of the plasmid
def extract_seq(gb_file):
    seqs = []
    for gb_record in SeqIO.parse(open(gb_file, 'r'), 'genbank'):
        seqs.append((gb_record.name, str(gb_record.seq)))

    return seqs


# creates gene database (and plasmids database) from given list of plasmids
def create_db(genes_file, from_accession, sources, blacklist, dereplicate = DEF_DEREPLICATE, keep_plasmids = DEF_KEEP_PLASMIDS,
              extend = DEF_EXTEND, min_length = DEF_MIN_LENGTH, max_length = DEF_MAX_LENGTH, min_gene_length = DEF_MIN_GENE_LENGTH,
              num_attempts = DEF_NUM_ATTEMPTS, verbose = DEF_VERBOSE):
    temp = genes_file + '_tmp_gb'

    # get genes from old database (and dereplicate them)
    if verbose:
        print('Reading genes from old database %s (and dereplicating them, if activated)...' % genes_file)
    old_genes = []
    gene_seqs = set()
    if extend:
        if os.path.isfile(genes_file):
            with open(genes_file, 'r') as in_genes:
                for record in SeqIO.parse(in_genes, 'fasta'):
                    seq = str(record.seq)
                    if ((not dereplicate) or (seq not in gene_seqs)) and len(seq) >= min_gene_length:
                        gene_seqs.add(seq)
                        old_genes.append((record.id, seq))
        else:
            print('WARNING: File %s with existing gene database does not exist. Starting with an empty gene database.' % genes_file)

    # clear existing plasmids database if one should be created but not extended
    if not extend and keep_plasmids != '':
        with open(keep_plasmids, 'w') as out_plasmids:
            out_plasmids.write('')

    # remove blacklisted accessions
    sources = [src for src in sources if src not in blacklist]
    if verbose:
        print('%i plasmid accessions after removing blacklisted plasmids.' % len(sources))

    # write old (dereplicated) contents and add genes from new plasmids
    if verbose:
        print('Adding genes (and plasmids) to database...')
    num_genes = len(old_genes)
    num_length_discarded = 0
    num_attempts_discarded = 0
    with open(genes_file, 'w') as out_genes:
        for id, seq in old_genes:
            out_genes.write('>%s\n%s\n' % (id, seq))
        for src in sources:
            plasmid_seqs = []
            if from_accession:
                res = -1
                cnt_failures = 0
                while res != 0 and cnt_failures < num_attempts:
                    res = call('curl -s "https://eutils.be-md.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=%s&rettype=gbwithparts&retmode=text" > %s' % (src, temp), shell = True)
                    if res != 0:
                        print('WARNING: Download of plasmid %s failed. Trying again.' % src)
                        cnt_failures += 1
                    else:
                        plasmid_seqs = extract_seq(temp)
                        if len(plasmid_seqs) == 0:
                            print('WARNING: No sequence for plasmid %s was found in the GenBank file. Download seems to be faulty. Repeating download.' % src)
                            res = -1
                            cnt_failures += 1

                if cnt_failures == num_attempts:
                    print('ERROR: Plasmid %s could not be downloaded properly (within %i attempts). Continuing with next plasmid.' % (src, num_attempts))
                    num_attempts_discarded += 1
                    plasmid_seqs = None
                    genes = []
                else:
                    genes = extract_all_genes(temp)

            else:
                plasmid_seqs = extract_seq(src)
                genes = extract_all_genes(src)

            if plasmid_seqs is not None:
                if len(plasmid_seqs) > 1:
                    print('WARNING: More than one plasmid sequence in the GenBank file. Using only the first one.')
                plasmid_name, plasmid_seq = plasmid_seqs[0]

                if min_length <= len(plasmid_seq) <= max_length:
                    for id, seq in genes:
                        if ((not dereplicate) or (seq not in gene_seqs)) and len(seq) >= min_gene_length:
                            gene_seqs.add(seq)
                            out_genes.write('>%s\n%s\n' % (id, seq))
                            num_genes += 1

                    if keep_plasmids != '':
                        with open(keep_plasmids, 'a') as out_plasmids:
                            out_plasmids.write('>%s\n%s\n' % (plasmid_name, plasmid_seq))
                else:
                    num_length_discarded += 1

    if verbose and num_length_discarded > 0:
        print('%i plasmids were discarded for their length.' % num_length_discarded)
    if verbose and num_attempts_discarded > 0:
        print('%i plasmids were discarded for reaching the download-attempt limit.' % num_attempts_discarded)

    print('Database comprises %i genes.' % num_genes)
    if from_accession:
        call('rm %s' % temp, shell = True)


# creates gene database (and plasmids database) from the NCBI plasmid table
def create_db_from_table(genes_file, table_file, blacklist, dereplicate = DEF_DEREPLICATE, keep_plasmids = DEF_KEEP_PLASMIDS,
                         extend = DEF_EXTEND, released_before = DEF_RELEASED_BEFORE, type = DEF_TYPE, min_length = DEF_MIN_LENGTH,
                         max_length = DEF_MAX_LENGTH, min_gene_length = DEF_MIN_GENE_LENGTH, num_attempts = DEF_NUM_ATTEMPTS,
                         verbose = DEF_VERBOSE):
    temp = genes_file + '_tmp_gb'

    # get genes from old database (and dereplicate them)
    if verbose:
        print('Reading genes from old database %s (and dereplicating them, if activated)...' % genes_file)
    old_genes = []
    gene_seqs = set()
    if extend:
        if os.path.isfile(genes_file):
            with open(genes_file, 'r') as in_genes:
                for record in SeqIO.parse(in_genes, 'fasta'):
                    seq = str(record.seq)
                    if ((not dereplicate) or (seq not in gene_seqs)) and len(seq) >= min_gene_length:
                        gene_seqs.add(seq)
                        old_genes.append((record.id, seq))
        else:
            print('WARNING: File %s with existing gene database does not exist. Starting with an empty gene database.' % genes_file)

    # clear existing plasmids database if one should be created but not extended
    if not extend and keep_plasmids != '':
        with open(keep_plasmids, 'w') as out_plasmids:
            out_plasmids.write('')

    plasmid_table = pd.read_csv(table_file, sep = ',', dtype = str)

    if released_before != '':
        plasmid_table = plasmid_table.loc[plasmid_table['Release Date'] < released_before]

    # Possible entry formats of the 'Replicon' column:
    # a) [<identifier>:]<RefSeq accession>/<GenBank accession>
    # b) [<identifier>:]<RefSeq accession>/
    # c) [<identifier>:]<GenBank accession>
    #
    # Rule 1: if the entry contains '/', the RefSeq accession is specified and is the first token of the split result
    # Rule 2: if the entry does not end in '/', the GenBank accession is specified and is the last (and maybe only) token of the split result
    accessions = []
    if type == 'both':
        accessions = [acc for entry in plasmid_table['Replicons'] for acc in entry.split(':')[-1].split('/') if acc != '']
    if type == 'RefSeq':
        accessions = [entry.split(':')[-1].split('/')[0] for entry in plasmid_table['Replicons'] if '/' in entry]
    if type == 'GenBank':
        accessions = [entry.split(':')[-1].split('/')[-1] for entry in plasmid_table['Replicons'] if entry[-1] != '/']
    if verbose:
        print('%i plasmid accessions extracted.' % len(accessions))

    # remove blacklisted accessions
    accessions = [acc for acc in accessions if acc not in blacklist]
    print('%i plasmid accessions after removing blacklisted plasmids.' % len(accessions))

    # write old (dereplicated) contents and add genes from new plasmids
    if verbose:
        print('Adding genes (and plasmids) to database...')
    num_genes = len(old_genes)
    num_length_discarded = 0
    num_attempts_discarded = 0
    with open(genes_file, 'w') as out_genes:
        for id, seq in old_genes:
            out_genes.write('>%s\n%s\n' % (id, seq))
        for i, acc in enumerate(accessions, start = 1):
            if verbose:
                print('%i / %i' % (i, len(accessions)))
            plasmid_seqs = []
            res = -1
            cnt_failures = 0
            while res != 0 and cnt_failures < num_attempts:
                res = call('curl -s "https://eutils.be-md.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=%s&rettype=gbwithparts&retmode=text" > %s' % (acc, temp), shell = True)
                if res != 0:
                    print('WARNING: Download of plasmid %s failed. Trying again.' % acc)
                    cnt_failures += 1
                else:
                    plasmid_seqs = extract_seq(temp)
                    if len(plasmid_seqs) == 0:
                        print('WARNING: No sequence for plasmid %s was found in the GenBank file. Download seems to be faulty. Repeating download.' % acc)
                        res = -1
                        cnt_failures += 1

            if cnt_failures == num_attempts:
                print('ERROR: Plasmid %s could not be downloaded properly (within %i attempts). Continuing with next plasmid.' % (acc, num_attempts))
                num_attempts_discarded += 1
                plasmid_seqs = None
                genes = []
            else:
                genes = extract_all_genes(temp)

            if plasmid_seqs is not None:
                if len(plasmid_seqs) > 1:
                    print('WARNING: More than one plasmid sequence in the GenBank file. Using only the first one.')
                plasmid_name, plasmid_seq = plasmid_seqs[0]

                if min_length <= len(plasmid_seq) <= max_length:
                    for id, seq in genes:
                        if ((not dereplicate) or (seq not in gene_seqs)) and len(seq) >= min_gene_length:
                            gene_seqs.add(seq)
                            out_genes.write('>%s\n%s\n' % (id, seq))
                            num_genes += 1

                    if keep_plasmids != '':
                        if verbose:
                            print('%s -> %s' % (acc, plasmid_name))
                        with open(keep_plasmids, 'a') as out_plasmids:
                            out_plasmids.write('>%s\n%s\n' % (plasmid_name, plasmid_seq))
                else:
                    num_length_discarded += 1

    if verbose and num_length_discarded > 0:
        print('%i plasmids were discarded for their length.' % num_length_discarded)
    if verbose and num_attempts_discarded > 0:
        print('%i plasmids were discarded for reaching the download-attempt limit.' % num_attempts_discarded)

    print('Database comprises %i genes.' % num_genes)
    if len(accessions) > 0:
        call('rm %s' % temp, shell = True)


# print configuration of database generation
def show_config(used_src, from_accession, from_genbank, from_plasmid_table, genes_file, keep_plasmids, dereplicate, extend,
                released_before, type, blacklist, min_length, max_length, min_gene_length, num_attempts, verbose):
    print('############################################')
    print('### Configuration of database generation ###\n')

    print('>>> Input / output')
    print('Plasmid source: %s' % used_src)
    if used_src == from_plasmid_table:
        print('Type: %s' % type)
    print('Genes file: %s' % genes_file)
    print('Plasmids file: %s' % (keep_plasmids if keep_plasmids != DEF_KEEP_PLASMIDS else '---'))

    print('\n>>> Filtering')
    print('Minimum plasmid length: %f' % min_length)
    print('Maximum plasmid length: %f' % max_length)
    print('Minimum gene length: %f' % min_gene_length)
    print('Released before: %s' % (released_before if released_before != DEF_RELEASED_BEFORE else '(filter not applied)'))
    print('Blacklist: %s' % (','.join(blacklist) if len(blacklist) > 0 else '(empty)'))

    print('\n>>> Other options')
    print('Dereplicate: %i' % dereplicate)
    print('Extend: %i' % extend)
    print('Maximum number of download attempts: %f' % num_attempts)
    print('Verbose: %i' % verbose)

    print('############################################\n')


# choose correct method depending on inputs and make sure that a valid selections of options is provided
def create(genes_file, from_accession = DEF_FROM_ACCESSION, from_genbank = DEF_FROM_GENBANK,
           from_plasmid_table = DEF_FROM_PLASMID_TABLE, keep_plasmids = DEF_KEEP_PLASMIDS, dereplicate = DEF_DEREPLICATE,
           from_command_line = DEF_FROM_COMMAND_LINE, extend = DEF_EXTEND, released_before = DEF_RELEASED_BEFORE,
           type = DEF_TYPE, blacklist = DEF_BLACKLIST, min_length = DEF_MIN_LENGTH, max_length = DEF_MAX_LENGTH,
           min_gene_length = DEF_MIN_GENE_LENGTH, num_attempts = DEF_NUM_ATTEMPTS, verbose = DEF_VERBOSE):

    if from_plasmid_table != '' and from_command_line:
        print('ERROR: A plasmid table (-p) cannot be read from command line (-c). Please specify the file path and remove option -c.')
        return

    if from_plasmid_table == '' and released_before != '':
        print('ERROR: Option -r is only available toegether with -p.')
        return

    if extend and blacklist != '':
        print('ERROR: Blacklisting (-b) cannot be combined with option -e.')
        return

    used_src = [src for src in [from_accession, from_genbank, from_plasmid_table] if src != ''] # exactly one will be != '' (form required but mutually exclusive group)
    if len(used_src) != 1:
        print('ERROR: Only one source (-a, -g or -p) can be used at once.')
        return
    used_src = used_src[0]

    blacklisted = []
    if blacklist != '':
        with open(blacklist, 'r') as in_file:
            for line in in_file:
                blacklisted.append(line.strip())

    show_config(used_src, from_accession, from_genbank, from_plasmid_table, genes_file, keep_plasmids, dereplicate,
                extend, released_before, type, blacklisted, min_length, max_length, min_gene_length, num_attempts, verbose)

    if used_src == from_plasmid_table:
        create_db_from_table(genes_file, from_plasmid_table, blacklisted, dereplicate, keep_plasmids, extend, released_before, type, min_length, max_length, min_gene_length, num_attempts, verbose)
    else:
        if from_command_line:
            sources = used_src.split(',')
        else:
            with open(used_src, 'r') as infile:
                sources = [line.strip() for line in infile if line.strip()]

        create_db(genes_file, from_accession != '', sources, blacklisted, dereplicate, keep_plasmids, extend, min_length, max_length, min_gene_length, num_attempts, verbose)
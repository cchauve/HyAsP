#!/usr/bin/env python

# Pipeline for finding plasmids in FASTQ read data.
#
# The pipeline takes FASTQ reads and a collection of (plasmid) genes as input.
# The reads are preprocessed (sickle, Trim Galore) and assembled (Unicycler).
# The preprocessing and an analysis of the read data using FastQC are optional.
# Subsequently, the genes are mapped to the assembly contigs using BLAST (blastn / megablast)
# and plasmids are predicted in the assembly graph using a greedy algorithm based on the gene occurrences
# as well as other characteristics such as read depth and GC content.
#
# Requirements:
#  - FastQC (fastqc / --fastqc; tested with v0.11.5)
#  - standard UNIX tools (mkdir)
#  # recursively
#  - sickle (sickle / --sickle; tested with v1.33)
#  - Trim Galore (trim_galore / --trim_galore; tested with v0.4.5_dev)
#  - Unicycler (unicycler-runner.py / --unicycler; tested with v0.4.5)
#  - Python (python / --python; tested with v3.5.4)
#  - SPAdes (spades.py / --spades; tested with v3.12.0)
#  - Racon (racon / --racon; tested with v1.3.0)
#  - Pilon (pilon-1.22.jar / --pilon; tested with v1.22)
#  - SAMtools (samtools / --samtools; tested with v1.5)
#  - Bowtie 2 (bowtie2 / --bowtie2 and bowtie2-build / --bowtie2_build; tested with v2.3.3.1)
#  - Java (java / --java; tested with v1.8.0_121)
#  - makeblastdb (makeblastdb / --makeblastdb; tested with BLAST+ v2.6.0)
#  - tblastn (tblastn / --tblastn; tested with BLAST+ v2.6.0)
#  - blastn (blastn / --blastn; tested with BLAST+ v2.6.0)
#  - standard UNIX tools (rm, cat)
#
# FastQC has to be in $PATH or explicitly specified using the respective path option.

import os.path

from subprocess import call

from HyAsP import preprocess_fastq as pf, assemble_unicycler as au, determine_plasmids as dp
from HyAsP import find_plasmids as fp # only for default values



# default values / constants
DEF_FASTQC_PATH = 'fastqc'
DEF_USE_FASTQC = False
DEF_NO_PREPROCESSING = False


# main method organising the different steps of the pipeline
def run_pipeline(analysis_dir, genes_file, first_short_reads = '', second_short_reads = '', single_short_reads = '',
                 long_reads = '', unicycler_mode = au.DEF_UNICYCLER_MODE, unicycler_verbosity = au.DEF_UNICYCLER_VERBOSITY,
                 unicycler_keep = au.DEF_UNICYCLER_KEEP, use_fastqc = DEF_USE_FASTQC, no_preprocessing = DEF_NO_PREPROCESSING,
                 read_qual_threshold = pf.DEF_QUAL_THRESHOLD, min_read_length = pf.DEF_MIN_LENGTH, min_gene_density = fp.DEF_MIN_GENE_DENSITY,
                 identity_threshold = dp.DEF_IDENTITY_THRESHOLD, length_threshold = dp.DEF_LENGTH_THRESHOLD,
                 min_seed_gene_density = fp.DEF_MIN_SEED_GENE_DENSITY, min_length = fp.DEF_MIN_LENGTH,
                 max_length = fp.DEF_MAX_LENGTH, min_read_depth = fp.DEF_MIN_READ_DEPTH, min_plasmid_read_depth = fp.DEF_MIN_PLASMID_READ_DEPTH,
                 max_gc_diff = fp.DEF_MAX_GC_DIFF, max_intermediate_contigs = fp.DEF_MAX_INTERMEDIATE_CONTIGS,
                 max_intermediate_nt = fp.DEF_MAX_INTERMEDIATE_NT, max_score = fp.DEF_MAX_SCORE, score_weights = fp.DEF_SCORE_WEIGHTS,
                 keep_subplasmids = fp.DEF_KEEP_SUBPLASMIDS, overlap_ends = fp.DEF_OVERLAP_ENDS, binning = fp.DEF_BINNING,
                 fanout = fp.DEF_FANOUT, probabilistic = fp.DEF_PROBABILISTIC, use_median = fp.DEF_USE_MEDIAN,
                 use_node_based = fp.DEF_USE_NODE_BASED, verbose = fp.DEF_VERBOSE, fastqc = DEF_FASTQC_PATH, sickle = pf.DEF_SICKLE_PATH,
                 trim_galore = pf.DEF_TRIM_GALORE_PATH, python = au.DEF_PYTHON_PATH, unicycler = au.DEF_UNICYCLER_PATH, spades = au.DEF_SPADES_PATH,
                 racon = au.DEF_RACON_PATH, pilon = au.DEF_PILON_PATH, samtools = au.DEF_SAMTOOLS_PATH, bowtie2 = au.DEF_BOWTIE2_PATH,
                 bowtie2_build = au.DEF_BOWTIE2_BUILD_PATH, java = au.DEF_JAVA_PATH, makeblastdb = dp.DEF_MAKEBLASTDB_PATH,
                 tblastn = au.DEF_TBLASTN_PATH, blastn = dp.DEF_BLASTN_PATH):

    data_dir = os.path.join(analysis_dir, 'data') # directory for (preprocessed) read data and quality-check results
    assembly_dir = os.path.join(analysis_dir, 'assembly') # directory for assembly outputs
    plasmids_dir = os.path.join(analysis_dir, 'plasmids') # directory for results of plasmid contruction

    paired_available = first_short_reads != '' and second_short_reads != ''
    unpaired_available = single_short_reads != ''
    long_available = long_reads != ''

    if not os.path.isdir(data_dir):
        call('mkdir -p ' + data_dir, shell = True)


    ### 1) Analysis of FASTQ reads with FastQC
    # Input:  FASTQ reads (first_short_reads, second_short_reads, single_short_reads, long_reads)
    #         name of data directory (data_dir)
    # Output: results of quality check
    if use_fastqc:
        print('===== Analysis of read quality (before preprocessing) =====\n')

        print('\nPerforming quality checks on FASTQ reads (with FastQC)...')
        if unpaired_available:
            print('Unpaired short reads: %s' % single_short_reads)
            call('%s -o %s %s' % (fastqc, data_dir, single_short_reads), shell = True)

        if paired_available:
            print('Paired short reads: %s, %s' % (first_short_reads, second_short_reads))
            call('%s -o %s %s %s' % (fastqc, data_dir, first_short_reads, second_short_reads), shell = True)

        if long_available:
            print('Long reads: %s' % (long_reads))
            call('%s -o %s %s' % (fastqc, data_dir, long_reads), shell = True)

        print('Results of checks available in %s/.' % data_dir)


    if no_preprocessing:
        # use FASTQ reads directly and skip second round of FastQC analysis
        first_prep = first_short_reads
        second_prep = second_short_reads
        single_short_prep = single_short_reads
        long_prep = long_reads

    else:
        ### 2) Preprocessing the read data
        # Input:  FASTQ read file(s) (first_short_reads, second_short_reads, single_short_reads, long_reads)
        #         name of data directory (data_dir)
        # Output: preprocessed FASTQ read file(s)
        print('\n===== Preprocessing of FASTQ read data =====\n')

        first_prep, second_prep, single_short_prep, long_prep = pf.preprocess(data_dir, first_short_reads, second_short_reads, single_short_reads, long_reads,
                                                                              qual_threshold = read_qual_threshold, min_length = min_read_length,
                                                                              sickle = sickle, trim_galore = trim_galore)


        ### 3) Analysis of preprocessed FASTQ reads with FastQC
        # Input:  FASTQ reads (first_prep, second_prep, single_short_prep, long_prep)
        #         name of data directory (data_dir)
        # Output: results of quality check
        if use_fastqc:
            print('===== Analysis of read quality (after preprocessing) =====\n')

            print('\nPerforming quality checks on FASTQ reads (with FastQC)...')
            if unpaired_available:
                print('Unpaired short reads: %s' % single_short_prep)
                call('%s -o %s %s' % (fastqc, data_dir, single_short_prep), shell = True)

            if paired_available:
                print('Paired short reads: %s, %s' % (first_prep, second_prep))
                call('%s -o %s %s %s' % (fastqc, data_dir, first_prep, second_prep), shell = True)

            if long_available:
                print('Long reads: %s' % (long_prep))
                call('%s -o %s %s' % (fastqc, data_dir, long_prep), shell = True)

            print('Results of checks available in %s/.' % data_dir)


    ### 4) Assembly of preprocessed FASTQ reads
    # Input:  preprocessed FASTQ read file(s) (first_prep, second_prep, single_short_prep, long_prep)
    #         assembly mode for unicycler (unicycler_mode)
    #         name of assembly directory (assembly_dir)
    # Output: Unicycler assembly
    print('\n===== Assemby of FASTQ reads =====\n')

    if not os.path.isdir(assembly_dir):
        call('mkdir -p ' + assembly_dir, shell = True)

    final_assembly_graph, final_assembly_fasta = au.assemble(assembly_dir, first_prep, second_prep, single_short_prep, long_prep,
                                                             unicycler_mode = unicycler_mode, unicycler_verbosity = unicycler_verbosity,
                                                             unicycler_keep = unicycler_keep, python = python, unicycler = unicycler,
                                                             spades = spades, racon = racon, pilon = pilon, samtools = samtools,
                                                             bowtie2 = bowtie2, bowtie2_build = bowtie2_build, java = java,
                                                             makeblastdb = makeblastdb, tblastn = tblastn)


    ### 5) Determination of putative plasmids
    # Input:  final assembly graph (final_assembly_graph)
    #         FASTA file of contigs of final assembly graph (final_assembly_fasta)
    #         FASTA files of genes to map (genes_file)
    # Output: results of greedy plasmid algorithm
    print('\n===== Step 5: Determination of putative plasmids via greedy algorithm =====\n')

    # map plasmid genes to assembly contigs and determine related contigs
    dp.map_and_find(plasmids_dir, genes_file, final_assembly_graph, final_assembly_fasta,
                    identity_threshold = identity_threshold, length_threshold = length_threshold,
                    min_gene_density = min_gene_density, min_seed_gene_density = min_seed_gene_density,
                    min_length = min_length, max_length = max_length, min_read_depth = min_read_depth,
                    min_plasmid_read_depth = min_plasmid_read_depth, max_gc_diff = max_gc_diff,
                    max_intermediate_contigs = max_intermediate_contigs, max_intermediate_nt = max_intermediate_nt,
                    max_score = max_score, score_weights = score_weights, keep_subplasmids = keep_subplasmids,
                    overlap_ends = overlap_ends, binning = binning, fanout = fanout, probabilistic = probabilistic,
                    use_median = use_median, use_node_based = use_node_based, verbose = verbose,
                    makeblastdb = makeblastdb, blastn = blastn)



def main():
    import argparse
    argparser = argparse.ArgumentParser()
    argparser.add_argument('analysis_dir', help = 'directory into which all data and outputs are stored')
    argparser.add_argument('genes_file', help = '(plasmid) genes to map to contigs')
    argparser.add_argument('-1', '--first_short_reads', default = '', help = 'FASTQ file of first paired reads (if any)')
    argparser.add_argument('-2', '--second_short_reads', default = '', help = 'FASTQ file of second paired reads (if any)')
    argparser.add_argument('-s', '--single_short_reads', default = '', help = 'FASTQ file of single / unpaired reads (if any)')
    argparser.add_argument('-l', '--long_reads', default = '', help = 'FASTQ file of long reads (if any)')
    argparser.add_argument('-u', '--unicycler_mode', default = au.DEF_UNICYCLER_MODE, help = 'Unicycler assembly mode (conservative, normal or bold)')
    argparser.add_argument('-f', '--use_fastqc', action = 'store_true', help = 'if set, FASTQ reads are analysed with FastQC before and after preprocessing')
    argparser.add_argument('-n', '--no_preprocessing', action = 'store_true', help = 'if set, the FASTQ reads are not preprocessed before the assembly')
    argparser.add_argument('--read_qual_threshold', type = float, default = pf.DEF_QUAL_THRESHOLD, help = 'threshold for trimming low-quality ends')
    argparser.add_argument('--min_read_length', type = float, default = pf.DEF_MIN_LENGTH, help = 'minimum length of reads after quality / adapter trimming')
    argparser.add_argument('--unicycler_verbosity', type = int, default = au.DEF_UNICYCLER_VERBOSITY, help = 'Unicycler level of stdout and log file information (0 - 3)')
    argparser.add_argument('--unicycler_keep', type = int, default = au.DEF_UNICYCLER_KEEP, help = 'Unicycler level of file retention (0 - 3)')
    argparser.add_argument('--identity_threshold', type = float, default = dp.DEF_IDENTITY_THRESHOLD, help = 'filtering the gene-contig mapping: minimum identity in a hit to keep it')
    argparser.add_argument('--length_threshold', type = float, default = dp.DEF_LENGTH_THRESHOLD, help = 'filtering the gene-contig mapping: minimum fraction of query that has to be matched to keep a hit')
    # options of greedy algorithm
    argparser.add_argument('--min_gene_density', type = float, default = fp.DEF_MIN_GENE_DENSITY, help = 'greedy algorithm: minimum gene density of predicted plasmids (are marked as questionable otherwise)')
    argparser.add_argument('--min_seed_gene_density', type = float, default = fp.DEF_MIN_SEED_GENE_DENSITY, help = 'greedy algorithm: minimum gene density of a contig to be considered as a seed')
    argparser.add_argument('--min_length', type = int, default = fp.DEF_MIN_LENGTH, help = 'greedy algorithm: minimum length of predicted plasmid (shorter ones are marked as questionable)')
    argparser.add_argument('--max_length', type = int, default = fp.DEF_MAX_LENGTH, help = 'greedy algorithm: maximum length of predicted plasmid (longer ones are marked as questionable)')
    argparser.add_argument('--min_read_depth', type = float, default = fp.DEF_MIN_READ_DEPTH, help = 'greedy algorithm: minimum read depth of contig to be used in plasmids')
    argparser.add_argument('--min_plasmid_read_depth', type = float, default = fp.DEF_MIN_PLASMID_READ_DEPTH, help = 'greedy algorithm: minimum read depth of predicted plasmids (are marked as questionable otherwise)')
    argparser.add_argument('--max_gc_diff', type = float, default = fp.DEF_MAX_GC_DIFF, help = 'greedy algorithm: maximum difference in GC content between contig to be added and current plasmid')
    argparser.add_argument('--max_intermediate_contigs', type = int, default = fp.DEF_MAX_INTERMEDIATE_CONTIGS, help = 'greedy algorithm: maximum number of contigs between gene-containing contigs in a plasmid')
    argparser.add_argument('--max_intermediate_nt', type = int, default = fp.DEF_MAX_INTERMEDIATE_NT, help = 'greedy algorithm: maximum sum of lengths of contigs between gene-containing contigs in a plasmid')
    argparser.add_argument('--max_score', type = float, default = fp.DEF_MAX_SCORE, help = 'greedy algorithm: maximum score for potential extensions')
    argparser.add_argument('--score_weights', type = str, default = fp.DEF_SCORE_WEIGHTS, help = 'greedy algorithm: weights of the score components (comma-separated list of elements <component>=<weight>)')
    argparser.add_argument('--keep_subplasmids', action = 'store_true', help = 'greedy algorithm: mark plasmids contained by others as questionable')
    argparser.add_argument('--overlap_ends', type = int, default = fp.DEF_OVERLAP_ENDS, help = 'greedy algorithm: minimum overlap between ends of plasmid sequence for circularisation')
    argparser.add_argument('--binning', type = float, default = fp.DEF_BINNING, help = 'greedy algorithm: setting a threshold value activates plasmid binning')
    argparser.add_argument('--fanout', type = int, default = fp.DEF_FANOUT, help = 'greedy algorithm: maximum degree of branching of contigs (per end)')
    argparser.add_argument('--probabilistic', action = 'store_true', help = 'greedy algorithm: choose extension probabilistically')
    argparser.add_argument('--use_median', action = 'store_true', help = 'greedy algorithm: determine average read depth of a plasmid using the median instead of the mean')
    argparser.add_argument('--use_node_based', action = 'store_true', help = 'greedy algorithm: activates node-based extensions')
    argparser.add_argument('--verbose', action = 'store_true', help = 'greedy algorithm: prints detailed information during the extension process')
    argparser.add_argument('--fastqc', default = DEF_FASTQC_PATH, help = 'path to FastQC executable')
    argparser.add_argument('--sickle', default = pf.DEF_SICKLE_PATH, help = 'path to sickle executable')
    argparser.add_argument('--trim_galore', default = pf.DEF_TRIM_GALORE_PATH, help = 'path to Trim Galore executable')
    argparser.add_argument('--python', default = au.DEF_PYTHON_PATH, help = 'path to python executable')
    argparser.add_argument('--unicycler', default = au.DEF_UNICYCLER_PATH, help = 'path to Unicycler executable')
    argparser.add_argument('--spades', default = au.DEF_SPADES_PATH, help = 'path to SPAdes executable')
    argparser.add_argument('--racon', default = au.DEF_RACON_PATH, help = 'path to Racon executable')
    argparser.add_argument('--pilon', default = au.DEF_PILON_PATH, help = 'path to Pilon executable')
    argparser.add_argument('--samtools', default = au.DEF_SAMTOOLS_PATH, help = 'path to SAMtools executable')
    argparser.add_argument('--bowtie2', default = au.DEF_BOWTIE2_PATH, help = 'path to Bowtie 2 executable')
    argparser.add_argument('--bowtie2_build', default = au.DEF_BOWTIE2_BUILD_PATH, help = 'path to bowtie2_build executable')
    argparser.add_argument('--java', default = au.DEF_JAVA_PATH, help = 'path to Java executable')
    argparser.add_argument('--makeblastdb', default = au.DEF_MAKEBLASTDB_PATH, help = 'path to makeblastdb executable')
    argparser.add_argument('--tblastn', default = au.DEF_TBLASTN_PATH, help = 'path to tblastn executable')
    argparser.add_argument('--blastn', default = dp.DEF_BLASTN_PATH, help = 'path to blastn executable')
    args = argparser.parse_args()

    if args.first_short_reads == '' and args.second_short_reads == '' and args.single_short_reads == '' and args.long_reads == '':
        print('ERROR: No read data is specified (at least -1 / -2 or -s or -l have to be used).')
    elif args.first_short_reads != '' and args.second_short_reads == '' or args.first_short_reads == '' and args.second_short_reads != '':
        print('ERROR: Specified paired read data is incomplete. Both options -1 and -2 have to be used when specifying paired read data.')
    else:
        run_pipeline(args.analysis_dir, args.genes_file, args.first_short_reads, second_short_reads = args.second_short_reads,
                     single_short_reads = args.single_short_reads, long_reads = args.long_reads,
                     unicycler_mode = args.unicycler_mode, unicycler_verbosity = args.unicycler_verbosity,
                     unicycler_keep = args.unicycler_keep, use_fastqc = args.use_fastqc, no_preprocessing = args.no_preprocessing,
                     read_qual_threshold = args.read_qual_threshold, min_read_length = args.min_read_length,
                     identity_threshold = args.identity_threshold, length_threshold = args.length_threshold, min_gene_density = args.min_gene_density,
                     min_seed_gene_density = args.min_seed_gene_density, min_length = args.min_length, max_length = args.max_length,
                     min_read_depth = args.min_read_depth, min_plasmid_read_depth = args.min_plasmid_read_depth, max_gc_diff = args.max_gc_diff,
                     max_intermediate_contigs = args.max_intermediate_contigs, max_intermediate_nt = args.max_intermediate_nt,
                     max_score = args.max_score, score_weights = args.score_weights, keep_subplasmids = args.keep_subplasmids,
                     overlap_ends = args.overlap_ends, binning = args.binning, fanout = args.fanout, probabilistic = args.probabilistic,
                     use_median = args.use_median, use_node_based = args.use_node_based, verbose = args.verbose,
                     fastqc = args.fastqc, sickle = args.sickle, trim_galore = args.trim_galore, python = args.python, unicycler = args.unicycler,
                     spades = args.spades, racon = args.racon, pilon = args.pilon, samtools = args.samtools,
                     bowtie2 = args.bowtie2, bowtie2_build = args.bowtie2_build, java = args.java,
                     makeblastdb = args.makeblastdb, tblastn = args.tblastn, blastn = args.blastn)


if __name__ == '__main__':
    main()
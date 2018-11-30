#!/usr/bin/env python

# Assembles the given FASTQ read data using Unicycler.
# Different combinations of (un)paired short and long reads are possible as input.
# Logging (--verbosity) and retention of intermediate files (--keep) are set to Unicycler's defaults.
#
# Requirements:
#  - Unicycler (as unicycler-runner.py in $PATH or specified via --unicycler; tested with v0.4.5)
#  - Python (as python in $PATH; tested with v3.5.4)
#
# Unicycler has a number of dependencies which have to be in $PATH or explicitly specified using the respective
# path option:
#  - SPAdes (spades.py / --spades; tested with v3.12.0)
#  - Racon (racon / --racon; tested with v1.3.0)
#  - Pilon (pilon-1.22.jar / --pilon; tested with v1.22)
#  - SAMtools (samtools / --samtools; tested with v1.5)
#  - Bowtie 2 (bowtie2 / --bowtie2 and bowtie2-build / --bowtie2_build; tested with v2.3.3.1)
#  - Java (java / --java; tested with v1.8.0_121)
#  - makeblastdb (makeblastdb / --makeblastdb; tested with BLAST+ v2.6.0)
#  - tblastn (tblastn / --tblastn; tested with BLAST+ v2.6.0)


import os

from subprocess import call


# default values / constants
DEF_UNICYCLER_MODE = 'normal'
DEF_UNICYCLER_VERBOSITY = 1
DEF_UNICYCLER_KEEP = 1
DEF_PYTHON_PATH = 'python'
DEF_UNICYCLER_PATH = 'unicycler-runner.py'
DEF_SPADES_PATH = 'spades.py'
DEF_RACON_PATH = 'racon'
DEF_PILON_PATH = 'pilon-1.22.jar'
DEF_SAMTOOLS_PATH = 'samtools'
DEF_BOWTIE2_PATH = 'bowtie2'
DEF_BOWTIE2_BUILD_PATH = 'bowtie2-build'
DEF_JAVA_PATH = 'java'
DEF_MAKEBLASTDB_PATH = 'makeblastdb'
DEF_TBLASTN_PATH = 'tblastn'
DEF_VERBOSE = False


# create assembly (GFA format) from the given FASTQ reads using Unicycler
def assemble(out_dir, first_short_reads = '', second_short_reads = '', single_short_reads = '', long_reads = '',
             unicycler_mode = DEF_UNICYCLER_MODE, unicycler_verbosity = DEF_UNICYCLER_VERBOSITY, unicycler_keep = DEF_UNICYCLER_KEEP,
             python = DEF_PYTHON_PATH, unicycler = DEF_UNICYCLER_PATH, spades = DEF_SPADES_PATH, racon = DEF_RACON_PATH,
             pilon = DEF_PILON_PATH, samtools = DEF_SAMTOOLS_PATH, bowtie2 = DEF_BOWTIE2_PATH, bowtie2_build = DEF_BOWTIE2_BUILD_PATH,
             java = DEF_JAVA_PATH, makeblastdb = DEF_MAKEBLASTDB_PATH, tblastn = DEF_TBLASTN_PATH, verbose = DEF_VERBOSE):

    gfa_output = os.path.join(out_dir, 'assembly.gfa')
    fasta_output = os.path.join(out_dir, 'assembly.fasta')

    # put together the input-files portion of the Unicycler command
    inputs = ''
    msg = 'Assembling'
    if first_short_reads != '' and second_short_reads != '':
        inputs += '--short1 %s --short2 %s' % (first_short_reads, second_short_reads)
        msg += ' paired short-read data from %s and %s' % (first_short_reads, second_short_reads)

    if single_short_reads != '':
        inputs += ' --unpaired %s' % single_short_reads
        msg += (' and' if inputs != '' else '') + ' unpaired short-read data from %s' % single_short_reads

    if long_reads != '':
        inputs += ' --long %s' % long_reads
        msg += (' and' if inputs != '' else '') + ' long-read data from %s' % long_reads

    msg += ' (mode: %s).' % unicycler_mode

    # put together the non-default dependencies
    dependencies = ''
    if spades != DEF_SPADES_PATH:
        dependencies += ' --spades_path %s' % spades
    if racon != DEF_RACON_PATH:
        dependencies += ' --racon_path %s' % racon
    if pilon != DEF_PILON_PATH:
        dependencies += ' --pilon_path %s' % pilon
    if samtools != DEF_SAMTOOLS_PATH:
        dependencies += ' --samtools_path %s' % samtools
    if bowtie2 != DEF_BOWTIE2_PATH:
        dependencies += ' --bowtie2_path %s' % bowtie2
    if bowtie2_build != DEF_BOWTIE2_BUILD_PATH:
        dependencies += ' --bowtie2_build_path %s' % bowtie2_build
    if java != DEF_JAVA_PATH:
        dependencies += ' --java_path %s' % java
    if makeblastdb != DEF_MAKEBLASTDB_PATH:
        dependencies += ' --makeblastdb_path %s' % makeblastdb
    if tblastn != DEF_TBLASTN_PATH:
        dependencies += ' --tblastn_path %s' % tblastn

    # perform assembly
    if verbose:
        print(msg)
        print('Assembly files: %s, %s' % (gfa_output, fasta_output))
    call('%s %s %s --out %s --mode %s %s --min_fasta_length 0 --verbosity %i --keep %i'
         % (python, unicycler, inputs, out_dir, unicycler_mode, dependencies, unicycler_verbosity, unicycler_keep), shell = True)

    return gfa_output, fasta_output


if __name__ == '__main__':
    import argparse
    argparser = argparse.ArgumentParser()
    argparser.add_argument('out_dir', help = 'output directory')
    argparser.add_argument('-1', '--first_short_reads', default = '', help = 'first reads of paired FASTQ read data (if any)')
    argparser.add_argument('-2', '--second_short_reads', default = '', help = 'second reads of paired FASTQ read data (if any)')
    argparser.add_argument('-s', '--single_short_reads', default = '', help = 'unpaired FASTQ read data (if any)')
    argparser.add_argument('-l', '--long_reads', default = '', help = 'long FASTQ read data (if any)')
    argparser.add_argument('-u', '--unicycler_mode', default = DEF_UNICYCLER_MODE, help = 'Unicycler assembly mode (conservative, normal or bold)')
    argparser.add_argument('-v', '--unicycler_verbosity', type = int, default = DEF_UNICYCLER_VERBOSITY, help = 'Unicycler level of stdout and log file information (0 - 3)')
    argparser.add_argument('-k', '--unicycler_keep', type = int, default = DEF_UNICYCLER_KEEP, help = 'Unicycler level of file retention (0 - 3)')
    argparser.add_argument('--verbose', action = 'store_true', help = 'print more information')
    argparser.add_argument('--python', default = DEF_PYTHON_PATH, help = 'path to python executable')
    argparser.add_argument('--unicycler', default = DEF_UNICYCLER_PATH, help = 'path to Unicycler executable')
    argparser.add_argument('--spades', default = DEF_SPADES_PATH, help = 'path to SPAdes executable')
    argparser.add_argument('--racon', default = DEF_RACON_PATH, help = 'path to Racon executable')
    argparser.add_argument('--pilon', default = DEF_PILON_PATH, help = 'path to Pilon executable')
    argparser.add_argument('--samtools', default = DEF_SAMTOOLS_PATH, help = 'path to SAMtools executable')
    argparser.add_argument('--bowtie2', default = DEF_BOWTIE2_PATH, help = 'path to Bowtie 2 executable')
    argparser.add_argument('--bowtie2_build', default = DEF_BOWTIE2_BUILD_PATH, help = 'path to bowtie2_build executable')
    argparser.add_argument('--java', default = DEF_JAVA_PATH, help = 'path to Java executable')
    argparser.add_argument('--makeblastdb', default = DEF_MAKEBLASTDB_PATH, help = 'path to makeblastdb executable')
    argparser.add_argument('--tblastn', default = DEF_TBLASTN_PATH, help = 'path to tblastn executable')
    args = argparser.parse_args()


    if args.first_short_reads == '' and args.second_short_reads == '' and args.single_short_reads == '' and args.long_reads == '':
        print('ERROR: No read data is specified (at least -1 / -2 or -s or -l have to be used).')
    elif args.first_short_reads != '' and args.second_short_reads == '' or args.first_short_reads == '' and args.second_short_reads != '':
        print('ERROR: Specified paired read data is incomplete. Both options -1 and -2 have to be used when specifying paired read data.')
    else:
        assemble(args.out_dir, args.first_short_reads, second_short_reads = args.second_short_reads,
                 single_short_reads = args.single_short_reads, long_reads = args.long_reads,
                 unicycler_mode = args.unicycler_mode, unicycler_verbosity = args.unicycler_verbosity,
                 unicycler_keep = args.unicycler_keep, python = args.python, unicycler = args.unicycler,
                 spades = args.spades, racon = args.racon, pilon = args.pilon, samtools = args.samtools,
                 bowtie2 = args.bowtie2, bowtie2_build = args.bowtie2_build, java = args.java,
                 makeblastdb = args.makeblastdb, tblastn = args.tblastn)
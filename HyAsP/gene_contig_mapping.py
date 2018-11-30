#!/usr/bin/env python

import pandas as pd

from Bio import SeqIO


# Representation of gene-hits on contigs (based on BLAST output file, default output format 6).
# Provides easy access to, e.g., number of genes found on a contig or the proportion of it covered by genes.
class GeneContigMapping:
    __slots__ = 'data_',\
                'gene_lengths_',\
                'contig_lengths_'

    def __init__(self, blast_file, genes_file, assembly_graph):
        col_names = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']  # outfmt 6
        self.data_ = pd.read_csv(blast_file, sep = '\t', names = col_names, dtype = str)

        self.gene_lengths_ = dict()
        with open(genes_file, 'r') as in_file:
            for record in SeqIO.parse(in_file, 'fasta'):
                self.gene_lengths_[record.id] = len(record.seq)

        self.contig_lengths_ = dict()
        with open(assembly_graph, 'r') as in_file:
            for line in in_file:
                if line.startswith('S'):
                    tokens = line.split('\t')
                    name = tokens[1]
                    seq = str(tokens[2])
                    self.contig_lengths_[name] = len(seq)

    # return identifiers of all genes in the mapping
    def genes(self):
        return sorted(self.data_.qseqid.unique())

    # return identifiers of all contigs in the mapping
    def contigs(self):
        return sorted(self.data_.sseqid.unique())

    # return identifiers of all genes mapping to the given contig
    def genes_of(self, contig):
        return sorted(self.data_.loc[self.data_.sseqid == contig].qseqid.unique())

    # return identifiers of all genes mapping to the given collection of contigs
    def genes_of_set(self, contigs):
        return sorted(self.data_.loc[self.data_.sseqid.isin(contigs)].qseqid.unique())

    # return number of hits for the given contig
    def num_gene_locations_of(self, contig):
        rows = self.data_.loc[self.data_.sseqid == contig]
        return len(set(zip(map(int, rows.sstart), map(int, rows.send))))

    # return number of hits for the given collection of contigs
    def num_gene_locations_of_set(self, contigs):
        rows = self.data_.loc[self.data_.sseqid.isin(contigs)]
        return len(set(zip(rows.sseqid, map(int, rows.sstart), map(int, rows.send))))

    # return identifiers of contigs to which the given gene maps
    def contigs_of(self, gene):
        return sorted(self.data_.loc[self.data_.qseqid == gene].sseqid.unique())

    # return identifiers of contigs to which the given genes map
    # a contig is reported as soon as one of the genes maps to it
    def contigs_of_set(self, genes):
        return sorted(self.data_.loc[self.data_.qseqid.isin(genes)].sseqid.unique())

    # return gene density (i.e. proportion of contig covered by genes) of the given contig
    def gene_density(self, contig):
        return self.num_gene_covered_bases(contig) / self.contig_lengths_[contig]

    # return gene frequency (i.e. number of gene hits per bases) of the given contig
    def gene_frequency(self, contig):
        return len(self.data_.loc[self.data_.sseqid == contig].qseqid) / self.contig_lengths_[contig]

    # return the number of bases covered by at least one gene
    def num_gene_covered_bases(self, contig):
        rows = self.data_.loc[self.data_.sseqid == contig]
        intervals = list(zip(map(int, rows.sstart), map(int, rows.send)))  # portions of contig covered by genes
        intervals = [(start, end) if start <= end else (end, start) for start, end in
                     intervals]  # swap start / end if gene was found on reverse complement
        intervals.sort(key=lambda x: x[0])  # intervals is now sorted by start position

        num_pos_covered = 0
        last_pos_covered = 0  # last (right-most) position of contig covered so far
        for start, end in intervals:
            if end <= last_pos_covered:
                pass  # contained in previous interval -> no new position covered
            else:
                num_pos_covered += end - max(last_pos_covered + 1, start) + 1
                last_pos_covered = end

        return num_pos_covered

    # return whether the given gene occurs in the mapping
    def contains_gene(self, gene):
        return gene in self.data_.qseqid.values

    # return whether the given contig occurs in the mapping
    def contains_contig(self, contig):
        return contig in self.data_.sseqid.values

    # check whether the gene-contig mapping is consistent with the used gene database (i.e. it contains no
    # gene that is not in the database)
    def check_consistency(self):
        mapping_genes = set(self.genes())
        db_genes = set(self.gene_lengths_)
        return mapping_genes.difference(db_genes)
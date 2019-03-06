#!/usr/bin/env python

# Implements command 'find' of HyAsP.
#
# Predicts plasmids based on an assembly graph and hits of plasmid genes on its contigs.
# Plasmids are contructed from gene-containig contigs (seeds) by greedily adding the best neighbour(s)
# in the assembly graph.
#
# See class SeedEnumerator for the order in which the seeds are chosen.
#
# See methods score_depth_and_gene_density(...) and pick_by_depth_and_genes(...) for how neighbours are scored
# and chosen as candidates for extension.
#
# See method extend(...) (resp. extend_probabilistic(...)) for how the greedy extension is performed.


import math
import os
import random

from Bio.Seq import Seq

from HyAsP import assembly_graph as ag_module
from HyAsP import gene_contig_mapping as gcm_module


# default values / constants
DEF_MIN_GENE_DENSITY = 0.3
DEF_MIN_SEED_GENE_DENSITY = math.nan
DEF_MIN_LENGTH = 1500
DEF_MAX_LENGTH = 1750000
DEF_MIN_READ_DEPTH = math.nan
DEF_MIN_PLASMID_READ_DEPTH = math.nan
DEF_MAX_GC_DIFF = 0.15
DEF_MAX_INTERMEDIATE_CONTIGS = 2
DEF_MAX_INTERMEDIATE_NT = 2000
DEF_MAX_SCORE = math.inf
DEF_SCORE_WEIGHTS = 'depth_diff=1,gene_density=1,gc_diff=1'
DEF_KEEP_SUBPLASMIDS = False
DEF_OVERLAP_ENDS = math.inf
DEF_BINNING = math.nan
DEF_FANOUT = 1
DEF_PROBABILISTIC = False
DEF_USE_MEDIAN = False
DEF_USE_NODE_BASED = False
DEF_VERBOSE = False


# Representation of a plasmid
#
# Ideally, the plasmid consists of a linear chain of contigs (with orientation) that can form a cycle.
# However, this representation can also handle 'degenerated' plasmids with branching,
# i.e. contigs can have multiple left / right extensions.
# Currently, each contig can occur only once in a particular plasmid instance (node-based)
# or each link can occur only once in a particular plasmid instance (link-based).
#
# Note: This representation assumes that the overlap between contigs is 0.
# Currently, computations of mean depth etc. for overlapping contigs have to be done externally.
# See generate_plasmid_infos(...).
# TODO redesign plasmid representation to handle both cases properly
class Plasmid:
    __slots__ = 'contigs_', \
                'orientations_', \
                'depths_', \
                'lengths_', \
                'gc_contents_', \
                'intermediate_contigs_', \
                'intermediate_nt_', \
                'scores_', \
                'left_extensions_', \
                'right_extensions_', \
                'extendabilities_', \
                'fan_out_', \
                'min_read_depth_', \
                'max_length_', \
                'max_gc_diff_', \
                'max_intermediate_contigs_', \
                'max_intermediate_nt_', \
                'use_median_', \
                'score_weights_', \
                'node_based_'

    def __init__(self, seed, depth, length, gc_content, gene_density, min_read_depth, max_length, max_gc_diff,
                 max_intermediate_contigs, max_intermediate_nt, weights, fan_out = 1, orientation = '+',
                 use_median = False, node_based = False):
        self.contigs_ = dict()
        self.orientations_ = dict()
        self.depths_ = dict()
        self.lengths_ = dict()
        self.gc_contents_ = dict()
        self.intermediate_contigs_ = dict()
        self.intermediate_nt_ = dict()
        self.scores_ = dict()
        self.left_extensions_ = dict()
        self.right_extensions_ = dict()

        self.contigs_[0] = seed  # 0 is always the id of the seed
        self.orientations_[0] = orientation
        self.depths_[0] = depth
        self.lengths_[0] = length
        self.gc_contents_[0] = gc_content
        self.intermediate_contigs_[0] = 0
        self.intermediate_nt_[0] = 0
        self.score_weights_ = weights
        self.scores_[0] = self.score_depth_and_gene_density(depth, gc_content, depth, gene_density, gc_content) # for seed, only gene density adds to score

        self.extendabilities_ = ['none', 'right', 'left', 'both']
        self.fan_out_ = fan_out
        self.min_read_depth_ = min_read_depth
        self.max_length_ = max_length
        self.max_gc_diff_ = max_gc_diff
        self.max_intermediate_contigs_ = max_intermediate_contigs
        self.max_intermediate_nt_ = max_intermediate_nt
        self.use_median_ = use_median
        self.node_based_ = node_based

    # return name of seed contig
    def seed(self):
        return self.contigs_[0]

    # return name of contig
    def contig_name(self, contig):
        return self.contigs_[contig]

    # return orientation of a given contig in this plasmid
    def orientation(self, contig):
        return self.orientations_[contig]

    # return read depth of a given contig in this plasmid
    def depth(self, contig):
        return self.depths_[contig]

    # return length of a given contig
    def length(self, contig):
        return self.lengths_[contig]

    # return length of a given contig
    def gc_content(self, contig):
        return self.gc_contents_[contig]

    # return number of contigs between the given contig and the nearest gene-containing contig
    def contigs_to_genes(self, contig):
        return self.intermediate_contigs_[contig]

    # return sum of lengths of contigs between the given contig and the nearest gene-containing contig
    def nt_to_genes(self, contig):
        return self.intermediate_nt_[contig]

    # return mean read depth of the plasmid
    def mean_depth(self):
        total_length = sum([self.length(c) for c in self.contigs_])
        return sum([self.depth(c) * self.length(c) for c in self.contigs_]) / (total_length if total_length > 0 else 1)

    # return median read depth of the plasmid
    def median_depth(self):
        sel = list(self.contigs_.keys())
        sel.sort(key = lambda s: self.depths_[s])
        total_num_values = sum([self.lengths_[seg] for seg in sel])
        index = total_num_values * 0.50
        fraction50 = index - math.floor(index)  # fraction for linear interpolation (median)
        rank50 = math.floor(index)  # rank of interest (median)

        count_sum = 0
        median_val = -1
        for i in range(0, len(sel)):
            count_sum += self.lengths_[sel[i]]

            if count_sum >= rank50:
                left = self.depths_[sel[i]]
                right = self.depths_[sel[i]] if count_sum >= (rank50 + 1) else self.depths_[sel[i + 1]]
                median_val = left + (right - left) * fraction50
                break

        return median_val

    # return read depth of plasmid (value depends on whether mean or median is used)
    def overall_depth(self):
        return self.median_depth() if self.use_median_ else self.mean_depth()

    # return overall GC content of plasmid
    def overall_gc_content(self):
        total_length = sum([self.length(c) for c in self.contigs_])
        return sum([self.gc_content(c) * self.length(c) for c in self.contigs_]) / (total_length if total_length > 0 else 1)

    # check whether contig can be extended to the left (i.e. does not already have left extensions)
    def is_left_extendable(self, contig):
        return contig not in self.left_extensions_

    # check whether contig can be extended to the right (i.e. does not already have right extensions)
    def is_right_extendable(self, contig):
        return contig not in self.right_extensions_

    # return integer indicating how the given contig can be extended
    # 0 = no extensions possible, 1 = only right extensions possible; 2 = only left extensions possible; 3 = left and right extensions possible
    def extendability(self, contig):
        type = 0
        type += 1 if self.is_right_extendable(contig) else 0
        type += 2 if self.is_left_extendable(contig) else 0
        return self.extendabilities_[type]

    # return list of contigs that can be extended
    def get_extendable_contigs(self):
        return sorted(set(self.left_extensions_) ^ set(self.right_extensions_)) if len(self.contigs_) > 1 else [0]

    # return list of all contigs (as names) underlying this plasmid
    def get_contig_names(self):
        return [self.contigs_[id] for id in self.get_contig_ids()]

    # return list of all contigs (as ids) underlying this plasmid
    def get_contig_ids(self):
        return sorted(self.contigs_.keys())

    # extend the plasmid by searching for possible extensions at its outer ends, scoring and filtering them and choosing the best
    def extend(self, ag, gcm, max_score, verbose):
        return self.extend_node_based(ag, gcm, max_score, verbose) if self.node_based_ else self.extend_link_based(ag, gcm, max_score, verbose)

    # extend the plasmid probabilistically by searching for possible extensions at its outer ends, filtering them
    # and choosing one randomly (probability of choosing a particular candidate depends on its share of the total read depth of all candidates)
    def extend_probabilistic(self, ag, gcm, max_score, verbose):
        if self.node_based_:
            self.extend_probabilistic_node_based(ag, gcm, max_score, verbose)
        else:
            self.extend_probabilistic_link_based(ag, gcm, max_score, verbose)

    # perform non-probabilistic extension(s), node-based mode
    def extend_node_based(self, ag, gcm, max_score, verbose):
        if verbose:
            print('Searching for next extension(s)...')

        extendable_contigs = self.get_extendable_contigs()
        if verbose and len(extendable_contigs) == 0:
            print('\tNo extendable contigs found.')

        candidates = [(contig, cand_contig_name, cand_orientation, cand_direction, cand_score) for \
                      contig in extendable_contigs for cand_contig_name, cand_orientation, cand_direction, cand_score in \
                      self.pick_by_depth_and_genes_node_based(ag, gcm, contig, self.orientation(contig), self.extendability(contig), verbose)]
        candidates = [c for c in candidates if c[1] is not None and c[4] <= max_score
                      and sum(self.lengths_.values()) + (ag.length(c[1]) if c[1] not in self.contigs_.values() else 0) <= self.max_length_]

        if verbose:
            print('%i candidates satisfy maximum score %f.' % (len(candidates), max_score))

        if len(candidates) > 0:
            candidates.sort(key = lambda c: c[4])
        candidates = candidates[:self.fan_out_]

        if verbose and len(candidates) > 0:
            print('Performing %i best extension(s)...' % len(candidates))
        for contig, cand_contig_name, cand_orientation, direction, score in candidates:
            prev_occ = [id for id in self.contigs_ if self.contigs_[id] == cand_contig_name] # not empty when another connection to a contig already included will be added (cycle is closed)
            cand_contig_id = len(self.contigs_)

            if verbose:
                print('\tFrom %s (id = %i) %s to %s%s (id = %i, score %f)' % (self.contig_name(contig), contig, direction, cand_contig_name, cand_orientation, cand_contig_id, score))

            if len(prev_occ) == 0:
                self.contigs_[cand_contig_id] = cand_contig_name
                self.scores_[cand_contig_id] = score
                self.orientations_[cand_contig_id] = cand_orientation
                self.depths_[cand_contig_id] = ag.depth(cand_contig_name)
                self.lengths_[cand_contig_id] = ag.length(cand_contig_name)
                self.gc_contents_[cand_contig_id] = ag.gc_content(cand_contig_name)

                if len(gcm.genes_of(cand_contig_name)) > 0:
                    self.intermediate_contigs_[cand_contig_id] = 0
                    self.intermediate_nt_[cand_contig_id] = 0
                else:
                    self.intermediate_contigs_[cand_contig_id] = self.intermediate_contigs_[contig] + 1
                    self.intermediate_nt_[cand_contig_id] = self.intermediate_nt_[contig] + ag.length(cand_contig_name)

                if direction == 'left':
                    if contig in self.left_extensions_:
                        self.left_extensions_[contig].append(cand_contig_id)
                    else:
                        self.left_extensions_[contig] = [cand_contig_id]
                    self.right_extensions_[cand_contig_id] = [contig]
                else:
                    if contig in self.right_extensions_:
                        self.right_extensions_[contig].append(cand_contig_id)
                    else:
                        self.right_extensions_[contig] = [cand_contig_id]
                    self.left_extensions_[cand_contig_id] = [contig]

            else:
                prev_id = prev_occ[0]
                self.scores_[prev_id] = min(self.scores_[prev_id], score)
                self.intermediate_contigs_[prev_id] = min(self.intermediate_contigs_[prev_id], self.intermediate_contigs_[contig] + 1)
                self.intermediate_nt_[prev_id] = min(self.intermediate_nt_[prev_id], self.intermediate_nt_[contig] + ag.length(cand_contig_name))

                if direction == 'left':
                    if contig in self.left_extensions_:
                        self.left_extensions_[contig].append(prev_id)
                    else:
                        self.left_extensions_[contig] = [prev_id]
                    if prev_id in self.right_extensions_:
                        self.right_extensions_[prev_id].append(contig)
                    else:
                        self.right_extensions_[prev_id] = [contig]
                else:
                    if contig in self.right_extensions_:
                        self.right_extensions_[contig].append(prev_id)
                    else:
                        self.right_extensions_[contig] = [prev_id]
                    if prev_id in self.left_extensions_:
                        self.left_extensions_[prev_id].append(contig)
                    else:
                        self.left_extensions_[prev_id] = [contig]

        if verbose:
            print('New total score: %f' % sum([self.scores_[c] for c in self.contigs_]))
            print('New average score: %f' % (sum([self.scores_[c] for c in self.contigs_]) / len(self.contigs_)))
            print()

        return len(candidates) > 0

    # perform probabilistic extension, node-based mode
    def extend_probabilistic_node_based(self, ag, gcm, max_score, verbose):
        if verbose:
            print('Searching for next extension(s)...')

        extendable_contigs = self.get_extendable_contigs()
        if verbose and len(extendable_contigs) == 0:
            print('\tNo extendable contigs found.')

        candidates = [(contig, cand_contig_name, cand_orientation, cand_direction, cand_score) for \
                      contig in extendable_contigs for cand_contig_name, cand_orientation, cand_direction, cand_score in \
                      self.pick_by_depth_and_genes_node_based(ag, gcm, contig, self.orientation(contig), self.extendability(contig), verbose)]
        candidates = [c for c in candidates if c[1] is not None and c[4] <= max_score
                      and sum(self.lengths_.values()) + (ag.length(c[1]) if c[1] not in self.contigs_.values() else 0) <= self.max_length_]

        if verbose:
            print('%i candidates satisfy maximum score %f.' % (len(candidates), max_score))

        if len(candidates) > 0:
            total_depth = sum([ag.depth(c[1]) for c in candidates])
            weights = [ag.depth(c[1]) / total_depth for c in candidates]
            last_sum = 0
            intervals = []
            for i in range(0, len(candidates)):
                intervals.append((last_sum, last_sum + weights[i], i))
                last_sum += weights[i]
            rnd = random.random()
            best_option = [index for left, right, index in intervals if left <= rnd < right][0]
            candidates = [candidates[best_option]]

        if verbose and len(candidates) > 0:
            print('Performing extension probabilistically...')
        for contig, cand_contig_name, cand_orientation, direction, score in candidates:
            prev_occ = [id for id in self.contigs_ if self.contigs_[id] == cand_contig_name] # not empty when another connection to a contig already included will be added (cycle is closed)
            cand_contig_id = len(self.contigs_)

            if verbose:
                print('\tFrom %s (id = %i) %s to %s%s (id = %i, score %f)' % (self.contig_name(contig), contig, direction, cand_contig_name, cand_orientation, cand_contig_id, score))

            if len(prev_occ) == 0:
                self.contigs_[cand_contig_id] = cand_contig_name
                self.scores_[cand_contig_id] = score
                self.orientations_[cand_contig_id] = cand_orientation
                self.depths_[cand_contig_id] = ag.depth(cand_contig_name)
                self.lengths_[cand_contig_id] = ag.length(cand_contig_name)
                self.gc_contents_[cand_contig_id] = ag.gc_content(cand_contig_name)

                if len(gcm.genes_of(cand_contig_name)) > 0:
                    self.intermediate_contigs_[cand_contig_id] = 0
                    self.intermediate_nt_[cand_contig_id] = 0
                else:
                    self.intermediate_contigs_[cand_contig_id] = self.intermediate_contigs_[contig] + 1
                    self.intermediate_nt_[cand_contig_id] = self.intermediate_nt_[contig] + ag.length(cand_contig_name)

                if direction == 'left':
                    if contig in self.left_extensions_:
                        self.left_extensions_[contig].append(cand_contig_id)
                    else:
                        self.left_extensions_[contig] = [cand_contig_id]
                    self.right_extensions_[cand_contig_id] = [contig]
                else:
                    if contig in self.right_extensions_:
                        self.right_extensions_[contig].append(cand_contig_id)
                    else:
                        self.right_extensions_[contig] = [cand_contig_id]
                    self.left_extensions_[cand_contig_id] = [contig]

            else:
                prev_id = prev_occ[0]
                self.scores_[prev_id] = min(self.scores_[prev_id], score)
                self.intermediate_contigs_[prev_id] = min(self.intermediate_contigs_[prev_id], self.intermediate_contigs_[contig] + 1)
                self.intermediate_nt_[prev_id] = min(self.intermediate_nt_[prev_id], self.intermediate_nt_[contig] + ag.length(cand_contig_name))

                if direction == 'left':
                    if contig in self.left_extensions_:
                        self.left_extensions_[contig].append(prev_id)
                    else:
                        self.left_extensions_[contig] = [prev_id]
                    if prev_id in self.right_extensions_:
                        self.right_extensions_[prev_id].append(contig)
                    else:
                        self.right_extensions_[prev_id] = [contig]
                else:
                    if contig in self.right_extensions_:
                        self.right_extensions_[contig].append(prev_id)
                    else:
                        self.right_extensions_[contig] = [prev_id]
                    if prev_id in self.left_extensions_:
                        self.left_extensions_[prev_id].append(contig)
                    else:
                        self.left_extensions_[prev_id] = [contig]

        if verbose:
            print('New total score: %f' % sum([self.scores_[c] for c in self.contigs_]))
            print('New average score: %f' % (sum([self.scores_[c] for c in self.contigs_]) / len(self.contigs_)))
            print()

        return len(candidates) > 0

    # perform non-probabilistic extension(s), link-based mode
    def extend_link_based(self, ag, gcm, max_score, verbose):
        if verbose:
            print('Searching for next extension(s)...')

        extendable_contigs = self.get_extendable_contigs()
        if verbose and len(extendable_contigs) == 0:
            print('\tNo extendable contigs found.')

        candidates = [(contig, cand_contig_name, cand_orientation, cand_direction, cand_score) for \
                      contig in extendable_contigs for cand_contig_name, cand_orientation, cand_direction, cand_score in \
                      self.pick_by_depth_and_genes_link_based(ag, gcm, contig, self.orientation(contig), self.extendability(contig), verbose)]
        candidates = [c for c in candidates if c[1] is not None and c[4] <= max_score
                      and sum(self.lengths_.values()) + ag.length(c[1]) <= self.max_length_]

        if verbose:
            print('%i candidates satisfy maximum score %f.' % (len(candidates), max_score))

        if len(candidates) > 0:
            candidates.sort(key = lambda c: c[4])
        candidates = candidates[:self.fan_out_]

        if verbose and len(candidates) > 0:
            print('Performing %i best extension(s)...' % len(candidates))
        for contig, cand_contig_name, cand_orientation, direction, score in candidates:
            cand_contig_id = len(self.contigs_)
            self.contigs_[cand_contig_id] = cand_contig_name
            occs = [id for id in self.contigs_ if self.contigs_[id] == cand_contig_name]
            share = ag.depth(cand_contig_name) / len(occs)
            for id in occs:
                self.depths_[id] = share

            if verbose:
                print('\tFrom %s (id = %i) %s to %s%s (id = %i, score %f)' % (self.contig_name(contig), contig, direction, cand_contig_name, cand_orientation, cand_contig_id, score))
            self.scores_[cand_contig_id] = score
            self.orientations_[cand_contig_id] = cand_orientation
            self.lengths_[cand_contig_id] = ag.length(cand_contig_name)
            self.gc_contents_[cand_contig_id] = ag.gc_content(cand_contig_name)

            if len(gcm.genes_of(cand_contig_name)) > 0:
                self.intermediate_contigs_[cand_contig_id] = 0
                self.intermediate_nt_[cand_contig_id] = 0
            else:
                self.intermediate_contigs_[cand_contig_id] = self.intermediate_contigs_[contig] + 1
                self.intermediate_nt_[cand_contig_id] = self.intermediate_nt_[contig] + ag.length(cand_contig_name)

            if direction == 'left':
                if contig in self.left_extensions_:
                    self.left_extensions_[contig].append(cand_contig_id)
                else:
                    self.left_extensions_[contig] = [cand_contig_id]
                if cand_contig_id in self.right_extensions_:
                    self.right_extensions_[cand_contig_id].append(contig)
                else:
                    self.right_extensions_[cand_contig_id] = [contig]
            else:
                if contig in self.right_extensions_:
                    self.right_extensions_[contig].append(cand_contig_id)
                else:
                    self.right_extensions_[contig] = [cand_contig_id]
                if cand_contig_id in self.left_extensions_:
                    self.left_extensions_[cand_contig_id].append(contig)
                else:
                    self.left_extensions_[cand_contig_id] = [contig]

        if verbose:
            print('New total score: %f' % sum([self.scores_[c] for c in self.contigs_]))
            print('New average score: %f' % (sum([self.scores_[c] for c in self.contigs_]) / len(self.contigs_)))
            print()

        return len(candidates) > 0

    # perform probabilistic extension, link-based mode
    def extend_probabilistic_link_based(self, ag, gcm, max_score, verbose):
        if verbose:
            print('Searching for next extension(s)...')

        extendable_contigs = self.get_extendable_contigs()
        if verbose and len(extendable_contigs) == 0:
            print('\tNo extendable contigs found.')

        candidates = [(contig, cand_contig_name, cand_orientation, cand_direction, cand_score) for \
                      contig in extendable_contigs for cand_contig_name, cand_orientation, cand_direction, cand_score in \
                      self.pick_by_depth_and_genes_link_based(ag, gcm, contig, self.orientation(contig), self.extendability(contig), verbose)]
        candidates = [c for c in candidates if c[1] is not None and c[4] <= max_score
                      and sum(self.lengths_.values()) + ag.length(c[1]) <= self.max_length_]

        if verbose:
            print('%i candidates satisfy maximum score %f.' % (len(candidates), max_score))

        if len(candidates) > 0:
            total_depth = sum([ag.depth(c[1]) for c in candidates])
            weights = [ag.depth(c[1]) / total_depth for c in candidates]
            last_sum = 0
            intervals = []
            for i in range(0, len(candidates)):
                intervals.append((last_sum, last_sum + weights[i], i))
                last_sum += weights[i]
            rnd = random.random()
            best_option = [index for left, right, index in intervals if left <= rnd < right][0]
            candidates = [candidates[best_option]]

        if verbose and len(candidates) > 0:
            print('Performing extension probabilistically...')
        for contig, cand_contig_name, cand_orientation, direction, score in candidates:
            cand_contig_id = len(self.contigs_)
            self.contigs_[cand_contig_id] = cand_contig_name
            occs = [id for id in self.contigs_ if self.contigs_[id] == cand_contig_name]
            share = ag.depth(cand_contig_name) / len(occs)
            for id in occs:
                self.depths_[id] = share

            if verbose:
                print('\tFrom %s (id = %i) %s to %s%s (id = %i, score %f)' % (self.contig_name(contig), contig, direction, cand_contig_name, cand_orientation, cand_contig_id, score))
            self.scores_[cand_contig_id] = score
            self.orientations_[cand_contig_id] = cand_orientation
            self.lengths_[cand_contig_id] = ag.length(cand_contig_name)
            self.gc_contents_[cand_contig_id] = ag.gc_content(cand_contig_name)

            if len(gcm.genes_of(cand_contig_name)) > 0:
                self.intermediate_contigs_[cand_contig_id] = 0
                self.intermediate_nt_[cand_contig_id] = 0
            else:
                self.intermediate_contigs_[cand_contig_id] = self.intermediate_contigs_[contig] + 1
                self.intermediate_nt_[cand_contig_id] = self.intermediate_nt_[contig] + ag.length(cand_contig_name)

            if direction == 'left':
                if contig in self.left_extensions_:
                    self.left_extensions_[contig].append(cand_contig_id)
                else:
                    self.left_extensions_[contig] = [cand_contig_id]
                if cand_contig_id in self.right_extensions_:
                    self.right_extensions_[cand_contig_id].append(contig)
                else:
                    self.right_extensions_[cand_contig_id] = [contig]
            else:
                if contig in self.right_extensions_:
                    self.right_extensions_[contig].append(cand_contig_id)
                else:
                    self.right_extensions_[contig] = [cand_contig_id]
                if cand_contig_id in self.left_extensions_:
                    self.left_extensions_[cand_contig_id].append(contig)
                else:
                    self.left_extensions_[cand_contig_id] = [contig]

        if verbose:
            print('New total score: %f' % sum([self.scores_[c] for c in self.contigs_]))
            print('New average score: %f' % (sum([self.scores_[c] for c in self.contigs_]) / len(self.contigs_)))
            print()

        return len(candidates) > 0

    # return list of contig ids in the order in which they form the plasmid
    # works only properly if each contig has at most one left and one right extension, extracts just one path for degenerated plasmids
    def get_contig_id_chain(self):
        cur = 0
        visited = set()
        while cur in self.left_extensions_ and self.left_extensions_[cur][0] not in visited:
            visited.add(cur)
            cur = self.left_extensions_[cur][0]

        res = [(cur, self.orientation(cur))]
        visited = {cur}
        while cur in self.right_extensions_ and self.right_extensions_[cur][0] not in visited:
            cur = self.right_extensions_[cur][0]
            res.append((cur, self.orientation(cur)))
            visited.add(cur)

        return res

    # return list of names of contigs in the order in which they form the plasmid
    # works only properly if each contig has at most one left and one right extension, extracts just one path for degenerated plasmids
    def get_contig_name_chain(self):
        return [(self.contig_name(id), orientation) for id, orientation in self.get_contig_id_chain()]

    # check whether this plasmid consists of a circular contig chain
    # works only properly if each contig has at most one left and one right extension, not suitable for degenerated plasmids
    def is_circular_contig_chain(self):
        cur = 0
        is_circular = False
        while not is_circular and cur in self.right_extensions_:
            next = self.right_extensions_[cur][0]
            is_circular = (next == 0)
            cur = next

        return is_circular

    # contigs can contribute to several plasmids by splitting their total read depth;
    # therefore, at the end of the contruction of each plasmid, the contribution of each contig is limited to the mean / median
    # read depth of the plasmid
    def finalise_depths(self):
        overall_depth = self.overall_depth()
        for c in self.contigs_:
            self.depths_[c] = min(self.depths_[c], overall_depth)

    # how potential extensions are scored
    # currently, the deviation from a reference read depth (e.g. current mean read depth of plasmid) and the gene density are used
    # NOTE: when changing this function, check whether the default value of parameter -w and determine_weights(...) have to be adapted
    def score_depth_and_gene_density(self, depth, gc_content, depth_partner, density_partner, gc_content_partner):
        return self.score_weights_['depth_diff'] * abs(1 - depth_partner / depth) \
               + self.score_weights_['gene_density'] * (1 - density_partner) \
               + self.score_weights_['gc_diff'] * abs(gc_content - gc_content_partner)

    # for node-based extensions; check whether contig is new or can be extended to the left (i.e. does not already have left extensions)
    def is_new_or_left_extendable(self, contig_name):
        occs = [id for id in self.contigs_ if self.contigs_[id] == contig_name]
        if len(occs) > 1:
            print('WARNING: Contig %s should not occur more than once per plasmid in node-based mode!' % contig_name)
        return len(occs) == 0 or occs[0] not in self.left_extensions_

    # for node-based extensions; check whether contig is new or can be extended to the right (i.e. does not already have right extensions)
    def is_new_or_right_extendable(self, contig_name):
        occs = [id for id in self.contigs_ if self.contigs_[id] == contig_name]
        if len(occs) > 1:
            print('WARNING: Contig %s should not occur more than once per plasmid in node-based mode!' % contig_name)
        return len(occs) == 0 or occs[0] not in self.right_extensions_

    # check whether the proposed link between to contigs has already been used in this plasmid (in either direction)
    # (backward = using reverse complements of contigs in specified orientations and going from to_contig to from_contig)
    def link_is_unused(self, from_contig, from_ori, direction, to_contig, to_ori):
        if to_contig not in self.contigs_.values():
            return True

        if direction == 'left':
            pred_contig = to_contig
            pred_ori = to_ori
            succ_contig = from_contig
            succ_ori = from_ori
        else:
            pred_contig = from_contig
            pred_ori = from_ori
            succ_contig = to_contig
            succ_ori = to_ori

        occs_pred = [id for id in self.contigs_ if self.contigs_[id] == pred_contig]
        neg_pred_ori = '-' if pred_ori == '+' else '+'
        neg_succ_ori = '-' if succ_ori == '+' else '+'

        forward_partners = [(self.contig_name(id), self.orientation(id)) for pred_id in occs_pred for id in (self.right_extensions_[pred_id] if pred_id in self.right_extensions_ else []) if self.orientation(pred_id) == pred_ori]
        backward_partners = [(self.contig_name(id), self.orientation(id)) for pred_id in occs_pred for id in (self.left_extensions_[pred_id] if pred_id in self.left_extensions_ else []) if self.orientation(pred_id) == neg_pred_ori]

        return (succ_contig, succ_ori) not in forward_partners and (succ_contig, neg_succ_ori) not in backward_partners

    # how potential extensions are determined
    # currently, all contigs that have sufficient read depth left and that are not to far away from gene-containing contigs
    # are considered as candidates
    def pick_by_depth_and_genes_node_based(self, ag, gcm, contig, orientation, direction, verbose):
        contig_name = self.contigs_[contig]

        left_candidates = []
        if direction in ['left', 'both']:
            preds = ag.predecessors_of_pos(contig_name) if orientation == '+' else ag.predecessors_of_neg(contig_name)
            for p, o in preds:
                if self.is_new_or_right_extendable(p) and (ag.depth(p) >= self.min_read_depth_) and (abs(self.overall_gc_content() - ag.gc_content(p)) <= self.max_gc_diff_) and \
                        (self.nt_to_genes(contig) + ag.length(p) <= self.max_intermediate_nt_ or self.contigs_to_genes(contig) < self.max_intermediate_contigs_ or gcm.contains_contig(p)):
                    left_candidates.append((p, o, 'left', self.score_depth_and_gene_density(self.overall_depth(), self.overall_gc_content(), ag.depth(p), gcm.gene_density(p), ag.gc_content(p))))

        right_candidates = []
        if direction in ['right', 'both']:
            succs = ag.successors_of_pos(contig_name) if orientation == '+' else ag.successors_of_neg(contig_name)
            for s, o in succs:
                if self.is_new_or_left_extendable(s) and (ag.depth(s) >= self.min_read_depth_) and (abs(self.overall_gc_content() - ag.gc_content(s)) <= self.max_gc_diff_) and \
                        (self.nt_to_genes(contig) + ag.length(s) <= self.max_intermediate_nt_ or self.contigs_to_genes(contig) < self.max_intermediate_contigs_ or gcm.contains_contig(s)):
                    right_candidates.append((s, o, 'right', self.score_depth_and_gene_density(self.overall_depth(), self.overall_gc_content(), ag.depth(s), gcm.gene_density(s), ag.gc_content(s))))

        candidates = left_candidates + right_candidates
        if len(candidates) > 0:
            candidates.sort(key = lambda c: c[3])

        if verbose:
            print('\tExtensions from %s (id = %i):' % (contig_name, contig))
            print('\tCandidates for left extension (%i):' % len(left_candidates), left_candidates)
            print('\tCandidates for right extension (%i):' % len(right_candidates), right_candidates)
            if len(left_candidates) > 0:
                sorted_lefts = [c for c in candidates if c[2] == 'left']
                print('\tBest left extension: %s%s, score %f' % (sorted_lefts[0][0], sorted_lefts[0][1], sorted_lefts[0][3]))
            else:
                print('\tBest left extension: none')
            if len(right_candidates) > 0:
                sorted_rights = [c for c in candidates if c[2] == 'right']
                print('\tBest right extension: %s%s, score %f' % (sorted_rights[0][0], sorted_rights[0][1], sorted_rights[0][3]))
            else:
                print('\tBest right extension: none')

        if len(candidates) == 0:
            return [(None, None, 'none', math.inf)]
        else:
            return candidates

    # how potential extensions are determined
    # currently, all contigs that have sufficient read depth left and that are not to far away from gene-containing contigs
    # are considered as candidates
    def pick_by_depth_and_genes_link_based(self, ag, gcm, contig, orientation, direction, verbose):
        contig_name = self.contigs_[contig]

        left_candidates = []
        if direction in ['left', 'both']:
            preds = ag.predecessors_of_pos(contig_name) if orientation == '+' else ag.predecessors_of_neg(contig_name)
            for p, o in preds:
                num = len([id for id in self.contigs_ if self.contigs_[id] == p]) + 1
                if self.link_is_unused(contig_name, orientation, 'left', p, o) and (ag.depth(p) / num >= self.min_read_depth_) and (abs(self.overall_gc_content() - ag.gc_content(p)) <= self.max_gc_diff_) and \
                        (self.nt_to_genes(contig) + ag.length(p) <= self.max_intermediate_nt_ or self.contigs_to_genes(contig) < self.max_intermediate_contigs_ or gcm.contains_contig(p)):
                    left_candidates.append((p, o, 'left', self.score_depth_and_gene_density(self.overall_depth(), self.overall_gc_content(), ag.depth(p) / num, gcm.gene_density(p), ag.gc_content(p))))

        right_candidates = []
        if direction in ['right', 'both']:
            succs = ag.successors_of_pos(contig_name) if orientation == '+' else ag.successors_of_neg(contig_name)
            for s, o in succs:
                num = len([id for id in self.contigs_ if self.contigs_[id] == s]) + 1
                if self.link_is_unused(contig_name, orientation, 'right', s, o) and (ag.depth(s) / num >= self.min_read_depth_) and (abs(self.overall_gc_content() - ag.gc_content(s)) <= self.max_gc_diff_) and \
                        (self.nt_to_genes(contig) + ag.length(s) <= self.max_intermediate_nt_ or self.contigs_to_genes(contig) < self.max_intermediate_contigs_ or gcm.contains_contig(s)):
                    right_candidates.append((s, o, 'right', self.score_depth_and_gene_density(self.overall_depth(), self.overall_gc_content(), ag.depth(s) / num, gcm.gene_density(s), ag.gc_content(s))))

        candidates = left_candidates + right_candidates
        if len(candidates) > 0:
            candidates.sort(key = lambda c: c[3])

        if verbose:
            print('\tExtensions from %s (id = %i):' % (contig_name, contig))
            print('\tCandidates for left extension (%i):' % len(left_candidates), left_candidates)
            print('\tCandidates for right extension (%i):' % len(right_candidates), right_candidates)
            if len(left_candidates) > 0:
                sorted_lefts = [c for c in candidates if c[2] == 'left']
                print('\tBest left extension: %s%s, score %f' % (sorted_lefts[0][0], sorted_lefts[0][1], sorted_lefts[0][3]))
            else:
                print('\tBest left extension: none')
            if len(right_candidates) > 0:
                sorted_rights = [c for c in candidates if c[2] == 'right']
                print('\tBest right extension: %s%s, score %f' % (sorted_rights[0][0], sorted_rights[0][1], sorted_rights[0][3]))
            else:
                print('\tBest right extension: none')

        if len(candidates) == 0:
            return [(None, None, 'none', math.inf)]
        else:
            return candidates

    # check whether the contig chain ends in the same contig in the same orientation
    # if this is the case, remove one of the duplicates and connect the two ends to obtain a circular contig chain
    # works only properly if each contig has at most one left and one right extension
    def merge_ends(self):
        # single-contig plasmid
        if not self.node_based_ and len(self.contigs_) == 2 and self.contig_name(0) == self.contig_name(1) and self.orientation(0) == self.orientation(1):
            del self.contigs_[1]
            del self.orientations_[1]
            self.depths_[0] += self.depths_[1]
            del self.depths_[1]
            del self.lengths_[1]
            del self.intermediate_contigs_[1]
            del self.intermediate_nt_[1]
            del self.scores_[1]
            self.left_extensions_[0] = [0]
            if 1 in self.left_extensions_:
                del self.left_extensions_[1]
            self.right_extensions_[0] = [0]
            if 1 in self.right_extensions_:
                del self.right_extensions_[1]

        # multi-contig plasmid whose outermost contigs are the same in the same orientation
        if not self.node_based_ and len(self.contigs_) > 2:
            left_end = [c for c in self.contigs_ if c not in self.left_extensions_][0]
            right_end = [c for c in self.contigs_ if c not in self.right_extensions_][0]

            if self.contig_name(left_end) == self.contig_name(right_end) and self.orientation(left_end) == self.orientation(right_end):
                kept_end = left_end if self.scores_[left_end] < self.scores_[right_end] else right_end
                removed_end = right_end if self.scores_[left_end] < self.scores_[right_end] else left_end
                contig_name = self.contig_name(kept_end)

                remaining_occs = [id for id in self.contigs_ if self.contigs_[id] == contig_name and id != removed_end]
                del self.contigs_[removed_end]
                total_depth = self.depth(kept_end) * (len(remaining_occs) + 1)

                del self.orientations_[removed_end]

                for id in remaining_occs:
                    self.depths_[id] = total_depth / len(remaining_occs)
                del self.depths_[removed_end]

                del self.lengths_[removed_end]

                self.intermediate_contigs_[kept_end] = min(self.intermediate_contigs_[kept_end], self.intermediate_contigs_[removed_end])
                del self.intermediate_contigs_[removed_end]

                self.intermediate_nt_[kept_end] = min(self.intermediate_nt_[kept_end], self.intermediate_nt_[removed_end])
                del self.intermediate_nt_[removed_end]

                del self.scores_[removed_end]

                if removed_end == left_end:
                    self.left_extensions_[self.right_extensions_[left_end][0]] = [right_end]
                    self.right_extensions_[right_end] = self.right_extensions_[left_end]
                else:
                    self.left_extensions_[left_end] = self.left_extensions_[right_end]
                    self.right_extensions_[self.left_extensions_[right_end][0]] = [left_end]

                if removed_end in self.left_extensions_:
                    del self.left_extensions_[removed_end]
                if removed_end in self.right_extensions_:
                    del self.right_extensions_[removed_end]


# Representation of seeds satisfying 'quality' criteria in order of preference
#
# Enumerates the seeds in order of decreasing quality: non-branching seeds are assumed better that branching seeds
# and in both groups the seeds are sorted by decreasing gene density and difference in GC content to overall
# GC content of assembly graph.
class SeedEnumerator:
    __slots__ = 'seed_stack_', \
                'depths_', \
                'depth_threshold_'

    def __init__(self, ag, gcm, depth_threshold, gene_density_threshold, max_length):
        seeds = [c for c in gcm.contigs() if
                 ag.depth(c) >= depth_threshold  # do not use contigs with very low depth (possibly a residue or sequencing fragment)
                 and gcm.gene_density(c) >= gene_density_threshold  # do not use contigs which barely contain genes as seeds
                 and ag.length(c) <= max_length]  # do not use overly long contigs as seeds (more likely to be chromosomal)
        overall_gc = ag.overall_gc_content()

        if len(seeds) > 0:
            total_length = sum([ag.length(s) for s in seeds])

            mean_dens = sum([gcm.gene_density(s) * ag.length(s) for s in seeds]) / total_length
            sd_dens = math.sqrt(sum([(gcm.gene_density(s) - mean_dens) ** 2 * ag.length(s) for s in seeds]) / total_length)

            mean_gc_diff = sum([abs(ag.gc_content(s) - overall_gc) * ag.length(s) for s in seeds]) / total_length
            sd_gc_diff = math.sqrt(sum([(abs(ag.gc_content(s) - overall_gc) - mean_gc_diff) ** 2 * ag.length(s) for s in seeds]) / total_length)

            if sd_dens == 0:
                sd_dens = 1
            if sd_gc_diff == 0:
                sd_gc_diff = 1

            # use normalised values to balance the influence of the two components (otherwise, the usually much higher gene density values dominate)
            seeds.sort(key = lambda s: (gcm.gene_density(s) - mean_dens) / sd_dens + (abs(ag.gc_content(s) - overall_gc) - mean_gc_diff) / sd_gc_diff)

        unambiguous_seeds = [s for s in seeds if len(ag.predecessors_of_pos(s)) <= 1 and len(ag.successors_of_pos(s)) <= 1]
        seeds = [s for s in seeds if s not in unambiguous_seeds]

        self.seed_stack_ = seeds + unambiguous_seeds # top of stack = right end of list
        self.depths_ = dict([(s, ag.depth(s)) for s in self.seed_stack_]) # depths of elements in seed_stack_
        self.depth_threshold_ = depth_threshold

    # return whether a contig is (still) a seed in the enumerator
    def contains(self, contig):
        return contig in self.seed_stack_

    # update the depth value of a contig in the enumerator (and remove it if necessary)
    def update_depth(self, contig, value, overwrite = False):
        # update values
        if contig in self.depths_:
            self.depths_[contig] = value if overwrite else (self.depths_[contig] + value)

        # remove contigs which have dropped below threshold
        self.seed_stack_ = [s for s in self.seed_stack_ if self.depths_[s] > self.depth_threshold_]

    # return currently best seed
    def get_next(self):
        return self.seed_stack_[-1]

    # check whether there is still a seed in the enumerator
    def has_next(self):
        return len(self.seed_stack_) > 0

    # return list of all seeds in the enumerator
    def get_all(self):
        return self.seed_stack_


# compute length of longest overlap between suffix and prefix of a sequence
def len_suffix_prefix_overlap(seq):
    overlaps = [0] * len(seq) # overlaps[i] = length of longest proper suffix-prefix overlap for S_{0, i - 1}
    len_overlap = 0 # length of previous longest proper suffix-prefix overlap

    i = 1
    while i < len(seq):
        if seq[i] == seq[len_overlap]: # current longest overlap can be extended
            len_overlap += 1
            overlaps[i] = len_overlap
            i += 1

        elif len_overlap != 0:
            len_overlap = overlaps[len_overlap - 1] # switch to shorter overlap, try to extend this in next iteration

        else:
            overlaps[i] = 0 # no extension possible with current letter, go to next one
            i += 1

    return overlaps[len(seq) - 1]


# determine reverse complement of nucleotide sequence
def reverse_complement(seq):
    return str(Seq(seq).reverse_complement())


# generate plasmid sequence from list of contigs with orientation
def generate_plasmid_sequence(ag, contig_list):
    seq = ''
    prev_contig = None
    prev_ori = None
    for contig, orientation in contig_list:
        if prev_contig is None:
            seq += ag.sequence(contig) if orientation == '+' else reverse_complement(ag.sequence(contig))
        else:
            ov_len = ag.overlap_length(prev_contig, prev_ori, contig, orientation)
            seq += (ag.sequence(contig) if orientation == '+' else reverse_complement(ag.sequence(contig)))[ov_len:]

        prev_contig = contig
        prev_ori = orientation

    return seq

### Helper methods for computing the characteristics of linear plasmids whose contigs overlap ###

# get depth and length values (while considering the overlaps)
def get_depth_length(ag, plasmid, contig_id_list):
    depth_per_contig = []
    length_per_contig = []
    prev_contig = None
    prev_ori = None
    for contig_id, orientation in contig_id_list:
        if prev_contig is None:
            depth_per_contig.append(plasmid.depth(contig_id))
            length_per_contig.append(plasmid.length(contig_id))
        else:
            ov_len = ag.overlap_length(prev_contig, prev_ori, plasmid.contig_name(contig_id), orientation)
            depth_per_contig.append(plasmid.depth(contig_id))
            length_per_contig.append(plasmid.length(contig_id) - ov_len)

        prev_contig = plasmid.contig_name(contig_id)
        prev_ori = orientation

    return depth_per_contig, length_per_contig


def mean_depth(depths, lengths):
    total_length = sum(lengths)
    return sum([depths[i] * lengths[i] for i in range(0, len(depths))]) / (total_length if total_length > 0 else 1)


def median_depth(depths, lengths):
    dep_len = sorted([(depths[i], lengths[i]) for i in range(0, len(depths))], key = lambda x: x[0])
    total_num_values = sum(lengths)
    index = total_num_values * 0.50
    fraction50 = index - math.floor(index)  # fraction for linear interpolation (median)
    rank50 = math.floor(index)  # rank of interest (median)

    count_sum = 0
    median_val = -1
    for i in range(len(dep_len)):
        count_sum += dep_len[i][1]

        if count_sum >= rank50:
            left = dep_len[i][0]
            right = dep_len[i][0] if count_sum >= (rank50 + 1) else dep_len[i + 1][0]
            median_val = left + (right - left) * fraction50
            break

    return median_val


# compute the number of gene-covered bases
# bases in an overlap that are covered by genes in both contigs still contribute only once
def num_gene_covered_bases(ag, gcm, contig_list, seq_length):
    offset = 0
    all_intervals = []
    prev_contig = None
    prev_ori = None
    for contig, orientation in contig_list:
        if prev_contig is None:
            all_intervals += gcm.get_gene_intervals(contig, orientation)
        else:
            offset += ag.length(prev_contig) - ag.overlap_length(prev_contig, prev_ori, contig, orientation)
            all_intervals += [(start + offset, end + offset) for start, end in gcm.get_gene_intervals(contig, orientation)]

        prev_contig = contig
        prev_ori = orientation

    # wrap around intervals of circular plasmids
    intervals = [(start, end) for start, end in all_intervals if end < seq_length]
    to_be_split = [(start, end) for start, end in all_intervals if start < seq_length and end >= seq_length]
    completely_wrapped = [(start, end) for start, end in all_intervals if start >= seq_length]

    for start, end in to_be_split:
        intervals.append((start, seq_length - 1))
        intervals.append((0, end % seq_length))

    for start, end in completely_wrapped:
        intervals.append((start % seq_length, end % seq_length))

    intervals.sort(key = lambda x: x[0])  # intervals is now sorted by start position

    # compute covered bases
    num_pos_covered = 0
    last_pos_covered = 0  # last (right-most) position of contig covered so far
    for start, end in intervals:
        if end <= last_pos_covered:
            pass  # contained in previous interval -> no new position covered
        else:
            num_pos_covered += end - max(last_pos_covered + 1, start) + 1
            last_pos_covered = end

    return num_pos_covered


# generate all plasmid information (while considering overlaps between contigs)
def generate_plasmid_infos(ag, gcm, plasmid_collection, self_overlaps, use_median):
    seqs = dict()
    average_read_depths = dict()
    gene_densities = dict()
    gc_contents = dict()

    for id, plasmid in plasmid_collection:
        contig_chain = plasmid.get_contig_name_chain()
        contig_id_chain = plasmid.get_contig_id_chain()

        depth_per_contig, length_per_contig = get_depth_length(ag, plasmid, contig_id_chain)

        if plasmid.is_circular_contig_chain():
            ov_len = ag.overlap_length(contig_chain[0][0], contig_chain[0][1], contig_chain[-1][0], contig_chain[-1][1])
            seqs[id] = generate_plasmid_sequence(ag, contig_chain)[:-ov_len]

            length_per_contig[-1] -= ov_len

        elif id in self_overlaps:
            seqs[id] = generate_plasmid_sequence(ag, contig_chain)[:-self_overlaps[id]]

            length_per_contig[-1] -= self_overlaps[id]

        else:
            seqs[id] = generate_plasmid_sequence(ag, contig_chain)

        gc_contents[id] = sum([1 for c in seqs[id] if c in ['g', 'G', 'c', 'C']]) / len(seqs[id])
        average_read_depths[id] = median_depth(depth_per_contig, length_per_contig) if use_median else mean_depth(depth_per_contig, length_per_contig)
        gene_densities[id] = num_gene_covered_bases(ag, gcm, contig_chain, len(seqs[id])) / len(seqs[id])

    return seqs, average_read_depths, gene_densities, gc_contents

### END of helper methods

# determine the weights of the components of the scoring function used to find plasmid extensions
# weights are data-dependent based on the seeds (the higher the variation (standard deviation) of an attribute, the higher the weight)
def determine_weights(argument, ag, gcm, seed_enumerator):
    seeds = seed_enumerator.get_all()
    weights = dict()
    components = argument.split(',')
    for comp in components:
        tokens = comp.split('=')
        if tokens[1] != '?':
            weights[tokens[0]] = float(tokens[1])
        else:
            if tokens[0] == 'depth_diff':
                mean_val = (sum([ag.depth(s) * ag.length(s) for s in seeds]) / sum([ag.length(s) for s in seeds])) if len(seeds) > 1 else 0.0
                sd_val = math.sqrt(sum([(ag.depth(s) - mean_val) ** 2 * ag.length(s) for s in seeds]) / sum([ag.length(s) for s in seeds])) if len(seeds) > 1 else 0.0

            elif tokens[0] == 'gene_density':
                mean_val = (sum([gcm.gene_density(s) * ag.length(s) for s in seeds]) / sum([ag.length(s) for s in seeds])) if len(seeds) > 1 else 0.0
                sd_val = math.sqrt(sum([(gcm.gene_density(s) - mean_val) ** 2 * ag.length(s) for s in seeds]) / sum([ag.length(s) for s in seeds])) if len(seeds) > 1 else 0.0

            elif tokens[0] == 'gc_diff':
                mean_val = (sum([ag.gc_content(s) * ag.length(s) for s in seeds]) / sum([ag.length(s) for s in seeds])) if len(seeds) > 1 else 0.0
                sd_val = math.sqrt(sum([(ag.gc_content(s) - mean_val) ** 2 * ag.length(s) for s in seeds]) / sum([ag.length(s) for s in seeds])) if len(seeds) > 1 else 0.0

            else:
                print('WARNING: Unknown weight component (%s). Check for typos and the needed components.' % tokens[0])
                sd_val = 0

            weights[tokens[0]] = sd_val

    return weights


# check the quality of plasmids
# plasmids with a gene density or read deapth that is too low and plasmids which are too long or too short are marked as questionable
# plasmids whose underlying contig collection is a subset of another plasmid are marked as questionable too
def check_plasmids(plasmids, plasmid_seqs, plasmid_gene_densities, plasmid_read_depths, min_gene_density, min_read_depth, min_length, max_length, keep_subplasmids, verbose):
    putative = []
    questionable = []

    for id, plasmid in plasmids:
        plasmid_seq = plasmid_seqs[id]
        gene_density = plasmid_gene_densities[id]
        read_depth = plasmid_read_depths[id]

        if gene_density >= min_gene_density and read_depth >= min_read_depth and min_length <= len(plasmid_seq) <= max_length:
            putative.append((id, plasmid))
        else:
            if verbose:
                if gene_density < min_gene_density:
                    print('Low gene density (< %f) detected for %s. Plasmid is marked as questionable.' % (min_gene_density, id))
                if read_depth < min_read_depth:
                    print('Low mean read depth (< %f) detected for %s. Plasmid is marked as questionable.' % (min_read_depth, id))
                if len(plasmid_seq) < min_length:
                    print('Plasmid %s is shorter than %i nt. It is marked as questionable.' % (id, min_length))
                if max_length < len(plasmid_seq):
                    print('Plasmid %s is longer than %i nt. It is marked as questionable.' % (id, min_length))

            questionable.append((id, plasmid))

    # plasmid A is a subplasmid of plasmid B if the set of contigs underlying A is a subset of the set of contigs underlying B (not multisets)
    if not keep_subplasmids:
        plasmid_collection = dict([(id, set(plasmid.get_contig_names())) for id, plasmid in putative])
        subplasmid_ids = []
        for id in plasmid_collection:
            superplasmids = [super_id for super_id in plasmid_collection if id != super_id and super_id not in subplasmid_ids # first condition: keep largest in chain of subplasmid relations, keep one of several equal collections
                             and len(plasmid_collection[id]) <= len(plasmid_collection[super_id])
                             and plasmid_collection[id].issubset(plasmid_collection[super_id])]
            if len(superplasmids) > 0:
                if verbose:
                    print('Plasmid %s is contained in %i plasmid(s), e.g. %s. It is marked as questionable.' % (id, len(superplasmids), superplasmids[0]))
                subplasmid_ids.append(id)

        questionable = questionable + [elem for elem in putative if elem[0] in subplasmid_ids]
        putative = [elem for elem in putative if elem[0] not in subplasmid_ids]

    return putative, questionable


# check quality of contig collections
# contig collections with a gene density or read deapth that is too low and plasmids which are too long or too short are marked as questionable
# contig collections which are a subset of another contig collection are marked as questionable too
# note: does not consider overlaps
def check_contig_collections(ag, gcm, contig_collections, min_gene_density, min_read_depth, min_length, max_length, keep_subcollection, verbose):
    putative = []
    questionable = []

    for id, plasmid in contig_collections:
        contigs = plasmid.get_contig_names()

        total_length = sum([ag.length(c) for c in contigs])
        gene_density = sum([gcm.num_gene_covered_bases(c) for c in contigs]) / total_length

        if gene_density >= min_gene_density and plasmid.overall_depth() >= min_read_depth and min_length <= total_length <= max_length:
            putative.append((id, plasmid))
        else:
            if verbose:
                if gene_density < min_gene_density:
                    print('Low gene density (< %f) detected for %s. Plasmid is marked as questionable.' % (min_gene_density, id))
                if plasmid.overall_depth() < min_read_depth:
                    print('Low mean read depth (< %f) detected for %s. Plasmid is marked as questionable.' % (min_read_depth, id))
                if total_length < min_length:
                    print('Plasmid %s is shorter than %i nt. It is marked as questionable.' % (id, min_length))
                if max_length < total_length:
                    print('Plasmid %s is longer than %i nt. It is marked as questionable.' % (id, min_length))

            questionable.append((id, plasmid))

    # plasmid A is a subplasmid of plasmid B if the set of contigs underlying A is a subset of the set of contigs underlying B (not multisets)
    if not keep_subcollection:
        plasmid_collection = dict([(id, set(plasmid.get_contig_names())) for id, plasmid in putative])
        subplasmid_ids = []
        for id in plasmid_collection:
            superplasmids = [super_id for super_id in plasmid_collection if plasmid_collection[id].issubset(plasmid_collection[super_id]) and \
                             id != super_id and super_id not in subplasmid_ids] # last condition: keep largest in chain of subplasmid relations, keep one of several equal collections
            if len(superplasmids) > 0:
                if verbose:
                    print('Plasmid %s is contained in %i plasmid(s), e.g. %s. It is marked as questionable.' % (id, len(superplasmids), superplasmids[0]))
                subplasmid_ids.append(id)

        questionable = questionable + [elem for elem in putative if elem[0] in subplasmid_ids]
        putative = [elem for elem in putative if elem[0] not in subplasmid_ids]

    return putative, questionable


# tries to find plasmids that actually belong together based on their read depth
# a bin is represented by a list of plasmid identifiers (without an orientation)
# does not concatenate the plasmids assumed to belong together
def bin_plasmids(ag, plasmids, plasmid_seqs, plasmid_gene_densities, plasmid_read_depths, plasmid_gc_contents, median_read_depth, threshold_factor):
    indices = list(range(0, len(plasmids)))
    overall_gc_content = ag.overall_gc_content()

    bins = []

    if len(plasmids) > 0:
        total_length = sum([len(plasmid_seqs[id]) for id, _ in plasmids])

        mean_depth_ratio = sum([plasmid_read_depths[id] / median_read_depth * len(plasmid_seqs[id]) for id, _ in plasmids]) / total_length
        sd_depth_ratio = math.sqrt(sum([(plasmid_read_depths[id] / median_read_depth - mean_depth_ratio) ** 2 * len(plasmid_seqs[id]) for id, _ in plasmids]) / total_length)

        mean_gene_density = sum([plasmid_gene_densities[id] * len(plasmid_seqs[id]) for id, _ in plasmids]) / total_length
        sd_gene_density = math.sqrt(sum([(plasmid_gene_densities[id] - mean_gene_density) ** 2 * len(plasmid_seqs[id]) for id, _ in plasmids]) / total_length)

        mean_gc_diff = sum([abs(plasmid_gc_contents[id] - overall_gc_content) * len(plasmid_seqs[id]) for id, _ in plasmids]) / total_length
        sd_gc_diff = math.sqrt(sum([(abs(plasmid_gc_contents[id] - overall_gc_content) - mean_gc_diff) ** 2 * len(plasmid_seqs[id]) for id, _ in plasmids]) / total_length)

        if sd_depth_ratio == 0:
            sd_depth_ratio = 1
        if sd_gene_density == 0:
            sd_gene_density = 1
        if sd_gc_diff == 0:
            sd_gc_diff = 1

        indices.sort(key = lambda i: ((plasmid_read_depths[plasmids[i][0]] / median_read_depth) - mean_depth_ratio) / sd_depth_ratio +
                                     (plasmid_gene_densities[plasmids[i][0]] - mean_gene_density) / sd_gene_density +
                                     (abs(plasmid_gc_contents[plasmids[i][0]] - overall_gc_content) - mean_gc_diff) / sd_gc_diff,
                     reverse = True)

        while len(indices) > 0:
            rep = indices[0]
            rep_depth = plasmid_read_depths[plasmids[rep][0]]
            rep_gc = plasmid_gc_contents[plasmids[rep][0]]

            # circular plasmids form their own bin
            if plasmids[rep][1].is_circular_contig_chain():
                bin = [rep]

            # add all plasmids not yet binned that fall within threshold_factor standard deviations of the current plasmid w.r.t. their read depth
            else:
                bin = [ind for ind in indices if abs(rep_depth - plasmid_read_depths[plasmids[ind][0]]) <= threshold_factor * sd_depth_ratio and abs(rep_gc - plasmid_gc_contents[plasmids[ind][0]]) <= threshold_factor * sd_gc_diff]

            bins.append(bin)
            indices = [ind for ind in indices if ind not in bin]

    bins = [[plasmids[i][0] for i in bin] for bin in bins]

    return bins


# tries to find contig collections (= plasmids with branching) that actually belong together based on their read depth
# a bin is represented by a list of contig-collection identifiers (without an orientation)
# note: does not consider overlaps
def bin_contig_collections(ag, gcm, contig_collections, median_read_depth, threshold_factor):
    indices = list(range(0, len(contig_collections)))
    plasmid_depths = dict([(id, plasmid.overall_depth()) for id, plasmid in contig_collections])
    plasmid_lengths = dict([(id, sum([ag.length(contig) for contig in plasmid.get_contig_names()])) for id, plasmid in contig_collections])
    plasmid_gene_densities = dict([(id, sum([gcm.num_gene_covered_bases(c) for c in plasmid.get_contig_names()]) / plasmid_lengths[id]) for id, plasmid in contig_collections])
    plasmid_gc_contents = dict([(id, sum([ag.gc_content(c) * ag.length(c) for c in plasmid.get_contig_names()]) / plasmid_lengths[id]) for id, plasmid in contig_collections])
    overall_gc_content = ag.overall_gc_content()

    bins = []

    if len(contig_collections) > 0:
        total_length = sum([plasmid_lengths[id] for id, _ in contig_collections])

        mean_depth_ratio = sum([plasmid_depths[id] / median_read_depth * plasmid_lengths[id] for id, _ in contig_collections]) / total_length
        sd_depth_ratio = math.sqrt(sum([(plasmid_depths[id] / median_read_depth - mean_depth_ratio) ** 2 * plasmid_lengths[id] for id, _ in contig_collections]) / total_length)

        mean_gene_density = sum([plasmid_gene_densities[id] * plasmid_lengths[id] for id, _ in contig_collections]) / total_length
        sd_gene_density = math.sqrt(sum([(plasmid_gene_densities[id] - mean_gene_density) ** 2 * plasmid_lengths[id] for id, _ in contig_collections]) / total_length)

        mean_gc_diff = sum([abs(plasmid_gc_contents[id] - overall_gc_content) * plasmid_lengths[id] for id, _ in contig_collections]) / total_length
        sd_gc_diff = math.sqrt(sum([(abs(plasmid_gc_contents[id] - overall_gc_content) - mean_gc_diff) ** 2 * plasmid_lengths[id] for id, _ in contig_collections]) / total_length)

        if sd_depth_ratio == 0:
            sd_depth_ratio = 1
        if sd_gene_density == 0:
            sd_gene_density = 1
        if sd_gc_diff == 0:
            sd_gc_diff = 1

        indices.sort(key = lambda i: ((plasmid_depths[contig_collections[i][0]] / median_read_depth) - mean_depth_ratio) / sd_depth_ratio +
                                     (plasmid_gene_densities[contig_collections[i][0]] - mean_gene_density) / sd_gene_density +
                                     (abs(plasmid_gc_contents[contig_collections[i][0]] - overall_gc_content) - mean_gc_diff) / sd_gc_diff,
                     reverse = True)

        while len(indices) > 0:
            rep = indices[0]
            rep_depth = plasmid_depths[contig_collections[rep][0]]
            rep_gc = plasmid_gc_contents[contig_collections[rep][0]]

            # add all plasmids not yet binned that fall within threshold_factor standard deviations of the current plasmid w.r.t. their read depth
            bin = [ind for ind in indices if abs(rep_depth - plasmid_depths[contig_collections[ind][0]]) <= threshold_factor * sd_depth_ratio and abs(rep_gc - plasmid_gc_contents[contig_collections[ind][0]]) <= threshold_factor * sd_gc_diff]
            bins.append(bin)
            indices = [ind for ind in indices if ind not in bin]

    bins = [[contig_collections[i][0] for i in bin] for bin in bins]

    return bins


# output contig chains to file (one plasmid per row)
# format: <plasmid id>;<comma-separated list of orientated contigs (= id followed by + or -)>
def output_chains(output_file, plasmids):
    with open(output_file, 'w') as out:
        for id_num, plasmid in plasmids:
            contig_chain = plasmid.get_contig_name_chain()
            out.write(id_num + ';' + ','.join([contig + orientation for contig, orientation in contig_chain]) + '\n')


# output contig collections to file (one plasmid per row)
# used for degenerated plasmids
# format: <plasmid id>;<comma-separated list of orientated contigs (= id followed by + or -)>
def output_collections(output_file, plasmids):
    with open(output_file, 'w') as out:
        for id_num, plasmid in plasmids:
            out.write(id_num + ';' + ','.join([plasmid.contig_name(contig) + plasmid.orientation(contig) for contig in plasmid.get_contig_ids()]) + '\n')


# output plasmid sequences in FASTA format (defline contains additional information, e.g. gene density)
def output_plasmids(ag, gcm, output_file, plasmids, plasmid_seqs, plasmid_gene_densities, plasmid_read_depths, plasmid_gc_contents, self_overlaps, use_median = False):
    with open(output_file, 'w') as out:
        for id, plasmid in plasmids:
            contig_chain = plasmid.get_contig_name_chain()
            contigs = [elem[0] for elem in contig_chain]

            plasmid_seq = plasmid_seqs[id]
            gene_density = plasmid_gene_densities[id]
            gc_content = plasmid_gc_contents[id]
            read_depth = plasmid_read_depths[id]
            out.write('>%s seed_contig=%s\tlength=%i\t%s_read_depth=%f\tgene_density=%f\tnum_cds=%i\tgc_content=%f\tcircular=%i\n%s\n' % (
                id, plasmid.seed(), len(plasmid_seq), 'median' if use_median else 'mean', read_depth, gene_density, gcm.num_gene_locations_of_set(contigs), gc_content, plasmid.is_circular_contig_chain() or id in self_overlaps, plasmid_seq))


# output the sequences of the contigs underlying the plasmids in FASTA format
def output_plasmid_contigs(ag, gcm, output_file, plasmids):
    with open(output_file, 'w') as out:
        for id, plasmid in plasmids:
            for contig in plasmid.get_contig_ids():
                contig_name = plasmid.contig_name(contig)
                out.write('>%s|%i_%s\n%s\n' % (contig_name, contig, id, ag.sequence(contig_name) if plasmid.orientation(contig) == '+' else reverse_complement(ag.sequence(contig_name))))


# output assembly graph in GFA format and add colour labels to the contigs indicating
# whether they belong to (questionable) plasmids
def output_tagged_gfa(ag, gcm, putative_plasmids, questionable_plasmids, input_gfa, output_gfa):
    gene_contigs = gcm.contigs()
    contig_labels = dict([(s, set()) for s in ag.segments()])
    putative_plasmid_contigs = set()
    questionable_plasmid_contigs = set()
    for name, plasmid in putative_plasmids:
        id = name.split('_')[-1]
        for contig in plasmid.get_contig_names():
            contig_labels[contig].add(id)
            putative_plasmid_contigs.add(contig)
    for name, plasmid in questionable_plasmids:
        id = name.split('_')[-1]
        for contig in plasmid.get_contig_names():
            contig_labels[contig].add(id)
            questionable_plasmid_contigs.add(contig)

    sep = '\t'

    # reduce graph to given contigs
    with open(input_gfa) as infile, open(output_gfa, 'w') as out:
        line = next(infile)

        while (line):
            if line.startswith('H'):  # header: keep
                out.write(line)

            if line.startswith('S'):  # segment: keep if segment is one of the plasmid contigs
                tokens = line.split(sep)
                if tokens[1] in putative_plasmid_contigs:
                    line = line.strip('\n').strip('\r')
                    line += sep + 'LB:z:' + ('*' if tokens[1] in gene_contigs else '') + ','.join(contig_labels[tokens[1]])
                    line += sep + 'CL:z:blue' + '\n'
                elif tokens[1] in questionable_plasmid_contigs:
                    line = line.strip('\n').strip('\r')
                    line += sep + 'LB:z:' + ('*' if tokens[1] in gene_contigs else '') + ','.join(contig_labels[tokens[1]])
                    line += sep + 'CL:z:lightblue' + '\n'
                out.write(line)

            if line.startswith('L'):  # link: keep if both segments of the link are plasmid contigs
                tokens = line.split(sep)
                out.write(line)

            if line.startswith('C'):  # containment: keep if both segments are plasmid contigs
                tokens = line.split(sep)
                out.write(line)

            if line.startswith('P'):  # path: keep if all segments are plasmid contigs
                tokens = line.split(sep)

                path_name = tokens[1]
                segments = [seg[:-1] for seg in tokens[2].split(',')]

                out.write(line)

            line = next(infile, None)


# output plasmids bins (comma-separated list of plasmid ids for each bin, one bin per line)
def output_bins(bins, output_file):
    with open(output_file, 'w') as out:
        for bin in bins:
            out.write(','.join(bin) + '\n')


# print configuration of greedy algorithm
def show_config(assembly_graph, gene_contig_mapping, output_dir, genes_file, min_gene_density,
                min_seed_gene_density, min_length, max_length, min_read_depth, min_plasmid_read_depth, max_gc_diff,
                max_intermediate_contigs, max_intermediate_nt, max_score, score_weights, keep_subplasmids, overlap_ends,
                binning, fanout, probabilistic, use_median, use_node_based, verbose):
    print('#########################################')
    print('### Configuration of greedy algorithm ###\n')

    print('>>> Input / output')
    print('Assembly graph: %s' % assembly_graph)
    print('Gene-contig mapping: %s' % gene_contig_mapping)
    print('Genes database: %s' % genes_file)
    print('Output directory: %s' % output_dir)

    print('\n>>> Thresholds')
    print('Minimum gene density: %f' % min_gene_density)
    print('Minimum seed gene density: %f' % min_seed_gene_density)
    print('Minimum plasmid length: %i' % min_length)
    print('Maximum plasmid length: %i' % max_length)
    print('Minimum contig read depth: %f' % min_read_depth)
    print('Minimum plasmid read depth: %f' % min_plasmid_read_depth)
    print('Maximum GC difference: %f' % max_gc_diff)
    print('Maximum number of intermediate contigs: %i' % max_intermediate_contigs)
    print('Maximum number of intermediate nucleotides: %i' % max_intermediate_nt)
    print('Maximum extension score: %f' % max_score)
    print('Score weights: %s' % score_weights)
    print('Minimum overlap for circularisation: %f' % overlap_ends)

    print('\n>>> Modes and other options')
    print('Probabilistic: %i' % probabilistic)
    print('Use median: %i' % use_median)
    print('Node-based: %i' % use_node_based)
    print('Keep subplasmids: %i' % keep_subplasmids)
    print('Binning: %s' % ('deactivated' if math.isnan(binning) else ('activated with factor %f' % binning)))
    print('Fanout: %i' % fanout)
    print('Verbose: %i' % verbose)

    print('#########################################\n')


# main method greedy algorithm
def greedy(assembly_graph, genes_file, gene_contig_mapping, output_dir,
           min_gene_density = DEF_MIN_GENE_DENSITY, min_seed_gene_density = DEF_MIN_SEED_GENE_DENSITY,
           min_length = DEF_MIN_LENGTH, max_length = DEF_MAX_LENGTH,
           min_read_depth = DEF_MIN_READ_DEPTH, min_plasmid_read_depth = DEF_MIN_PLASMID_READ_DEPTH,
           max_gc_diff = DEF_MAX_GC_DIFF, max_intermediate_contigs = DEF_MAX_INTERMEDIATE_CONTIGS, max_intermediate_nt = DEF_MAX_INTERMEDIATE_NT,
           max_score = DEF_MAX_SCORE, score_weights = DEF_SCORE_WEIGHTS,
           keep_subplasmids = DEF_KEEP_SUBPLASMIDS, overlap_ends = DEF_OVERLAP_ENDS,
           binning = DEF_BINNING, fanout = DEF_FANOUT, probabilistic = DEF_PROBABILISTIC,
           use_median = DEF_USE_MEDIAN, use_node_based = DEF_USE_NODE_BASED, verbose = DEF_VERBOSE):

    if probabilistic and fanout != 1:
        print('WARNING: Probabilistic mode resets fanout to 1!')
        fanout = 1

    ## ----- Initialise data structures from inputs and set up parameters -----

    ag = ag_module.AssemblyGraph(assembly_graph)
    ag_median_read_depth =  ag.median_depth()
    gcm = gcm_module.GeneContigMapping(gene_contig_mapping, genes_file, assembly_graph)
    contig_sequences = dict()
    with open(assembly_graph, 'r') as in_file:
        for line in in_file:
            if line.startswith('S'):
                tokens = line.split('\t')
                name = tokens[1]
                seq = str(tokens[2])
                contig_sequences[name] = seq

    # check consistency of input data
    inconsistent_genes = gcm.check_consistency()
    if len(inconsistent_genes) > 0:
        raise ValueError('ERROR: The used gene-contig mapping and gene database are not consistent '
                         + '(%i genes are in the mapping but not in the database). ' % len(inconsistent_genes)
                         + 'Please check whether mapping %s was created using database %s.'
                         % (gene_contig_mapping, genes_file))
    inconsistent_contigs = set(gcm.contigs()).difference(contig_sequences)
    if len(inconsistent_contigs) > 0:
        raise ValueError('ERROR: The used gene-contig mapping and assembly graph are not consistent '
                         + '(%i contigs are in the mapping but not in the assembly graph). ' % len(inconsistent_contigs)
                         + 'Please check whether mapping %s was created using assembly graph %s.'
                         % (gene_contig_mapping, assembly_graph))

    # default value of min_seed_gene_density is 50 % higher than min_gene_density (of plasmid)
    if min_seed_gene_density is math.nan:
        min_seed_gene_density = 1.5 * min_gene_density

    # default value for min_plasmid_read_depth is intended to be similar to read depth of chromosomes
    # set to median of largest n contigs, which are assumed to be chromosomal for this computation due to their length
    if min_plasmid_read_depth is math.nan:
        p = 0.4
        n = len(ag.segments())
        min_plasmid_read_depth = p * ag.median_depth(sorted(ag.segments(), key = lambda s: ag.length(s), reverse = True)[:n])

    # min_read_depth per contig is less strict in order to take some fluctuation into consideration
    if min_read_depth is math.nan:
        min_read_depth = 0.75 * min_plasmid_read_depth


    # determine contig seeds and establish order of preference
    seed_enum = SeedEnumerator(ag, gcm, min_read_depth, min_seed_gene_density, max_length)

    # set weights of components of scoring function
    weights = determine_weights(score_weights, ag, gcm, seed_enum)

    show_config(assembly_graph, gene_contig_mapping, output_dir, genes_file, min_gene_density,
                min_seed_gene_density, min_length, max_length, min_read_depth, min_plasmid_read_depth, max_gc_diff,
                max_intermediate_contigs, max_intermediate_nt, max_score, score_weights, keep_subplasmids, overlap_ends,
                binning, fanout, probabilistic, use_median, use_node_based, verbose)


    ## ----- Determine new plasmid while there is a seed available -----

    plasmid_collection = []
    id_num = 0

    while (seed_enum.has_next()):
        if verbose:
            print('\n--- plasmid_%i ---' % id_num)

        # get seed and initialise new plasmid
        seed = seed_enum.get_next()
        plasmid = Plasmid(seed, ag.depth(seed), ag.length(seed), ag.gc_content(seed), gcm.gene_density(seed), min_read_depth,
                          max_length, max_gc_diff, max_intermediate_contigs, max_intermediate_nt, weights,
                          fan_out = fanout, use_median = use_median, node_based = use_node_based)
        if verbose:
            print('Seed: %s\n' % seed)

        # try to extend plasmid at its ends until no more suitable contig is found
        cont = True
        while (cont):
            if probabilistic:
                cont = plasmid.extend_probabilistic(ag, gcm, max_score, verbose)
            else:
                cont = plasmid.extend(ag, gcm, max_score, verbose)

        # if contig chain starts and ends in the same contig with the same orientation, remove one duplicate
        if fanout <= 1:
            plasmid.merge_ends()

        # peel off plasmid
        plasmid.finalise_depths()
        for c in plasmid.get_contig_ids():
            ag.update_depth(plasmid.contig_name(c), -plasmid.depth(c))
            seed_enum.update_depth(plasmid.contig_name(c), -plasmid.depth(c))

        plasmid_collection.append(('plasmid_' + str(id_num), plasmid))
        id_num += 1


    ## ----- Statistical analyses and postprocessing -----

    # overall statistics for putative chromosomal contigs
    all_plasmid_contigs = [contig for id, plasmid in plasmid_collection for contig in plasmid.get_contig_names()]
    chromosomal_contigs = set(ag.segments()).difference(all_plasmid_contigs)
    if verbose:
        print('\nStatistics of putative chromosomal contigs:')
        print(' - number: %i' % len(chromosomal_contigs))
        print(' - total length: %i' % sum([ag.length(c) for c in chromosomal_contigs]))
        print(' - mean read depth: %f' % ag.mean_depth(chromosomal_contigs))
        print(' - median read depth: %f' % ag.median_depth(chromosomal_contigs))
        print(' - GC content: %f' % ag.overall_gc_content(chromosomal_contigs))
        print()

    # mark plasmid as circular when there is sufficient suffix-prefix overlap or when the contig chain is a cycle
    if fanout <= 1:
        if verbose:
            print('Searching for circular plasmids...')
        self_overlaps = dict()
        circular_found = False
        for id, plasmid in plasmid_collection:
            if plasmid.is_circular_contig_chain():
                if verbose:
                    print('Plasmid %s consists of a circular contig chain.' % id)
                circular_found = True
            else:
                len_overlap = len_suffix_prefix_overlap(generate_plasmid_sequence(ag, plasmid.get_contig_name_chain()))

                if len_overlap >= overlap_ends:
                    self_overlaps[id] = len_overlap
                    if verbose:
                        print('Plasmid %s has sufficient suffix-prefix overlap (%i nt).' % (id, len_overlap))
                    circular_found = True
        if verbose and not circular_found:
            print('No circular plasmid found.')
        print()

        # generate plasmid infos (considering the circularity information)
        plasmid_seqs, average_read_depths, gene_densities, gc_contents = generate_plasmid_infos(ag, gcm, plasmid_collection, self_overlaps, use_median)

    # filter putative plasmids
    if verbose:
        print('Inspecting putative plasmids...')
    putative_plasmids, questionable_plasmids = \
        check_plasmids(plasmid_collection, plasmid_seqs, gene_densities, average_read_depths, min_gene_density, min_plasmid_read_depth, min_length, max_length, keep_subplasmids, verbose) if fanout <= 1 \
        else check_contig_collections(ag, gcm, plasmid_collection, min_gene_density, min_plasmid_read_depth, min_length, max_length, keep_subplasmids, verbose)
    if verbose and len(questionable_plasmids) == 0:
        print('No putative plasmid was deemed questionable.')

    # plasmid binning
    if binning is not math.nan:
        plasmid_bins_putative = bin_plasmids(ag, putative_plasmids, plasmid_seqs, gene_densities, average_read_depths, gc_contents, ag_median_read_depth, binning) if fanout <= 1 else bin_contig_collections(ag, gcm, putative_plasmids, ag_median_read_depth, binning)
        plasmid_bins_all = bin_plasmids(ag, plasmid_collection, plasmid_seqs, gene_densities, average_read_depths, gc_contents, ag_median_read_depth, binning) if fanout <= 1 else bin_contig_collections(ag, gcm, plasmid_collection, ag_median_read_depth, binning)

    # ----- Output generation -----

    file_putative_contig_lists = 'putative_plasmid_contigs_list.csv' # contig lists of all putative plasmids
    file_questionable_contig_lists = 'questionable_plasmid_contigs_list.csv' # contig lists of all questionable plasmids
    file_chains = 'contig_chains.csv' # contig chains of all plasmids
    file_putative = 'putative_plasmids.fasta' # sequences of all putative plasmids
    file_putative_contigs = 'putative_plasmid_contigs.fasta' # sequences of contigs of all putative plasmids
    file_questionable = 'questionable_plasmids.fasta' # sequences of all questionable plasmids
    file_questionable_contigs = 'questionable_plasmid_contigs.fasta' # sequences of contigs of all questionable plasmids
    file_binned_putative = 'plasmid_bins_putative.csv' # bins of putative 'plasmids' potentially belonging to each other
    file_binned_all = 'plasmid_bins_all.csv' # bins of all 'plasmids' potentially belonging to each other
    file_tagged_gfa = 'tagged_assembly.gfa' # plasmid contigs are labelled with plasmid ids and coloured

    if fanout <= 1:
        output_chains(os.path.join(output_dir, file_chains), plasmid_collection)
        output_plasmids(ag, gcm, os.path.join(output_dir, file_putative), putative_plasmids, plasmid_seqs, gene_densities, average_read_depths, gc_contents, self_overlaps, use_median)
        output_plasmids(ag, gcm, os.path.join(output_dir, file_questionable), questionable_plasmids, plasmid_seqs, gene_densities, average_read_depths, gc_contents, self_overlaps, use_median)
    else:
        output_collections(os.path.join(output_dir, file_putative_contig_lists), putative_plasmids)
        output_collections(os.path.join(output_dir, file_questionable_contig_lists), questionable_plasmids)

    output_plasmid_contigs(ag, gcm, os.path.join(output_dir, file_putative_contigs), putative_plasmids)
    output_plasmid_contigs(ag, gcm, os.path.join(output_dir, file_questionable_contigs), questionable_plasmids)
    output_tagged_gfa(ag, gcm, putative_plasmids, questionable_plasmids, assembly_graph, os.path.join(output_dir, file_tagged_gfa))

    if binning is not math.nan:
        output_bins(plasmid_bins_putative, os.path.join(output_dir, file_binned_putative))
        output_bins(plasmid_bins_all, os.path.join(output_dir, file_binned_all))
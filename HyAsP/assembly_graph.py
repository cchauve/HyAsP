#!/usr/bin/env python

import math


# Representation of an assembly graph with analytical and navigational capabilities.
# The connected components of the graph as well as the neighbours and characteristics (read depth, length, GC content)
# of individual contigs can be determined.
# Read-depth values of contigs can be changed.
# The mean / median read depth and the overall GC content of a collection of contigs can be computed.
#
# The constructor assumes an assembly graph in GFA format (as it is created by Unicycler).
class AssemblyGraph:
    __slots__ = 'version_', \
                'sequences_', \
                'lengths_', \
                'read_counts_', \
                'frag_counts_', \
                'kmer_counts_', \
                'depths_', \
                'gc_contents_', \
                'predecessors_', \
                'predecessors_rc_', \
                'successors_', \
                'successors_rc_', \
                'overlaps_', \
                'map_qual_', \
                'num_diffs_', \
                'read_counts_links_', \
                'frag_counts_links_', \
                'kmer_counts_links_', \
                'containers_', \
                'containers_rc_', \
                'containees_', \
                'containees_rc_', \
                'containments_', \
                'read_coverage_', \
                'num_diffs_containments_', \
                'path_segments_', \
                'path_overlaps_', \
                'connected_components_'

    def __init__(self, gfa_path):
        sep = '\t'

        version_ = 'unspecified'

        # nucleotide sequences of contigs / segments
        self.sequences_ = {}

        # segments - optional fields (specification)
        self.lengths_ = {}
        self.read_counts_ = {}
        self.frag_counts_ = {}
        self.kmer_counts_ = {}

        # segments - optional fields (Unicycler)
        self.depths_ = {}

        # segments - additional information
        self.gc_contents_ = {}

        # links - connectivity information
        self.predecessors_ = {}
        self.predecessors_rc_ = {}
        self.successors_ = {}
        self.successors_rc_ = {}
        self.overlaps_ = {}

        # links - optional fields (specification)
        self.map_qual_ = {}
        self.num_diffs_ = {}
        self.read_counts_links_ = {}
        self.frag_counts_links_ = {}
        self.kmer_counts_links_ = {}

        # containments - main / required information
        self.containers_ = {}  # for each segment 'seg', segments containing seg+
        self.containers_rc_ = {}  # for each segment 'seg', segments containing seg-
        self.containees_ = {}  # for each segment 'seg', segments contained in seg+
        self.containees_rc_ = {}  # for each segment 'seg', segments contained in seg-
        self.containments_ = {}  # 'pos' and 'overlap' fields of containment lines

        # containments - optional fields (specification)
        self.read_coverage_ = {}
        self.num_diffs_containments_ = {}

        # paths - main / required information
        self.path_segments_ = {}
        self.path_overlaps_ = {}

        with open(gfa_path, 'r') as assembly_file:
            line = next(assembly_file)

            while (line):
                if line.startswith('#'):  # comment: skip
                    continue

                if line.startswith('H'):  # header: extract version number
                    self.version_ = line.split(sep)[1]

                if line.startswith('S'):  # segment: associate segment name with sequence length and read depth
                    tokens = line.split(sep)

                    name = tokens[1]
                    seq = tokens[2]

                    # sequence
                    self.sequences_[name] = seq

                    # segment length
                    tmp = [field for field in tokens[3:] if field.startswith('LN:')]
                    ln = -1  # -1 means information not available
                    if tmp:
                        ln = int(tmp[0].split(':')[2])
                    if ln == -1 and seq != '*':
                        ln = len(seq)

                    # read count
                    tmp = [field for field in tokens[3:] if field.startswith('RC:')]
                    rc = -1  # -1 means information not available
                    if tmp:
                        rc = int(tmp[0].split(':')[2])

                    # fragment count
                    tmp = [field for field in tokens[3:] if field.startswith('FC:')]
                    fc = -1  # -1 means information not available
                    if tmp:
                        fc = int(tmp[0].split(':')[2])

                    # k-mer count
                    tmp = [field for field in tokens[3:] if field.startswith('KC:')]
                    kc = -1  # -1 means information not available
                    if tmp:
                        kc = int(tmp[0].split(':')[2])

                    # read depth
                    tmp = [field for field in tokens[3:] if field.startswith('dp:')]
                    dp = -1  # depth of sequence, -1 means information not available
                    if tmp:
                        dp = float(tmp[0].split(':')[2])

                    self.lengths_[name] = ln
                    self.read_counts_[name] = rc
                    self.frag_counts_[name] = fc
                    self.kmer_counts_[name] = kc
                    self.depths_[name] = dp

                    self.successors_[name] = []
                    self.successors_rc_[name] = []
                    self.predecessors_[name] = []
                    self.predecessors_rc_[name] = []
                    self.containers_[name] = []
                    self.containers_rc_[name] = []
                    self.containees_[name] = []
                    self.containees_rc_[name] = []

                    # GC content
                    gc = sum([1 for c in seq if c in ['g', 'G', 'c', 'C']])
                    self.gc_contents_[name] = (gc / len(seq)) if len(seq) > 0 else 0

                if line.startswith('L'):  # link: extract connectivity information for each segment (distinguishing as-is and reverse complement)
                    tokens = line.split(sep)

                    from_seg = tokens[1]
                    from_ori = tokens[2]
                    to_seg = tokens[3]
                    to_ori = tokens[4]
                    overlap = tokens[5].rstrip()

                    if from_ori == '+':
                        self.successors_[from_seg].append((to_seg, to_ori))
                    else:
                        self.successors_rc_[from_seg].append((to_seg, to_ori))

                    if to_ori == '+':
                        self.predecessors_[to_seg].append((from_seg, from_ori))
                    else:
                        self.predecessors_rc_[to_seg].append((from_seg, from_ori))

                    # BEGIN - link between reverse complements (e.g. if A+ -> B-, then also B+ -> A-)
                    if from_ori == '+':
                        self.predecessors_rc_[from_seg].append((to_seg, '-' if to_ori == '+' else '+'))
                    else:
                        self.predecessors_[from_seg].append((to_seg, '-' if to_ori == '+' else '+'))

                    if to_ori == '-':
                        self.successors_[to_seg].append((from_seg, '-' if from_ori == '+' else '+'))
                    else:
                        self.successors_rc_[to_seg].append((from_seg, '-' if from_ori == '+' else '+'))
                    # END - link between reverse complements

                    self.overlaps_[(from_seg, to_seg)] = overlap

                    # mapping quality
                    tmp = [field for field in tokens[6:] if field.startswith('MQ:')]
                    mq = -1  # -1 means information not available
                    if tmp:
                        mq = int(tmp[0].split(':')[2])

                    # number of mismatches / gaps
                    tmp = [field for field in tokens[6:] if field.startswith('NM:')]
                    nm = -1  # -1 means information not available
                    if tmp:
                        nm = int(tmp[0].split(':')[2])

                    # read count
                    tmp = [field for field in tokens[6:] if field.startswith('RC:')]
                    rc = -1  # -1 means information not available
                    if tmp:
                        rc = int(tmp[0].split(':')[2])

                    # fragment count
                    tmp = [field for field in tokens[6:] if field.startswith('FC:')]
                    fc = -1  # -1 means information not available
                    if tmp:
                        fc = int(tmp[0].split(':')[2])

                    # k-mer count
                    tmp = [field for field in tokens[6:] if field.startswith('KC:')]
                    kc = -1  # -1 means information not available
                    if tmp:
                        kc = int(tmp[0].split(':')[2])

                    self.map_qual_[(from_seg, to_seg)] = mq
                    self.num_diffs_[(from_seg, to_seg)] = nm
                    self.read_counts_links_[(from_seg, to_seg)] = rc
                    self.frag_counts_links_[(from_seg, to_seg)] = fc
                    self.kmer_counts_links_[(from_seg, to_seg)] = kc

                if line.startswith('C'):  # containment
                    tokens = line.split(sep)

                    container_seg = tokens[1]
                    container_ori = tokens[2]
                    containee_seg = tokens[3]
                    containee_ori = tokens[4]
                    pos = int(tokens[5])
                    overlap = tokens[6]

                    if containee_ori == '+':
                        self.containers_[containee_seg].append((container_seg, container_ori))
                    else:
                        self.containers_rc_[containee_seg].append((container_seg, container_ori))

                    if container_ori == '+':
                        self.containees_[container_seg].append((containee_seg, containee_ori))
                    else:
                        self.containees_rc_[container_seg].append((containee_seg, containee_ori))

                        self.containments_[(container_seg, containee_seg)] = (pos, overlap)

                    # read coverage
                    tmp = [field for field in tokens[6:] if field.startswith('RC:')]
                    rc = -1  # -1 means information not available
                    if tmp:
                        rc = int(tmp[0].split(':')[2])

                    # number of mismatches / gaps
                    tmp = [field for field in tokens[6:] if field.startswith('NM:')]
                    nm = -1  # -1 means information not available
                    if tmp:
                        nm = int(tmp[0].split(':')[2])

                    self.read_coverage_[(container_seg, containee_seg)] = rc
                    self.num_diffs_containments_[(container_seg, containee_seg)] = nm

                if line.startswith('P'):  # path
                    tokens = line.split(sep)

                    path_name = tokens[1]
                    segments = [(seg[:-1], seg[-1]) for seg in tokens[2].split(',')]
                    overlap = tokens[3].split(',')

                    self.path_segments_[path_name] = segments
                    self.path_overlaps_[path_name] = overlap

                line = next(assembly_file, None)

        self.connected_components_ = self.connected_components()

    # list of segments in assembly graph
    def segments(self):
        return sorted(self.lengths_.keys())

    # list of segments (with orientation) succeeding seg+
    def successors_of_pos(self, seg):
        return self.successors_[seg]

    # list of segments (with orientation) succeeding seg- (reverse complement)
    def successors_of_neg(self, seg):
        return self.successors_rc_[seg]

    # list of segments (with orientation) preceding seg+
    def predecessors_of_pos(self, seg):
        return self.predecessors_[seg]

    # list of segments (with orientation) preceding seg- (reverse complement)
    def predecessors_of_neg(self, seg):
        return self.predecessors_rc_[seg]

    # list of segments (with orientation) succeeding seg (+ or -)
    def successors(self, seg):
        return list(set(self.successors_of_pos(seg) + self.successors_of_neg(seg)))

    # list of segments (with orientation) preceding seg (+ or -)
    def predecessors(self, seg):
        return list(set(self.predecessors_of_pos(seg) + self.predecessors_of_neg(seg)))

    # list of segments (with orientation) preceding or succeeding seg (+ or -)
    # note that a successor is also a predecessor (considering reverse complements)
    def neighbours(self, seg):
        return list(set(self.successors(seg) + self.predecessors(seg)))

    # count number of incoming and outgoing edges
    #  - edge connecting segment to itself contributes 2
    #  - multiple edges with the same target are counted separately
    def degree(self, seg):
        return (len(self.predecessors_[seg]) + len(self.predecessors_rc_[seg]) + len(self.successors_[seg]) + len(self.successors_rc_[seg])) / 2

    # determine connected components
    def connected_components(self):
        num_cc = 0
        visited = dict([(s, -1) for s in self.segments()])
        queue = []

        for s in self.segments():
            if visited[s] == -1:
                queue.append(s)

                while queue:
                    cur = queue.pop(0)

                    if visited[cur] == -1:
                        visited[cur] = num_cc

                        neighbours = [n[0] for n in self.predecessors_[cur]]
                        neighbours.extend([n[0] for n in self.predecessors_rc_[cur]])
                        neighbours.extend([n[0] for n in self.successors_[cur]])
                        neighbours.extend([n[0] for n in self.successors_rc_[cur]])

                        queue.extend(neighbours)


                num_cc += 1

        seg_comp_rel = [(k, visited[k]) for k in visited]
        return dict([(cid, [elem[0] for elem in seg_comp_rel if elem[1] == cid]) for cid in range(0, num_cc)])

    # determine connected components in subgraph defined by the given segments
    def connected_components_restricted(self, segs):
        num_cc = 0
        visited = dict([(s, -1) for s in self.segments()])
        queue = []

        for s in self.segments():
            if visited[s] == -1 and s in segs:
                queue.append(s)

                while queue:
                    cur = queue.pop(0)

                    if visited[cur] == -1:
                        visited[cur] = num_cc

                        neighbours = [n[0] for n in self.predecessors_[cur] if n[0] in segs]
                        neighbours.extend([n[0] for n in self.predecessors_rc_[cur] if n[0] in segs])
                        neighbours.extend([n[0] for n in self.successors_[cur] if n[0] in segs])
                        neighbours.extend([n[0] for n in self.successors_rc_[cur] if n[0] in segs])

                        queue.extend(neighbours)


                num_cc += 1

        seg_comp_rel = [(k, visited[k]) for k in visited if k != -1]
        return dict([(cid, [elem[0] for elem in seg_comp_rel if elem[1] == cid]) for cid in range(0, num_cc)])

    # get list of segments of connected component (by integer id)
    def get_connected_component(self, c):
        return self.connected_components_[c]

    # get 'representative' segment (= longest) of connected component (by integer id)
    def rep_of_component(self, c):
        segs = self.connected_components_[c]
        return max(segs, key = lambda s: self.lengths_[s])

    # get integer id of connected component (by segment)
    def component_of_seg(self, seg):
        return [key for key, segs in self.connected_components_.items() if seg in segs][0]

    # get nucleotide sequence of segment
    def sequence(self, seg):
        return self.sequences_[seg]

    # get length of segment
    def length(self, seg):
        return self.lengths_[seg]

    # get read depth of segment
    def depth(self, seg):
        return self.depths_[seg]

    # get GC content of segment
    def gc_content(self, seg):
        return self.gc_contents_[seg]

    # update depth value of segment by setting it to or modifying it by a given value
    def update_depth(self, seg, value, overwrite = False):
        self.depths_[seg] = value if overwrite else (self.depths_[seg] + value)

    # get mean read depth of selection of segments (or all segments if no selection is specified)
    def mean_depth(self, sel = None):
        if sel is None:
            sel = self.segments()
        total_length = sum([self.lengths_[seg] for seg in sel])
        return sum([self.depths_[seg] * self.lengths_[seg] for seg in sel]) / (total_length if total_length > 0 else 1)

    # get median read depth of selection of segments (or all segments if no selection is specified)
    def median_depth(self, sel = None):
        if sel is None:
            sel = self.segments()

        sel = list(sel)
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

    # get GC content of selection of segments (or all segments if no selection is specified)
    def overall_gc_content(self, sel = None):
        if sel is None:
            sel = self.segments()
        total_length = sum([self.lengths_[seg] for seg in sel])
        return sum([self.gc_contents_[seg] * self.lengths_[seg] for seg in sel]) / (total_length if total_length > 0 else 1)
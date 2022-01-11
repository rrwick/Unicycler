"""
Copyright 2017 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Unicycler

This module contains functionality related to minimap, which Unicycler uses to seed long read
alignments.

This file is part of Unicycler. Unicycler is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Unicycler is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Unicycler. If
not, see <http://www.gnu.org/licenses/>.
"""

import os
import sys
from collections import defaultdict
from .misc import get_nice_header, dim, line_iterator, range_overlap, range_is_contained, \
    range_overlap_size, simplify_ranges
from . import log
from . import settings

try:
    from .cpp_wrappers import minimap_align_reads
except AttributeError as e:
    sys.exit('Error when importing C++ library: ' + str(e) + '\n'
             'Have you successfully built the library file using make?')


class MinimapAlignment(object):

    def __init__(self, paf_line=None):
        if paf_line is None:
            self.paf_line = ''
            self.read_name = ''
            self.read_length = 0
            self.read_start = 0
            self.read_end = 0
            self.read_strand = '+'
            self.ref_name = ''
            self.ref_length = 0
            self.ref_start = 0
            self.ref_end = 0
            self.matching_bases = 0
            self.num_bases = 0
            self.minimiser_count = 0
            self.read_end_gap = 0
            self.ref_end_gap = 0

        else:
            self.paf_line = paf_line.strip()
            line_parts = self.paf_line.split('\t')

            self.read_name = line_parts[0]
            self.read_length = int(line_parts[1])
            self.read_start = int(line_parts[2])
            self.read_end = int(line_parts[3])
            self.read_strand = line_parts[4]

            self.ref_name = get_nice_header(line_parts[5])
            self.ref_length = int(line_parts[6])
            self.ref_start = int(line_parts[7])
            self.ref_end = int(line_parts[8])

            self.matching_bases = int(line_parts[9])
            self.num_bases = int(line_parts[10])
            # Mapping quality is part 11, not currently used
            self.minimiser_count = int(line_parts[12].split('cm:i:')[-1])

            self.read_end_gap = self.read_length - self.read_end
            self.ref_end_gap = self.ref_length - self.ref_end


    def get_read_length(self):
        return self.read_length


    def get_concise_string(self):
        return ','.join([str(x) for x in [self.read_start, self.read_end, self.read_strand,
                                          self.ref_name, self.ref_start, self.ref_end]])

    def get_string_for_cpp_scrub(self):
        return '\t'.join([str(x) for x in [self.read_length, self.read_start, self.read_end,
                                           self.ref_length, self.ref_start, self.ref_end]])

    def __repr__(self):
        return str(self.read_start) + '-' + str(self.read_end) + '(' + self.read_strand + '):' + \
            self.ref_name + ':' + str(self.ref_start) + '-' + str(self.ref_end) + \
            '(' + str(self.matching_bases) + '/' + str(self.num_bases) + ',' + \
            str(self.minimiser_count) + ')'

    def get_signed_ref_name(self):
        """
        Returns the ref name preceded by a '-' if the read strand is '-'.
        """
        return ('-' if self.read_strand == '-' else '') + self.ref_name

    def overlaps_reference(self):
        """
        Returns true if the alignment overlaps either end of the reference.
        """
        adjusted_contig_start = self.ref_start - self.read_start
        adjusted_contig_end = self.ref_end + self.read_end_gap
        return adjusted_contig_start < 0 or adjusted_contig_end >= self.ref_length

    def fraction_read_aligned(self):
        try:
            return (self.read_end - self.read_start) / self.read_length
        except ZeroDivisionError:
            return 0.0

    def fraction_ref_aligned(self):
        try:
            return (self.ref_end - self.ref_start) / self.ref_length
        except ZeroDivisionError:
            return 0.0

    def get_start_overhang(self):
        # return min(self.read_start, self.ref_start)
        return self.read_start if self.read_start < self.ref_start else self.ref_start

    def get_end_overhang(self):
        # return min(self.read_end_gap, self.ref_end_gap)
        return self.read_end_gap if self.read_end_gap < self.ref_end_gap else self.ref_end_gap

    def get_total_overhang(self):
        return self.get_start_overhang() + self.get_end_overhang()

    def get_smallest_overhang(self):
        start_overhang = self.get_start_overhang()
        end_overhang = self.get_end_overhang()
        return start_overhang if start_overhang < end_overhang else end_overhang


def load_minimap_alignments_basic(minimap_alignments_str):
    """
    This simple function just loads the minimap alignments in a list of MinimapAlignment objects.
    It doesn't filter, group by reads, or any of the other fancy stuff that the
    load_minimap_alignments function does.
    """
    alignments = []
    for line in line_iterator(minimap_alignments_str):
        alignments.append(MinimapAlignment(line))
    return alignments


def load_minimap_alignments(minimap_alignments_str, filter_by_minimisers=False,
                            minimiser_ratio=10, filter_overlaps=False, allowed_overlap=0):
    """
    Loads minimap's output string into MinimapAlignment objects, grouped by read.
    If filter_by_minimisers is True, it will remove low minimiser count hits.
    If filter_overlaps is True, it will exclude hits which overlap better hits.
    """
    alignments = defaultdict(list)
    for line in line_iterator(minimap_alignments_str):
        try:
            log.log(dim(line), 3)
            alignment = MinimapAlignment(line)
            read_alignments = alignments[alignment.read_name]
            read_alignments.append(alignment)
            read_alignments = sorted(read_alignments, key=lambda x: x.minimiser_count, reverse=True)
            if filter_by_minimisers:
                best_minimiser_count = read_alignments[0].minimiser_count
                min_minimiser_count = best_minimiser_count / minimiser_ratio
                read_alignments = [x for x in read_alignments
                                   if x.minimiser_count >= min_minimiser_count]
            if filter_overlaps:
                kept_alignments = []
                for alignment in read_alignments:
                    if not alignments_overlap(alignment, kept_alignments, allowed_overlap):
                        kept_alignments.append(alignment)
                read_alignments = kept_alignments
            alignments[alignment.read_name] = sorted(read_alignments, key=lambda x: x.read_start)
        except (IndexError, ValueError):
            pass
    return alignments


def alignments_overlap(a, other, allowed_overlap):
    adjusted_start = a.read_start + allowed_overlap
    return any(range_overlap(adjusted_start, a.read_end, x.read_start, x.read_end) for x in other)


def align_long_reads_to_assembly_graph(graph, long_read_filename, working_dir, threads):
    """
    Aligns all long reads to all graph segments and returns a dictionary of alignments (key =
    read name, value = list of MinimapAlignment objects).
    """
    segments_fasta = os.path.join(working_dir, 'all_segments.fasta')
    log.log('Aligning long reads to graph using minimap', 1)
    graph.save_to_fasta(segments_fasta, verbosity=2)
    minimap_alignments_str = minimap_align_reads(segments_fasta, long_read_filename, threads, 3,
                                                 'default')
    minimap_alignments = \
        load_minimap_alignments(minimap_alignments_str, filter_overlaps=True,
                                allowed_overlap=settings.ALLOWED_MINIMAP_OVERLAP,
                                filter_by_minimisers=True,
                                minimiser_ratio=settings.MAX_TO_MIN_MINIMISER_RATIO)
    log.log('Number of minimap alignments: ' + str(len(minimap_alignments)), 2)
    log.log('', 1)
    return minimap_alignments


def build_start_end_overlap_sets(minimap_alignments):
    """
    Build indices of start and end contig overlaps so we can quickly determine which reads
    overlap which end of a particular contig. These dictionaries have a key of a signed segment
    number and a value of a set of read names.
    """
    start_overlap_reads = defaultdict(set)
    end_overlap_reads = defaultdict(set)
    min_overlap_amount = 100
    for read_name, alignments in minimap_alignments.items():
        for a in alignments:
            seg_num = int(a.ref_name)
            if a.read_strand == '+':
                seg_start = a.ref_start
                seg_end = a.ref_end
            else:  # a.read_strand == '-'
                seg_num *= -1
                seg_start = a.ref_length - a.ref_end
                seg_end = a.ref_length - a.ref_start
            adjusted_seg_start = seg_start - a.read_start
            adjusted_seg_end = seg_end + a.read_end_gap
            if adjusted_seg_start < -min_overlap_amount:
                start_overlap_reads[seg_num].add(read_name)
            if adjusted_seg_end > a.ref_length + min_overlap_amount:
                end_overlap_reads[seg_num].add(read_name)
    return start_overlap_reads, end_overlap_reads


def combine_close_hits(alignments, min_read_ref_ratio, max_read_ref_ratio):
    """
    This function takes a list of alignments and it combines alignments if doing so would stay
    within the allowed ratio range.
    """
    # If there's just one alignment, no grouping is necessary.
    if len(alignments) <= 1:
        return alignments

    # Group alignments by read/ref/strand.
    alignment_groups = defaultdict(list)
    for a in alignments:
        alignment_groups[(a.read_name, a.ref_name, a.read_strand)].append(a)

    # For each potentially mergeable group of hits...
    merged_alignments = []
    for read_ref_strand, grouped_alignments in alignment_groups.items():
        _, seg_name, strand = read_ref_strand
        grouped_alignments = sorted(grouped_alignments, key=lambda x: x.read_start)

        current = grouped_alignments[0]
        for i in range(1, len(grouped_alignments)):
            a = grouped_alignments[i]
            a_ref_start = a.ref_start
            a_ref_end = a.ref_end

            # 'read' in this context is the contig were are trying to place.
            potential_merge_read_range = a.read_end - current.read_start
            if strand == '+':
                potential_merge_ref_range = a_ref_end - current.ref_start
            else:  # strand == '-'
                potential_merge_ref_range = current.ref_end - a_ref_start

            merge_would_extend = (current.read_start < a.read_start and
                                  current.read_end < a.read_end)
            read_ref_ratio = abs(potential_merge_read_range) / abs(potential_merge_ref_range)
            ratio_okay = (min_read_ref_ratio <= read_ref_ratio <= max_read_ref_ratio)

            if merge_would_extend and ratio_okay:
                current.read_start = min(current.read_start, a.read_start)
                current.read_end = max(current.read_end, a.read_end)
                current.ref_start = min(current.ref_start, a_ref_start)
                current.ref_end = max(current.ref_end, a_ref_end)
                current.matching_bases += a.matching_bases
                current.minimiser_count += a.minimiser_count
                current.read_end_gap = current.read_length - current.read_end
                current.ref_end_gap = current.ref_length - current.ref_end
                current.num_bases = max(current.ref_end - current.ref_start,
                                        current.read_end - current.read_start)
            # If not, make a new alignment.
            else:
                merged_alignments.append(current)
                current = a
        merged_alignments.append(current)

    return merged_alignments


def remove_conflicting_alignments(alignments, allowed_overlap):
    """
    This function takes alignments for one read, removes conflicting and/or redundant alignments
    and returns the clean-up list of alignments.
    """
    alignments = sorted(alignments, reverse=True,
                        key=lambda x: (x.matching_bases, x.minimiser_count, x.ref_name))
    kept_alignments = []
    kept_alignment_ranges = []
    for a in alignments:
        this_range = (a.read_start, a.read_end)

        # Don't keep alignments for which their part of the read is already aligned.
        if range_is_contained(this_range, kept_alignment_ranges):
            continue

        # Don't keep alignments which overlap too much with existing alignments.
        if range_overlap_size(this_range, kept_alignment_ranges) > allowed_overlap:
            continue

        kept_alignments.append(a)
        kept_alignment_ranges = simplify_ranges(kept_alignment_ranges + [this_range])

    return sorted(kept_alignments, key=lambda x: x.read_start)


def get_opposite_alignment(alignment):
    opposite_alignment = MinimapAlignment()

    opposite_alignment.read_name = alignment.ref_name
    opposite_alignment.read_length = alignment.ref_length
    opposite_alignment.read_start = alignment.ref_start
    opposite_alignment.read_end = alignment.ref_end

    opposite_alignment.ref_name = alignment.read_name
    opposite_alignment.ref_length = alignment.read_length
    opposite_alignment.ref_start = alignment.read_start
    opposite_alignment.ref_end = alignment.read_end

    opposite_alignment.read_strand = alignment.read_strand
    opposite_alignment.matching_bases = alignment.matching_bases
    opposite_alignment.num_bases = alignment.num_bases
    opposite_alignment.minimiser_count = alignment.minimiser_count

    opposite_alignment.read_end_gap = opposite_alignment.read_length - opposite_alignment.read_end
    opposite_alignment.ref_end_gap = opposite_alignment.ref_length - opposite_alignment.ref_end

    return opposite_alignment

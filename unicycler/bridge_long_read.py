"""
Copyright 2017 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Unicycler

Long read bridges are the big important type of bridge for a hybrid Unicycler assembly. They are
made using long reads which align to multiple segments in the graph.

This file is part of Unicycler. Unicycler is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Unicycler is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Unicycler. If
not, see <http://www.gnu.org/licenses/>.
"""

from multiprocessing.dummy import Pool as ThreadPool
import time
import math
import statistics
import sys
from collections import defaultdict
from .bridge_common import get_bridge_str, get_mean_depth, get_depth_agreement_factor, \
    get_bridge_table_parameters, print_bridge_table_header, print_bridge_table_row
from .misc import float_to_str, reverse_complement, flip_number_order, score_function
from . import settings
from .path_finding import get_best_paths_for_seq
from . import log

try:
    from .cpp_wrappers import consensus_alignment
except AttributeError as e:
    sys.exit('Error when importing C++ library: ' + str(e) + '\n'
             'Have you successfully built the library file using make?')


class LongReadBridge(object):
    """
    This class describes a bridge created from long read alignments.
    """
    def __init__(self, graph, start, end):

        # The numbers of the two single copy segments which are being bridged.
        self.start_segment = start
        self.end_segment = end

        # The individual reads contributing to the bridge. The sequences/qualities are not for the
        # entire read, just the part in the bridge. In the case of overlapping alignments, this has
        # the overlap size instead of the sequence.
        self.reads = []

        # The consensus of all read sequences. If there is only one read, then this is the same as
        # that read's sequence. If there are multiple reads, this is hopefully a more accurate
        # sequence. If this is a case where the start and end segments overlap, then there will not
        # be a consensus sequence and consensus_length will be the mean of read_lengths (a negative
        # value).
        self.consensus_sequence = ''

        # The path through the unbridged graph, if one was found.
        self.graph_path = []
        self.all_paths = []

        # The bridge sequence, gotten from the graph path if a good path was found. Otherwise it's
        # from the consensus sequence.
        self.bridge_sequence = ''

        # The bridge depth, a weighted mean of the start and end depths.
        self.depth = get_mean_depth(graph.segments[abs(self.start_segment)],
                                    graph.segments[abs(self.end_segment)], graph)

        # A score used to determine the order of bridge application.
        self.quality = 1.0

        # When a bridge is applied, the segments in the bridge may have their depth reduced
        # accordingly. This member stores which segments have had their depth reduced and by how
        # much due to this bridge's application. It is stored so if this bridge is later deleted,
        # we can restore the depth to the segments.
        self.segments_reduced_depth = []

        self.graph = graph

    def __repr__(self):
        return 'long read bridge: ' + get_bridge_str(self) + \
               ' (quality = ' + float_to_str(self.quality, 2) + ')'

    def predicted_time_to_finalise(self):
        """
        This function very roughly predicts how long the bridge will take to finalise. It's not
        meant to be particularly accurate, but can hopefully be used to roughly order the bridges
        from slow to fast.
        """
        total_seq_length = 0
        seq_count = 0
        for read in self.reads:
            if not isinstance(read[0], int):
                total_seq_length += len(read[0])
                seq_count += 1
        if not seq_count:
            mean_seq_length = 0.0
        else:
            mean_seq_length = total_seq_length / seq_count

        if seq_count > 1:
            predicted_consensus_time = ((1.34e-9 * (total_seq_length ** 2)) +
                                        (2.76e-5 * total_seq_length))
        else:
            predicted_consensus_time = 0.0
        predicted_path_time = ((1.78e-7 * (mean_seq_length ** 2)) + (3.75e-3 * mean_seq_length))

        return predicted_consensus_time + predicted_path_time

    def finalise(self, scoring_scheme, min_alignment_length, read_lengths, estimated_genome_size,
                 expected_linear_seqs):
        """
        Determines the consensus sequence for the bridge, attempts to find it in the graph and
        assigns a quality score to the bridge. This is the big performance-intensive step of long
        read bridging!
        """
        start_seg = self.graph.segments[abs(self.start_segment)]
        end_seg = self.graph.segments[abs(self.end_segment)]

        output = [str(self.start_segment), str(self.end_segment), str(len(self.reads))]

        start_alignment_scaled_scores = [x[2].scaled_score for x in self.reads]
        end_alignment_scaled_scores = [x[3].scaled_score for x in self.reads]
        best_overall_scaled_score = min(max(start_alignment_scaled_scores),
                                        max(end_alignment_scaled_scores))
        alignment_scaled_scores = start_alignment_scaled_scores + end_alignment_scaled_scores
        mean_alignment_scaled_score = statistics.mean(alignment_scaled_scores)
        read_to_ref_ratios = [x[2].get_read_to_ref_ratio() for x in self.reads] + \
                             [x[3].get_read_to_ref_ratio() for x in self.reads]
        mean_read_to_ref_ratio = statistics.mean(read_to_ref_ratios)

        # Partition the full span reads into two groups: those with negative numbers (implying that
        # the two segments overlap) and those with actual sequences.
        reads_without_seq = []
        reads_with_seq = []
        for read in self.reads:
            if isinstance(read[0], int):
                reads_without_seq.append(read)
            else:
                reads_with_seq.append(read)

        # There shouldn't usually be both full spans with sequence and without. If there are some
        # of each, we'll throw out the minority group.
        if reads_with_seq and reads_without_seq:
            if len(reads_without_seq) > len(reads_with_seq):
                reads_with_seq = []
            else:
                reads_without_seq = []

        # For reads with sequence, we perform a MSA and get a consensus sequence.
        if reads_with_seq:

            self.consensus_sequence = get_consensus_sequence(reads_with_seq, scoring_scheme,
                                                             output)

            # We now make an expected scaled score for an alignment between the consensus and a
            # graph path. I.e. when we find a path in the graph for this consensus, this is about
            # how well it should align. This goes up with alignment scores (better alignment scores
            # imply better read sequences) and read count (more reads can give a better consensus).
            #     https://www.desmos.com/calculator
            #     y=100\cdot \left(\left(1-\frac{a}{100}\right)
            #       \left(1-\frac{3}{2+x}\right)+\frac{a}{100}\right)
            #     y = expected_scaled_score, x = num_span_reads, a = mean_alignment_scaled_score
            num_span_reads = len(self.reads)
            expected_scaled_score = 100.0 * ((1.0 - mean_alignment_scaled_score / 100.0) *
                                             (1.0 - (3.0 / (2.0 + num_span_reads))) +
                                             mean_alignment_scaled_score / 100.0)
            expected_scaled_score = max(expected_scaled_score, best_overall_scaled_score)

            # We can also predict the ratio in length between the consensus sequence and the
            # graph path. For a low number of input reads, this is close to the mean ratio for
            # the read alignments, but it approaches 1 as we have more reads and expect our
            # consensus to be better.
            #     https://www.desmos.com/calculator
            #     y=\left(a-1\right)\left(\frac{b}{x+b-1}\right)+1
            expected_consensus_to_ref_ratio = 1.0 + (mean_read_to_ref_ratio - 1.0) * \
                                                    (4 / (4 + num_span_reads - 1))
            target_path_length = int(round((len(self.consensus_sequence) /
                                            expected_consensus_to_ref_ratio)))

        # For reads without sequence, we simply need a mean distance.
        else:
            self.consensus_sequence = ''
            mean_overlap = int(round(sum(x[0] for x in reads_without_seq) / len(reads_without_seq)))
            output.append(str(mean_overlap))
            output.append('')
            target_path_length = 0
            expected_scaled_score = 100.0

        output.append(str(target_path_length))

        path_start_time = time.time()
        self.all_paths, progressive_path_search = \
            get_best_paths_for_seq(self.graph, self.start_segment, self.end_segment,
                                   target_path_length, self.consensus_sequence, scoring_scheme,
                                   expected_scaled_score)
        path_time = time.time() - path_start_time

        output.append(str(len(self.all_paths)))
        output.append('progressive' if progressive_path_search else 'exhaustive')
        output.append(float_to_str(path_time, 1))

        # If paths were found, use a path sequence for the bridge.
        if self.all_paths:
            if self.all_paths[0][0]:
                output.append(', '.join(str(x) for x in self.all_paths[0][0]))
            else:
                output.append('direct connection')
            best_path_len = self.graph.get_bridge_path_length(self.all_paths[0][0])
            output.append(str(best_path_len))
            output.append(float_to_str(self.all_paths[0][1], 1))
            output.append(float_to_str(self.all_paths[0][3], 2))
            output.append(str(self.all_paths[0][2]))

            self.graph_path = self.all_paths[0][0]
            self.bridge_sequence = self.graph.get_path_sequence(self.graph_path)

            # We start this bridge's quality using a function that takes into account the
            # actual, expected and minimum acceptable scores. If the actual scaled score is 100,
            # this function gives 1. If it is the expected value, this function gives a number
            # around 0.7. Then as the actual score approaches the minimum acceptable score, this
            # function approaches 0.
            #     https://www.desmos.com/calculator
            #     y=\left(\frac{1}{1+2^{a-x}}\right)^{0.5}
            #     y = self.quality, x = actual_scaled_score, a = expected_scaled_score,
            #     b = min_acceptable_scaled_score
            actual_scaled_score = self.all_paths[0][3]
            self.quality = math.sqrt(1.0 /
                                     (1.0 + 2.0 ** (expected_scaled_score - actual_scaled_score)))

        # If a path wasn't found, the consensus sequence is the bridge.
        else:
            self.graph_path = []
            output += ['', '', '', '', '']

            if self.consensus_sequence:
                self.bridge_sequence = self.consensus_sequence
            else:
                self.bridge_sequence = ''

            # The quality of non-graph-path-supported bridges depends on the number of dead ends
            # and whether or not linear sequences are expected.
            dead_end_count = 0
            if self.graph.ends_with_dead_end(self.start_segment):
                dead_end_count += 1
            if self.graph.starts_with_dead_end(self.end_segment):
                dead_end_count += 1

            if expected_linear_seqs:
                if dead_end_count == 2:
                    self.quality = settings.PATHLESS_BRIDGE_QUAL_TWO_DEAD_ENDS_WITH_LINEAR_SEQS
                elif dead_end_count == 1:
                    self.quality = settings.PATHLESS_BRIDGE_QUAL_ONE_DEAD_END_WITH_LINEAR_SEQS
                else:  # dead_end_count == 0
                    self.quality = settings.PATHLESS_BRIDGE_QUAL_NO_DEAD_ENDS_WITH_LINEAR_SEQS
            else:
                if dead_end_count == 2:
                    self.quality = settings.PATHLESS_BRIDGE_QUAL_TWO_DEAD_ENDS
                elif dead_end_count == 1:
                    self.quality = settings.PATHLESS_BRIDGE_QUAL_ONE_DEAD_END
                else:  # dead_end_count == 0
                    self.quality = settings.PATHLESS_BRIDGE_QUAL_NO_DEAD_ENDS

            # Bridge length affects quality too: short bridges are better.
            bridge_len = max(0, len(self.bridge_sequence))
            half_qual_len = settings.LONG_READ_BRIDGE_HALF_QUAL_LENGTH
            self.quality *= half_qual_len / (bridge_len + half_qual_len)

        # Expected read count is determined using the read lengths and bridge size. For a given
        # read length and bridge, there are an estimable number of positions where a read of that
        # length would be able to contribute to the bridge. This is used to get the probability
        # that any read would create a bridge, and totalling those up gives us our estimated count.
        min_read_len = (2 * min_alignment_length) + len(self.bridge_sequence)
        total_possible_placements = 0
        for read_len, count in read_lengths.items():
            if read_len < min_read_len:
                continue
            possible_read_placements = read_len - min_read_len + 1
            possible_read_placements *= count
            possible_read_placements *= max(self.depth, 1)
            total_possible_placements += possible_read_placements
        expected_read_count = total_possible_placements / estimated_genome_size
        actual_read_count = len(self.reads)

        # Adjust the expected read count down, especially for higher values.
        # TO DO: reevaluate this step - is it necessary?
        expected_read_count = reduce_expected_count(expected_read_count, 30, 0.5)

        # The start segment and end segment should agree in depth. If they don't, that's very bad,
        # as it implies that they aren't actually single copy or on the same piece of DNA.
        depth_agreement_factor = get_depth_agreement_factor(start_seg.depth, end_seg.depth)
        self.quality *= depth_agreement_factor

        # The number of reads which contribute to a bridge is a big deal, so the read count factor
        # scales linearly. This value is capped at 1, which means that bridges with too few reads
        # are punished but bridges with excess reads are not rewarded.
        try:
            read_count_factor = min(1.0, actual_read_count / expected_read_count)
            self.quality *= read_count_factor
        except ZeroDivisionError:
            pass

        # The length of alignments to the start/end segments is positively correlated with quality
        # to reward bridges with long alignments. Specifically, we want there to be at least one
        # spanning read with a long alignment to the start segment and at least one spanning read
        # with a long alignment to the end segment.
        longest_start_alignment = max(x[2].get_aligned_ref_length() for x in self.reads)
        longest_end_alignment = max(x[3].get_aligned_ref_length() for x in self.reads)
        alignment_length = min(longest_start_alignment, longest_end_alignment)
        align_length_factor = score_function(alignment_length, min_alignment_length * 4)
        self.quality *= align_length_factor

        # The mean alignment score to the start/end segments is positively correlated with quality,
        # so bridges with high quality alignments are rewarded. Specifically, we want there to be at
        # least one spanning read with a high quality alignment to the start segment and at least
        # one spanning read with a high quality alignment to the end segment.
        best_start_alignment = max(x[2].scaled_score for x in self.reads)
        best_end_alignment = max(x[3].scaled_score for x in self.reads)
        alignment_quality = min(best_start_alignment, best_end_alignment)
        align_score_factor = alignment_quality / 100.0
        self.quality *= align_score_factor

        # Bridges between long start/end segments are rewarded, as they are more likely to actually
        # be single copy. We apply a length factor for both the start and the end segments,
        # and then apply the smaller of two again. This is to punish cases where both segments
        # are not long.
        start_length_factor = score_function(start_seg.get_length(), min_alignment_length * 4)
        self.quality *= start_length_factor
        end_length_factor = score_function(end_seg.get_length(), min_alignment_length * 4)
        self.quality *= end_length_factor
        smaller_length_factor = min(start_length_factor, end_length_factor)
        self.quality *= smaller_length_factor

        # We finalise the quality to a range of 0 to 100. We also use the sqrt function to pull
        # the scores up a bit (otherwise they tend to hang near the bottom of the range).
        self.quality = 100.0 * math.sqrt(self.quality)

        # noinspection PyTypeChecker
        output.append(self.quality)

        return output

    def set_path_based_on_availability(self, graph, unbridged_graph):
        """
        This function will change a bridge's graph path based on what's currently available. This
        is to handle the case where a bridge has multiple possible graph paths, but its first
        choice isn't available anymore (because other bridges used up those segments).
        It has to balance the path quality with the path availability to make a choice.
        """
        best_path = self.all_paths[0][0]
        best_sequence = unbridged_graph.get_path_sequence(best_path)
        best_scaled_score = self.all_paths[0][3]
        best_availability = graph.get_path_availability(best_path)
        for i in range(1, len(self.all_paths)):
            potential_path = self.all_paths[i][0]
            potential_scaled_score = self.all_paths[i][3]
            potential_availability = graph.get_path_availability(potential_path)

            # relative_score measures how much worse this path aligned than the current best.
            # Differences matter more close to scores of 100. E.g. 99 to 97 is a big drop in
            # score, but 79 to 77 is not as large of a drop.
            if potential_scaled_score == 100.0:
                relative_score = 1.0
            else:
                relative_score = (100.0 - best_scaled_score) / (100.0 - potential_scaled_score)
                relative_score = min(1.0, relative_score)

            # relative_availability measures how much more available this path is than the current
            # best. We use 1.1 (instead of 1.0) in the equation to attenuate the affect of
            # availability a bit (because score is more important).
            relative_availability = (1.1 - best_availability) / (1.1 - potential_availability)
            relative_availability = min(2.0, relative_availability)

            # If this path looks better than our current best (considering both score and
            # availability), then it becomes the new best.
            if relative_score * relative_availability > 1.0:
                best_path = potential_path
                best_sequence = unbridged_graph.get_path_sequence(potential_path)
                best_scaled_score = potential_scaled_score
                best_availability = potential_availability

        self.graph_path = best_path
        self.bridge_sequence = best_sequence

    @staticmethod
    def get_type_score():
        """
        Returns a score indicating the relative importance of the bridge types:
        LongReadBridge = 2, SpadesContigBridge = 1, LoopUnrollingBridge = 0
        """
        return 2

    @staticmethod
    def get_type_name():
        """
        Returns the name of the bridge type.
        """
        return 'long read'


def create_long_read_bridges(graph, read_dict, read_names, anchor_segments, verbosity,
                             min_scaled_score, threads, scoring_scheme, min_alignment_length,
                             expected_linear_seqs, min_bridge_qual):
    """
    Makes bridges between single copy segments using the alignments in the long reads.
    """
    log.log_section_header('Building long read bridges')
    log.log_explanation('Unicycler uses the long read alignments to produce bridges between '
                        'anchor segments. These bridges can be formed using as few as one long '
                        'read, giving Unicycler the ability to bridge the graph even when '
                        'long-read depth is low.')

    anchor_seg_nums = set(x.number for x in anchor_segments)

    # This dictionary will collect the read sequences which span between two single copy segments.
    # Key = tuple of signed segment numbers (the segments being bridged)
    # Value = list of tuples containing the bridging sequence and the single copy segment
    #         alignments.
    spanning_read_seqs = defaultdict(list)

    for read_name in read_names:
        read = read_dict[read_name]
        alignments = get_single_copy_alignments(read, anchor_seg_nums, min_scaled_score)
        if len(alignments) < 2:
            continue

        # If the code got here, then we have some alignments to anchor segments. We grab
        # neighbouring pairs of alignments, starting with the highest scoring ones and work our
        # way down. This means that we should have a pair for each neighbouring alignment, but
        # potentially also more distant pairs if the alignments are strong.
        already_added = set()
        sorted_alignments = sorted(alignments, key=lambda x: x.raw_score, reverse=True)
        available_alignments = []
        for alignment in sorted_alignments:

            # If the alignment being added is to a reference that has already been added but in the
            # opposite direction, then we don't include it. E.g. we don't add an alignment for 10
            # if we already have an alignment for -10. This is because there's no legitimate way
            # for a single copy segment to appear in the same read in two different directions. The
            # same direction is okay, as that can happen with a circular piece of DNA, but opposite
            # directions implies multi-copy.
            opposite_num = -alignment.get_signed_ref_num()
            if opposite_num in set(x.get_signed_ref_num() for x in available_alignments):
                continue

            available_alignments.append(alignment)
            available_alignments = sorted(available_alignments,
                                          key=lambda x: x.read_start_positive_strand())
            if len(available_alignments) < 2:
                continue

            for i in range(len(available_alignments)):
                if i < len(available_alignments) - 1:
                    alignment_1 = available_alignments[i]
                    alignment_2 = available_alignments[i + 1]

                # Special case: when the first and last alignments are to the same graph segment,
                # make a bridge for them, even if they aren't a particularly high scoring pair of
                # alignments. This can help to circularise plasmids which are very tied up with
                # other, similar plasmids.
                elif available_alignments[0].ref.name == available_alignments[-1].ref.name:
                    alignment_1 = available_alignments[0]
                    alignment_2 = available_alignments[-1]
                else:
                    continue

                # Standardise the order so we don't end up with both directions (e.g. 5 to -6 and
                # 6 to -5) in spanning_read_seqs.
                seg_nums, flipped = flip_number_order(alignment_1.get_signed_ref_num(),
                                                      alignment_2.get_signed_ref_num())
                if seg_nums not in already_added:
                    bridge_start = alignment_1.read_end_positive_strand()
                    bridge_end = alignment_2.read_start_positive_strand()

                    if bridge_end > bridge_start:
                        bridge_seq = read.sequence[bridge_start:bridge_end]
                        bridge_qual = read.qualities[bridge_start:bridge_end]
                        if flipped:
                            bridge_seq = reverse_complement(bridge_seq)
                            bridge_qual = bridge_qual[::-1]
                    else:
                        bridge_seq = bridge_end - bridge_start  # 0 or a negative number
                        bridge_qual = ''

                    spanning_read_seqs[seg_nums].append((bridge_seq, bridge_qual, alignment_1,
                                                         alignment_2))
                    already_added.add(seg_nums)

    # If a bridge already exists for a spanning sequence, we add the sequence to the bridge. If
    # not, we create a new bridge and add it.
    new_bridges = []
    for seg_nums, span in spanning_read_seqs.items():
        start, end = seg_nums

        # If the start and end are the same and already exclusively connect (i.e. if this segment
        # is already circular), then skip - there's no need to bridge.
        if start == end and graph.get_downstream_seg_nums(start) == [start] and \
                graph.get_upstream_seg_nums(start) == [start]:
            continue

        new_bridge = LongReadBridge(graph, start, end)
        new_bridge.reads += span
        new_bridges.append(new_bridge)

    new_bridges = sorted(new_bridges, key=lambda x: (x.start_segment, x.end_segment))

    # During finalisation, we will compare the expected read count to the actual read count for
    # each bridge. To do this, we'll need the lengths of all reads (excluding those with no
    # alignments). We also need an estimate of the genome size.
    read_lengths = defaultdict(int)
    for read_name in read_names:
        read = read_dict[read_name]
        if read.alignments:
            read_lengths[read.get_length()] += 1
    estimated_genome_size = graph.get_estimated_sequence_len()

    # Now we need to finalise the bridges. This is the intensive step, as it involves creating a
    # consensus sequence, finding graph paths and doing alignments between the consensus and the
    # graph paths. We can therefore use threads to make this faster.
    num_long_read_bridges = len(new_bridges)

    # We want to display this table one row at a time, so we have to fix all of the column widths
    # at the start.
    alignments, col_widths = get_bridge_table_parameters(graph, num_long_read_bridges, verbosity,
                                                         'LongReadBridge')
    print_bridge_table_header(alignments, col_widths, verbosity, 'LongReadBridge')
    completed_count = 0

    # Use a simple loop if we only have one thread.
    if threads == 1:
        for bridge in new_bridges:
            output = bridge.finalise(scoring_scheme, min_alignment_length, read_lengths,
                                     estimated_genome_size, expected_linear_seqs)
            completed_count += 1
            print_bridge_table_row(alignments, col_widths, output, completed_count,
                                   num_long_read_bridges, min_bridge_qual, verbosity,
                                   'LongReadBridge')

    # Use a thread pool if we have more than one thread.
    else:
        pool = ThreadPool(threads)
        arg_list = []

        # Sort the bridges based on how long they're predicted to take to finalise. This will make
        # the big ones runs first which helps to more efficiently use the CPU cores.
        # E.g. if the biggest bridge was at the end, we'd be left waiting for it to finish with
        # only one core (bad), but if it was at the start, other work could be done in parallel.
        long_read_bridges = sorted(new_bridges, reverse=True,
                                   key=lambda x: x.predicted_time_to_finalise())
        for bridge in long_read_bridges:
            arg_list.append((bridge, scoring_scheme, min_alignment_length, read_lengths,
                             estimated_genome_size, expected_linear_seqs))
        for output in pool.imap_unordered(finalise_bridge, arg_list):
            completed_count += 1
            print_bridge_table_row(alignments, col_widths, output, completed_count,
                                   num_long_read_bridges, min_bridge_qual, verbosity,
                                   'LongReadBridge')

    # Now that the bridges are finalised, we split bridges that contain anchor segments in their
    # path such that all bridges start and end on an anchor segment but contain no anchor segments
    # in their path.
    split_bridges = []
    for bridge in new_bridges:
        if not bridge.graph_path:
            split_bridges.append(bridge)
        else:  # bridge has a path
            if not any(abs(x) in anchor_seg_nums for x in bridge.graph_path):  # already good
                split_bridges.append(bridge)

            # If the bridge path contains one or more anchor segments, it must be split!
            else:
                full_path = [bridge.start_segment] + bridge.graph_path + [bridge.end_segment]
                anchor_indices = []
                for i, seg_num in enumerate(full_path):
                    if abs(seg_num) in anchor_seg_nums:
                        anchor_indices.append(i)
                anchor_indices = sorted(anchor_indices)
                for i in range(len(anchor_indices) - 1):
                    start_i = anchor_indices[i]
                    end_i = anchor_indices[i+1]
                    start_seg_num = full_path[start_i]
                    end_seg_num = full_path[end_i]
                    new_path = full_path[start_i+1:end_i]
                    split_bridge = LongReadBridge(graph, start_seg_num, end_seg_num)
                    split_bridge.graph_path = new_path
                    split_bridge.all_paths = [new_path]
                    split_bridge.bridge_sequence = graph.get_path_sequence(new_path)
                    split_bridge.quality = bridge.quality
                    split_bridges.append(split_bridge)

    return split_bridges


def get_single_copy_alignments(read, single_copy_num_set, min_scaled_score):
    """
    Returns a list of single copy segment alignments for the read.
    """
    sc_alignments = []
    for alignment in read.alignments:
        if alignment.ref.number in single_copy_num_set and \
                        alignment.scaled_score >= min_scaled_score:
            sc_alignments.append(alignment)
    return sc_alignments


def finalise_bridge(all_args):
    """
    Just a one-argument version of bridge.finalise, for pool.imap.
    """
    bridge, scoring_scheme, min_alignment_length, read_lengths, estimated_genome_size,\
        expected_linear_seqs = all_args
    return bridge.finalise(scoring_scheme, min_alignment_length, read_lengths,
                           estimated_genome_size, expected_linear_seqs)


def reduce_expected_count(expected_count, a, b):
    r"""
    This function reduces the expected read count. It reduces by a factor which is a function of
    the read count, so low expected values aren't reduced much, but high expected values are
    reduced more. This is to help with high read depth cases where expected counts get quite high.

    https://www.desmos.com/calculator
    y=x\cdot \left(\left(\frac{a}{a+x}\right)\cdot \left(1-b\right)+b\right)
    """
    return expected_count * ((a / (a + expected_count)) * (1.0 - b) + b)


def get_consensus_sequence(reads, scoring_scheme, output):
    consensus_start_time = time.time()

    # Sort the reads from best to worst, as judged by their scaled scores (specifically,
    # whichever scaled score is smaller of their two).
    reads = sorted(reads, reverse=True, key=lambda x: min(x[2].scaled_score, x[3].scaled_score))

    # We toss out any reads which are too much worse than the best, as they are unlikely to
    # improve the consensus much.
    best_scaled_score = min(reads[0][2].scaled_score, reads[0][3].scaled_score)
    min_allowed_score = best_scaled_score - 10.0
    reads = [x for x in reads if min(x[2].scaled_score, x[3].scaled_score) >= min_allowed_score]

    # If we have exactly two sequences, we throw out the worse one unless it is very nearly as good
    # as the best one. This is because two-sequence consensuses are difficult and if there's a
    # major difference in quality, just using the best read will probably be better.
    if len(reads) == 2:
        score_diff = min(reads[0][2].scaled_score, reads[0][3].scaled_score) - \
                     min(reads[1][2].scaled_score, reads[1][3].scaled_score)
        if score_diff > 2.0:
            reads = reads[0:1]

    # Too many reads contributing to the consensus will probably be slow and not worth it, so we
    # set an upper limit.
    if len(reads) > settings.MAX_READS_FOR_CONSENSUS:
        reads = reads[:settings.MAX_READS_FOR_CONSENSUS]

    # If there's only one read, there's no consensus to be done.
    if len(reads) == 1:
        consensus_sequence = reads[0][0]

    # If there's more than one, we make a consensus!
    else:
        read_seqs = [x[0] for x in reads]
        read_quals = [x[1] for x in reads]
        consensus_sequence = consensus_alignment(read_seqs, read_quals, scoring_scheme)[0]

    consensus_time = time.time() - consensus_start_time
    output.append(str(len(consensus_sequence)))
    output.append(float_to_str(consensus_time, 1))
    return consensus_sequence

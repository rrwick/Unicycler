"""
Copyright 2017 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Unicycler

Bridges are links between two single copy segments in an assembly graph. Bridges can come from
multiple sources, each described in a separate module. This module has a few functions that are
common to multiple bridge types.

This file is part of Unicycler. Unicycler is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Unicycler is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Unicycler. If
not, see <http://www.gnu.org/licenses/>.
"""

import math
from .misc import weighted_average, print_table, get_right_arrow, float_to_str


def get_mean_depth(seg_1, seg_2, graph):
    """
    Returns the mean depth of the two segments, weighted by their length.
    """
    return weighted_average(seg_1.depth, seg_2.depth,
                            seg_1.get_length_no_overlap(graph.overlap),
                            seg_2.get_length_no_overlap(graph.overlap))


def get_bridge_str(bridge):
    """
    Returns a bridge sequence in human-readable form.
    """
    bridge_str = str(bridge.start_segment) + ' -> '
    if bridge.graph_path:
        bridge_str += ', '.join([str(x) for x in bridge.graph_path]) + ' -> '
    bridge_str += str(bridge.end_segment)
    return bridge_str


def get_depth_agreement_factor(start_seg_depth, end_seg_depth):
    r"""
    This function is set up such that:
      * equal depths return 1.0
      * similar depths return a value near 1.0
      * more divergent depths return a much lower value:
          a ratio of 1.35 return a value of about 0.5
          a ratio of 2.06 return a value of about 0.1
      * very different depths return a value near 0.0
    https://www.desmos.com/calculator
        y=\frac{1}{1+10^{2\left(\log \left(x-1\right)+0.45\right)}}
        y=\frac{1}{1+10^{2\left(\log \left(\frac{1}{x}-1\right)+0.45\right)}}
    """
    larger_depth = max(start_seg_depth, end_seg_depth)
    smaller_depth = min(start_seg_depth, end_seg_depth)
    if larger_depth == 0.0 or smaller_depth == 0.0:
        return 0.0
    elif larger_depth == smaller_depth:
        return 1.0
    else:
        ratio = larger_depth / smaller_depth
        return 1.0 / (1.0 + 10.0 ** (2 * (math.log10(ratio - 1.0) + 0.45)))


def get_bridge_table_parameters(graph, num_bridges, verbosity, bridge_type):
    """
    Used for LongReadBridge and MiniasmBridge objects.
    """
    assert bridge_type == 'LongReadBridge' or bridge_type == 'MiniasmBridge'
    max_seg_num_len = len(str(max(graph.segments.keys()))) + 1

    table_alignments = 'RL'
    table_col_widths = [2 * len(str(num_bridges)) + 1,  # Number
                        max_seg_num_len + 10]           # Start to end

    if verbosity > 1 and bridge_type == 'LongReadBridge':
        table_alignments += 'R'
        table_col_widths += [5]   # Read count

    if verbosity > 1:
        table_alignments += 'R'
        table_col_widths += [9]  # Consensus length

    if verbosity > 2 and bridge_type == 'LongReadBridge':
        table_alignments += 'RR'
        table_col_widths += [9, 8]   # Consensus time, Target length

    if verbosity > 1:
        table_alignments += 'LRR'
        table_col_widths += [11, 8, 5]   # Search type, Search time, Path count

    table_alignments += 'L'
    table_col_widths += [40]  # Best path

    if verbosity > 2:
        table_alignments += 'RRRR'
        table_col_widths += [9, 9, 12, 11]  # Best path length, raw score, scaled score, len disc

    table_alignments += 'R'
    table_col_widths += [7]  # Quality

    return table_alignments, table_col_widths


def print_bridge_table_header(alignments, col_widths, verbosity, bridge_type):
    """
    Used for LongReadBridge and MiniasmBridge objects.
    """
    assert bridge_type == 'LongReadBridge' or bridge_type == 'MiniasmBridge'
    header_line_1 = ['', '']
    header_line_2 = ['', 'Start ' + get_right_arrow() + ' end']

    if verbosity > 1 and bridge_type == 'LongReadBridge':
        header_line_1 += ['']
        header_line_2 += ['Reads']
    if verbosity > 1:
        header_line_1 += ['Consensus']
        header_line_2 += ['len (bp)']
    if verbosity > 2 and bridge_type == 'LongReadBridge':
        header_line_1 += ['Consensus', 'Target']
        header_line_2 += ['time (s)',  'len (bp)']
    if verbosity > 1:
        header_line_1 += ['',            'Search',   'Path']
        header_line_2 += ['Search type', 'time (s)', 'count']

    header_line_1 += ['']
    header_line_2 += ['Best path']

    if verbosity > 2:
        header_line_1 += ['Best path', 'Best path', 'Best path',    'Best path']
        header_line_2 += ['len (bp)',  'raw score', 'scaled score', 'length disc']

    header_line_1 += ['']
    header_line_2 += ['Quality']

    if any(x for x in header_line_1):
        print_table([header_line_1], col_separation=2, alignments=alignments,
                    header_format='normal', fixed_col_widths=col_widths, indent=0)
    print_table([header_line_2], col_separation=2, alignments=alignments,
                header_format='underline', fixed_col_widths=col_widths, indent=0)


def print_bridge_table_row(alignments, col_widths, output, completed_count, num_bridges,
                           min_bridge_qual, verbosity, bridge_type):
    """
    Used for LongReadBridge and MiniasmBridge objects.
    """
    assert bridge_type == 'LongReadBridge' or bridge_type == 'MiniasmBridge'
    fraction = str(completed_count) + '/' + str(num_bridges)

    start, end, read_count, consensus_length, consensus_time, target_length, path_count, \
        search_type, search_time, best_path, best_path_len, best_path_raw_score, \
        best_path_scaled_score, best_path_length_discrepancy, quality = output

    start_to_end = (start + ' ' + get_right_arrow()).rjust(7) + ' ' + end
    quality_str = float_to_str(quality, 3)

    table_row = [fraction, start_to_end]

    if verbosity > 1 and bridge_type == 'LongReadBridge':
        table_row.append(read_count)
    if verbosity > 1:
        table_row.append(consensus_length)
    if verbosity > 2 and bridge_type == 'LongReadBridge':
        table_row += [consensus_time, target_length]
    if verbosity > 1:
        table_row += [search_type, search_time, path_count]

    table_row += [best_path]

    if verbosity > 2:
        table_row += [best_path_len, best_path_raw_score, best_path_scaled_score,
                      best_path_length_discrepancy]

    table_row += [quality_str]

    sub_colour = {}
    if quality < min_bridge_qual:
        sub_colour[quality_str] = 'red'
    print_table([table_row], col_separation=2, header_format='normal', indent=0,
                left_align_header=False, alignments=alignments, fixed_col_widths=col_widths,
                sub_colour=sub_colour, bottom_align_header=False)

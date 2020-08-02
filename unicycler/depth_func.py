import csv
import sys
import os
import glob
import subprocess
from .unicycler_align import load_references
from .cpp_wrappers import minimap_align_reads
from .minimap_alignment import load_minimap_alignments
from . import log


class CannotCalculateDepth(Exception):
    def __init__(self, message):
        self.message = message

    def __str__(self):
        return repr(self.message)


def get_read_depth(final_fasta, args):
    log.log_section_header('Calculating read depth', single_newline=True)

    short_reads_available = bool(args.short1) or bool(args.unpaired)
    long_reads_available = bool(args.long)

    if short_reads_available:
        bam_filename = get_short_read_alignments(final_fasta, args)
        short_read_depth = get_short_read_depth(bam_filename, args)
        short_read_depth = round(short_read_depth, 2)
        log.log('Short read depth:   ' + str(short_read_depth))

    if long_reads_available:
        aligned_long_reads = get_long_read_alignments(final_fasta, args.long, args.threads)
        long_read_depth = get_long_read_depth(aligned_long_reads, final_fasta)
        long_read_depth = round(long_read_depth, 2)
        log.log('Long read depth:    ' + str(long_read_depth))

    if short_reads_available and long_reads_available:
        total_read_depth = long_read_depth + short_read_depth
        log.log('Total read depth:   ' + str(total_read_depth))


def get_short_read_depth(bam_filename, args):
    """
    this command will generate a 3 column tsv: contig, contig_pos, depth
    """
    depth_tsv_file = os.path.join(args.out, 'depth.tsv')
    samtools_depth_command = ['samtools', 'depth', '-a', bam_filename]

    try:
        f = open(depth_tsv_file, 'w')
        subprocess.call(samtools_depth_command, stdout=f, stderr=subprocess.STDOUT)
        f.close()
    except subprocess.CalledProcessError as e:
        raise CannotCalculateDepth('samtools depth encountered an error:\n' + e.output.decode())

    depth_tsv = csv.reader(open(depth_tsv_file), delimiter='\t')

    depth_at_bases_sum = 0
    total_ref_length = 0
    for row in depth_tsv:
        depth_at_bases_sum += int(row[2])
        total_ref_length += 1

    # rm depth.tsv file
    try:
        os.remove(depth_tsv_file)
    except FileNotFoundError:
        pass

    # rm short read alignments bam
    try:
        os.remove(bam_filename)
    except FileNotFoundError:
        pass

    # rm short read alignments bam index
    try:
        os.remove(bam_filename + ".bai")
    except FileNotFoundError:
        pass

    return depth_at_bases_sum / total_ref_length


def get_short_read_alignments(final_fasta, args):
    final_fasta = os.path.join(args.out, final_fasta)
    using_paired_reads = bool(args.short1) and bool(args.short2)
    using_unpaired_reads = bool(args.unpaired)
    assert using_paired_reads or using_unpaired_reads # just in case --no-pilon

    # build the bowtie reference index
    bowtie2_build_command = [args.bowtie2_build_path, final_fasta, final_fasta]
    try:
        subprocess.check_output(bowtie2_build_command, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        raise CannotCalculateDepth('bowtie2-build encountered an error:\n' + e.output.decode())

    # run the alignments
    sam_filename = os.path.join(args.out, 'depth_alignments.sam')
    if using_paired_reads:
        bowtie2_command = [args.bowtie2_path, '-S', sam_filename, '-1',
                           args.short1, '-2', args.short2, '-x', final_fasta,
                           '--threads', str(args.threads),
                           '--local', '--fast-local']

    if using_unpaired_reads:
        bowtie2_command = [args.bowtie2_path, '-S', sam_filename,
                           '-U', args.unpaired, '-x', final_fasta,
                           '--threads', str(args.threads),
                           '--local', '--fast-local']
    try:
        subprocess.check_output(bowtie2_command, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        raise CannotCalculateDepth('Bowtie2 encountered an error:\n' + e.output.decode())

    # sort the alignments
    bam_filename = os.path.join(args.out, 'depth_alignments.bam')
    samtools_sort_command = [args.samtools_path, 'sort', '-@', str(args.threads),
                             '-o', bam_filename, '-O', 'bam', '-T', 'temp', sam_filename]

    try:
        subprocess.check_output(samtools_sort_command, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        raise CannotCalculateDepth('Samtools encountered an error:\n' + e.output.decode())

    # rm remaining sam file
    try:
        os.remove(sam_filename)
    except FileNotFoundError:
        pass

    # rm .bt2 files
    bowtie2_old_files = glob.glob(os.path.join(args.out, "assembly.fasta.*.bt2"))
    for file in bowtie2_old_files:
        try:
            os.remove(file)
        except FileNotFoundError:
            pass

    # index the alignments
    samtools_index_command = [args.samtools_path, 'index', bam_filename]
    try:
        subprocess.check_output(samtools_index_command, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        raise CannotCalculateDepth('Samtools encountered an error:\n' + e.output.decode())

    return bam_filename


def get_long_read_alignments(final_fasta, reads_fastq, threads):
    log.log(str(final_fasta) + "  " + str(reads_fastq))
    minimap_alignments_str = minimap_align_reads(final_fasta, reads_fastq, threads, 0, 'default')
    minimap_alignments = load_minimap_alignments(minimap_alignments_str)
    return minimap_alignments


def get_long_read_depth(minimap_alignments, final_fasta):
    total_ref_length = 0
    total_read_length = 0
    reference_contigs = load_references(final_fasta, show_progress=False, section_header='')
    for contig in reference_contigs:
        total_ref_length += contig.get_length()
    for key in minimap_alignments.keys():
        if len(minimap_alignments[key]) > 0:
            total_read_length += minimap_alignments[key][0].get_read_length()
    return total_read_length / total_ref_length

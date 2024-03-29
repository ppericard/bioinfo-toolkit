#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import gzip


def open_fastq(fastq_filepath, read_write):
    """
    Opens a FASTQ file for reading or writing. Determines if the file should be
    treated as gzipped based on its extension. If opening for writing and the file
    does not exist, it will be created.
    """
    is_gzip = fastq_filepath.endswith(".gz")
    mode = f"{read_write}t"

    try:
        if is_gzip:
            file_handle = gzip.open(fastq_filepath, mode)
        else:
            file_handle = open(fastq_filepath, mode)
    except IOError as e:
        raise IOError(f"Failed to open '{fastq_filepath}' in mode '{mode}': {e}")

    return file_handle


def parse_fastq(fastq_filehandle):
    """
    Generator function to parse fastq files. Yields each read as a tuple
    containing the header, sequence, and quality string.
    """
    line_number = 0
    for line in fastq_filehandle:
        stripped_line = line.strip()
        if isinstance(line, bytes):  # Handle binary mode for gzip
            stripped_line = stripped_line.decode("utf-8")
        if stripped_line:
            line_number += 1
            if line_number % 4 == 1:
                read_header = stripped_line[1:]
            elif line_number % 4 == 2:
                read_seq = stripped_line
            elif line_number % 4 == 3:
                continue  # Skip the '+' line
            elif line_number % 4 == 0:
                read_qual = stripped_line
                yield (read_header, read_seq, read_qual)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Merge a reads fastq file with the corresponding UMI fastq file. The UMI will be placed in 5' of the read sequence."
    )
    parser.add_argument(
        "-i",
        "--input",
        metavar="INFQ",
        type=str,
        required=True,
        help="input fastq file",
    )
    parser.add_argument(
        "-u", 
        "--umi", 
        metavar="UMIFQ", 
        type=str, 
        required=True, 
        help="UMIs fastq file"
    )
    parser.add_argument(
        "-o",
        "--output",
        metavar="OUTFQ",
        type=str,
        required=True,
        help="output fastq file",
    )
    parser.add_argument(
        "--header",
        action="store_true",
        help="append the UMI to the read header instead of merging it with the sequence",
    )

    args = parser.parse_args()

    with open_fastq(args.input, "r") as input_fastq_filehandle, \
         open_fastq(args.umi, "r") as umi_fastq_filehandle, \
         open_fastq(args.output, "w") as output_fastq_filehandle:
        
        if args.header:
            for (input_header, input_seq, input_qual), (umi_header, umi_seq, umi_qual) \
            in zip(parse_fastq(input_fastq_filehandle), parse_fastq(umi_fastq_filehandle)):
                split_input_header = input_header.split()
                output_header = f"{split_input_header[0]}_{umi_seq} {split_input_header[1:].join(' ')}"
                output_fastq_filehandle.write(
                    f"@{output_header}\n{input_seq}\n+\n{input_qual}\n"
                )
        else:
            for (input_header, input_seq, input_qual), (umi_header, umi_seq, umi_qual) \
            in zip(parse_fastq(input_fastq_filehandle), parse_fastq(umi_fastq_filehandle)):
                merged_seq = umi_seq + input_seq
                merged_qual = umi_qual + input_qual
                output_fastq_filehandle.write(
                    f"@{input_header}\n{merged_seq}\n+\n{merged_qual}\n"
                )

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import gzip
import time

def open_fastq(fastq_filepath, mode):
    """
    Opens a FASTQ file, which can be optionally compressed with gzip.
    
    Parameters:
    - fastq_filepath: Path to the FASTQ file.
    - mode: Mode to open the file in ('r' for read, 'w' for write, etc.).
    
    Returns:
    - File handle to the opened file.
    
    Raises:
    - FileNotFoundError if the file does not exist.
    - PermissionError if there is no permission to access the file.
    - IOError for other IO related errors.
    """
    try:
        if fastq_filepath.endswith(".gz"):
            # Handle gzip compressed files
            return gzip.open(fastq_filepath, mode + "t")
        # Handle uncompressed files
        return open(fastq_filepath, mode)
    except FileNotFoundError:
        raise FileNotFoundError(f"File '{fastq_filepath}' not found.")
    except PermissionError:
        raise PermissionError(f"Permission denied when accessing '{fastq_filepath}'.")
    except Exception as e:
        raise IOError(f"Failed to open '{fastq_filepath}' in mode '{mode}': {e}")

def parse_fastq(fastq_filehandle):
    """
    Generator function to parse fastq files, yielding each read as a tuple
    containing the header, sequence, and quality string.
    
    Parameters:
    - fastq_filehandle: An open file handle to a FASTQ file.
    
    Yields:
    - Tuple of (header, sequence, quality) for each read in the FASTQ file.
    """
    while True:
        try:
            # Attempt to read 4 lines at a time (FASTQ format standard)
            lines = [next(fastq_filehandle).strip() for _ in range(4)]
        except StopIteration:
            break  # End of file reached
        # Yield the parsed components of each read
        yield lines[0][1:], lines[1], lines[3]

def merge_sequences(args, input_fastq, umi_fastq, output_fastq):
    """
    Merge sequences or headers from input and UMI FASTQ files based on user arguments.
    
    Parameters:
    - args: Command line arguments specifying user preferences.
    - input_fastq: File handle to the input FASTQ file.
    - umi_fastq: File handle to the UMI FASTQ file.
    - output_fastq: File handle for the output FASTQ file where the results are written.
    """
    read_count = 0
    start_time = time.time()
    for (input_header, input_seq, input_qual), (_, umi_seq, _) in zip(parse_fastq(input_fastq), parse_fastq(umi_fastq)):
        if args.header:
            # Append the UMI to the read header
            output_header = f"{input_header.split()[0]}{args.delimiter}{umi_seq} {' '.join(input_header.split()[1:])}"
            output_fastq.write(f"@{output_header}\n{input_seq}\n+\n{input_qual}\n")
        else:
            # Merge the UMI sequence with the input sequence
            merged_seq = umi_seq + input_seq
            output_fastq.write(f"@{input_header}\n{merged_seq}\n+\n{input_qual}\n")
        
        read_count += 1
        if read_count % 10000 == 0:  # Update for every 10000 reads processed
            elapsed_time = time.time() - start_time
             # Work in kilo-reads for display
            kreads_per_minute = (read_count / elapsed_time) * 60 / 1000
            kreads_processed = read_count / 1000
            print(f"Processed {kreads_processed:.2f} kilo-reads at {kreads_per_minute:.2f} kilo-reads/minute", end='\n')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merge FASTQ sequences with UMIs.")
    # Definitions for command-line arguments with descriptions
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
    parser.add_argument(
        "-d",
        "--delimiter",
        type=str,
        default='_',
        help="delimiter between the read name and the UMI",
    )
    args = parser.parse_args()

    try:
        # Open the input, UMI, and output files, handling them within a context manager
        with open_fastq(args.input, "r") as input_fastq, \
             open_fastq(args.umi, "r") as umi_fastq, \
             open_fastq(args.output, "w") as output_fastq:
            
            # Merge sequences or headers based on user arguments
            merge_sequences(args, input_fastq, umi_fastq, output_fastq)
    except IOError as e:
        print(f"Error: {e}")

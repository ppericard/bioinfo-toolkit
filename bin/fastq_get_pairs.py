#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Get Pairs

Description : Get separately paired reads and singletons 
                from two fastq files (left and right)

  get_pairs.py file1.fastq file2.fastq

-----------------------------------------------------------------------

Author: This software is written and maintained by Pierre Pericard
(pierre.pericard@ed.univ-lille1.fr)
Created: 2012-11-09
Last Modified: 2016-01-11
Licence: GNU GPL 3.0

Copyright 2013-2016 Pierre Pericard

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import argparse
import sys
import gzip

def open_fastq_gz(fastq_filepath, read_write):
    fastq_extension = fastq_filepath.split('.')[-1]
    fastq_filehandle = None
    if fastq_extension == 'gz':
        fastq_filehandle = gzip.open(fastq_filepath, f'{read_write}t')
    else:
        fastq_filehandle = open(fastq_filepath, f'{read_write}')
    return fastq_filehandle

def parse_fastq(fastq_filehandle):
    line_number = 0
    read_header = ""
    read_seq = ""
    read_qual = ""
    for line in fastq_filehandle:
        l = line.strip()
        if l:
            line_number += 1
            if line_number%4==1:
                read_header = l[1:]
            elif line_number%4==2:
                read_seq = l
            elif line_number%4==0:
                read_qual = l
                yield (read_header, read_seq, read_qual)


if __name__ == '__main__':

    # Arguments
    parser = argparse.ArgumentParser(description='Get separately paired reads and singletons from two fastq files (left and right)')
    parser.add_argument('leftreads', metavar='leftreads', type=argparse.FileType('r'), help='left reads fastq')
    parser.add_argument('rightreads', metavar='rightreads', type=argparse.FileType('r'), help='right reads fastq')

    args = parser.parse_args()
    

    left_fastq_filepath = args.leftreads.name
    right_fastq_filepath = args.rightreads.name

    args.leftreads.close()
    args.rightreads.close()

    print(left_fastq_filepath)
    print(right_fastq_filepath)
    
    left_fastq_filehandle = open_fastq_gz(left_fastq_filepath, 'r')
    right_fastq_filehandle = open_fastq_gz(right_fastq_filepath, 'r')

    left_fastq_reads_id_list = list()
    right_fastq_reads_id_list = list()

    for fastq_filehandle, fastq_reads_id_list in ((left_fastq_filehandle, left_fastq_reads_id_list), (right_fastq_filehandle, right_fastq_reads_id_list)):
        print(f'Reading {fastq_filehandle.name}')
        for (header, seq, qual) in parse_fastq(fastq_filehandle):
            read_id = header.split()[0]
            fastq_reads_id_list.append(read_id)
        fastq_filehandle.close()

    not_common_read_ids = set(left_fastq_reads_id_list) ^ set(right_fastq_reads_id_list)

    print(f'There are {len(not_common_read_ids)} singletons')

    for fastq_filepath in (left_fastq_filepath, right_fastq_filepath):
        basefilepath = '.'.join(fastq_filepath.split('.')[:-1])
        extension = fastq_filepath.split('.')[-1]
        if extension == "gz":
            extension = '.'.join(fastq_filepath.split('.')[-2:])
            basefilepath = '.'.join(fastq_filepath.split('.')[:-2])

        paired_fastq_filehandle = open_fastq_gz(f'{basefilepath}.paired.{extension}', 'w')
        unpaired_fastq_filehandle = open_fastq_gz(f'{basefilepath}.unpaired.{extension}', 'w')

        fastq_filehandle = open_fastq_gz(fastq_filepath, 'r')

        for (header, seq, qual) in parse_fastq(fastq_filehandle):
            read_id = header.split()[0]
            paired = read_id not in not_common_read_ids
            output_fastq_filehandle = unpaired_fastq_filehandle
            if paired:
                output_fastq_filehandle = paired_fastq_filehandle

            output_fastq_filehandle.write(f'@{header}\n')
            output_fastq_filehandle.write(f'{seq}\n')
            output_fastq_filehandle.write(f'+\n')
            output_fastq_filehandle.write(f'{qual}\n')

        paired_fastq_filehandle.close()
        unpaired_fastq_filehandle.close()
        fastq_filehandle.close()

        

    # for f in (leftreads, rightreads):
        
    #     basefilename = '.'.join(f.split('.')[:-1])
    #     extension = f.split('.')[-1]
        
    #     if args.galaxy is True :
    #         extension = "fastq"
    #         if f == leftreads:
    #             basefilename = "left"
    #         else :
    #             basefilename = "right"
        
    #     pfh = open(basefilename+'.paired.'+extension, 'w')
    #     ufh = open(basefilename+'.unpaired.'+extension, 'w')
    #     with open(f, 'r') as fh:
    #         c = 0
    #         paired = False
    #         for line in fh:
    #             l = line.strip()
    #             if l:
    #                 c += 1
    #                 if c%4==1:
    #                     paired = l.split()[0][1:].split('/')[0] not in notcommon
    #                     if c%40000==1:
    #                         sys.stdout.write("\r%.2f M reads writen" % (c/4000000.0))
    #                 if paired:
    #                     pfh.write("%s\n" % l)
    #                 else:
    #                     ufh.write("%s\n" % l)
    #         sys.stdout.write("\r%.2f M reads writen\n" % (c/4000000.0))
    #     pfh.close()
    #     ufh.close()

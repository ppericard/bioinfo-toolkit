#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
FastqNameFilter

Description: Filter a fastq file based on a string to find in the
               sequences headers

  fastqNameFilter.py -i input.fq -o output.fq -s "stringtofind"

-----------------------------------------------------------------------

Author: This software is written and maintained by Pierre Pericard
(pierre.pericard@ed.univ-lille1.fr)
Created: 2014
Last Modified: 2016-01-13
Licence: GNU GPL 3.0

Copyright 2014-2016 Pierre Pericard

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


import sys
import os
import argparse
import string
import re

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Filter a fastq file based on sequence name.')
    parser.add_argument('-i', '--input_fastq', metavar='input', 
                        type=argparse.FileType('r', 0), default='-',
                        help='input fastq file')
    parser.add_argument('-o', '--output_fastq', metavar='output', 
                        type=argparse.FileType('w', 0), default='-',
                        help='ouput fastq file')
    parser.add_argument('-s', '--stringtofind', metavar='string',
                        required=True, help='String to filter on')
    args = parser.parse_args()
    
    tofind = re.compile(args.stringtofind, flags=re.IGNORECASE)
    
    count = 0
    buff = ''
    to_write = False
    
    for line in args.input_fastq:
        line = line.strip()
        if line:
            count += 1
            if count % 4 == 1:
                if to_write:
                    args.output_fastq.write(buff)
                to_write = False
                buff = ''
                if tofind.search(line):
                    to_write = True
            buff += "{0}\n".format(line)
    # Last line
    if to_write:
        args.output_fastq.write(buff)

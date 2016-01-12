#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Get Pairs

Description : Get separately paired reads and singletons 
                from two fastq files (left and right)

  get_pairs.py file1.fastq file2.fastq

---------------------------------------------------------

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


if __name__ == '__main__':

    # Arguments
    parser = argparse.ArgumentParser(description='Get separately paired reads and singletons from two fastq files (left and right)')
    parser.add_argument('leftreads', metavar='leftreads', type=argparse.FileType('r'), help='left reads fastq')
    parser.add_argument('rightreads', metavar='rightreads', type=argparse.FileType('r'), help='right reads fastq')
    parser.add_argument('--galaxy', dest='galaxy', action="store_true", default=False, help="Enable the galaxy mode which create 4 static outputs")

    args = parser.parse_args()
    
    leftreads = args.leftreads.name
    rightreads = args.rightreads.name
    
    #(leftreads, righttreads) = sys.argv[1:3]
    (n1, n2) = (list(), list())

    for f, n in ((leftreads, n1), (rightreads, n2)):
        with open(f, 'r') as fh:
            c = 0
            for line in fh:
                l = line.strip()
                if l:
                    c += 1
                    if c%4==1:
                        n.append(l.split()[0][1:].split('/')[0])
                        if c%40000==1:
                            sys.stdout.write("\r%.2f M reads read" % (c/4000000.0))
            sys.stdout.write("\r%.2f M reads read\n" % (c/4000000.0))

    notcommon = set(n1) ^ set(n2)

    for f in (leftreads, rightreads):
        
        basefilename = '.'.join(f.split('.')[:-1])
        extension = f.split('.')[-1]
        
        if args.galaxy is True :
            extension = "fastq"
            if f == leftreads:
                basefilename = "left"
            else :
                basefilename = "right"
        
        pfh = open(basefilename+'.paired.'+extension, 'w')
        ufh = open(basefilename+'.unpaired.'+extension, 'w')
        with open(f, 'r') as fh:
            c = 0
            paired = False
            for line in fh:
                l = line.strip()
                if l:
                    c += 1
                    if c%4==1:
                        paired = l.split()[0][1:].split('/')[0] not in notcommon
                        if c%40000==1:
                            sys.stdout.write("\r%.2f M reads writen" % (c/4000000.0))
                    if paired:
                        pfh.write("%s\n" % l)
                    else:
                        ufh.write("%s\n" % l)
            sys.stdout.write("\r%.2f M reads writen\n" % (c/4000000.0))
        pfh.close()
        ufh.close()

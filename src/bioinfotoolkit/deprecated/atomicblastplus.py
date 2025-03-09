#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
DEPRECATED: This script has been permanently deprecated and archived.
Please consider using alternative BLAST tools or contact the maintainer for recommendations.

Atomic Blast+ v5.0

Description: Submit a massively parallel Blast+ job-array to a
               computer cluster running on Oracle Grid Engine
               (previously Sun Grid Engine)

  atomicblastplus.py -p blastp -i input.fa -d nr -o input_vs_nr

-----------------------------------------------------------------------

Based upon a Wilfrid Carre original idea ^^

Author: This software is written and maintained by Pierre Pericard
(pierre.pericard@ed.univ-lille1.fr)
Created: 2013
Last Modified: 2016-01-11
Archived: 2024
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

import os
import sys
import argparse
import getpass
import re
import subprocess

# Add deprecation warning that prints when the script is run
print("\n⚠️ DEPRECATION WARNING ⚠️")
print("This script has been permanently deprecated and archived.")
print("It is kept here for historical reference only and should not be used.")
print("The script will likely not work with modern systems and may contain security vulnerabilities.")
print("Please consider using alternative BLAST tools or contact the maintainer for recommendations.\n")

# Exit immediately with warning message
sys.exit("Script execution stopped due to deprecation. This script is archived for reference only.")

class DefaultHelpParser(argparse.ArgumentParser):     
    """
    This is a slightly modified argparse parser to display the full help 
    on parser error instead of only usage 
    """                              
    def error(self, message):                                                         
        sys.stderr.write('\nerror: %s\n\n' % message)                                       
        self.print_help()                                                               
        sys.exit(2)

def get_args():
    """
    Return pre-processed command-line arguments
    """
    parser = DefaultHelpParser(description='Atomic Blast+ can split an input multifasta file, generate a SGE script and submit an array job of Blast+ jobs to the cluster.',
                               epilog='atomicblastplus.py -p blastp -i input.fa -d /db/blast/all/nr -o myproject/tmp/blast/input_vs_nr --options " -seg yes"',
                               # to precisely format help display
                               formatter_class=lambda prog: argparse.HelpFormatter(prog, width=175, max_help_position=80))
    parser.add_argument('-p', '--program', metavar='PROGRAM',
                        default='blastp', type=str, choices=['blastn', 'blastp', 'blastx', 'tblastn', 'tblastx', 'megablast'],
                        help='Blast+ algorithm to use from: blastn, blastp, blastx, tblastn, tblastx, megablast (default: blastp)')
    parser.add_argument('-i', '-query', '--input_file', metavar='INPUT',
                        required=True, type=str, help='Fasta input file with query sequences')
    parser.add_argument('-d', '-db', '--database', metavar='DATABASE',
                        default='/db/blast/all/nr', type=str, help='Blast formated database (default: nr)')
    parser.add_argument('-o', '-out', '--output', metavar='OUTPUT',
                        type=str, help='Output file basename (default: $input_file_basename.atomic_$program_vs_$database_name)')
    parser.add_argument('-n', '--nb_seq', metavar='NBSEQ', default='50',
                        type=int, help='Number of sequences per subfile (default: 50)')
    parser.add_argument('-f', '--nb_file', metavar='NBFILE',
                        type=int, help='Number of subfiles to be created. Incompatible with the -n/--nb_seq argument')
    parser.add_argument('-m', '--email', metavar='EMAIL',
                        type=str, help='Email address for the SGE reports (default: username@sb-roscoff.fr)')
    parser.add_argument('-q', '--queue', metavar='QUEUE',
                        default='short.q', type=str, help=argparse.SUPPRESS)
    parser.add_argument('-s', '--steps', metavar='STEPS', default='1,2,3',
                        type=str, help='Steps to run 1:split, 2:shell, 3:submit (default: 1,2,3)')
    parser.add_argument('-b', '--batch', metavar='BATCH', default='all',
                        type=str, help='Job batches to submit, eg 1-50 or all (default: all)')
    parser.add_argument('-outfmt', '--output_format', metavar='OUTFMT',
                        default='tabular', type=str, choices=['tabular', 'pairwise', 'xml', 'extended'],
                        help='Output format to choose from: tabular, pairwise, xml, extended (tabular + additional columns: qlen, slen) (default: tabular)')
    parser.add_argument('-e', '-evalue', '--evalue', metavar='EVALUE',
                        default='10', type=float, help='Expectation value (E) threshold for saving hits (default: 10)')
    parser.add_argument('-c', '--cpu', metavar='CPU',
                        default='100', type=int, help=argparse.SUPPRESS)
    parser.add_argument('-max_target_seqs', metavar='MAXTARGET', default='500',
                        type=int, help='Maximum number of aligned sequences to keep (default: 500)')
    parser.add_argument('--options', metavar='OPTIONS', type=str, default='',
                        help='Additional options to pass to the program (Blast+ format) [WARNING: the string between quote should always begin with a space], eg. " -seg yes -max_target_seqs 10"')
    parser.add_argument('--dont_wait', action='store_true',
                        help='Do not wait for the SGE array job to complete (equivalent to the "qsub -sync no" option')
    parser.add_argument('--no_cleanup', action='store_true',
                        help='Do not clean all intermediary files and directory after output files concatenation')
    parser.add_argument('--send_mails', action='store_true', help='Send begin, end, and abort mails (Default: only send abort emails)')
    parser.add_argument('--verbose', action='store_true', help='Display additional informations about the run')
    parser.add_argument('-v', '--version', action='version', version='Atomic Blast+ v5.0.3')
    args = parser.parse_args()
    # Additional arguments checking
    if args.input_file == '-' and args.nb_file:
        sys.stdout.write('ERROR: [Arguments] Can\'t count sequences number on STDIN. Please use -n/--nb_seq instead.\n\n')
        parser.print_help()
        sys.exit(1)
    return args

# NOTE: The rest of the code has been truncated as this script is deprecated and archived.
# The full source code is maintained in the repository history for reference purposes only. 
#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
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

def read_fasta_file_handle(fasta_file_handle):
    """
    Parse a fasta file and return a generator
    """
    # Variables initialization
    header = ''
    seqlines = list()
    sequence_nb = 0
    # Reading input file
    for line in fasta_file_handle:
        if line[0] == '>':
            # Yield the last read header and sequence
            if sequence_nb:
                yield (header, ''.join(seqlines))
                del seqlines[:]
            # Get header
            header = line[1:].rstrip()
            sequence_nb += 1
        else:
            # Concatenate sequence
            seqlines.append(line.strip())
    # Yield the input file last sequence
    yield (header, ''.join(seqlines))
    # Close input file
    fasta_file_handle.close()

def format_seq(seq, linereturn=80):
    """
    Format an input sequence
    """
    buff = list()
    for i in range(0, len(seq), linereturn):
        buff.append("{0}\n".format(seq[i:(i + linereturn)]))
    return ''.join(buff)

def split_multifasta_by_increment(input_file_handle, directory_path, output_files_basename, max_nb_seq, verbose):
    """
    Split a multifasta file in subfiles of nb_seq_increment sequences
    """
    # Variable initialization
    subfile_nb = 0
    sequence_nb = 0
    subfile_fh = None
    # Reading input multifasta file
    for header, sequence in read_fasta_file_handle(input_file_handle):
        # Open a new subfile at the beginning and when reaching the maximum sequence number
        if sequence_nb % max_nb_seq == 0:
            # Close the last opened subfile
            if subfile_nb:
                subfile_fh.close()
            # Open a new subfile
            subfile_nb += 1
            subfile_name = "{0}.{1}.fasta".format(subfile_nb, output_files_basename)
            if verbose:
                sys.stdout.write("\nFILE: writing {0}".format(subfile_name))
            subfile_path = "{0}/{1}".format(directory_path, subfile_name)
            try:
                subfile_fh = open(subfile_path, 'w')
            except:
                sys.stderr.write("\nERROR: [Subfile] {0} cannot be created\n\n".format(subfile_path))
                raise
        # Write sequences in the subfiles
        subfile_fh.write(">{0}\n".format(header))
        subfile_fh.write("{0}\n".format(format_seq(sequence, linereturn=80)))
        sequence_nb += 1
    if verbose:
        sys.stdout.write('\n')
    # Return the total subfiles number
    return subfile_nb

def run_shell(cmd_line, verbose=False):
    """
    Run a shell command and return STDOUT, STDERR, and return code
    """
    if verbose:
        sys.stdout.write("\nCMD: {0}\n".format(cmd_line))
    process = subprocess.Popen(cmd_line, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    # Retrieve STDOUT and STDERR
    process_out, process_err = process.communicate()
    # Get return code
    process_errcode = process.returncode
    # For debug and/or print
    if verbose:
        if process_out.strip():
            # If there is something in STDOUT
            sys.stdout.write("\n{0}\n".format(process_out.strip()))
        if process_err.strip():
            # If there is something in STDERR
            sys.stderr.write("\n{0}\n".format(process_err.strip()))
    # Return STDOUT, STDERR, and return code
    return process_out, process_err, process_errcode
    
def main(args):
    """
    Main program
    """
    # Initialise variables
    global_error_code = 0
    
    # Start message
    sys.stdout.write("\n!!! This is Atomic Blast+ !!!\n\n")
    
    # Megablast is blastn+ default algorithm
    # real blastn algorithm has to be set
    blast_program = args.program
    if args.program == 'blastn':
        blast_program = 'blastn -task blastn'
    elif args.program == 'megablast':
        blast_program = 'blastn -task megablast'
    
    # Get Blast version
    blast_version = subprocess.check_output("{0} -version".format(blast_program), shell=True).strip()
    sys.stdout.write("PROGRAM:\n{0}\n\n".format(blast_version))
    
    # Open the input file
    try:
        sys.stdout.write("QUERY: {input}\n".format(input=args.input_file))
        input_fh = open(args.input_file, 'r')
    except IOError:
        sys.stderr.write("\nERROR: [Input] {0} cannot be open\n\n".format(args.input_file))
        raise
    
    # Get database name
    sys.stdout.write("DB: {database}\n".format(database=args.database))
    database_path = args.database
    database_name = database_path.split('/')[-1]
    
    # Generate the subfiles basename
    # get the input basename and remove the fasta extension
    subfiles_basename = re.sub(r'(\.fasta|\.fa|\.fsa)$', '', args.input_file.strip().split('/')[-1])
    # complete the basename
    subfiles_basename += ".atomic_{0}_vs_{1}".format(args.program, database_name)  
    
    # Define default tmp output directory and add the pid
    output_directory = "{0}".format(args.output)
    if not args.output:
        args.output = subfiles_basename
        output_directory = "{0}.pid{1}".format(subfiles_basename, os.getpid())
    # Remove the last backslash if any
    if output_directory[-1] == '/':
        output_directory = output_directory[:-1]
    sys.stdout.write("OUTDIR: {outdir}\n".format(outdir=output_directory))
    
    # Create the output directory if it doesn't exists
    try:
        if not os.path.exists(output_directory):
            os.makedirs(output_directory)
    except OSError:
        sys.stderr.write("\nERROR: [Outdir] {0} cannot be created\n\n".format(output_directory))
        raise
    
    ##### Step n°1: Split multifasta file #####
    # Compute subfile sequence number if a subfile number is defined
    if args.nb_file:
        # Count sequences number
        input_file_nb_seq = int(subprocess.check_output("grep -c '>' {0}".format(args.input_file), shell=True).strip())
        sys.stdout.write("\nINFO: {0} sequences in the query\n".format(input_file_nb_seq))
        # Compute subfiles sequence number
        args.nb_seq = (input_file_nb_seq / float(args.nb_file))
        # 10 seq / 2 files --> 5 seq/subfile || 10 seq / 3 files --> 3.333 --> 4 seq/subfiles
        if args.nb_seq % 1 != 0:
            args.nb_seq += 1
        args.nb_seq = int(args.nb_seq)
    
    total_subfiles_nb = args.nb_file
    # If step 1 has to be run
    if args.steps.count('1'):
        sys.stdout.write('\nINFO: Splitting input file...\n')
        total_subfiles_nb = split_multifasta_by_increment(input_fh, output_directory, subfiles_basename, args.nb_seq, args.verbose)
        sys.stdout.write("\nINFO: The input file was splitted into {0} subfiles\n".format(total_subfiles_nb))
    
    ##### Step n°2: Generate SGE qsub script #####
    # Define default email
    sge_email = args.email
    if not args.email:
        current_user_name = getpass.getuser()
        sge_email = "{0}@sb-roscoff.fr".format(current_user_name)
    
    # Define SGE mail verbosity
    sge_email_options = 'a'
    if args.send_mails:
        sge_email_options = 'bea'
    
    # Define output files extension
    output_file_extension = 'tab'
    if args.output_format == 'xml':
        output_file_extension = 'xml'
    elif args.output_format == 'pairwise':
        output_file_extension = 'txt'
    
    # Define Blast output format    
    if args.output_format == 'tabular':
        blast_output_format = '6'
    elif args.output_format == 'pairwise':
        blast_output_format = '0'
    elif args.output_format == 'xml':
        blast_output_format = '5'
    elif args.output_format == 'extended':
        # Custom fields can be added to the standard tabular format
        blast_output_format = '\'6 std qlen slen\''
    
    # Open SGE qsub script
    sge_qsub_script_path = "{0}/qsub.{1}.sh".format(output_directory, subfiles_basename)
    try:
        sge_qsub_script_fh = open(sge_qsub_script_path, 'w')
    except:
        sys.stderr.write("\nERROR: {0} cannot be created\n\n".format(sge_qsub_script_path))
        raise
    
    # Build the Blast command line
    blast_command_line = "{program}".format(program=blast_program)
    blast_command_line += " -query {outdir}/$SGE_TASK_ID.{subname}.fasta".format(outdir=output_directory,
                                                                                 subname=subfiles_basename)
#     blast_command_line += " -out {outdir}/$JOB_ID.$SGE_TASK_ID.{subname}.{extension}".format(outdir=output_directory,
#                                                                                              subname=subfiles_basename,
#                                                                                              extension=output_file_extension)
    blast_command_line += " -out {outdir}/$JOB_ID.$SGE_TASK_ID.{subname}.{extension}".format(outdir=output_directory,
                                                                                             subname=subfiles_basename,
                                                                                             extension=output_file_extension)
    blast_command_line += " -db {database}".format(database=database_path)
    blast_command_line += " -evalue {evalue}".format(evalue=args.evalue)
    blast_command_line += " -outfmt {outformat}".format(outformat=blast_output_format)
    if args.output_format == 'pairwise':
        blast_command_line += " -num_descriptions {maxtarget} -num_alignments {maxtarget}".format(maxtarget=args.max_target_seqs)
    else:
        blast_command_line += " -max_target_seqs {maxtarget}".format(maxtarget=args.max_target_seqs)
    blast_command_line += " -num_threads 1 {options}".format(options=args.options.strip())
    
    # Write SGE qsub script, if step 2 has to be run
    if args.steps.count('2'):
        sge_qsub_script_fh.write("""
#!/bin/bash

# Shell to use for the job execution
#$ -S /bin/bash

#$ -o /dev/null
#$ -e /dev/null

# Run from the current working directory
#$ -cwd

# User to inform by email
#$ -M {email}
 
# Email at job (b)egin, (e)nd, (a)bort and (s)uspend
#$ -m {mailoptions}
 
# Export de toutes les variables d'environnement
#$ -V

{commandline} > {outdir}/$JOB_ID.$SGE_TASK_ID.{subname}.out 2> {outdir}/$JOB_ID.$SGE_TASK_ID.{subname}.err
        """.format(outdir=output_directory, subname=subfiles_basename, email=sge_email, 
                   mailoptions=sge_email_options, commandline=blast_command_line))
        sge_qsub_script_fh.close()
        if args.verbose:
            sys.stdout.write("\nINFO: SGE qsub script was written to {0}\n".format(sge_qsub_script_path))
    
    ##### Step n°3: Submit SGE array-job #####
    # Build the SGE qsub command line
    sge_qsub_command_line = "qsub -q {queue}".format(queue=args.queue)
    # batch option
    if args.batch == 'all':
        args.batch = "1-{subfilenb}".format(subfilenb=total_subfiles_nb)
    else:
        args.no_cleanup = True
    sge_qsub_command_line += " -t {batch}".format(batch=args.batch)
    # cpu number
    if args.cpu != 'all':
        sge_qsub_command_line += " -tc {cpu}".format(cpu=args.cpu)
    # sync option
    if not args.dont_wait:
        sge_qsub_command_line += ' -sync yes'
    # job name
    sge_qsub_command_line += " -N at_{program}_{subname}".format(program=args.program, subname=subfiles_basename)
    # SGE qsub script file
    sge_qsub_command_line += " {qsub_script}".format(qsub_script=sge_qsub_script_path)
    
    # Submit the SGE qsub, if step 3 has to be run
    if args.steps.count('3'):
        sys.stdout.write('\nINFO: Running job-array on SGE...\n')
        sge_qsub_out, sge_qsub_err, sge_qsub_errcode = run_shell(sge_qsub_command_line, args.verbose)
        global_error_code += sge_qsub_errcode
        # Check SGE qsub error code
        if sge_qsub_errcode > 0:
            sys.stderr.write('\nWARNING: SGE job-array exited with a positive error code. There might be an error, please check your outputs. Reverting to no_cleanup\n')
            args.no_cleanup = True
        # Get the sge job-array id
        sge_job_id = re.findall('Your job-array (\d+).\d+-\d+:.*', sge_qsub_out.split('\n')[0])[0]        
        if not args.dont_wait:
            # For pairwise output, change concatenation command to insert empty lines between subfiles
            cat_program = 'cat'
            if args.output_format == 'pairwise':
                cat_program = "awk 'FNR==1{print \"\\n---------------------\\n\"}1'"
            # Concatenate output, STDOUT, and STDERR subfiles
            for extension in (output_file_extension, 'err', 'out'):
                cat_command_line= "{program} {outdir}/{jobid}.*.{extension} > {output}.{jobid}.{extension}".format(program=cat_program,
                                                                                                                   outdir=output_directory,
                                                                                                                   jobid=sge_job_id,
                                                                                                                   extension=extension,
                                                                                                                   output=args.output)
                cat_out, cat_err, cat_errcode = run_shell(cat_command_line, args.verbose)
                # Check cat error code
                if cat_errcode > 0:
                    sys.stderr.write('\nWARNING: Concatenation might not have been successfull, please check output subfiles. Reverting to no_cleanup\n')
                    args.no_cleanup = True
                global_error_code += cat_errcode
            sys.stdout.write('\nINFO: Subfiles were concatenated\n')
            sys.stdout.write("\nOUTPUT: {output}.{jobid}.{extension}\n".format(output=args.output,
                                                                               jobid=sge_job_id,
                                                                               extension=output_file_extension))
            if not args.no_cleanup:
                # Delete temporary directory
                rm_command_line = "rm -rf {outdir}".format(outdir=output_directory)
                rm_out, rm_err, rm_errcode = run_shell(rm_command_line, args.verbose)
                # Check rm error code
                if rm_errcode > 0:
                    sys.stderr.write('\nWARNING: Temporary directory deletion might not have been completed\n')
                global_error_code += rm_errcode
                sys.stdout.write('\nINFO: Temporary files were cleaned\n')
    
    # Exit the program
    if global_error_code == 0:
        sys.stdout.write('\nINFO: Atomic Blast+ was completed without incidents. Now it\'s time to work ^^\n')
    else:
        sys.stderr.write('\nWARNING: At least one command did not end well. Sorry for the inconvenience, but thanks anyway for all the fish\n')
        
    sys.stdout.write("\n<3 <3 <3 <3 <3 <3 <3 <3 <3 <3 <3\nThank you for using Atomic Blast+\n\n")
    
    return 0

if __name__ == '__main__':
    args = get_args()
    main(args)

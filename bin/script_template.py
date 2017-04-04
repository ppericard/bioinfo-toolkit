#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

"""

import os
import sys
import argparse
import subprocess
import time
import logging
import multiprocessing
#~ import humanfriendly

# Create logger
logger = logging.getLogger(__name__)

# Get program filename, filepath and directory
program_filepath = os.path.realpath(sys.argv[0])
program_dirpath, program_filename = os.path.split(program_filepath)
program_name, program_extension = os.path.splitext(program_filename)


class DefaultHelpParser(argparse.ArgumentParser):
    """
    This is a slightly modified argparse parser to display the full help
    on parser error instead of only usage
    """
    def error(self, message):
        sys.stderr.write('\nError: %s\n\n' % message)
        self.print_help()
        sys.exit(2)


def parse_arguments():
    """
    Parse the command line, and check if arguments are correct
    """
    # Initiate argument parser
    parser = DefaultHelpParser(description='Program description',
                               # to precisely format help display
                               formatter_class=lambda prog: argparse.HelpFormatter(prog, width=120, max_help_position=80))

    # Main parameters
    group_main = parser.add_argument_group('Main parameters')
    # -ifh / --input_file_handle
    group_main.add_argument('-ifh', '--input_file_handle',
                            action = 'store',
                            metavar = 'INFH',
                            type = argparse.FileType('r'),
                            required = True,
                            help = 'Input file, opened and managed by argparse')
    # -ifs / --input_file_string
    group_main.add_argument('-ifs', '--input_file_string',
                            action = 'store',
                            metavar = 'INFS',
                            type = str,
                            default = '/dev/null',
                            help = 'Input file, as a string. '
                                   'Has to be dealt with manually. '
                                   'Default is $(default)s')
    # -ofh / --output_file_handle
    group_main.add_argument('-ofh', '--output_file_handle',
                            action = 'store',
                            metavar = 'OUTFH',
                            type = argparse.FileType('w'),
                            required = True,
                            help = 'Output file, opened and managed by argparse')
    # -d / --output_directory
    group_main.add_argument('-d', '--output_directory',
                            action = 'store',
                            metavar = 'OUTDIR',
                            type = str,
                            default = 'default_outdir',
                            help = 'Output directory, as a string. '
                                   'Has to be dealt with manually. '
                                   'Default is $(default)s')
    # -v / --verbose
    group_main.add_argument('-v', '--verbose',
                            action = 'store_true',
                            help = 'Increase verbosity')

    # Performance parameters
    group_perf = parser.add_argument_group('Performance')
    # --cpu
    group_perf.add_argument('--cpu',
                            action = 'store',
                            metavar = 'CPU',
                            type = int,
                            default = 1,
                            help = 'Max number of CPU to use. '
                                   'Default is %(default)s cpu')
    #~ # --max_memory
    #~ group_perf.add_argument('--max_memory',
                            #~ action = 'store',
                            #~ metavar = 'MAXMEM',
                            #~ type = str,
                            #~ default = '4 GB',
                            #~ help = 'Maximum memory to use. '
                                   #~ 'Default is %(default)s. '
                                   #~ 'Most size notations can be used.')

    # Advanced parameters
    group_adv = parser.add_argument_group('Advanced parameters')
    # --keep_tmp
    group_adv.add_argument('--keep_tmp',
                            action = 'store_true',
                            help = 'Do not remove tmp files')

    # Debug
    group_debug = parser.add_argument_group('Debug parameters')
    # --debug
    group_debug.add_argument('--debug',
                             action = 'store_true',
                             help = 'Output debug infos')

    #
    args = parser.parse_args()

    # Set debug parameters
    if args.debug:
        args.verbose = True
        args.keep_tmp = True

    # Get absolute path for all arguments
    args.input_file_string = os.path.abspath(args.input_file_string)
    args.output_directory = os.path.abspath(args.output_directory)

    #
    return args


def print_intro(args):
    """
    Print the introduction
    """

    sys.stderr.write("""
#################################
         Program name
#################################\n\n""")

    # Retrieve complete cmd line
    cmd_line = '{binpath} '.format(binpath = program_filename)

    # Verbose
    if args.verbose:
        cmd_line += '--verbose '

    # Advanced parameters
    if args.debug:
        cmd_line += '--debug '

    if args.keep_tmp:
        cmd_line += '--keep_tmp '

    # Performance
    cmd_line += '--cpu {0} '.format(args.cpu)
    #~ cmd_line += '--max_memory {0} '.format(args.max_memory)

    # Main parameters
    cmd_line += '--out_dir {0}'.format(args.out_dir)
    cmd_line += '--ref_db {0} '.format(args.ref_db)
    cmd_line += '--input_fastx {0}'.format(args.input_fastx)

    # Print cmd line
    sys.stderr.write('CMD: {0}\n\n'.format(cmd_line))

    return 0


if __name__ == '__main__':

    # Set global t0
    global_t0_wall = time.time()

    # Init global error code
    global_error_code = 0

    # Arguments parsing
    args = parse_arguments()

    # Set logging
    # create console handler
    ch = logging.StreamHandler()
    #
    if args.debug:
        logger.setLevel(logging.DEBUG)
        # create formatter for debug level
        formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    else:
        if args.verbose:
            logger.setLevel(logging.INFO)
        else:
            logger.setLevel(logging.WARNING)
        # create default formatter
        formatter = logging.Formatter('%(levelname)s - %(message)s')
    # add the formatter to the console handler
    ch.setFormatter(formatter)
    # add the handler to logger
    logger.addHandler(ch)

    # Print intro infos
    if args.verbose:
        print_intro(args)

    #
    sys.stderr.write('\n')

    #######
    # Exit

    logger.info('-- Program complete --')

    if global_error_code > 0:
        logger.warning('Problems might have happened during program execution. Please check log above')
    else:
        logger.debug('Execution completed in {0:.2f} seconds'.format(time.time() - global_t0_wall))

    sys.stderr.write('\n')

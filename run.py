#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import os
import sys
import tempfile

import glacon


def main():
    # filename = '/home/mblum/work/BingRen/AdrenalGland.gla'
    #
    # with glacon.Glacon(filename, 'r') as fh:
    #     result = fh.query(('chrun_gl000217', 0, 3100000), resolution=1000000)
    #     for chrom1, chrom2, rows, cols, values in result:
    #         for i, j, v in zip(rows, cols, values):
    #             print(chrom1, chrom2, i, j, v)
    #
    #
    #
    #     # # chr1:0-chrM:16,571 & chr4:7,862,541-chr15:94,379,315 [offset 0,4306017:0,0]
    #     # result = fh.query(('chr1', 0, 'chrM', 0), ('chr4', 7000000, 'chr15', 94000000), resolution=500000)
    #     # #result = fh.query(('chr1', 0, 3100000), resolution=1000000)
    #     #
    #     # for chrom1, chrom2, rows, cols, values in result:
    #     #     for i, j, v in zip(rows, cols, values):
    #     #         print(chrom1, chrom2, i, j, v)
    #
    #
    #
    # return


    parser = argparse.ArgumentParser(description='Storage format for Hi-C data')
    subparsers = parser.add_subparsers(dest='command', help='command')
    subparsers.required = True

    parser_create = subparsers.add_parser('store', help='Store Hi-C contact maps')
    parser_create.add_argument('assembly', metavar='assembly', help='UCSC assembly (e.g. hg19) or chrom sizes file')
    parser_create.add_argument('output', help='output file')
    parser_create.add_argument('-f', dest='bed_files', help='fragment BED files', nargs='+', required=True)
    parser_create.add_argument('-m', dest='mat_files', help='matrix files', nargs='+', required=True)

    parser_create.add_argument('--buffersize', type=int, metavar='INT', default=0, help='buffer size (default: 0)')
    parser_create.add_argument('-p', dest='processes', type=int, metavar='INT', default=1, help='process threads (default: 1)')
    parser_create.add_argument('-q', '--quiet', action='store_true', default=False,
                               help='do not print progress messages (default: off)')
    parser_create.add_argument('--tmp', default=tempfile.gettempdir(),
                               help='temporary directory (default: {})'.format(tempfile.gettempdir()))
    parser_create.add_argument('--no-agg', dest='no_agg', action='store_true', default=False,
                               help='do not aggregate maps (default: off)')

    args = parser.parse_args()

    if args.command == 'store':
        for filename in args.bed_files + args.mat_files:
            if not os.path.isfile(filename):
                sys.stderr.write('cannot access {}: no such file or directory\n'.format(filename))
                exit(1)

        if len(args.bed_files) != len(args.mat_files):
            sys.stderr.write('different number of fragment files ({}) '
                             'and matrix files ({})\n'.format(len(args.bed_files), len(args.mat_files)))
            exit(1)

        glacon.create(args.assembly, args.bed_files, args.mat_files, args.output,
                      buffersize=args.buffersize,
                      processes=args.processes,
                      verbose=(not args.quiet),
                      aggregate=(not args.no_agg),
                      tmp=args.tmp)


if __name__ == '__main__':
    main()
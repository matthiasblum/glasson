#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import os
import sys

import glacon


def main():
    # url = 'http://hglab.org/data/SRX2179260.gla'
    # filename = '/home/mblum/Dropbox/Perso/labs/igbmc_2013_strasbourg/HiC/Ren/SRX2179260/SRX2179260.gla'
    # if not os.path.isfile(filename):
    #     filename = '/home/matthias/Dropbox/Perso/labs/igbmc_2013_strasbourg/HiC/Ren/SRX2179260/SRX2179260.gla'
    #
    # fh = glasson.open(filename)
    #
    # # chr1:0-chrM:16,571 & chr4:7,862,541-chr15:94,379,315 [offset 0,4306017:0,0]
    # result = fh.query(('chr1', 0, 'chrM', 0), ('chr4', 7000000, 'chr15', 94000000), resolution=500000)
    # #result = fh.query(('chr1', 0, 3100000), resolution=1000000)
    #
    # for chrom1, chrom2, rows, cols, values in result:
    #     for i, j, v in zip(rows, cols, values):
    #         print(chrom1, chrom2, i, j, v)
    #
    #
    #
    # fh.close()
    # return


    parser = argparse.ArgumentParser(description='Storage format for Hi-C data')
    subparsers = parser.add_subparsers(dest='command', help='command')
    subparsers.required = True

    parser_create = subparsers.add_parser('store', help='Store Hi-C contact maps')
    parser_create.add_argument('chrom_sizes', metavar='chrom.sizes', help='UCSC chromosome sizes file')
    parser_create.add_argument('output', help='output file')
    parser_create.add_argument('-f', dest='bed_files', help='fragment BED files', nargs='+', required=True)
    parser_create.add_argument('-m', dest='mat_files', help='matrix files', nargs='+', required=True)

    parser_create.add_argument('--buffersize', type=int, metavar='INT', default=0, help='buffer size')
    parser_create.add_argument('-p', dest='processes', type=int, metavar='INT', default=1, help='process threads')
    parser_create.add_argument('-q', dest='quiet', action='store_true', default=False, help='do not print progress messages')

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

        glacon.create(args.chrom_sizes, args.bed_files, args.mat_files, args.output,
                      buffersize=args.buffersize,
                      processes=args.processes,
                      verbose=(not args.quiet),
                      tmp='tmp')


if __name__ == '__main__':
    main()
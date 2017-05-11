#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import os
import sys


import gla





def create_glasson(chrom_sizes_file, bed_files, mat_files, labels, output, buffersize=0):
    fh = gla.open(output, mode='w')
    chrom_sizes = gla.load_chrom_sizes(chrom_sizes_file)
    fh.set_chrom_sizes(chrom_sizes)

    for bed_file, mat_file, label in zip(bed_files, mat_files, labels):
        cmap = gla.ContactMap(bed_file, mat_file, chrom_sizes, verbose=False, buffersize=buffersize, tmpdir='.')
        fh.add(cmap, label)


    fh.write()



def main():
    parser = argparse.ArgumentParser(description='Storage format for Hi-C data')
    subparsers = parser.add_subparsers(dest='command', help='command')
    subparsers.required = True

    parser_create = subparsers.add_parser('store', help='Store Hi-C contact maps')
    parser_create.add_argument('chrom_sizes', metavar='chrom.sizes', help='UCSC chromosome sizes file')
    parser_create.add_argument('output', help='output file')
    parser_create.add_argument('-f', dest='bed_files', help='fragment BED files', nargs='+', required=True)
    parser_create.add_argument('-m', dest='mat_files', help='matrix files', nargs='+', required=True)
    parser_create.add_argument('-l', dest='labels', help='matrix labels', nargs='*')

    parser_create.add_argument('--buffersize', type=int, metavar='INT', default=10000000, help='buffer size')

    args = parser.parse_args()

    if args.command == 'store':
        for filename in [args.chrom_sizes] + args.bed_files + args.mat_files:
            if not os.path.isfile(filename):
                sys.stderr.write('cannot access {}: no such file or directory\n'.format(filename))
                exit(1)

        if len(args.bed_files) != len(args.mat_files):
            sys.stderr.write('different number of fragment files ({}) '
                             'and matrix files ({})\n'.format(len(args.bed_files), len(args.mat_files)))
            exit(1)

        labels = [None] * len(args.bed_files)
        if args.labels:
            if len(args.labels) > len(args.bed_files):
                sys.stderr.write('too many labels: '
                                 'expected {}, got {}\n'.format(len(args.labels), len(args.bed_files)))

            for i, e in enumerate(args.labels):
                try:
                    x, label = e.split(':', 1)
                    x = int(x) - 1
                except ValueError:
                    if len(args.labels) == len(args.bed_files):
                        labels[i] = e
                    else:
                        sys.stderr.write('invalid labels\n')  # todo: better error message
                        exit(1)
                else:
                    labels[x] = label

        create_glasson(args.chrom_sizes , args.bed_files, args.mat_files, labels, args.output, buffersize=args.buffersize)





if __name__ == '__main__':
    main()
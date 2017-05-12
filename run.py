#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import os
import sys


import gla as glasson


def create_glasson(chrom_sizes_file, bed_files, mat_files, labels, output, buffersize=0):
    fh = glasson.open(output, mode='w')
    chrom_sizes = glasson.load_chrom_sizes(chrom_sizes_file)
    fh.chrom_sizes(chrom_sizes)

    for bed_file, mat_file, label in zip(bed_files, mat_files, labels):
        cmap = glasson.ContactMap(bed_file, mat_file, chrom_sizes)
        fh.add(cmap, label)

    fh.write(verbose=True, buffersize=buffersize, tmpdir='.')
    fh.close()


def create(chrom_sizes_file, bed_files, mat_files, mat_labels, output):
    chrom_sizes = glasson.load_chrom_sizes(chrom_sizes_file)

    gla = glasson.Glasson(output, mode='w', chrom_sizes=chrom_sizes)

    for bed, mat, label in zip(bed_files, mat_files, mat_labels):
        gla.add_mat(bed, mat, label)
        break

    gla.freeze(processes=4, verbose=True, buffersize=10000000, tmpdir='.')

    gla.close()


def main():
    # url = 'http://hglab.org/data/SRX2179260.gla'
    # filename = '/home/mblum/Dropbox/Perso/labs/igbmc_2013_strasbourg/HiC/Ren/SRX2179260/SRX2179260.gla'
    # if not os.path.isfile(filename):
    #     filename = '/home/matthias/Dropbox/Perso/labs/igbmc_2013_strasbourg/HiC/Ren/SRX2179260/SRX2179260.gla'
    #
    # fh = gla.open(filename)
    #
    # # chr1:0-chrM:16,571 & chr4:7,862,541-chr15:94,379,315 [offset 0,4306017:0,0]
    # fh.query(('chr1', 0, 'chrM', 0), ('chr4', 7e6, 'chr15', 94e6), resolution=1000000)
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

        #create_glasson(args.chrom_sizes , args.bed_files, args.mat_files, labels, args.output, buffersize=args.buffersize)
        create(args.chrom_sizes, args.bed_files, args.mat_files, labels, args.output)


if __name__ == '__main__':
    main()
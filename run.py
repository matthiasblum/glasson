#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import bisect
import gzip
import os
import sys
import tempfile

try:
    import cPickle as pickle
except ImportError:
    import pickle

import gla


def _open_text(filename, mode='rt'):
    if filename.lower().endswith('.gz'):
        return gzip.open(filename, mode=mode)
    else:
        return open(filename, mode=mode)


def load_chrom_sizes(filename):
    chrom_sizes = {}

    fh = _open_text(filename, mode='rt')

    for i, line in enumerate(fh):
        try:
            chrom, size = line.rstrip().split()
            size = int(size)
        except ValueError:
            sys.stderr.write('invalid format at line {} in file {}\n'.format(i + 1, filename))
            fh.close()
            exit(1)
        else:
            chrom_sizes[chrom] = size

    fh.close()
    return chrom_sizes


def load_bed(filename, chrom_sizes):
    fragments = {}
    offsets = {}
    _chrom = None
    bin_size = 0
    c = 0

    fh = _open_text(filename, mode='rt')

    for i, line in enumerate(fh):
        try:
            # file is expected to be sorted by chrom + frag position
            chrom, start, end, frag_id = line.rstrip().split()
            start = int(start)
            end = int(end)
            frag_id = int(frag_id)
            chrom_size = chrom_sizes[chrom]
        except ValueError:
            sys.stderr.write('invalid format at line {} in file {}\n'.format(i + 1, filename))
            fh.close()
            exit(1)
        except KeyError:
            sys.stderr.write('unknown chromosome \'{}\' at line {} in file {}\n'.format(chrom, i + 1, filename))
            fh.close()
            exit(1)
        else:
            fragments[frag_id] = (chrom, start)
            if chrom != _chrom:
                offsets[chrom] = i
            _chrom = chrom

            if end < chrom_size:
                bin_size += end - start
                c += 1

    fh.close()

    return fragments, offsets, bin_size // c if bin_size // c == float(bin_size) / c else None


def _clean_pickles(files):
    for chrom1 in files:
        for chrom2 in files[chrom1]:
            os.unlink(files[chrom1][chrom2])


def load_mat(filename, fragments, offsets, buffersize=None, tmpdir=tempfile.gettempdir()):
    pk_files = {}
    matrices = {}
    cnt = 0

    fh = _open_text(filename, mode='rt')

    for x, line in enumerate(fh):
        try:
            fields = line.rstrip().split()
            i = int(fields[0])
            j = int(fields[1])
            v = float(fields[2])
        except (IndexError, ValueError):
            sys.stderr.write('invalid format at line {} in file {}'.format(x + 1, line))
            fh.close()
            _clean_pickles(pk_files)
            exit(1)

        try:
            chrom1 = fragments[i][0]
        except KeyError:
            sys.stderr.write('invalid fragment ID ({}) at line {} in file {}\n'.format(i, x + 1, filename))
            fh.close()
            _clean_pickles(pk_files)
            exit(1)

        try:
            chrom2 = fragments[j][0]
        except KeyError:
            sys.stderr.write('invalid fragment ID ({}) at line {} in file {}\n'.format(j, x + 1, filename))
            fh.close()
            _clean_pickles(pk_files)
            exit(1)

        i -= offsets[chrom1]
        j -= offsets[chrom2]

        if chrom1 not in matrices:
            matrices[chrom1] = {chrom2: [(i, j, v)]}
        elif chrom2 not in matrices[chrom1]:
            matrices[chrom1][chrom2] = [(i, j, v)]
        else:
            matrices[chrom1][chrom2].append((i, j, v))

        cnt += 1
        if not cnt % 1000000:
            sys.stderr.write('{} contacts loaded\n'.format(cnt))

        if buffersize and not cnt % buffersize:
            for chrom1 in matrices:
                for chrom2 in matrices[chrom1]:
                    data = matrices[chrom1][chrom2]

                    if chrom1 in pk_files and chrom2 in pk_files[chrom1]:
                        with open(pk_files[chrom1][chrom2], 'rb') as pfh:
                            data += pickle.load(pfh)
                    else:
                        fd, path = tempfile.mkstemp(dir=tmpdir)
                        os.close(fd)

                        if chrom1 not in pk_files:
                            pk_files[chrom1] = {chrom2: path}
                        else:
                            pk_files[chrom1][chrom2] = path

                    with open(pk_files[chrom1][chrom2], 'wb') as pfh:
                        pickle.dump(data, pfh)

                    matrices[chrom1][chrom2] = None  # Can I haz memory back?

            matrices = {}

    fh.close()

    for chrom1 in matrices:
        for chrom2 in matrices[chrom1]:
            data = matrices[chrom1][chrom2]

            if chrom1 in pk_files and chrom2 in pk_files[chrom1]:
                with open(pk_files[chrom1][chrom2], 'rb') as pfh:
                    data += pickle.load(pfh)
            else:
                fd, path = tempfile.mkstemp(dir=tmpdir)
                os.close(fd)

                if chrom1 not in pk_files:
                    pk_files[chrom1] = {chrom2: path}
                else:
                    pk_files[chrom1][chrom2] = path

                with open(pk_files[chrom1][chrom2], 'wb') as pfh:
                    pickle.dump(data, pfh)

                matrices[chrom1][chrom2] = None

    return pk_files


def create_glasson(chrom_sizes_file, bed_files, mat_files, labels, output, buffersize=None, tmpdir=tempfile.gettempdir()):
    fh = gla.open(output, mode='w')
    chrom_sizes = load_chrom_sizes(chrom_sizes_file)
    fh.set_chrom_sizes(chrom_sizes)

    for bed_file, mat_file, label in zip(bed_files, mat_files, labels):
        fragments, offsets, bin_size = load_bed(bed_file, chrom_sizes)
        pk_files = load_mat(mat_file, fragments, offsets, buffersize=buffersize, tmpdir=tmpdir)

        break


def main():
    parser = argparse.ArgumentParser(description='Storage format for Hi-C data')
    subparsers = parser.add_subparsers(dest='command', help='command')
    subparsers.required = True

    parser_create = subparsers.add_parser('store', help='Store Hi-C contact maps')
    parser_create.add_argument('-g', dest='chrom_sizes', metavar='chrom.sizes', help='UCSC chromosome sizes file')
    parser_create.add_argument('-f', dest='bed_files', help='fragment BED files', nargs='+', required=True)
    parser_create.add_argument('-m', dest='mat_files', help='matrix files', nargs='+', required=True)
    parser_create.add_argument('-l', dest='labels', help='matrix labels', nargs='*')
    parser_create.add_argument('output', help='output file')
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

        create_glasson(args.chrom_sizes , args.bed_files, args.mat_files, labels, args.output, buffersize=args.buffersize, tmpdir='.')





if __name__ == '__main__':
    main()
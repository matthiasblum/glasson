#!/usr/bin/env python
# -*- coding: utf-8 -*-

import gzip
import os
import struct
import sys
import zlib

try:
    import urllib.request as request
except ImportError:
    import urllib2 as request


_SIGNATURE = 'GLA'
_VERSION = '0.1'


def load_index(filename):
    index = {}

    with open(filename, 'rt') as fh:
        for i, line in enumerate(fh):
            try:
                fields = line.rstrip().split()
                chrom = fields[0]
                resolution = int(fields[1])
                matrix = fields[2]
            except (IndexError, ValueError):
                sys.stderr.write('glasson: {}: invalid format at line {}\n'.format(filename, i + 1))
                exit(1)

            if not os.path.isfile(matrix):
                sys.stderr.write('glasson: {}: no such file or directory\n'.format(matrix))
                exit(1)
            elif chrom not in index:
                index[chrom] = {resolution: matrix}
            else:
                index[chrom][resolution] = matrix

    return index


def load_chrom_lens(filename):
    chrom_lens = {}

    with open(filename) as fh:
        for i, line in enumerate(fh):
            try:
                chrom, size = line.rstrip().split('\t')
                chrom_lens[chrom] = int(size)
            except ValueError:
                sys.stderr.write('glasson: {}: invalid format at line {}\n'.format(filename, i + 1))
                exit(1)

    return chrom_lens


class Matrix:
    def __init__(self, filename, chromsize, resolution):
        self.filename = filename
        self.chromsize = chromsize
        self.resolution = resolution
        self.size = (chromsize + resolution - 1) // resolution
        self.col = []
        self.data = []
        self.nzrow = []

    def load(self):
        self.col = []
        self.data = []
        self.nzrow = []

        with gzip.open(self.filename, 'rt') if self.filename.lower().endswith('.gz') else open(self.filename,
                                                                                               'rt') as fh:
            for i, row in enumerate(fh):
                if i == self.size:
                    sys.stderr.write('glasson: {}: expected {} rows\n'.format(self.filename, self.size))
                    exit(1)

                col = row.rstrip().split('\t')
                if len(col) != self.size:
                    sys.stderr.write(
                        'glasson: {}: expected {} fields at line {}, found {}\n'.format(self.filename, self.size, i + 1,
                                                                                        len(col)))
                    exit(1)

                try:
                    col = [float(val) for val in col[i:]]
                except ValueError:
                    sys.stderr.write('glasson: {}: invalid float value at line {}\n'.format(self.filename, i + 1))
                    exit(1)
                else:
                    nnz = 0
                    for x in range(self.size - i):
                        if col[x]:
                            self.col.append(i + x)
                            self.data.append(col[x])
                            nnz += 1

                    self.nzrow.append(nnz)

    def iter(self):
        x = 0

        for i, nnz in enumerate(self.nzrow):
            if nnz:
                yield i, self.col[x:x + nnz], self.data[x:x + nnz]
                x += nnz

    def resize(self, new_resolution):
        new_size = (self.chromsize + new_resolution - 1) // new_resolution
        new_nzrows = [0] * new_size
        new_cols = []
        new_data = []
        rowdata = []

        x = 0
        prev_i = None
        for i, nnz in enumerate(self.nzrow):
            _i = i * self.resolution // new_resolution

            if _i != prev_i:
                if rowdata:
                    for xx, val in enumerate(rowdata):
                        if val:
                            new_cols.append(prev_i + xx)
                            new_data.append(val)
                            new_nzrows[prev_i] += 1

                rowdata = [0] * (new_size - _i)

            for xx in range(x, x + nnz):
                _j = self.col[xx] * self.resolution // new_resolution
                rowdata[_j - _i] += self.data[xx]

            prev_i = _i
            x += nnz

        for xx, val in enumerate(rowdata):
            if val:
                new_cols.append(prev_i + xx)
                new_data.append(val)
                new_nzrows[prev_i] += 1

        self.nzrow = new_nzrows
        self.col = new_cols
        self.data = new_data
        self.size = new_size
        self.resolution = new_resolution

        x = 0
        for i, nnz in enumerate(self.nzrow):
            for xx in range(x, x + nnz):
                j = self.col[xx]

                if i > j:
                    print(self.resolution, new_resolution, i, j)
                    exit()

            x += nnz


class File:
    def __init__(self, filename, mode='r'):
        if mode not in ('r', 'w'):
            raise ValueError("invalid mode: '{}'".format(mode))

        self.filename = filename
        self.mode = mode
        self.remote = False
        self.fh = None
        self.index = {}
        self.offset = None

    def __enter__(self):
        if self.mode == 'w':
            self.fh = open(self.filename, 'wb')
            self.fh.write(struct.pack('>3sQ', _SIGNATURE.encode('utf-8'), 0))
            self.offset = 11
        elif os.path.isfile(self.filename):
            self.fh = open(self.filename, 'rb')
            self.remote = False
            self._parseindex()
        else:
            self.remote = True
            self._parseindex()

        return self

    def _readindex(self):
        signature, offset, = struct.unpack('>3sQ', self.fh.read(11))
        self.fh.seek(offset)
        return self.fh.read()

    def _fetchheader(self):
        req = request.Request(self.filename, headers={'Range': 'bytes=0-10'})
        signature, offset = struct.unpack('>3sQ', request.urlopen(req).read())

        req = request.Request(self.filename, headers={'Range': 'bytes={}-'.format(offset)})
        return request.urlopen(req).read()

    def _parseindex(self):
        if self.remote:
            data = self._fetchheader()
        else:
            data = self._readindex()

        x = 0
        while True:
            if x == len(data):
                break

            chrom_len, l, n = struct.unpack('>3I', data[x:x + 12])
            x += 12
            chrom = data[x:x + l].decode('utf-8')
            x += l

            self.index[chrom] = {'size': chrom_len, 'dsets': {}}

            for i in range(n):
                resolution, n = struct.unpack('>2I', data[x:x + 8])
                x += 8

                rows = struct.unpack('>{}I'.format(n), data[x:x + 4 * n])
                x += 4 * n

                offsets = struct.unpack('>{}Q'.format(n), data[x:x + 8 * n])
                x += 8 * n

                self.index[chrom]['dsets'][resolution] = (rows, offsets)

    def _writeindex(self):
        for chrom, chrom_info in self.index.items():
            l = len(chrom)
            n = len(chrom_info['dsets'])
            self.fh.write(struct.pack('>3I{}s'.format(l), chrom_info['size'], l, n, chrom.encode('utf-8')))

            for resolution, row_offsets in chrom_info['dsets'].items():
                n = len(row_offsets) // 2

                rows = [row_offsets[i * 2] for i in range(n)]
                offsets = [row_offsets[i * 2 + 1] for i in range(n)]
                self.fh.write(struct.pack('>2I{0}I{0}Q'.format(n), resolution, n, *rows + offsets))

        self.fh.seek(3)
        self.fh.write(struct.pack('>Q', self.offset))

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.mode == 'w':
            self._writeindex()
            self.fh.close()
        elif not self.remote:
            self.fh.close()

    def _extract(self, start, end):
        if self.remote:
            req = request.Request(self.filename, headers={'Range': 'bytes={}-{}'.format(start, end - 1)})
            return request.urlopen(req).read()
        else:
            self.fh.seek(start)
            return self.fh.read(end - start)

    def query(self, chrom, start=0, end=0, resolution=None):
        chrom_info = self.index.get(chrom)

        if chrom_info is None:
            sys.stderr.write("glasson: {}: no chromosome '{}'\n".format(self.filename, chrom))
            exit(1)

        chrom_len = chrom_info['size']
        dsets = chrom_info['dsets']

        if not end:
            end = chrom_len - 1

        region_size = end - start
        for optimal_res in sorted(dsets):
            if optimal_res * 100 >= region_size:
                break

        if not resolution or resolution not in dsets:
            resolution = optimal_res

        start //= resolution
        end = (end + resolution - 1) // resolution

        rows, offsets = dsets[resolution]
        off_start = None
        off_end = None
        i = None

        for x, _row in enumerate(rows):
            if off_start is None or _row < start:
                i = _row
                off_start = offsets[x]
                off_end = offsets[x]
            elif _row >= end:
                off_end = offsets[x]
                break

        body = self._extract(off_start, off_end)
        l = len(body)
        x = 0
        row = []
        col = []
        data = []

        while x < l:
            n, zl = struct.unpack('>2I', body[x:x + 8])
            x += 8

            _data = struct.unpack('>{0}I{0}f'.format(n), zlib.decompress(body[x:x + zl]))
            x += zl

            for j in range(n):
                if _data[j] >= end:
                    break
                row.append(i)
                col.append(_data[j])
                data.append(_data[n + j])

            i += 1

        return row, col, data, resolution

    def write(self, chrom, chrom_len, resolution, mat):
        if chrom not in self.index:
            self.index[chrom] = {'size': chrom_len, 'dsets': {resolution: []}}
        else:
            self.index[chrom]['dsets'][resolution] = []

        for row, cols, data in mat.iter():
            n = len(cols)
            s = zlib.compress(struct.pack('>{0}I{0}f'.format(n), *cols + data))
            l = len(s)

            self.fh.write(struct.pack('>2I{}s'.format(l), n, l, s))
            self.index[chrom]['dsets'][resolution] += [row, self.offset]
            self.offset += 8 + l


def create(index_file, chrom_lens_file, output_file):
    index_matrices = load_index(index_file)
    chrom_lens = load_chrom_lens(chrom_lens_file)

    with File(output_file, 'w') as fh:
        for chrom, resolutions in index_matrices.items():
            chrom_len = chrom_lens.get(chrom)

            if chrom_len is None:
                sys.stderr.write("glasson: {} does not contain '{}'\n".format(chrom_lens_file, chrom))
                exit(1)

            for resolution, matrix_file in resolutions.items():
                sys.stderr.write('{:<8}\t{:>8}\n'.format(chrom, resolution))
                mat = Matrix(matrix_file, chrom_len, resolution)
                mat.load()
                fh.write(chrom, chrom_len, resolution, mat)

            while resolution < 1000000:
                new_res = resolution * 4

                if resolution < 100000 < new_res:
                    new_res = 100000
                elif new_res > 100000:
                    new_res = 1000000

                sys.stderr.write('{:<8}\t{:>8}\n'.format('', new_res))
                mat_size = (chrom_len + resolution - 1) // resolution

                if mat_size == 1:
                    break

                mat.resize(new_res)
                resolution = new_res
                fh.write(chrom, chrom_len, resolution, mat)


def read(filename, chrom, start, end, resolution=None):
    with File(filename, 'r') as fh:
        row, col, data, res = fh.query(chrom, start, end, resolution=resolution)

    for r, c, d in zip(row, col, data):
        print('{}\t{}\t{}'.format(r, c, d))


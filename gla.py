import builtins
import gzip
import os
import struct
import sys
import tempfile
import zlib

try:
    import cPickle as pickle
except ImportError:
    import pickle

try:
    range = xrange
except NameError:
    pass


_SIGNATURE = 'GLA'
_VERSION = 1


def _open_text(filename, mode='rt'):
    if filename.lower().endswith('.gz'):
        return gzip.open(filename, mode=mode)
    else:
        return builtins.open(filename, mode=mode)


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


def load_fragments(filename, chrom_sizes):
    fragments = []

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
            fragments.append((chrom, start, end, frag_id))

    fh.close()

    return fragments


class ContactMap:
    def __init__(self, bed_file, mat_file, chrom_sizes, verbose=False, buffersize=0, tmpdir=tempfile.gettempdir()):
        self._chrom_sizes = chrom_sizes
        self._frag_ids = {}
        self._chrom_frags = {}
        self._chrom_offsets = {}
        self._contacts = {}
        self._pk_files = {}
        self._bin_size = 0

        self._load_fragments(bed_file)
        self._load_contacts(mat_file, verbose=verbose, buffersize=buffersize, tmpdir=tmpdir)

    def _clean(self):
        for chrom1 in self._pk_files:
            for chrom2 in self._pk_files[chrom1]:
                os.unlink(self._pk_files[chrom1][chrom2])

    def __del__(self):
        self._clean()

    def _load_fragments(self, filename):
        frag_ids = {}
        chrom_frags = {}
        chrom_offsets = {}
        bin_size = 0
        cnt = 0

        fh = _open_text(filename, mode='rt')

        i = 0  # just in case the file is empty
        for i, line in enumerate(fh):
            try:
                # file is expected to be sorted by chrom + frag position
                chrom, start, end, frag_id = line.rstrip().split()
                start = int(start)
                end = int(end)
                frag_id = int(frag_id)
                chrom_size = self._chrom_sizes[chrom]
            except ValueError:
                sys.stderr.write('invalid format at line {} in file {}\n'.format(i + 1, filename))
                fh.close()
                exit(1)
            except KeyError:
                sys.stderr.write('unknown chromosome \'{}\' at line {} in file {}\n'.format(chrom, i + 1, filename))
                fh.close()
                exit(1)
            else:
                frag_ids[frag_id] = chrom

                if chrom not in chrom_frags:
                    chrom_frags[chrom] = [start]
                    chrom_offsets[chrom] = i
                else:
                    chrom_frags[chrom].append(start)

                if end < chrom_size:
                    cnt += 1
                    bin_size += end - start

        fh.close()

        if i + 1 != len(frag_ids):
            sys.stderr.write('expected {} unique fragment IDs, got {}\n'.format(i + 1, len(frag_ids)))
            exit(1)

        if float(bin_size) / cnt == bin_size // cnt:
            # fixed fragment size
            self._bin_size = bin_size // cnt
        else:
            # variable (most likely RE fragments)
            self._bin_size = 0

        self._frag_ids = frag_ids
        self._chrom_offsets = chrom_offsets
        self._chrom_frags = chrom_frags

    def _load_contacts(self, filename, verbose=False, buffersize=0, tmpdir=tempfile.gettempdir()):
        self._clean()
        self._pk_files = {}
        self._contacts = {}
        cnt = 0
        fh = _open_text(filename, mode='rt')

        for x, line in enumerate(fh):
            try:
                fields = line.rstrip().split()
                i = int(fields[0])
                j = int(fields[1])
                v = float(fields[2])
            except (IndexError, ValueError):
                sys.stderr.write('invalid format at line {} in file {}\n'.format(x + 1, line))
                fh.close()
                exit(1)

            if i > j or not v:  # todo: should we raise an error here?
                continue

            try:
                chrom1 = self._frag_ids[i]
                offset = self._chrom_offsets[chrom1]
            except KeyError:
                sys.stderr.write('invalid fragment ID ({}) at line {} in file {}\n'.format(i, x + 1, filename))
                fh.close()
                exit(1)
            else:
                i -= offset

            try:
                chrom2 = self._frag_ids[j]
                offset = self._chrom_offsets[chrom2]
            except KeyError:
                sys.stderr.write('invalid fragment ID ({}) at line {} in file {}\n'.format(i, x + 1, filename))
                fh.close()
                exit(1)
            else:
                j -= offset

            if chrom1 not in self._contacts:
                self._contacts[chrom1] = {chrom2: [(i, j, v)]}
            elif chrom2 not in self._contacts[chrom1]:
                self._contacts[chrom1][chrom2] = [(i, j, v)]
            else:
                self._contacts[chrom1][chrom2].append((i, j, v))

            cnt += 1
            if verbose and not cnt % 1000000:
                sys.stderr.write('{} contacts loaded\n'.format(cnt))

            if buffersize and not cnt % buffersize:
                for chrom1 in self._contacts:
                    for chrom2 in self._contacts[chrom1]:
                        data = self._contacts[chrom1][chrom2]

                        if chrom1 in self._pk_files and chrom2 in self._pk_files[chrom1]:
                            with builtins.open(self._pk_files[chrom1][chrom2], 'rb') as pfh:
                                data += pickle.load(pfh)
                        else:
                            fd, path = tempfile.mkstemp(dir=tmpdir)
                            os.close(fd)

                            if chrom1 not in self._pk_files:
                                self._pk_files[chrom1] = {chrom2: path}
                            else:
                                self._pk_files[chrom1][chrom2] = path

                        with builtins.open(self._pk_files[chrom1][chrom2], 'wb') as pfh:
                            pickle.dump(data, pfh)

                        self._contacts[chrom1][chrom2] = []

        fh.close()

        if buffersize:
            for chrom1 in self._contacts:
                for chrom2 in self._contacts[chrom1]:
                    data = self._contacts[chrom1][chrom2]

                    if chrom1 in self._pk_files and chrom2 in self._pk_files[chrom1]:
                        with builtins.open(self._pk_files[chrom1][chrom2], 'rb') as pfh:
                            data += pickle.load(pfh)
                    else:
                        fd, path = tempfile.mkstemp(dir=tmpdir)
                        os.close(fd)

                        if chrom1 not in self._pk_files:
                            self._pk_files[chrom1] = {chrom2: path}
                        else:
                            self._pk_files[chrom1][chrom2] = path

                        with builtins.open(self._pk_files[chrom1][chrom2], 'wb') as pfh:
                            pickle.dump(data, pfh)

                        self._contacts[chrom1][chrom2] = []

    def iter_contacts(self):
        for chrom1 in self._contacts:
            for chrom2 in self._contacts[chrom1]:
                if self._contacts[chrom1][chrom2]:
                    data = self._contacts[chrom1][chrom2]
                else:
                    with builtins.open(self._pk_files[chrom1][chrom2], 'rb') as fh:
                        data = pickle.load(fh)

                yield chrom1, chrom2, data

    def get_bin_size(self):
        return self._bin_size

    def iter_frags(self):
        for chrom in self._chrom_frags:
            yield chrom, self._chrom_frags[chrom]

    def get_chrom_frags(self, chrom):
        return self._chrom_frags[chrom]


class File:
    def __init__(self, filename, mode='r'):
        if mode not in ('r', 'w'):
            raise RuntimeError('invalid mode: \'{}\'\n'.format(mode))

        self._filename = filename
        self._mode = mode
        self._index = {}
        self._offset = None
        self._chrom_sizes = {}
        self._maps = []

        if self._mode == 'w':
            self._fh = builtins.open(self._filename, 'wb')
            self._remote = False
        elif os.path.isfile(self._filename):
            self._fh = builtins.open(self._filename, 'rb')
            self._remote = False
            # self._parseindex()
        elif self._filename.lower().startswith(('http://', 'https://')):
            self._fh = None
            self.remote = True
            # self._parseindex()
        else:
            raise RuntimeError('{} is not an existing file nor a valid HTTP URL'.format(self._filename))

    def add(self, cmap, name=None):
        self._maps.append((cmap, name))

    def set_chrom_sizes(self, chrom_sizes):
        self._chrom_sizes = chrom_sizes

    def write(self):
        offset = 0
        fh = self._fh

        # Header
        fh.write(struct.pack('>3s2H', _SIGNATURE.encode('utf-8'), _VERSION, len(self._chrom_sizes)))
        offset += 3 + 2 + 2

        chroms = list(self._chrom_sizes.keys())

        for chrom_name in chroms:
            l = len(chrom_name)
            fh.write(struct.pack('<HI{}s'.format(l), l, self._chrom_sizes[chrom_name], chrom_name.encode('utf-8')))
            offset += 2 + 4 + l

        self._maps.sort(key=lambda tp: tp[0].get_bin_size())
        for cmap, name in self._maps:
            bin_size = cmap.get_bin_size()
            if name:
                l = len(name)
                fh.write(struct.pack('<I{}sI'.format(l), l, name.encode('utf-8'), bin_size))
            else:
                l = 0
                fh.write(struct.pack('<2I'.format(l), l, bin_size))
            offset += 2 + l + 2

            if bin_size:
                for chrom, frags in cmap.iter_frags():
                    chrom_idx = chroms.index(chrom)
                    n_frags = len(frags)
                    fh.write(struct.pack('<H{}I'.format(n_frags+1), chrom_idx, n_frags, *frags))
                    offset = 2 + 4 + n_frags * 4

        _offset = offset
        fh.write(struct.pack('<Q', 0))

        cmap_offsets = []

        # Body
        for cmap, name in self._maps:
            _cmap_offsets = []
            # bin_size = cmap.get_bin_size()

            for chrom1, chrom2, contacts in cmap.iter_contacts():
                chrom1_idx = chroms.index(chrom1)
                chrom2_idx = chroms.index(chrom2)

                #
                # if bin_size:
                #     n_rows = (self._chrom_sizes[chrom1] + bin_size - 1) // bin_size
                # else:
                #     n_rows = len(cmap.get_chrom_frags(chrom1))

                x = 0
                row_offsets = [offset]
                cols = []
                data = []
                for i, j, v in contacts:
                    if i == x:
                        cols.append(j)
                        data.append(v)
                    elif i > x:
                        n = len(cols)
                        s = zlib.compress(struct.pack('<I{0}I{0}d'.format(n), n, *(cols + data)), 6)
                        l = len(s)
                        fh.write(struct.pack('<I{}s'.format(l), l, s))
                        offset += 4 + l

                        while x + 1 < i:
                            # Empty rows
                            x += 1
                            row_offsets.append(offset)
                            fh.write(struct.pack('<I', 0))
                            offset += 4

                        x += 1
                        row_offsets.append(offset)
                        cols = []
                        data = []

                _cmap_offsets.append(
                    (chrom1_idx, chrom2_idx, row_offsets)
                )

            cmap_offsets.append(_cmap_offsets)

        # Index
        for i, (cmap, name) in enumerate(self._maps):
            for chrom1_idx, chrom2_idx, row_offsets in cmap_offsets[i]:
                s = zlib.compress(struct.pack('<{}Q'.format(len(row_offsets))))



    def _close(self):
        try:
            self._fh.close()
        except AttributeError:
            pass
        finally:
            self._fh = None

    def __del__(self):
        self._close()


def open(filename, mode='r'):
    try:
        fh = File(filename, mode)
    except RuntimeError as e:
        sys.stderr.write('{}\n'.format(e.args[0]))
        fh = None
    finally:
        return fh

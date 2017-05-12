import builtins
import gzip
import logging
import os
import struct
import sys
import tempfile
from urllib.error import HTTPError
import zlib

if sys.version_info[0] == 3:
    import pickle
    import urllib.request as request
else:
    try:
        import cPickle as pickle
    except ImportError:
        import pickle

    import urllib2 as request

    range = xrange


logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s: %(message)s',
    datefmt='%y-%m-%d %H:%M:%S',
    stream=sys.stdout
)


_SIGNATURE = 'GLA'
_VERSION = 1
_COMPRESS_LEVEL = 6


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
            logging.critical('invalid format at line {} in file {}'.format(i + 1, filename))
            fh.close()
            exit(1)
        else:
            chrom_sizes[chrom] = size

    fh.close()
    return chrom_sizes


class ContactMap:
    def __init__(self, bed_file, mat_file, chrom_sizes):
        self._bed_file = bed_file
        self._mat_file = mat_file
        self._chrom_sizes = chrom_sizes
        self._frag_ids = {}
        self._chrom_frags = {}
        self._chrom_offsets = {}
        self._contacts = {}
        self._pk_files = {}
        self._bin_size = 0

        self._load_fragments()

    def clean(self):
        for chrom1 in self._pk_files:
            for chrom2 in self._pk_files[chrom1]:
                os.unlink(self._pk_files[chrom1][chrom2])
        self._pk_files = {}
        self._contacts = {}

    def __del__(self):
        self.clean()

    def _load_fragments(self):
        frag_ids = {}
        chrom_frags = {}
        chrom_offsets = {}
        bin_size = 0
        cnt = 0

        fh = _open_text(self._bed_file, mode='rt')

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
                logging.critical('invalid format at line {} in file {}'.format(i + 1, self._bed_file))
                fh.close()
                exit(1)
            except KeyError:
                logging.critical('unknown chromosome \'{}\' at line {} in file {}'.format(chrom, i + 1, self._bed_file))
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
            logging.critical('expected {} unique fragment IDs, got {}'.format(i + 1, len(frag_ids)))
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

    def load_contacts(self, verbose=False, buffersize=0, tmpdir=tempfile.gettempdir()):
        self.clean()
        cnt = 0
        fh = _open_text(self._mat_file, mode='rt')

        for x, line in enumerate(fh):
            try:
                fields = line.rstrip().split('\t')
                i = int(fields[0])
                j = int(fields[1])
                v = float(fields[2])
            except (IndexError, ValueError):
                logging.critical('invalid format at line {} in file {}'.format(x + 1, line))
                fh.close()
                exit(1)

            if i > j or not v:  # todo: should we raise an error here?
                continue

            try:
                chrom1 = self._frag_ids[i]
                offset = self._chrom_offsets[chrom1]
            except KeyError:
                logging.critical('invalid fragment ID ({}) at line {} in file {}'.format(i, x + 1, self._mat_file))
                fh.close()
                exit(1)
            else:
                i -= offset

            try:
                chrom2 = self._frag_ids[j]
                offset = self._chrom_offsets[chrom2]
            except KeyError:
                logging.critical('invalid fragment ID ({}) at line {} in file {}'.format(i, x + 1, self._mat_file))
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
            if verbose and not (cnt % 100000):
                logging.info('{} contacts loaded'.format(cnt))

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
            self._parseindex()
        elif self._filename.lower().startswith(('http://', 'https://')):
            self._fh = None
            self._remote = True
            self._parseindex()
        else:
            raise RuntimeError('{} is not an existing file nor a valid HTTP URL'.format(self._filename))

    def _fetch(self, offset, size=None):
        if size:
            range = 'bytes={}-{}'.format(offset, size-1)
        else:
            range = 'bytes={}-'.format(offset)

        req = request.Request(self._filename, headers={'Range': range})

        try:
            res = request.urlopen(req)
        except HTTPError:
            self.close()
            raise RuntimeError('glasson: file not found')
        else:
            return res.read()

    def _parseindex(self):
        if self._remote:
            data = self._fetch(0, 17)
        else:
            data = self._fh.read(17)

        signature, version, offset, blocksize = struct.unpack('<3sHQI', data)

        try:
            is_gla = signature.decode('utf-8') == 'GLA'
        except UnicodeDecodeError:
            is_gla = False

        if not is_gla:
            self.close()
            raise RuntimeError('glasson: not a valid glasson file')

        if self._remote:
            data = self._fetch(17, blocksize)
        else:
            data = self._fh.read(blocksize)

        n_chroms, = struct.unpack('<H', data[:2])
        x = 2

        chroms = []

        for _ in range(n_chroms):
            l, chrom_size = struct.unpack('<HI', data[x:x+6])
            x += 6

            chrom_name = data[x:x+l].decode('utf-8')
            x += l

            chroms.append({
                'name': chrom_name,
                'size': chrom_size
            })

        n_maps, = struct.unpack('<H', data[x:x+2])
        x += 2

        maps = []

        for _ in range(n_maps):
            l, = struct.unpack('<I', data[x:x+4])
            x += 4

            if l:
                name = data[x:x+l].decode('utf-8')
                x += l
            else:
                name = None

            bin_size, = struct.unpack('<I', data[x:x+4])
            x += 4

            m = {
                'name': name,
                'binsize': bin_size,
                'chrom_frags': []
            }

            if not bin_size:
                for i in range(n_chroms):
                    n_frags, = struct.unpack('<I', data[x:x+4])
                    x += 4

                    frags = struct.unpack('<{}I'.format(n_frags), data[x:x+4*n_frags])
                    x += 4 * n_frags

                    m['chrom_frags'].append(frags)

            maps.append(m)

        if self._remote:
            data = self._fetch(offset)
        else:
            self._fh.seek(offset)
            data = self._fh.read()

        x = 0
        for i in range(n_maps):
            chrom1_idx, chrom2_idx, l = struct.unpack('<2HI', data[x:x+8])
            x += 2 + 2 + 4

            if maps[i]['binsize']:
                chrom_size = chroms[chrom1_idx]['size']
                n_rows = (chrom_size + bin_size - 1) // bin_size
            else:
                n_rows = len(maps[i]['chrom_frags'][chrom1_idx])

            s = zlib.decompress(data[x:x+l])
            x += l

            row_offsets = struct.unpack('<{}Q'.format(n_rows), s)

    def add(self, cmap, name=None):
        self._maps.append((cmap, name))

    def set_chrom_sizes(self, chrom_sizes):
        self._chrom_sizes = chrom_sizes

    def write(self, verbose=False, buffersize=0, tmpdir=tempfile.gettempdir()):
        offset = 0
        fh = self._fh

        # Header
        fh.write(struct.pack('<3sH', _SIGNATURE.encode('utf-8'), _VERSION))
        offset += 3 + 2

        _offset = offset

        # Index file position and chrom/res block size
        fh.write(struct.pack('<QI', 0, 0))
        offset += 8 + 4
        blocksize = 0

        chroms = list(self._chrom_sizes.keys())
        fh.write(struct.pack('<H', len(chroms)))
        blocksize += 2

        for chrom_name in chroms:
            l = len(chrom_name)
            fh.write(struct.pack('<HI{}s'.format(l), l, self._chrom_sizes[chrom_name], chrom_name.encode('utf-8')))
            blocksize += 2 + 4 + l

        self._maps.sort(key=lambda tp: tp[0].get_bin_size())
        fh.write(struct.pack('<H', len(self._maps)))
        blocksize += 2
        for cmap, name in self._maps:
            bin_size = cmap.get_bin_size()
            if name:
                l = len(name)
                fh.write(struct.pack('<I{}sI'.format(l), l, name.encode('utf-8'), bin_size))
                blocksize += 4 + l + 4
            else:
                l = 0
                fh.write(struct.pack('<2I'.format(l), l, bin_size))
                blocksize += 4 + 4

            if not bin_size:
                for chrom in chroms:
                    frags = cmap.get_chrom_frags(chrom)
                    n_frags = len(frags)
                    fh.write(struct.pack('<{}I'.format(n_frags+1), n_frags, *frags))
                    blocksize += 4 + n_frags * 4

        # current file position
        offset += blocksize

        cmap_offsets = []

        # Body
        for cmap, name in self._maps:
            bin_size = cmap.get_bin_size()
            cmap.load_contacts(verbose=verbose, buffersize=buffersize, tmpdir=tmpdir)
            sub_cmap_offsets = []

            for chrom1, chrom2, contacts in cmap.iter_contacts():
                if bin_size:
                    n_rows = (self._chrom_sizes[chrom1] + bin_size - 1) // bin_size
                else:
                    n_rows = len(cmap.get_chrom_frags(chrom1))

                row_offsets = []
                it = iter(contacts)
                row, col, val = next(it)
                for i in range(n_rows):
                    row_offsets.append(offset)

                    while row < i:
                        try:
                            row, col, val = next(it)
                        except StopIteration:
                            break

                    cols = []
                    data = []
                    while row == i:
                        cols.append(row)
                        data.append(val)
                        try:
                            row, col, val = next(it)
                        except StopIteration:
                            break

                    n = len(cols)
                    if n:
                        s = zlib.compress(struct.pack('<I{0}I{0}d'.format(n), n, *(cols + data)), _COMPRESS_LEVEL)
                        l = len(s)
                        fh.write(struct.pack('<I{}s'.format(l), l, s))
                        offset += 4 + l
                    else:
                        fh.write(struct.pack('<I', 0))
                        offset += 4

                sub_cmap_offsets.append(
                    (chrom1, chrom2, row_offsets)
                )

            cmap.clean()
            cmap_offsets.append(sub_cmap_offsets)

        # Index
        for sub_cmap_offsets in cmap_offsets:
            for chrom1, chrom2, row_offsets in sub_cmap_offsets:
                chrom1_idx = chroms.index(chrom1)
                chrom2_idx = chroms.index(chrom2)
                s = zlib.compress(struct.pack('<{}Q'.format(len(row_offsets)), *row_offsets), _COMPRESS_LEVEL)
                l = len(s)
                fh.write(struct.pack('<2HI{}s'.format(l), chrom1_idx, chrom2_idx, l, s))

        fh.seek(_offset)
        fh.write(struct.pack('<QI', offset, blocksize))

    def close(self):
        try:
            self._fh.close()
        except AttributeError:
            pass
        finally:
            self._fh = None

    def __del__(self):
        self.close()


def open(filename, mode='r'):
    try:
        fh = File(filename, mode)
    except RuntimeError as e:
        logging.critical(format(e.args[0]))
        fh = None
    finally:
        return fh

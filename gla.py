import builtins
import gzip
import logging
import os
import struct
import sys
import tempfile
import zlib

from multiprocessing import Pool
from urllib.error import HTTPError

try:
    import h5py
except ImportError:
    _SUPPORT_COOLER = False
else:
    _SUPPORT_COOLER = True

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
    format='glasson: %(asctime)s: %(message)s',
    datefmt='%y-%m-%d %H:%M:%S',
    stream=sys.stdout
)


_SIGNATURE = 'GLA'
_VERSION = 1
_COMPRESS_LEVEL = 6


def _open_text(filename, mode='rt'):
    """
    Open a text file and return a file handler.

    Parameters
    ----------
    filename : str
        File path. May be gzip-compressed.
    mode : str
        Mode in which the file is opened.

    Returns
    -------
        File handler

    """
    if filename.lower().endswith('.gz'):
        return gzip.open(filename, mode=mode)
    else:
        return builtins.open(filename, mode=mode)


def load_chrom_sizes(filename):
    """
    Parse a UCSC chrom.sizes file.

    Parameters
    ----------
    filename : str
        File path. May be gzip-compressed.

    Returns
    -------
    List of tuples ``(chrom_name, chrom_size)``.

    """
    chrom_sizes = []

    fh = _open_text(filename, mode='rt')

    for i, line in enumerate(fh):
        try:
            name, size = line.rstrip().split()
            size = int(size)
        except ValueError:
            logging.critical('invalid format at line {} in file {}'.format(i + 1, filename))
            fh.close()
            exit(1)
        else:
            chrom_sizes.append((name, size))

    fh.close()
    return chrom_sizes


class ContactMap:
    def __init__(self, bed_file, mat_file, chrom_sizes):
        self._bed_file = bed_file
        self._mat_file = mat_file
        self._chrom_sizes = {chrom.lower(): size for chrom, size in chrom_sizes}
        self._format = self._detect_format()

        # frag ID -> chrom
        self._frag_ids = {}

        # chrom -> [frag1_start, frag2_start, ...]
        self._frags = {}

        # chrom name -> offset (e.g. chr1 = 0, chr2 = n_frags of chr1, chr3 = n_frags of chr1 + n_frags of chr2, ...)
        self._offsets = {}

        # bin size: > 0 for fixed-size (i.e. bins), 0 otherwise (i.e. RE fragments)
        self._bs = 0

        # chrom1 -> chrom2 -> [(i, j, v), (i, j, v), (i, j, v), ...]
        self._contacts = {}

        # chrom1 -> chrom2 -> path to pickle file
        self._files = {}

    @property
    def format(self):
        return self._format

    @property
    def bin_size(self):
        return self._bs

    def _detect_format(self):
        if self._is_coo(self._mat_file):
            return 'coo'
        elif self._is_dense(self._mat_file):
            return 'dense'
        elif _SUPPORT_COOLER:
            # todo: implement
            pass

        return None

    def load_frags(self):
        """
        Parse the fragment file and load the fragments in memory.

        Returns
        -------
        bool

        """
        status = True
        self._frag_ids = {}
        self._frags = {}
        self._offsets = {}
        self._bs = 0
        cnt = 0

        fh = _open_text(self._bed_file, mode='rt')

        i = 0
        for i, line in enumerate(fh):
            fields = line.rstrip().split()

            try:
                chrom = fields[0].lower()
                start = int(fields[1])
                end = int(fields[2])
                frag_id = int(fields[3])
                chrom_size = self._chrom_sizes[chrom]
            except (IndexError, KeyError, ValueError):
                logging.critical('invalid format at line {} in file {}'.format(i + 1, self._bed_file))
                status = False
                break

            self._frag_ids[frag_id] = chrom

            if chrom not in self._frags:
                self._frags[chrom] = [start]
                self._offsets[chrom] = i
            else:
                self._frags[chrom].append(start)

            if end < chrom_size:
                self._bs += end - start
                cnt += 1

        fh.close()

        if i + 1 != len(self._frag_ids):
            status = False

        if cnt and float(self._bs) / cnt == self._bs // cnt:
            self._bs //= cnt
        else:
            self._bs = 0

        return status

    def load_contacts(self, thread_id=None, **kwargs):
        if self._format == 'coo':
            return self._load_coo(thread_id, **kwargs)
        elif self._format == 'dense':
            pass
        elif _SUPPORT_COOLER:
            pass

        return False

    def _load_coo(self, thread_id=None, **kwargs):
        buffersize = kwargs.get('buffersize', 0)
        tmpdir = kwargs.get('tmpdir', tempfile.gettempdir())
        verbose = kwargs.get('verbose', False)

        self._clean()
        self._contacts = {}
        status = True

        fh = _open_text(self._mat_file, mode='rt')

        for i, line in enumerate(fh):
            if verbose and not (i + 1) % 1000000:
                logging.info('{}{} contacts parsed'.format(
                    'thread #' + str(thread_id) + ': ' if thread_id else '', i + 1)
                )

            fields = line.rstrip().split()

            try:
                frag1_id = int(fields[0])
                frag2_id = int(fields[1])
                value = float(fields[2])
            except ValueError:
                logging.critical('invalid format at line {} in file {}'.format(i + 1, self._mat_file))
                status = False
                break

            if frag1_id > frag2_id or not value:
                # not upper triangular (or null value)
                continue

            try:
                chrom1 = self._frag_ids[frag1_id]
                offset = self._offsets[chrom1]
            except KeyError:
                logging.critical('invalid fragment ID ({}) at line {} in file {}'.format(
                    frag1_id, i + 1, self._mat_file)
                )
                status = False
                break
            else:
                frag1_id -= offset

            try:
                chrom2 = self._frag_ids[frag2_id]
                offset = self._offsets[chrom2]
            except KeyError:
                logging.critical('invalid fragment ID ({}) at line {} in file {}'.format(
                    frag1_id, i + 1, self._mat_file)
                )
                status = False
                break
            else:
                frag2_id -= offset

            if chrom1 not in self._contacts:
                self._contacts[chrom1] = {chrom2: [(frag1_id, frag2_id, value)]}
            elif chrom2 not in self._contacts[chrom1]:
                self._contacts[chrom1][chrom2] = [(frag1_id, frag2_id, value)]
            else:
                self._contacts[chrom1][chrom2].append((frag1_id, frag2_id, value))

            if buffersize and not (i + 1) % buffersize:
                for chrom1 in self._contacts:
                    for chrom2 in self._contacts[chrom1]:
                        data = self._contacts[chrom1][chrom2]

                        if chrom1 in self._files and chrom2 in self._files[chrom1]:
                            with builtins.open(self._files[chrom1][chrom2], 'rb') as pfh:
                                data += pickle.load(pfh)
                        else:
                            fd, path = tempfile.mkstemp(dir=tmpdir)
                            os.close(fd)

                            if chrom1 not in self._files:
                                self._files[chrom1] = {chrom2: path}
                            else:
                                self._files[chrom1][chrom2] = path

                        with builtins.open(self._files[chrom1][chrom2], 'wb') as pfh:
                            pickle.dump(data, pfh)

                        self._contacts[chrom1][chrom2] = []

        fh.close()

        if status and buffersize:
            for chrom1 in self._contacts:
                for chrom2 in self._contacts[chrom1]:
                    data = self._contacts[chrom1][chrom2]

                    if chrom1 in self._files and chrom2 in self._files[chrom1]:
                        with builtins.open(self._files[chrom1][chrom2], 'rb') as pfh:
                            data += pickle.load(pfh)
                    else:
                        fd, path = tempfile.mkstemp(dir=tmpdir)
                        os.close(fd)

                        if chrom1 not in self._files:
                            self._files[chrom1] = {chrom2: path}
                        else:
                            self._files[chrom1][chrom2] = path

                    with builtins.open(self._files[chrom1][chrom2], 'wb') as pfh:
                        pickle.dump(data, pfh)

                    self._contacts[chrom1][chrom2] = []

        return status

    @staticmethod
    def _is_coo(filename, limit=10000):
        b = True

        fh = _open_text(filename, mode='rt')

        for i, line in enumerate(fh):
            if i == limit:
                break

            try:
                frag1_id, c, value = line.rstrip().split()
                int(frag1_id)
                int(frag1_id)
                float(value)
            except ValueError:
                b = False
                break

        fh.close()

        return b

    @staticmethod
    def _is_dense(filename, limit=10000):
        b = True
        n = None

        fh = _open_text(filename, mode='rt')

        for i, line in enumerate(fh):
            if i == limit:
                break

            fields = line.rstrip().split()
            if n is None:
                n = len(fields)
            elif len(fields) != n:
                b = False
                break

            try:
                values = [float(e) for e in fields]
            except ValueError:
                b = False
                break

        fh.close()

        return b

    def _clean(self):
        print(self._mat_file)
        for chrom1 in self._files:
            for chrom2 in self._files[chrom1]:
                os.unlink(self._files[chrom1][chrom2])

        self._files = {}

    # def __del__(self):
    #     self._clean()

class Glasson:
    def __init__(self, path, mode='r', chrom_sizes=list()):
        """
        Parameters
        ----------
        path : str
            Path to a GLA file, or URL
        mode : str
            Mode in which the file is opened: ``'r'`` which means open for reading,
            ``'w'`` for writing (overwrite the file if it already exists).
        chrom_sizes : list
            List of chromosome size tuples ``(name, size)``.
        """

        self._path = path
        self._chrom_sizes = chrom_sizes
        self._maps = []
        self._labels = []
        self._fh = None
        self._remote = False

        if mode == 'r':
            if os.path.isfile(path):
                self._fh = builtins.open(path, 'rb')
            elif path.lower().startswith(('http://', 'https://')):
                self._remote = True
            else:
                raise RuntimeError('{} is not an existing file nor a valid HTTP URL'.format(path))

            self._sniff()
        elif mode == 'w':
            self._fh = builtins.open(path, 'wb')
        else:
            raise RuntimeError("invalid mode: '{}'\n".format(mode))

    @property
    def chrom_sizes(self):
        return self._chrom_sizes

    def add_mat(self, bed_file, mat_file, mat_label):
        cmap = ContactMap(bed_file, mat_file, self._chrom_sizes)

        if not cmap.format:
            raise RuntimeError('cannot detect format for file {}'.format(mat_file))

        self._maps.append(cmap)
        self._labels.append(mat_label)

    def freeze(self, aggregate=0, buffersize=0, processes=1, tmpdir=tempfile.gettempdir(), verbose=False):
        if verbose:
            logging.info('loading fragments')

        with Pool(processes) as pool:
            result = pool.map(self._load_fragments, self._maps)

        if not all(result):
            return

        self._maps = result

        if verbose:
            logging.info('loading contacts')

        with Pool(processes) as pool:
            kwargs = dict(buffersize=buffersize, verbose=verbose, tmpdir=tmpdir)
            items = [(i + 1, cmap, kwargs) for i, cmap in enumerate(self._maps)]

            pool.map(self._load_contacts, items)

    @staticmethod
    def _load_contacts(args):
        thread_id, cmap, kwargs = args
        return cmap if cmap.load_contacts(thread_id, **kwargs) else None


    @staticmethod
    def _load_fragments(cmap):
        return cmap if cmap.load_frags() else None

    def _sniff(self):
        pass

    def close(self):
        try:
            self._fh.close()
        except AttributeError:
            pass
        finally:
            self._fh = None

    def __del__(self):
        self.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()























class _ContactMap:
    def __init__(self, bed_file, mat_file, chrom_sizes):
        """

        Parameters
        ----------
        bed_file : str
            Path to a BED-like file with columns `chrom`, `start`, `end`, `frag_id`.
        mat_file
        chrom_sizes
        """
        self._bed_file = bed_file
        self._mat_file = mat_file
        self._chrom_sizes = {chrom.lower(): size for chrom, size in chrom_sizes}
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
                # logging.critical('invalid format at line {} in file {}'.format(i + 1, self._bed_file))
                fh.close()
                exit(1)
            except KeyError:
                # logging.critical('unknown chromosome \'{}\' at line {} in file {}'.format(chrom, i + 1, self._bed_file))
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
            # logging.critical('expected {} unique fragment IDs, got {}'.format(i + 1, len(frag_ids)))
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
                # logging.critical('invalid format at line {} in file {}'.format(x + 1, line))
                fh.close()
                exit(1)

            if i > j or not v:  # todo: should we raise an error here?
                continue

            try:
                chrom1 = self._frag_ids[i]
                offset = self._chrom_offsets[chrom1]
            except KeyError:
                # logging.critical('invalid fragment ID ({}) at line {} in file {}'.format(i, x + 1, self._mat_file))
                fh.close()
                exit(1)
            else:
                i -= offset

            try:
                chrom2 = self._frag_ids[j]
                offset = self._chrom_offsets[chrom2]
            except KeyError:
                # logging.critical('invalid fragment ID ({}) at line {} in file {}'.format(i, x + 1, self._mat_file))
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
                # logging.info('{} contacts loaded'.format(cnt))
                pass

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

    @property
    def bin_size(self):
        return self._bin_size

    @property
    def chroms(self):
        return [chrom for chrom, offset in sorted(self._chrom_offsets.items(), key=lambda tp: tp[1])]

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
        self._chrom_sizes = []
        self._maps = []  # list of ContactMap in write mode, list of dict in read mode

        if self._mode == 'w':
            self._fh = builtins.open(self._filename, 'wb')
            self._remote = False
        elif os.path.isfile(self._filename):
            self._fh = builtins.open(self._filename, 'rb')
            self._remote = False
            self._sniff()
        elif self._filename.lower().startswith(('http://', 'https://')):
            self._fh = None
            self._remote = True
            self._sniff()
        else:
            raise RuntimeError('{} is not an existing file nor a valid HTTP URL'.format(self._filename))

    def _fetch(self, offset, size=None):
        if size:
            range = 'bytes={}-{}'.format(offset, offset+size-1)
        else:
            range = 'bytes={}-'.format(offset)

        req = request.Request(self._filename, headers={'Range': range})

        try:
            res = request.urlopen(req)
        except HTTPError:
            self.close()
            raise RuntimeError('glasson: file not found')
        else:
            data = res.read()
            return data

    def query(self, x_range, y_range=None, resolution=None):
        if len(x_range) == 3:
            x_chrom1, x_start, x_end = x_range
            x_chrom2 = x_chrom1
        else:
            x_chrom1, x_start, x_chrom2, x_end = x_range

        if not y_range:
            y_chrom1, y_start, y_chrom2, y_end = x_chrom1, x_start, x_chrom2, x_end
        elif len(y_range) == 3:
            y_chrom1, y_start, y_end = y_range
            y_chrom2 = y_chrom1
        else:
            y_chrom1, y_start, y_chrom2, y_end = y_range

        # x_chrom1 = x_chrom1.lower()
        # x_chrom2 = x_chrom2.lower()
        # y_chrom1 = y_chrom1.lower()
        # y_chrom2 = y_chrom2.lower()

        for chrom_name in (x_chrom1, x_chrom2, y_chrom1, y_chrom2):
            try:
                chrom_size = self._chrom_sizes[chrom_name]
            except KeyError:
                # logging.critical('unknown chromosome {}'.format(chrom_name))
                exit(1)

        _map = None
        for m in self._maps:
            # todo: handle RE frags
            bin_size = m['binsize']
            if bin_size and bin_size == resolution:
                _map = m
                break

        if _map is None:
            # todo handle error
            return

        x_i = self._chroms.index(x_chrom1)
        x_j = self._chroms.index(x_chrom2)
        y_i = self._chroms.index(y_chrom1)
        y_j = self._chroms.index(y_chrom2)

        print(x_i)
        print(x_j)
        print(x_chrom1)
        print(x_chrom2)

        # for chrom1, chrom2, row_offsets in _map['submaps']:
        #      if chrom1 == x_chrom1 and chrom2 == x_chrom2:
        #          print(chrom1, chrom2, len(row_offsets))
        #          break

    def _sniff(self):
        self._chrom_sizes = []
        self._maps = []

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

        for _ in range(n_chroms):
            l, chrom_size = struct.unpack('<HI', data[x:x+6])
            x += 6

            chrom_name = data[x:x+l].decode('utf-8')
            x += l

            self._chrom_sizes[chrom_name] = chrom_size

        n_maps, = struct.unpack('<H', data[x:x+2])
        x += 2

        for i in range(n_maps):
            l, = struct.unpack('<I', data[x:x+4])
            x += 4

            if l:
                name = data[x:x+l].decode('utf-8')
                x += l
            else:
                name = None

            bin_size, = struct.unpack('<I', data[x:x+4])
            x += 4

            _map = {
                'name': name,
                'binsize': bin_size,
                'chrom_frags': [],
                'submaps': []
            }

            if not bin_size:
                for j in range(n_chroms):
                    n_frags, = struct.unpack('<I', data[x:x+4])
                    x += 4

                    frags = struct.unpack('<{}I'.format(n_frags), data[x:x+4*n_frags])
                    x += 4 * n_frags

                    _map['chrom_frags'].append(frags)

            self._maps.append(_map)

        if self._remote:
            data = self._fetch(offset)
        else:
            self._fh.seek(offset)
            data = self._fh.read()

        data = zlib.decompress(data)

        x = 0
        for i in range(n_maps):
            n_combinations, = struct.unpack('<I', data[x:x+4])
            x += 4

            for j in range(n_combinations):
                chrom1_idx, chrom2_idx = struct.unpack('<2H', data[x:x+4])
                x += 2 + 2

                chrom1 = self._chroms[chrom1_idx]
                chrom2 = self._chroms[chrom2_idx]

                if self._maps[i]['binsize']:
                    chrom_size = self._chrom_sizes[chrom1]
                    n_rows = (chrom_size + self._maps[i]['binsize'] - 1) // self._maps[i]['binsize']
                else:
                    n_rows = len(self._maps[i]['chrom_frags'][chrom1_idx])

                row_offsets = struct.unpack('<{}Q'.format(n_rows), data[x:x+n_rows*8])
                x += n_rows * 8

                self._maps[i]['submaps'].append((chrom1, chrom2, row_offsets))

    def add(self, cmap, name=None):
        self._maps.append((cmap, name))

    @property
    def chrom_sizes(self):
        return self._chrom_sizes

    @chrom_sizes.setter
    def chrom_sizes(self, value):
        self._chrom_sizes = value

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
            fh.write(struct.pack('<HI{}s'.format(l), l, self._chrom_sizes[chrom_name], chrom_name.lower().encode('utf-8')))
            blocksize += 2 + 4 + l

        self._maps.sort(key=lambda tp: tp[0].bin_size)
        fh.write(struct.pack('<H', len(self._maps)))
        blocksize += 2
        for cmap, name in self._maps:
            bin_size = cmap.bin_size
            if name:
                l = len(name)
                fh.write(struct.pack('<I{}sI'.format(l), l, name.encode('utf-8'), bin_size))
                blocksize += 4 + l + 4
            else:
                l = 0
                fh.write(struct.pack('<2I'.format(l), l, bin_size))
                blocksize += 4 + 4

            if not bin_size:
                for chrom in cmap.chroms:
                    frags = cmap.get_chrom_frags(chrom)
                    n_frags = len(frags)
                    fh.write(struct.pack('<{}I'.format(n_frags+1), n_frags, *frags))
                    blocksize += 4 + n_frags * 4

        # current file position
        offset += blocksize

        cmap_offsets = []

        # Body
        for cmap, name in self._maps:
            bin_size = cmap.bin_size
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
        index = b''
        for i, sub_cmap_offsets in enumerate(cmap_offsets):
            cmap = self._maps[i][0]
            chroms = cmap.chroms

            index += struct.pack('<I', len(sub_cmap_offsets))
            for chrom1, chrom2, row_offsets in sub_cmap_offsets:
                chrom1_idx = chroms.index(chrom1)
                chrom2_idx = chroms.index(chrom2)
                index += struct.pack('<2H{}Q'.format(len(row_offsets)), chrom1_idx, chrom2_idx, *row_offsets)

        fh.write(zlib.compress(index, _COMPRESS_LEVEL))
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
        # logging.critical(format(e.args[0]))
        fh = None
    finally:
        return fh

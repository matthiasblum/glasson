import gzip
import logging
import os
import struct
import sys
import tempfile
import zlib

from multiprocessing import Pool

try:
    import h5py
except ImportError:
    _SUPPORT_COOLER = False
else:
    _SUPPORT_COOLER = True

if sys.version_info[0] == 3:
    import pickle
    import urllib.request as request
    from urllib.error import HTTPError, URLError
else:
    try:
        import cPickle as pickle
    except ImportError:
        import pickle

    import urllib2 as request
    from urllib2 import HTTPError, URLError
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


def _open(filename, mode='rt'):
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
        File handler or None (if the file can not be opened)

    """
    try:
        if filename.lower().endswith('.gz'):
            fh = gzip.open(filename, mode=mode)
        else:
            fh = open(filename, mode=mode)
    except IOError:
        return None
    else:
        return fh


def _read_chromsizes(path_or_db):
    """
    Read a UCSC-style chrom.sizes file.
    
    Parameters
    ----------
    path : str
        File path or database (e.g.: hg38, hg18, mm9, ...)

    Returns
    -------
    List of chromosomes (chrom name, chrom size in bp)

    """
    chrom_sizes = []

    fh = _open(path_or_db, mode='rt')
    if fh:
        content = fh.read().strip().decode()
        fh.close()
    else:
        try:
            res = request.urlopen('http://hgdownload.cse.ucsc.edu/goldenPath/{0}/bigZips/{0}.chrom.sizes'.format(path_or_db))
        except (HTTPError, URLError):
            logging.critical('invalid UCSC database: {}'.format(path_or_db))
            exit(1)
        else:
            content = res.read().strip().decode()

    for i, line in enumerate(content.split('\n')):
        try:
            name, size = line.rstrip().split()
            size = int(size)
        except ValueError:
            print(line)
            logging.critical('invalid format at line {} in file/database {}'.format(i + 1, path_or_db))
            exit(1)
        else:
            chrom_sizes.append((name, size))

    return chrom_sizes


def _fetch(url, offset, size=None):
    """
    Get data from a remote file.
    
    Parameters
    ----------
    url : str
        HTTP[s] URL
    offset : int
        Position at which we start to read.
    size : int, optional
        Number of bytes to read.

    Returns
    -------
        Byte string

    """
    if size:
        range = 'bytes={}-{}'.format(offset, offset+size-1)
    else:
        range = 'bytes={}-'.format(offset)

    req = request.Request(url, headers={'Range': range})

    try:
        res = request.urlopen(req)
    except (HTTPError, URLError):
        raise RuntimeError('glasson: file not found')
    else:
        return res.read()


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

        # chrom -> offset (e.g. chr1 = 0, chr2 = n_frags of chr1, chr3 = n_frags of chr1 + n_frags of chr2, ...)
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

    def get_frags(self, chrom):
        return self._frags.get(chrom.lower(), [])

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

        fh = _open(self._bed_file, mode='rt')

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
        self.empty()
        self._contacts = {}

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

        status = True

        fh = _open(self._mat_file, mode='rt')

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
                            with open(self._files[chrom1][chrom2], 'rb') as pfh:
                                data += pickle.load(pfh)
                        else:
                            fd, path = tempfile.mkstemp(dir=tmpdir)
                            os.close(fd)

                            if chrom1 not in self._files:
                                self._files[chrom1] = {chrom2: path}
                            else:
                                self._files[chrom1][chrom2] = path

                        with open(self._files[chrom1][chrom2], 'wb') as pfh:
                            pickle.dump(data, pfh)

                        self._contacts[chrom1][chrom2] = []

        fh.close()
        return status

    @staticmethod
    def _is_coo(filename, limit=10000):
        b = True

        fh = _open(filename, mode='rt')

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

        fh = _open(filename, mode='rt')

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

    def empty(self):
        """
        Delete temporary pickle files. 
        """
        for chrom1 in self._files:
            for chrom2 in self._files[chrom1]:
                os.unlink(self._files[chrom1][chrom2])

        self._files = {}

    def get_contacts(self):
        """
        Generator that returns all the chromosome-chromosome combinations of the contact map with their contacts.
        
        Returns
        -------
        tuple : (``'chrom1'``, ``'chrom1'``, ``'contacts'``). ``'contacts'`` is an iterator of tuples (``'i'``, ``'j'``, ``'v'``).

        """
        for chrom1 in self._contacts:
            for chrom2 in self._contacts[chrom1]:
                path = self._files.get(chrom1, {}).get(chrom2)
                if path:
                    with open(self._files[chrom1][chrom2], 'rb') as fh:
                        yield chrom1, chrom2, iter(pickle.load(fh) + self._contacts[chrom1][chrom2])
                else:
                    yield chrom1, chrom2, iter(self._contacts[chrom1][chrom2])


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
            The list order defines the order which will be followed when querying inter-chromosomal regions.
        """

        self._path = path
        self._chrom_sizes = chrom_sizes
        self._maps = []
        self._labels = []
        self._fh = None
        self._remote = False

        if mode == 'r':
            if os.path.isfile(path):
                self._fh = open(path, 'rb')
            elif path.lower().startswith(('http://', 'https://')):
                self._remote = True
            else:
                raise RuntimeError('{} is not an existing file nor a valid HTTP URL'.format(path))

            self._sniff()
        elif mode == 'w':
            self._fh = open(path, 'wb')
        else:
            raise RuntimeError("invalid mode: '{}'\n".format(mode))

    @property
    def chrom_sizes(self):
        return self._chrom_sizes

    def add_mat(self, bed_file, mat_file, mat_label):
        m = ContactMap(bed_file, mat_file, self._chrom_sizes)

        if not m.format:
            raise RuntimeError('cannot detect format for file {}'.format(mat_file))

        self._maps.append(m)
        self._labels.append(mat_label)

    def freeze(self, aggregate=0, buffersize=0, processes=1, tmpdir=tempfile.gettempdir(), verbose=False):
        processes = min([processes, len(self._maps)])
        if buffersize:
            buffersize //= processes

        # Load fragments
        with Pool(processes) as pool:
            if verbose:
                logging.info('loading fragments')

            result = pool.map(self._load_fragments, self._maps)

        if not all(result):
            return

        # Reassign result, as ContactMap instances processed by map() are copies
        self._maps = result

        # Load contacts
        # todo: load asynchronously and write header while waiting the first map to be ready
        with Pool(processes) as pool:
            if verbose:
                logging.info('loading contacts')
            kwargs = dict(buffersize=buffersize, verbose=verbose, tmpdir=tmpdir)
            items = [(i + 1, cmap, kwargs) for i, cmap in enumerate(self._maps)]
            result = pool.map(self._load_contacts, items)

        if not all(result):
            return

        # Reassign again
        self._maps = result

        if verbose:
            logging.info('writing output')

        # For convenience only #lazyass
        fh = self._fh

        # Write header
        fh.write(struct.pack('<3sHQI', _SIGNATURE.encode(), _VERSION, 0, 0))
        offset = 17
        blocksize = 0

        # Write chromosomes info
        fh.write(struct.pack('<H', len(self._chrom_sizes)))
        blocksize += 2
        for chrom_name, chrom_size in self._chrom_sizes:
            l = len(chrom_name)
            fh.write(struct.pack('<HI{}s'.format(l), l, chrom_size, chrom_name.encode()))
            blocksize += 6 + l

        # Write maps info
        fh.write(struct.pack('<H', len(self._maps)))
        blocksize += 2

        # Get the indices that would sort the maps by bin_size
        # Thus maps with bin_size = 0 (RE frags) are first
        indices = sorted(range(len(self._maps)), key=lambda i: self._maps[i].bin_size)

        # Sort maps and labels
        self._maps = [self._maps[i] for i in indices]
        self._labels = [self._labels[i] for i in indices]

        for m, label in zip(self._maps, self._labels):
            l = len(label) if label else 0
            bin_size = m.bin_size
            fh.write(struct.pack('<2I', l, bin_size))

            if l:
                fh.write(label.encode())

            blocksize += 8 + l

            if not bin_size:
                # variable-size windows (RE fragments): store start positions
                for chrom_name, chrom_size in self._chrom_sizes:
                    chrom_frags = m.get_frags(chrom_name)
                    n_frags = len(chrom_frags)
                    if n_frags:
                        fh.write(struct.pack('<{}I'.format(n_frags+1), n_frags, *chrom_frags))
                        blocksize += (n_frags + 1) * 4

        # current file position
        offset += blocksize

        # Body
        map_offsets = []
        # Chrom names in lower case (required, as maps store the names in lower case)
        chrom_sizes = {chrom_name.lower(): chrom_size for chrom_name, chrom_size in self._chrom_sizes}
        for m in self._maps:
            bin_size = m.bin_size

            submap_offsets = []

            for chrom1, chrom2, contacts in m.get_contacts():
                if bin_size:
                    n_rows = (chrom_sizes[chrom1] + bin_size - 1) // bin_size
                else:
                    n_rows = len(m.get_frags(chrom1))

                row_offsets = []
                row, col, val = next(contacts)
                for i in range(n_rows):
                    row_offsets.append(offset)

                    while row < i:
                        try:
                            row, col, val = next(contacts)
                        except StopIteration:
                            break

                    cols = []
                    values = []
                    while row == i:
                        cols.append(col)
                        values.append(val)
                        try:
                            row, col, val = next(contacts)
                        except StopIteration:
                            break

                    n = len(cols)
                    if n:
                        s = zlib.compress(struct.pack('<I{0}I{0}d'.format(n), n, *(cols + values)), _COMPRESS_LEVEL)
                        l = len(s)
                        fh.write(struct.pack('<I{}s'.format(l), l, s))
                        offset += 4 + l
                    else:
                        fh.write(struct.pack('<I', 0))
                        offset += 4

                submap_offsets.append((chrom1, chrom2, row_offsets))

            m.empty()
            map_offsets.append(submap_offsets)

        # Index/footer
        index = b''
        chrom_sizes = [chrom_name.lower() for chrom_name, chrom_size in self._chrom_sizes]
        for m, submap_offsets in zip(self._maps, map_offsets):
            index += struct.pack('<I', len(submap_offsets))
            for chrom1, chrom2, row_offsets in submap_offsets:
                chrom1_idx = chrom_sizes.index(chrom1)
                chrom2_idx = chrom_sizes.index(chrom2)
                index += struct.pack('<2H{}Q'.format(len(row_offsets)), chrom1_idx, chrom2_idx, *row_offsets)

        fh.write(zlib.compress(index, _COMPRESS_LEVEL))
        fh.seek(5)
        fh.write(struct.pack('<QI', offset, blocksize))

    @staticmethod
    def _load_contacts(args):
        thread_id, cmap, kwargs = args
        return cmap if cmap.load_contacts(thread_id, **kwargs) else None

    @staticmethod
    def _load_fragments(cmap):
        return cmap if cmap.load_frags() else None

    def _sniff(self):
        """
        Verify that the file passed to the constructor is a valid GLA file, and parse its header/index if it is.
        """
        self._chrom_sizes = []
        self._maps = []
        self._labels = []

        if self._remote:
            data = _fetch(self._path, 0, 17)
        else:
            data = self._fh.read(17)

        signature, version, offset, blocksize = struct.unpack('<3sHQI', data)

        try:
            is_gla = signature.decode() == 'GLA'
        except UnicodeDecodeError:
            is_gla = False

        if not is_gla:
            raise RuntimeError('glasson: not a valid GLA file')

        if self._remote:
            data = _fetch(self._path, 17, blocksize)
        else:
            data = self._fh.read(blocksize)

        n_chroms, = struct.unpack('<H', data[:2])
        x = 2

        for _ in range(n_chroms):
            l, chrom_size = struct.unpack('<HI', data[x:x+6])
            x += 6

            chrom_name = data[x:x+l].decode()
            x += l

            self._chrom_sizes.append((chrom_name, chrom_size))

        n_maps, = struct.unpack('<H', data[x:x + 2])
        x += 2

        for i in range(n_maps):
            l, bin_size = struct.unpack('<2I', data[x:x+8])
            x += 8

            if l:
                label = data[x:x+l].decode()
                x += l
            else:
                label = None

            m = {
                'name': label,
                'binsize': bin_size,
                'chrom_frags': [],
                'submaps': {}
            }

            if not bin_size:
                for j in range(n_chroms):
                    n_frags, = struct.unpack('<I', data[x:x+4])
                    x += 4

                    frags = struct.unpack('<{}I'.format(n_frags), data[x:x+4*n_frags])
                    x += 4 * n_frags

                    m['chrom_frags'].append(frags)

            self._maps.append(m)

        if self._remote:
            data = _fetch(offset)
        else:
            self._fh.seek(offset)
            data = self._fh.read()

        data = zlib.decompress(data)
        x = 0
        for i in range(n_maps):
            n_combinations, = struct.unpack('<I', data[x:x + 4])
            x += 4

            for j in range(n_combinations):
                chrom1_idx, chrom2_idx = struct.unpack('<2H', data[x:x + 4])
                x += 4

                chrom1_name, chrom1_size = self._chrom_sizes[chrom1_idx]
                chrom2_name, chrom2_size = self._chrom_sizes[chrom2_idx]

                if self._maps[i]['binsize']:
                    n_rows = (chrom1_size + self._maps[i]['binsize'] - 1) // self._maps[i]['binsize']
                else:
                    n_rows = len(self._maps[i]['chrom_frags'][chrom1_idx])

                row_offsets = struct.unpack('<{}Q'.format(n_rows), data[x:x + n_rows * 8])
                x += n_rows * 8

                chrom1_name = chrom1_name.lower()
                chrom2_name = chrom2_name.lower()

                # todo: really ugly with dictionaries. Replace by object.
                if chrom1_name not in self._maps[i]['submaps']:
                    self._maps[i]['submaps'][chrom1_name] = {chrom2_name: row_offsets}
                else:
                    self._maps[i]['submaps'][chrom1_name][chrom2_name] = row_offsets

    @staticmethod
    def _get_chrom_range(chrom1, chrom2, chroms):
        """
        Return the list of chromosome indices between two chromosomes.
        
        Parameters
        ----------
        chrom1 : str
            chromosome name
        chrom2 : str
            chromosome name
        chrom_names : list
            list of chromosome names

        Returns
        -------
        tuple of ints

        """
        i = 0
        j = 0

        try:
            i = chroms.index(chrom1.lower())
        except ValueError:
            logging.critical('unknown chromosome {}'.format(chrom1))
            exit(1)

        try:
            j = chroms.index(chrom2.lower())
        except ValueError:
            logging.critical('unknown chromosome {}'.format(chrom2))
            exit(1)

        if i > j:
            i, j = j, i

        return i, j

    def query(self, range1, range2=None, resolution=None):
        if len(range1) == 3:
            x_chrom1, x_start, x_end = range1
            x_chrom2 = x_chrom1
        else:
            x_chrom1, x_start, x_chrom2, x_end = range1

        if not range2:
            y_chrom1, y_start, y_chrom2, y_end = x_chrom1, x_start, x_chrom2, x_end
        elif len(range2) == 3:
            y_chrom1, y_start, y_end = range2
            y_chrom2 = y_chrom1
        else:
            y_chrom1, y_start, y_chrom2, y_end = range2

        chrom_names = [chrom_name.lower() for chrom_name, chrom_size in self._chrom_sizes]
        x_i, x_j = self._get_chrom_range(x_chrom1, x_chrom2, chrom_names)
        y_i, y_j = self._get_chrom_range(y_chrom1, y_chrom2, chrom_names)

        m = None
        for _m in self._maps:
            # todo: handle RE frags
            bin_size = _m['binsize']
            if bin_size and bin_size == resolution:
                m = _m
                break

        if m is None:
            # todo handle error
            return

        bin_size = m['binsize']

        result = []

        for j in range(y_i, y_j + 1):
            chrom2, chrom2_size = self._chrom_sizes[j]

            if j == y_i:
                _ys = y_start // bin_size
                _ye = (y_end + bin_size - 1) // bin_size if j == y_j else 0
            elif j < y_j:
                # from the firs row
                _ys = 0
                _ye = 0
            else:
                # not all the rows
                _ys = 0
                _ye = (y_end + bin_size - 1) // bin_size

            for i in range(x_i, x_j + 1):
                if i > j:
                    continue

                chrom1, chrom1_size = self._chrom_sizes[i]
                row_offsets = m['submaps'].get(chrom1.lower(), {}).get(chrom2.lower(), [])

                if i == x_i:
                    _xs = x_start // bin_size
                    _xe = (x_end + bin_size - 1) // bin_size if i == x_j else 0
                elif i < x_j:
                    # all the rows
                    _xs = 0
                    _xe = 0
                else:
                    # trim end of the rows
                    _xs = 0
                    _xe = (x_end + bin_size - 1) // bin_size

                _rows = []
                rows = []
                cols = []
                values = []

                _os = None
                _oe = None

                for k, o in enumerate(row_offsets):
                    if _ye and k >= _ye:
                        _oe = o
                        break
                    elif k >= _ys:
                        if _os is None:
                            _os = o
                        _oe = o
                        _rows.append(k)

                if self._remote:
                    data = _fetch(self._path, _os, _oe - _os)
                else:
                    self._fh.seek(_os)
                    data = self._fh.read(_oe - _os)

                l = len(data)
                o = 0
                k = 0
                while o < l:
                    ls, = struct.unpack('<I', data[o:o+4])
                    o += 4

                    if ls:
                        s = zlib.decompress(data[o:o+ls])
                        o += ls

                        n, = struct.unpack('<I', s[:4])
                        _values = list(struct.unpack('<{0}I{0}d'.format(n), s[4:]))

                        for x in range(n):
                            if _xe and _values[x] >= _xe:
                                break
                            elif _values[x] >= _xs:
                                rows.append(_rows[k])
                                cols.append(_values[x])
                                values.append(_values[x+n])

                    k += 1

                result.append((chrom1, chrom2, rows, cols, values))

        return result

    def close(self):
        for m in self._maps:
            try:
                m.empty()
            except AttributeError:
                pass

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


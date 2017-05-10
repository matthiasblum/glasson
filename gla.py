import builtins
import os
import sys


class File:
    def __init__(self, filename, mode='r'):
        if mode not in ('r', 'w'):
            raise RuntimeError('invalid mode: \'{}\'\n'.format(mode))

        self._filename = filename
        self._mode = mode
        self._remote = False
        self._fh = None
        self._index = {}
        self._offset = None
        self._chrom_sizes = {}

        if self._mode == 'w':
            self._fh = builtins.open(self._filename, 'wb')
        elif os.path.isfile(self._filename):
            self._fh = builtins.open(self._filename, 'rb')
            self._remote = False
            # self._parseindex()
        elif self._filename.lower().startswith(('http://', 'https://')):
            self.remote = True
            # self._parseindex()
        else:
            raise RuntimeError('{} is not an existing file nor a valid HTTP URL'.format(self._filename))

    def add_resolution(self, matrices, bin_size=0, fragments=list(), name=None):
        pass

    def set_chrom_sizes(self, chrom_sizes):
        self._chrom_sizes = chrom_sizes

    def write(self):
        pass


def open(filename, mode='r'):
    try:
        fh = File(filename, mode)
    except RuntimeError as e:
        sys.stderr.write('{}\n'.format(e.args[0]))
        fh = None
    finally:
        return fh

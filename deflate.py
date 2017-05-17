import struct
import zlib

# with open('test.dat', 'wb') as fh:
#     s = zlib.compress('16 May 2017'.encode(), 6)
#     l = len(s)
#     fh.write(struct.pack('<I', l))
#     fh.write(s)
#
# with open('test.dat', 'rb') as fh:
#     l, = struct.unpack('<I', fh.read(4))
#     s = zlib.decompress(fh.read(l))
#     print(l)
#     print(s)

with open('test.dat', 'rb') as fh:
    l, = struct.unpack('<I', fh.read(4))
    print(zlib.decompress(fh.read(l)))

    l, = struct.unpack('<I', fh.read(4))
    s = zlib.decompress(fh.read(l))
    print(struct.unpack('<10B', s))

    print(struct.unpack('<10I', fh.read(10 * 4)))

    l, = struct.unpack('<I', fh.read(4))
    s = zlib.decompress(fh.read(l))
    print(struct.unpack('<10I', s))




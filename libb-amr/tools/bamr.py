"""
A .bamr file (libb-amr's tree_output (), tree-output.h) read in Python, with
the standard library alone (no numpy):

    import bamr
    g = bamr.read ('step-120.bamr')
    g['dimension'], g['m'], g['origin'], g['size'], g['names']
    g['leaves']           # [(level, i, j, k), ...] along the Morton curve
    g['data']['s0']       # array of every cell's value, leaf by leaf, a
                          # leaf's cells i fastest, then j, then k
    for x, h, v in bamr.cells (g, 's0'):   # cell centres, sizes, values
        ...

As a script:  python3 bamr.py FILE [--verify]
prints the header and each scalar's range; --verify also checks the file's
size, that the leaves run strictly along the curve and cover the domain
(their volumes sum to size^D), and exits 1 if not.
"""

import array
import struct
import sys

HEAD = 80


def read (path):
    with open (path, 'rb') as fp:
        b = fp.read ()
    if b[:8] != b'BAMRGRID':
        raise ValueError ('%s: not a .bamr file' % path)
    version, order, dim, m, rs, n = struct.unpack_from ('<6I', b, 8)
    if order != 0x01020304:
        raise ValueError ('%s: byte order %#x' % (path, order))
    if version != 1:
        raise ValueError ('%s: version %d' % (path, version))
    (N,) = struct.unpack_from ('<Q', b, 32)
    origin = struct.unpack_from ('<3d', b, 40)
    (size,) = struct.unpack_from ('<d', b, 64)
    codec, _ = struct.unpack_from ('<2I', b, 72)
    if codec != 0:
        raise ValueError ('%s: codec %d' % (path, codec))
    at, names = HEAD, []
    for _ in range (n):
        (l,) = struct.unpack_from ('<I', b, at)
        names.append (b[at + 4:at + 4 + l].decode ('utf-8'))
        at += 4 + l
    H = (at + 7) // 8 * 8
    leaves = [struct.unpack_from ('<4i', b, H + 16 * q) for q in range (N)]
    mm = m ** dim
    data, start = {}, H + 16 * N
    for k, name in enumerate (names):
        a = array.array ('d' if rs == 8 else 'f')
        a.frombytes (b[start + k * N * mm * rs:start + (k + 1) * N * mm * rs])
        if sys.byteorder != 'little':
            a.byteswap ()
        data[name] = a
    return {'dimension': dim, 'm': m, 'real': rs, 'origin': origin[:dim],
            'size': size, 'names': names, 'leaves': leaves, 'data': data,
            'bytes': len (b), 'expected_bytes': H + 16 * N + n * N * mm * rs}


def cells (g, name):
    """Every cell of scalar name: (centre, size, value), in the file's
    order."""
    dim, m, o, size = g['dimension'], g['m'], g['origin'], g['size']
    v, c = g['data'][name], 0
    for level, *ijk in g['leaves']:
        h = size / 2 ** level
        dx = h / m
        for e in range (m ** dim):
            x, r = [], e
            for d in range (dim):
                x.append (o[d] + ijk[d] * h + (r % m + 0.5) * dx)
                r //= m
            yield tuple (x), dx, v[c]
            c += 1


def key (dim, leaf, depth = 21):
    """A leaf's place on the curve: its first finest cell's Morton key."""
    level, *ijk = leaf
    q = [ijk[d] << (depth - level) for d in range (dim)]
    k = 0
    for bit in range (depth):
        for d in range (dim):
            k |= (q[d] >> bit & 1) << (dim * bit + d)
    return k


def verify (g):
    """The problems found, [] if none."""
    bad = []
    if g['bytes'] != g['expected_bytes']:
        bad.append ('%d bytes, %d expected' % (g['bytes'], g['expected_bytes']))
    keys = [key (g['dimension'], L) for L in g['leaves']]
    if any (b <= a for a, b in zip (keys, keys[1:])):
        bad.append ('the leaves not along the curve')
    vol = sum ((g['size'] / 2 ** L[0]) ** g['dimension'] for L in g['leaves'])
    if abs (vol - g['size'] ** g['dimension']) > 1e-12 * g['size'] ** g['dimension']:
        bad.append ('the leaves cover %r of %r' % (vol, g['size'] ** g['dimension']))
    return bad


def main (argv):
    if len (argv) < 2:
        print (__doc__)
        return 2
    g = read (argv[1])
    levels = sorted ({L[0] for L in g['leaves']})
    print ('%s: %dD, %d leaves of %d^%d cells (levels %s), %s, origin %s, '
           'size %r' % (argv[1], g['dimension'], len (g['leaves']), g['m'],
                        g['dimension'], levels,
                        'double' if g['real'] == 8 else 'float',
                        list (g['origin']), g['size']))
    for name in g['names']:
        v = [x for x in g['data'][name] if x == x]
        print ('  %-12s min %-12.6g max %-12.6g' % (name, min (v, default = 0),
                                                    max (v, default = 0)))
    if '--verify' in argv:
        bad = verify (g)
        for b in bad:
            print ('  FAIL', b)
        if bad:
            return 1
        print ('  verified')
    return 0


if __name__ == '__main__':
    sys.exit (main (sys.argv))

# Copyright (C) 2021 OliverBScott
"""
pyalignit
---------
Python wrappers for the Align-it™ tool from Silicos-it.

"""
import warnings

from .cpyalignit import *
from .draw import PHARM_COLORS, draw_pharmacophore

__version__ = GetVersion()


class PharmacophoreSupplier(object):

    def __init__(self, filename):
        self._handle = open(filename, 'r')
        self._done = False

    @property
    def at_end(self):
        return self._done

    def _get_next_block(self):
        block, curr = '', None
        while curr != '$$$$\n':
            curr = self._handle.readline()
            if curr == '':
                return None
            if curr.startswith('#') or curr.startswith('NAME'):
                continue
            if curr != '$$$$\n':
                block += curr
        return block[:-1]

    @staticmethod
    def _process_block(block):
        pharm = Pharmacophore()
        for line in block.split('\n'):
            tokens = line.strip().split('\t')
            if len(tokens) != 9:
                warnings.warn(f'Invalid line format: {line}')
                return None
            p = PharmacophoreSupplier._tokens_to_point(tokens)
            pharm.append(p)
        return pharm

    @staticmethod
    def _tokens_to_point(tokens):
        point = PharmacophorePoint()
        point.func = FuncGroup.names.get(
            tokens[0], FuncGroup.UNDEF)
        point.point.x = float(tokens[1])
        point.point.y = float(tokens[2])
        point.point.z = float(tokens[3])
        point.alpha = float(tokens[4])
        point.hasNormal = bool(int(tokens[5]))
        point.normal.x = float(tokens[6])
        point.normal.y = float(tokens[7])
        point.normal.z = float(tokens[8])
        if point.func == FuncGroup.UNDEF:
            warnings.warn(
                f'Unknown functional group {tokens[0]},'
                ' reverting to UNDEF'
            )
        return point

    def __iter__(self):
        return self

    def __next__(self):
        if self._handle.closed:
            raise IOError('File closed unexpectedly')
        block = self._get_next_block()
        if block is None:
            self._done = True
            self._handle.close()
            raise StopIteration
        return self._process_block(block)


class PharmacophoreWriter(object):

    def __init__(self, filename):
        self._handle = open(filename, 'w')
        self._num = 0

    @property
    def num_pharmacophores(self):
        return self._num

    @property
    def closed(self):
        return self._handle.closed

    def get_text(self, pharmacophore, name=None):
        name = name if name else ''
        text = f'NAME\t{name}\n'
        for point in pharmacophore:
            text += self._point_to_text(point)
        text += '$$$$\n'
        return text

    @staticmethod
    def _point_to_text(p):
        text = p.func.name + '\t'
        text += f"{p.point.x:.6f}\t{p.point.y:.6f}\t{p.point.z:.6f}\t"
        text += f"{p.alpha}\t{int(p.hasNormal)}\t"
        text += f"{p.normal.x:.6f}\t{p.normal.y:.6f}\t{p.normal.z:.6f}\n"
        return text

    def write(self, pharmacophore, name=None):
        if self._handle.closed:
            raise IOError('Cannot write to a closed file')
        self._handle.write(self.get_text(pharmacophore, name))
        self._handle.flush()
        self._num += 1

    def flush(self):
        self._handle.flush()

    def close(self):
        self._handle.close()

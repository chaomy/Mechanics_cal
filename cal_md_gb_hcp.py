#!/usr/bin/env python
# encoding: utf-8
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-07-05 08:12:30
# @Last Modified by:   chaomy
# @Last Modified time: 2017-11-01 14:24:42

from ase import Atoms
import ase.lattice.hexagonal as Hexagonal
import ase.io
from ase.lattice.orthorhombic import SimpleOrthorhombicFactory
from math import sqrt

class othoHCP(SimpleOrthorhombicFactory):
    bravais_basis = [[0.0, 0.0, 0.0],
                     [0.5, 0.5, 0.0],
                     [0.0, 1. / 3., 1. / 2.],
                     [1 / 2., 5. / 6., 1. / 2.]]


class md_gb_hcp(object):

    def build_hcp_gb(self):
        # atoms = Hexagonal.HexagonalClosedPacked(
        #     latticeconstant={'a': self.pot['ahcp'],
        #                      'c': self.pot['chcp']},
        #     size=(10, 40, 4),
        #     symbol=self.pot['element'],
        #     pbc=(1, 1, 1))
        atoms = othoHCP(latticeconstant={'a': self.pot['ahcp'], 
        								 'b': self.pot['ahcp'] * sqrt(3.),
        								 'c': self.pot['chcp']})
        print atoms.get_cell()
        ase.io.write("poscar", images=atoms, format='vasp')
        self.write_lmp_config_data(atoms, 'lmp_init.txt') 
        return

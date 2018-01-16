#!/usr/bin/env python
# encoding: utf-8
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-07-05 08:12:30
# @Last Modified by:   chaomy
# @Last Modified time: 2017-12-03 21:15:53

from ase import Atoms
import ase.lattice.orthorhombic as otho
import ase.io
from ase.lattice.orthorhombic import SimpleOrthorhombicFactory
from math import sqrt

class othoHCPFractory(otho.SimpleOrthorhombicFactory):
    bravais_basis = [[0.0, 0.0, 0.0],
                     [0.5, 0.5, 0.0],
                     [0.0, 1. / 3., 1. / 2.],
                     [1 / 2., 5. / 6., 1. / 2.]]
othoHCP = othoHCPFractory()

class md_gb_hcp(object):

    def build_hcp_gb(self):
        # atoms = Hexagonal.HexagonalClosedPacked(
        #     latticeconstant={'a': self.pot['ahcp'],
        #                      'c': self.pot['chcp']},
        #     size=(10, 40, 4),
        #     symbol=self.pot['element'],
        #     pbc=(1, 1, 1))
        atoms = othoHCP(latticeconstant=(self.pot['ahcp'], 
        								 self.pot['ahcp'] * sqrt(3.), 
                                         self.pot['chcp']),
                                         size = (4, 4, 2),
                                         symbol=self.pot['element']) 

        print atoms.get_cell()

        # pmat = self.hcp_til_mtx(30)

        ase.io.write("perf_hcp", images=atoms, format='vasp')

        ang = 30
        atoms.rotate(30, 'z', center='COU', rotate_cell=True)

        # rotate 
        # cell = atoms.get_cell()
        # cell = (pmat * cell.transpose()).transpose()
        # atoms.set_cell(cell)

        # pos = atoms.get_positions()
        # pos = (pmat * pos.transpose()).transpose()
        # atoms.set_positions(pos) 

        # atoms.wrap()
        ase.io.write("rot30", images=atoms, format='vasp')
        print atoms.get_cell() 
        # self.write_lmp_config_data(atoms, 'lmp_init.txt') 
        return

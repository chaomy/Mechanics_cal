# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-12-03 11:07:29
# @Last Modified by:   chaomy
# @Last Modified time: 2018-06-08 00:16:22

import ase.lattice.orthorhombic as otho
import ase.io
import os
import numpy as np
import atomman as am
from utils import stroh_solve
from ase import Atoms
from ase.lattice.orthorhombic import SimpleOrthorhombicFactory
from numpy import sqrt


class othoHCPFractory(otho.SimpleOrthorhombicFactory):
    bravais_basis = [[0.0, 0.0, 0.0],
                     [0.5, 0.0, 0.5],
                     [1. / 3., 0.5, 0.0],
                     [5. / 6., 0.5, 0.5]]
othoHCP = othoHCPFractory()


class othoHCPFractoryB(otho.SimpleOrthorhombicFactory):
    bravais_basis = [[0.0, 0.0, 0.0],
                     [0.5, 0.0, 0.5],
                     [2. / 3., 0.5, 0.0],
                     [1. / 6., 0.5, 0.5]]
othoHCPB = othoHCPFractoryB()


class md_gb_ase_1210(object):

    def build_hcp_ase_1210_small(self):  # to examine the GB structures
        # self.find_angles_1210(il=[[1], [1]], jl=[1])    # 43.138
        self.find_angles_1210(il=[[], [1]], jl=[2])    # 61.91
        print(self.ag[0])

        ux = self.pot['ahcp'] * sqrt(3.)
        uy = self.pot['chcp']
        uz = self.pot['ahcp']

        #angle, length, i, j
        print(self.ag)
        atoms = othoHCP(latticeconstant=(ux, uy, uz), size=(
            80, 80, 2), symbol=self.pot['element'])
        atoms.rotate(self.ag[0][0], 'z')

        cell = atoms.get_cell()
        cell[0, 0], cell[1, 1] = 2 * self.ag[0][1], 160
        atoms.translate(np.array([cell[0, 0], -3 * uy, 0]))

        lob = np.array([0.0, 10, 0.0])
        # for 43.13
        # hib = np.array([cell[0, 0], 0.5 * cell[1, 1], cell[2, 2]])
        hib = np.array([cell[0, 0], 0.5 * cell[1, 1] - 1.0, cell[2, 2]])
        atoms = self.make_cubic('out', atoms, lob, hib)

        # the other grain
        atoms2 = othoHCP(latticeconstant=(ux, uy, uz), size=(
            80, 80, 2), symbol=self.pot['element'])

        # for 43.13
        lob = np.array([0.0, 0.5 * cell[1, 1], 0.0])
        hib = np.array([cell[0, 0], cell[1, 1] - 10, cell[2, 2]])

        atoms2.rotate(-self.ag[0][0], 'z')
        atoms2.translate(np.array([-33 * ux, cell[1, 1], 0]))  # for 43.1383
        atoms2 = self.make_cubic('out', atoms2, lob, hib)

        # atoms2.translate(np.array([0.0, 0.0, 0.5 * uz]))   # for 43.1383
        atoms2.translate(np.array([0.0, -1.0, 0.5 * uz]))   # for 43.1383
        atoms.extend(atoms2)

        # assign gb regio
        lob = np.array([0.0, 40, -1.])
        hib = np.array([cell[0, 0], cell[1, 1] - 40, cell[2, 2] + 1.0])
        atoms = self.assign_cubic(atoms, 'out', 'Mo', lob, hib)
        lob = np.array([0.0, 20, -1.])
        hib = np.array([cell[0, 0], cell[1, 1] - 20, cell[2, 2] + 1.0])
        atoms = self.assign_cubic(atoms, 'out', 'W', lob, hib)

        atoms.set_cell(cell)
        self.write_lmp_config_data(atoms, "lmp_init.txt")

    def intro_edge(self):
        self.find_angles_1100(il=[[1], [1]], jl=[1])    # 58.361
        ux = self.pot['ahcp']
        uy = self.pot['chcp']
        uz = self.pot['ahcp'] * sqrt(3.)
        atoms = ase.io.read("dump/dump.00017", format='lammps-dump')

        atoms.rotate(self.ag[0][0], 'z')  # rotate
        atoms.translate(np.array([80 * ux, 0.0, 0]))
        # # introduce dislocation
        axes = np.array([[1, 0, 0], [0, 0, 1], [0, -1, 0]])
        burgers = self.pot['lattice'] * np.array([1, 0, 0])
        c = am.ElasticConstants()
        c.hexagonal(C11=self.pot['C11'], C12=self.pot['C12'],
                    C33=self.pot['C33'], C13=self.pot['C13'], C44=self.pot['C44'])
        stroh = stroh_solve.Stroh(c, burgers, axes=axes)
        pos = atoms.get_positions()
        cx = 60 * ux + 0.01
        cy = 60 * uy + 0.01
        shift = np.ones(pos.shape) * np.array([cx, cy, 0.0])
        disp = stroh.displacement(pos - shift)
        atoms.set_positions(pos + np.real(disp))

        atoms.translate(np.array([-80 * ux, 0.0, 0]))
        atoms.rotate(-self.ag[0][0], 'z')  # rotate back

        # try with non periodict boundary conditions
        cell = atoms.get_cell()
        cell[0, 0] += 30
        cell[1, 1] += 30
        atoms.translate(np.array([15, 15, 0.0]))
        atoms.set_cell(cell)
        # ase.io.write("pos02", images=atoms, format='cfg')
        self.write_lmp_config_data(atoms, "pos02")

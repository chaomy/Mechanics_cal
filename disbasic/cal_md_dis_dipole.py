#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-06-25 14:28:58
# @Last Modified by:   chaomy
# @Last Modified time: 2018-05-01 22:50:11

from numpy import cos, sin, sqrt, mat
from utils import stroh_solve
from crack import cal_md_crack_ini
import numpy as np
import ase
import ase.io
import tool_elastic_constants
import ase.lattice
import atomman as am

strain = mat([[1.0, 0.0, 0.0],
              [0.5, 1.0, 0.5],
              [0.0, 0.0, 1.0]])

axes = np.array([[1, 1, -2],
                 [-1, 1, 0],
                 [1, 1, 1]])


class cal_dis_dipole(object):

    def set_dipole_box(self, sizen=1):  # 7 x 11 x 1
        n = 7 * sizen
        m = 11 * sizen
        t = 1 * sizen
        atoms = ase.lattice.cubic.BodyCenteredCubic(
            directions=[[1., 1., -2.], [-1., 1., 0], [0.5, 0.5, 0.5]],
            latticeconstant=self.pot['lattice'],
            size=(n, m, t),
            symbol=self.pot['element'])
        atoms = self.cut_half_atoms_new(atoms, "cuty")
        supercell = atoms.get_cell()
        supercell = strain * supercell
        atoms.set_cell(supercell)
        atoms.wrap(pbc=[1, 1, 1])
        return atoms

    def bcc_screw_dipole_configs_alongz(self, sizen=1):
        c = tool_elastic_constants.elastic_constants(
            C11=self.pot['c11'], C12=self.pot['c12'], C44=self.pot['c44'])
        burgers = self.pot['lattice'] / 2 * np.array([1., 1., 1.])
        stroh = stroh_solve.Stroh(c, burgers, axes=axes)

        atoms = self.set_dipole_box()
        atoms_perf = atoms.copy()
        pos = atoms.get_positions()

        unitx = sqrt(6) / 3. * self.pot['lattice']
        unity = sqrt(2) / 2. * self.pot['lattice']
        unitz = sqrt(3) / 2. * self.pot['lattice']

        sx = 10.0 * sizen
        sy = 5 * sizen
        ix = 10.5 * sizen

        # c1 = 1. / 3. * np.sum(self.pot['core1'], axis=0)
        # c2 = 1. / 3. * np.sum(self.pot['core2'], axis=0)
        # shiftc1 = \
        # np.ones(np.shape(pos)) * np.array([c1[0, 0], c1[0, 1], 0.0])
        # shiftc2 = \
        # np.ones(np.shape(pos)) * np.array([c2[0, 0], c2[0, 1], 0.0])

        opt = 'original'
        if opt in ['split']:
            c1 = self.pot['posleft'] + \
                np.array([0.0, 0.21 * self.pot['yunit']])
            c2 = self.pot['posrigh'] + \
                np.array([0.0, -0.21 * self.pot['yunit']])
        elif opt in ['move']:
            c1 = [(sx + 0.5) * unitx, (sy + 1. / 3. - 0.95 * 1. / 3.) * unity]
            c2 = [(sx + ix + 0.5) * unitx,
                  (sy + 2. / 3. + 0.95 * 1. / 3.) * unity]
        else:
            c1 = [(sx) * unitx, (sy + 1. / 3.) * unity]
            c2 = [(sx + ix) * unitx, (sy + 2. / 3.) * unity]

        shiftc1 = np.ones(np.shape(pos)) * np.array([c1[0], c1[1], 0.0])
        shiftc2 = np.ones(np.shape(pos)) * np.array([c2[0], c2[1], 0.0])

        disp1 = stroh.displacement(pos - shiftc1)
        disp2 = stroh.displacement(pos - shiftc2)

        if opt in ['pull']:
            radius = 2.0  # find the atoms near the center
            for ps, dp in zip(pos, disp1):
                dis = np.linalg.norm(ps[:2] - c1)
                if (dis < radius):
                    print(dis)
                    dp[2] += 1. / 6. * unitz
                    # add shirt
            for ps, dp in zip(pos, disp2):
                dis = np.linalg.norm(ps[:2] - c2)
                if (dis < radius):
                    print(dis)
                    dp[2] -= 1. / 6. * unitz

        atoms.set_positions(pos + np.real(disp1) - np.real(disp2))

        # periodic boundary conditions ???
        # c2l = [(sx - ix + 0.5) * unitx, (sy + 2. / 3. + 0.95 * 1. / 3.) * unity]
        # shft = np.ones(np.shape(pos)) * np.array([c2l[0], c2l[1], 0.0])
        # disp3 = stroh.displacement(pos - shft)
        # atoms.set_positions(atoms.get_positions() - np.real(disp3))

        atoms.wrap(pbc=[1, 1, 1])
        ase.io.write("perf_poscar", atoms_perf, format='vasp')
        ase.io.write('POSCAR', atoms, format='vasp')
        return (atoms, atoms_perf)

    def bcc_screw_dipole_configs_alongz(self, sizen=1):
        atoms, atoms_perf = self.bcc_screw_dipole_configs_alongz()

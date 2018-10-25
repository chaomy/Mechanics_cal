#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-06-25 14:28:58
# @Last Modified by:   chaomy
# @Last Modified time: 2018-10-20 15:39:16

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

        (layerList, distance) = ase.utils.geometry.get_layers(
            atoms, [0, 0, 1], tolerance=0.05)  # the normal of z direction
        layer1_index, layer2_index, layer3_index = [], [], []

        for i in range(len(atoms)):
            if layerList[i] in [1, 4]:
                layer1_index.append(i)
                atoms[i].symbol = 'W'

            if layerList[i] in [2]:
                layer2_index.append(i)
                atoms[i].symbol = 'Mo'

            if layerList[i] in [3, 0]:
                layer3_index.append(i)
                atoms[i].symbol = 'Ta'

        atoms = self.cut_half_atoms_new(atoms, "cuty")
        supercell = atoms.get_cell()
        supercell = strain * supercell
        atoms.set_cell(supercell)
        atoms.wrap(pbc=[1, 1, 1])
        ase.io.write("perf_poscar", images=self.convert_alongz_to_alongy(
            atoms.copy()), format='vasp')
        # change symbol back
        for atom in atoms:
            atom.symbol = 'Nb'
        return atoms

    def set_dipole_box_alongy(self, sizen=1):
        e1 = [1., 1., -2.]
        e2 = [0.5, 0.5, 0.5]
        e3 = [1, -1, 0]

        sizen = 1
        n = 7 * sizen
        m = 11 * sizen
        t = 1

        atoms = self.set_bcc_convention([e1, e2, e3], (n, t, m))
        (layerList, distance) = ase.utils.geometry.get_layers(
            atoms, [0, 1, 0], tolerance=0.05)  # the normal of z direction
        layer1_index, layer2_index, layer3_index = [], [], []

        for i in range(len(atoms)):
            if layerList[i] in [1, 4]:
                layer1_index.append(i)
                atoms[i].symbol = 'W'

            if layerList[i] in [2]:
                layer2_index.append(i)
                atoms[i].symbol = 'Mo'

            if layerList[i] in [3, 0]:
                layer3_index.append(i)
                atoms[i].symbol = 'Ta'

        # add shiftment to the supercell #
        atoms = self.cut_half_atoms_new(atoms, "cutz")
        supercell = atoms.get_cell()
        strain = np.mat([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.5, 0.5, 1.0]])
        supercell = strain * supercell
        atoms.set_cell(supercell)
        ase.io.write("perf_poscar", images=atoms, format='vasp')
        # change symbol back
        for atom in atoms:
            atom.symbol = 'Nb'
        return atoms

    def bcc_screw_dipole_alongz_atoms(self, atoms, c1, c2):
        c = tool_elastic_constants.elastic_constants(
            C11=self.pot['c11'], C12=self.pot['c12'], C44=self.pot['c44'])
        burgers = self.pot['lattice'] / 2 * np.array([1., 1., 1.])
        stroh = stroh_solve.Stroh(c, burgers, axes=axes)
        pos = atoms.get_positions()
        atoms = self.bcc_screw_dipole_tools(stroh, pos, atoms, c1, 1)
        atoms = self.bcc_screw_dipole_tools(stroh, pos, atoms, c2, -1)
        return atoms

    def bcc_screw_dipole_tools(self, stroh, pos, atoms, c, sign):
        shiftc = np.ones(np.shape(pos)) * np.array([c[0], c[1], 0.0])
        disp = stroh.displacement(pos - shiftc)
        atoms.set_positions(atoms.get_positions() + sign * np.real(disp))
        return atoms

    def bcc_screw_dipole_alongz_with_image(self, atoms, center):
        # x direction  +10.5 unitx,  +21
        #              -10.5 unitx,  -21
        # y direction  +(11 + 1./3.) unity  + 22
        #              -(11 - 1./3.) unity  - 22

        c = tool_elastic_constants.elastic_constants(
            C11=self.pot['c11'], C12=self.pot['c12'], C44=self.pot['c44'])
        burgers = self.pot['lattice'] / 2 * np.array([1., 1., 1.])

        stroh = stroh_solve.Stroh(c, burgers, axes=axes)

        px = 10.5
        py = 11

        # ydisp = [4 * py, 3 * py + 1. / 3., 2 * py, py + 1. / 3.,
        #          0., -py + 1. / 3., -2 * py, -3 * py + 1. / 3., -4 * py]
        # xdisp = [-4 * px, -3 * px, -2 * px, -1 *
        #          px, 0.0, px, 2 * px, 3 * px, 4 * px]

        ydisp = [2 * py, py + 1. / 3., 0., -py + 1. / 3., -2 * py]
        xdisp = [-2 * px, -1 * px, 0.0, px, 2 * px]

        pos = atoms.get_positions()

        ux = np.sqrt(6) / 3. * self.pot['lattice']
        uy = np.sqrt(2) / 2. * self.pot['lattice']

        sign = 1.0
        for dy in ydisp:
            for dx in xdisp:
                ctmp = np.array([center[0] + dx * ux, center[1] + dy * uy])
                atoms = self.bcc_screw_dipole_tools(
                    stroh, pos, atoms, ctmp, sign)
                sign *= -1
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

        ux = sqrt(6) / 3. * self.pot['lattice']
        uy = sqrt(2) / 2. * self.pot['lattice']
        uz = sqrt(3) / 2. * self.pot['lattice']

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
            c1 = [(sx + 0.5) * ux, (sy + 1. / 3. - 0.95 * 1. / 3.) * uy]
            c2 = [(sx + ix + 0.5) * ux,
                  (sy + 2. / 3. + 0.95 * 1. / 3.) * uy]
        else:
            c1 = [(sx) * ux, (sy + 1. / 3.) * uy]
            c2 = [(sx + ix) * ux, (sy + 2. / 3.) * uy]

        shiftc1 = np.ones(np.shape(pos)) * np.array([c1[0], c1[1], 0.0])
        shiftc2 = np.ones(np.shape(pos)) * np.array([c2[0], c2[1], 0.0])

        disp1 = stroh.displacement(pos - shiftc1)
        disp2 = stroh.displacement(pos - shiftc2)

        if opt in ['pull']:
            radius = 2.0  # find the atoms near the center
            for ps, dp in zip(pos, disp1):
                dis = np.linalg.norm(ps[:2] - c1)
                if (dis < radius):
                    dp[2] += 1. / 6. * uz
                    # add shirt
            for ps, dp in zip(pos, disp2):
                dis = np.linalg.norm(ps[:2] - c2)
                if (dis < radius):
                    dp[2] -= 1. / 6. * uz

        atoms.set_positions(pos + np.real(disp1) - np.real(disp2))

        # periodic boundary conditions ???
        # c2l = [(sx - ix + 0.5) * ux, (sy + 2. / 3. + 0.95 * 1. / 3.) * uy]
        # shft = np.ones(np.shape(pos)) * np.array([c2l[0], c2l[1], 0.0])
        # disp3 = stroh.displacement(pos - shft)
        # atoms.set_positions(atoms.get_positions() - np.real(disp3))

        atoms.wrap(pbc=[1, 1, 1])
        ase.io.write("perf_poscar", atoms_perf, format='vasp')
        ase.io.write('POSCAR', atoms, format='vasp')
        return (atoms, atoms_perf)

    def convert_alongz_to_alongy(self, atoms):
        cell, pos = atoms.get_cell(), atoms.get_positions()
        tm = np.copy(cell[1])
        cell[1] = cell[2]
        cell[2] = tm

        tm = np.copy(cell[:, 1])
        cell[:, 1] = cell[:, 2]
        cell[:, 2] = tm

        tm = np.copy(pos[:, 1])
        pos[:, 1] = pos[:, 2]
        pos[:, 2] = tm

        atoms.set_cell(cell)
        atoms.set_positions(pos)
        return atoms

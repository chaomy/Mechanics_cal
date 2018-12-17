#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-07-05 08:12:30
# @Last Modified by:   chaomy
# @Last Modified time: 2018-12-04 18:26:58

import numpy as np
import atomman as am
from utils import stroh_solve
from numpy import sqrt
import ase.lattice.orthorhombic as otho
import ase.io


class othoHCPFractory(otho.SimpleOrthorhombicFactory):
    bravais_basis = [[0.0, 0.0, 0.0],
                     [0.5, 0.0, 0.5],
                     [0.0, 1. / 2., 1. / 3.],
                     [1. / 2., 1. / 2., 5. / 6.]]


class othoHCPFractoryB(otho.SimpleOrthorhombicFactory):
    bravais_basis = [[0.0, 0.0, 0.0],
                     [0.5, 0.0, 0.5],
                     [1. / 3., 1. / 2., 0.0],
                     [5. / 6., 1. / 2., 1. / 2.]]


othoHCP = othoHCPFractory()
othoHCPB = othoHCPFractoryB()


class md_dislocation_hcp(object):
    # x [1  1  -2 0], y [0, 0, 0, 1], z [1, -1, 0, 0]

    def hcp_edge_dislocation(self):
        sz = (40, 40, 10)
        atoms = othoHCP(latticeconstant=(self.pot['ahcp'],
                                         self.pot['chcp'],
                                         self.pot['ahcp'] * sqrt(3.)),
                        size=sz, symbol=self.pot['element'])

        # let the axes consistent with the elastic constants
        axes = np.array([[1, 0, 0], [0, 0, 1], [0, -1, 0]])
        burgers = self.pot['lattice'] * np.array([1, 0, 0])
        c = am.ElasticConstants()

        # the lattice constant has conventional direction as
        # x - [1 -2 1 0]  y[1, 0, -1, 0] z [0, 0, 0, 1]
        c.hexagonal(C11=61.8, C33=67.5, C12=25.9, C13=21.9, C44=18.2)  # kim
        stroh = stroh_solve.Stroh(c, burgers, axes=axes)
        pos = atoms.get_positions()

        cx = 20 * self.pot["ahcp"] + 0.01
        cy = 20 * self.pot["chcp"] + 0.01
        s1 = np.ones(pos.shape) * np.array([cx, cy, 0.0])
        d1 = stroh.displacement(pos - s1)

        # cy = 60 * self.pot["chcp"] - 0.01
        # s2 = np.ones(pos.shape) * np.array([cx, cy, 0.0])
        # d2 = stroh.displacement(pos - s2)
        atoms.set_positions(pos + np.real(d1))

        # cut a layer normal the burger direction
        # atoms = self.cut_x_normal_atoms(atoms, lata, 1, sqrt(3) / 4.0)
        # atoms = self.cut_x_normal_atoms(atoms)
        atoms = self.cut_y_normal_atoms(atoms)
        atoms = self.cut_x_normal_atoms(atoms)
        self.write_lmp_config_data(atoms)

    def build_screw_basal_hcp_atoms_dipole(self, atoms):
        axes = np.array([[0, 1, 0], [0, 0, 1], [1, 0, 0]])
        burgers = self.pot['lattice'] * np.array([1, 0, 0])

        c = am.ElasticConstants()
        c.hexagonal(C11=64.3218889159844, C33=70.9452244231304, C12=25.4175328222502,
                    C13=20.306421440903, C44=18.0690056385527)  # curtin

        stroh1 = stroh_solve.Stroh(c, burgers, axes=axes)
        # stroh2 = stroh_solve.Stroh(c, -burgers, axes=axes)

        cell = atoms.get_cell()
        pos = atoms.get_positions()

        cx = 0.25 * cell[0, 0] + 0.01
        cy1 = 0.25 * cell[1, 1] + 0.01
        cy2 = 0.75 * cell[1, 1] + 0.01

        pos = atoms.get_positions()
        print(cy1, cy2)

        shift1 = np.ones(pos.shape) * np.array([cx, cy1, 0.0])
        shift2 = np.ones(pos.shape) * np.array([cx, cy2, 0.0])

        disp1 = stroh1.displacement(pos - shift1)
        disp2 = stroh1.displacement(pos - shift2)

        # atoms.set_positions(pos + np.real(disp1))
        atoms.set_positions(pos + np.real(disp1) - np.real(disp2))
        return atoms

    def build_screw_basal_hcp_atoms(self, atoms):
        axes = np.array([[0, 1, 0], [0, 0, 1], [1, 0, 0]])
        burgers = self.pot['lattice'] * np.array([1, 0, 0])
        c = am.ElasticConstants()

        # the lattice constant has conventional direction as
        # x - [1 -2 1 0]  y[1, 0, -1, 0] z [0, 0, 0, 1]
        # c.hexagonal(C11=61.8, C33=67.5, C12=25.9, C13=21.9, C44=18.2)  # kim
        c.hexagonal(C11=64.3218889159844, C33=70.9452244231304, C12=25.4175328222502,
                    C13=20.306421440903, C44=18.0690056385527)  # curtin

        stroh = stroh_solve.Stroh(c, burgers, axes=axes)

        cell = atoms.get_cell()
        cx = 0.25 * cell[0, 0] + 0.01
        cy = 0.50 * cell[1, 1] + 0.01

        pos = atoms.get_positions()
        shift = np.ones(pos.shape) * np.array([cx, cy, 0.0])
        disp = stroh.displacement(pos - shift)
        atoms.set_positions(pos + np.real(disp))

        # pos = atoms.get_positions()
        # shiftN = np.ones(pos.shape) * np.array([cx + 0.5 * cell[0, 0],
        #                                         cy, 0.0])
        # dispN = stroh.displacement(pos - shiftN)
        # atoms.set_positions(pos + np.real(dispN))
        return atoms

    def build_edge_basal_hcp_atoms(self, atoms, center, sign=1):
        # let the axes consistent with the elastic constants
        axes = np.array([[1, 0, 0], [0, 0, 1], [0, -1, 0]])
        burgers = self.pot['lattice'] * np.array([1, 0, 0]) * sign
        c = am.ElasticConstants()

        print("build edge basial hcp")
        # x - [1 -2 1 0]  y[1, 0, -1, 0] z [0, 0, 0, 1]
        # c.hexagonal(C11=61.8, C33=67.5, C12=25.9, C13=21.9, C44=18.2)  # kim
        c.hexagonal(C11=64.3218889159844, C33=70.9452244231304, C12=25.4175328222502,
                    C13=20.306421440903, C44=18.0690056385527)  # curtin
        stroh = stroh_solve.Stroh(c, burgers, axes=axes)
        pos = atoms.get_positions()
        cx = center[0] + 0.01
        cy = center[1] + 0.01
        print(cx, cy)
        shift = np.ones(pos.shape) * np.array([cx, cy, 0.0])
        disp = stroh.displacement(pos - shift)
        atoms.set_positions(pos + np.real(disp))
        return atoms

    def build_screw_basal_hcp(self):
        sz = (100, 100, 5)
        ux = self.pot["ahcp"] * sqrt(3)
        uy = self.pot["chcp"]
        uz = self.pot['ahcp']

        atoms = othoHCPB(latticeconstant=(ux, uy, uz),
                         size=sz, symbol=self.pot['element'])

        # let the axes consistent with the elastic constants
        axes = np.array([[0, 1, 0], [0, 0, 1], [1, 0, 0]])
        burgers = self.pot['lattice'] * np.array([1, 0, 0])
        c = am.ElasticConstants()

        # the lattice constant has conventional direction as
        # x - [1 -2 1 0]  y[1, 0, -1, 0] z [0, 0, 0, 1]
        # c.hexagonal(C11=61.8, C33=67.5, C12=25.9, C13=21.9, C44=18.2)
        # c.hexagonal(C11=62.807, C33=69.615, C12=25.974,
        #             C13=21.184, C44=17.138)                       # Kim
        c.hexagonal(C11=64.3218889159844, C33=70.9452244231304, C12=25.4175328222502,
                    C13=20.306421440903, C44=18.0690056385527)      # curtin

        stroh = stroh_solve.Stroh(c, burgers, axes=axes)
        pos = atoms.get_positions()

        cx = 0.5 * sz[0] * ux + 0.01
        cy = 0.5 * sz[1] * uy + 0.01

        print(cx, cy)
        self.write_lmp_config_data(atoms, "before.txt")
        shift = np.ones(pos.shape) * np.array([cx, cy, 0.0])
        disp = stroh.displacement(pos - shift)
        atoms.set_positions(pos + np.real(disp))
        self.write_lmp_config_data(atoms)
        ase.io.write("SCREW.cfg", atoms, format="cfg")
        return atoms

    def build_edge_basal_hcp(self):
        sz = (200, 100, 3)
        atoms = othoHCP(latticeconstant=(self.pot['ahcp'], self.pot['chcp'],
                                         self.pot['ahcp'] * sqrt(3.)),
                        size=sz, symbol=self.pot['element'])

        # let the axes consistent with the elastic constants
        axes = np.array([[1, 0, 0], [0, 0, 1], [0, -1, 0]])
        burgers = self.pot['lattice'] * np.array([1, 0, 0])
        c = am.ElasticConstants()

        # the lattice constant has conventional direction as
        # x - [1 -2 1 0]  y[1, 0, -1, 0] z [0, 0, 0, 1]
        # c.hexagonal(C11=61.8, C33=67.5, C12=25.9, C13=21.9, C44=18.2)  # kim
        # c.hexagonal(C11=62.369, C33=67.788, C12=26.252,
        #             C13=22.113, C44=18.274)  # Coco
        # c.hexagonal(C11=62.807, C33=69.615, C12=25.974,
        #             C13=21.184, C44=17.138)                       # Kim
        c.hexagonal(C11=64.3556, C33=70.9849, C12=25.460289,
                    C13=20.333408, C44=18.06331)                       # Curtin

        stroh = stroh_solve.Stroh(c, burgers, axes=axes)
        pos = atoms.get_positions()

        cx = 0.5 * sz[0] * self.pot["ahcp"] + 0.01
        cy = 0.5 * sz[1] * self.pot["chcp"] + 0.01

        print(cx, cy)
        shift = np.ones(pos.shape) * np.array([cx, cy, 0.0])
        disp = stroh.displacement(pos - shift)
        atoms.set_positions(pos + np.real(disp))
        self.write_lmp_config_data(atoms)
        return atoms

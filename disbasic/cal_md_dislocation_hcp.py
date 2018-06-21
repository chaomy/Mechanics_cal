#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-07-05 08:12:30
# @Last Modified by:   chaomy
# @Last Modified time: 2018-06-06 15:58:37

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
        atoms = othoHCP(latticeconstant=(self.pot['ahcp'], self.pot['chcp'],
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

    def build_screw_basal_hcp_atoms(self, atoms):
        axes = np.array([[0, 1, 0], [0, 0, 1], [1, 0, 0]])
        burgers = self.pot['lattice'] * np.array([1, 0, 0])
        c = am.ElasticConstants()

        # the lattice constant has conventional direction as
        # x - [1 -2 1 0]  y[1, 0, -1, 0] z [0, 0, 0, 1]
        c.hexagonal(C11=61.8, C33=67.5, C12=25.9, C13=21.9, C44=18.2)  # kim

        stroh = stroh_solve.Stroh(c, burgers, axes=axes)
        strohN = stroh_solve.Stroh(c, -burgers, axes=axes)

        cell = atoms.get_cell()
        cx = 0.25 * cell[0, 0] + 0.01
        cy = 0.50 * cell[1, 1] + 0.01

        pos = atoms.get_positions()
        shift = np.ones(pos.shape) * np.array([cx, cy, 0.0])
        disp = stroh.displacement(pos - shift)
        atoms.set_positions(pos + np.real(disp))

        pos = atoms.get_positions()
        shiftN = np.ones(pos.shape) * np.array([cx + 0.5 * cell[0, 0],
                                                cy, 0.0])
        dispN = stroh.displacement(pos - shiftN)
        atoms.set_positions(pos + np.real(dispN))
        return atoms

    def build_edge_basal_hcp_atoms(self, atoms, center):
        # let the axes consistent with the elastic constants
        axes = np.array([[1, 0, 0], [0, 0, 1], [0, -1, 0]])
        burgers = self.pot['lattice'] * np.array([1, 0, 0])
        c = am.ElasticConstants()

        # x - [1 -2 1 0]  y[1, 0, -1, 0] z [0, 0, 0, 1]
        c.hexagonal(C11=61.8, C33=67.5, C12=25.9, C13=21.9, C44=18.2)  # kim
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
        sz = (30, 30, 10)

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
        c.hexagonal(C11=61.8, C33=67.5, C12=25.9, C13=21.9, C44=18.2)  # kim

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
        sz = (30, 30, 10)
        atoms = othoHCP(latticeconstant=(self.pot['ahcp'], self.pot['chcp'],
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

        cx = 0.5 * sz[0] * self.pot["ahcp"] + 0.01
        cy = 0.5 * sz[1] * self.pot["chcp"] + 0.01

        print(cx, cy)
        shift = np.ones(pos.shape) * np.array([cx, cy, 0.0])
        disp = stroh.displacement(pos - shift)
        atoms.set_positions(pos + np.real(disp))

        self.write_lmp_config_data(atoms)
        ase.io.write("EDGE.cfg", atoms, format="cfg")
        return atoms

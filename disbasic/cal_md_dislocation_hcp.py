#!/usr/bin/env python
# encoding: utf-8
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-07-05 08:12:30
# @Last Modified by:   chaomy
# @Last Modified time: 2018-03-27 16:10:27

import numpy as np
import atomman as am
from utils import stroh_solve
from numpy import sqrt
import ase.lattice.orthorhombic as otho


class othoHCPFractory(otho.SimpleOrthorhombicFactory):
    bravais_basis = [[0.0, 0.0, 0.0],
                     [0.5, 0.0, 0.5],
                     [0.0, 1. / 2., 1. / 3.],
                     [1. / 2., 1. / 2., 5. / 6.]]

    # bravais_basis = [[0.0, 0.0, 0.0],
    #                  [0.5, 0.5, 0.0],
    #                  [0.0, 1. / 3., 1. / 2.],
    #                  [1 / 2., 5. / 6., 1. / 2.]]

othoHCP = othoHCPFractory()


class md_dislocation_hcp(object):

    def hcp_edge_dislocation(self):
        sz = (40, 40, 10)

        # x [1  1  -2 0], y [0, 0, 0, 1], z [1, -1, 0, 0]
        atoms = othoHCP(latticeconstant=(self.pot['ahcp'],
                                         self.pot['chcp'],
                                         self.pot['ahcp'] * sqrt(3.)),
                        size=sz,
                        symbol=self.pot['element'])

        lata, latc = self.pot["ahcp"], self.pot["chcp"]
        self.burger = self.pot["lattice"]

        # cut a layers to generate free surface
        cell = atoms.get_cell()

        c1 = [1. / 2. * cell[0, 0], 1. / 4. * cell[1, 1]]
        c2 = [1. / 2. * cell[0, 0], 3. / 4. * cell[1, 1]]

        atoms = self.intro_dipole_edge_atoms(atoms, c1, c2)

        # cut a layer normal the burger direction
        # atoms = self.cut_x_normal_atoms(atoms, lata, 1, sqrt(3) / 4.0)
        # atoms = self.cut_x_normal_atoms(atoms)

        self.write_lmp_config_data(atoms)

    def buildhcp(self):
        sz = (30, 30, 10)
        atoms = othoHCP(latticeconstant=(self.pot['ahcp'],
                                         self.pot['chcp'],
                                         self.pot['ahcp'] * sqrt(3.)),
                        size=sz,
                        symbol=self.pot['element'])

        axes = np.array([[1, 0, 0],
                         [0, 1, 0],
                         [0, 0, 1]])

        burgers = 0.25 * self.pot['lattice'] * np.array([1, 0, 0])

        c = am.ElasticConstants()
        c.hexagonal(C11=326.08, C33=357.50, C12=129.56, C13=119.48, C44=92.54)

        stroh = stroh_solve.Stroh(c, burgers, axes=axes)
        pos = atoms.get_positions()

        cx = 0.5 * sz[0] * self.pot["ahcp"]
        cy = 0.5 * sz[1] * self.pot["ahcp"] * sqrt(3)
        print(cx, cy)
        shift = np.ones(pos.shape) * np.array([cx, cy, 0.0])

        disp = stroh.displacement(pos - shift)
        atoms.set_positions(pos + np.real(disp))

        self.write_lmp_config_data(atoms)
        return atoms

# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2018-02-20 14:11:07
# @Last Modified by:   chaomy
# @Last Modified time: 2018-02-21 03:51:39

import ase.lattice.orthorhombic as otho
import ase.lattice.cubic as cubic
from numpy import sqrt
import numpy as np


class othoHCPFractory(otho.SimpleOrthorhombicFactory):
    bravais_basis = [[0.0, 0.0, 0.0],
                     [0.5, 0.0, 0.5],
                     [0.0, 1. / 2., 1. / 3.],
                     [1. / 2., 1. / 2., 5. / 6.]]

othoHCP = othoHCPFractory()


class gb_hcp_dis(object):

    def make_gb(self):
        ux, uy, uz = self.pot['ahcp'], self.pot[
            'chcp'], self.pot['ahcp'] * sqrt(3.)

        sz = (160, 70, 4)
        atoms = othoHCP(latticeconstant=(ux, uy, uz),
                        size=sz,
                        symbol=self.pot['element'])

        lata, latc = self.pot["ahcp"], self.pot["chcp"]
        self.burger = self.pot["lattice"]

        cell = atoms.get_cell()

        lob = np.array([ux * 130,  0.0 + 20, 0])
        hib = np.array([ux * 150, uy * 70 - 20, uz * 5])

        atoms = self.intro_single_edge_atoms(
            atoms, center=[ux * 40, 35 * uy, 25 * uz])
        atoms = self.make_cubic("in", atoms, lob, hib)

        atoms2 = othoHCP(latticeconstant=(ux, uy, uz),
                         size=(140, 150, 5),
                         symbol="Nd")

        atoms2.rotate(23.0, 'z')
        pos = atoms2.get_positions()
        pos += np.array([ux * 120,  -40, 0])
        atoms2.set_positions(pos)
        atoms2 = self.make_cubic("out", atoms2, lob - 1.0, hib + 1.0)

        atoms = self.cut_y_normal_atoms(atoms)
        # atoms = self.cut_x_normal_atoms(atoms)

        atoms.extend(atoms2)
        atoms = self.assign_ynormal_fixatoms(atoms)
        self.write_lmp_config_data(atoms, "lmp_init.txt")

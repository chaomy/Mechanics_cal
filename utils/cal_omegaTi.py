#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2018-03-28 21:31:56
# @Last Modified by:   chaomy
# @Last Modified time: 2018-03-29 00:52:01

import ase
import gn_config
import ase.lattice.orthorhombic as otho
import md_pot_data


class omegaTiFactory(otho.SimpleOrthorhombicFactory):
    bravais_basis = [[0.00000, 0.00000, 0.00000],
                     [2. / 3., 1. / 3., 1. / 2.],
                     [1. / 3., 2. / 3., 1. / 2.]]

omegaTi = omegaTiFactory()


class cal_omegaTi(gn_config.gnStructure):

    def __init__(self):
        self.pot = md_pot_data.va_pot.Nb_pbe
        gn_config.gnStructure.__init__(self, self.pot)

    def build_omegaTi(self):
        atoms = omegaTi(latticeconstant=(4.887, 4.887, 2.678),
                       size=(1, 1, 1), symbol=('Nb'))
        ase.io.write("POSCAR_omegaTi", atoms, "vasp")
        self.write_lmp_config_data(atoms)

if __name__ == '__main__':
    drv = cal_omegaTi()
    drv.build_omegaTi()

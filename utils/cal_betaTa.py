#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2018-03-28 21:31:56
# @Last Modified by:   chaomy
# @Last Modified time: 2018-03-28 23:23:05

import ase
import gn_config
import ase.lattice.orthorhombic as otho
import md_pot_data


class betaTaFactory(otho.SimpleOrthorhombicFactory):
    bravais_basis = [[0.50000, 0.00000, 0.22800],
                     [0.00000, 0.50000, 0.77200],

                     [0.81420, 0.31420, 0.00260],
                     [0.18580, 0.68580, 0.00260],
                     [0.31420, 0.18580, 0.99740],
                     [0.68580, 0.81420, 0.99740],

                     [0.03430, 0.12670, 0.25460],
                     [0.12670, 0.96570, 0.74540],
                     [0.96570, 0.87330, 0.25460],
                     [0.87330, 0.03430, 0.74540],
                     [0.37330, 0.46570, 0.25460],
                     [0.46570, 0.62670, 0.74540],
                     [0.62670, 0.53430, 0.25460],
                     [0.53430, 0.37330, 0.74540],

                     [0.60330, 0.10330, 0.76400],
                     [0.10330, 0.39670, 0.23600],
                     [0.39670, 0.89670, 0.76400],
                     [0.89670, 0.60330, 0.23600],

                     [0.75980, 0.06770, 0.23500],
                     [0.06770, 0.24020, 0.76500],
                     [0.24020, 0.93230, 0.23500],
                     [0.93230, 0.75980, 0.76500],
                     [0.56770, 0.25980, 0.23500],
                     [0.25980, 0.43230, 0.76500],
                     [0.43230, 0.74020, 0.23500],
                     [0.74020, 0.56770, 0.76500],

                     [0.31960, 0.18040, 0.49120],
                     [0.18040, 0.68040, 0.50880],
                     [0.68040, 0.81960, 0.49120],
                     [0.81960, 0.31960, 0.50880]]

betaTa = betaTaFactory()


class cal_betaTa(gn_config.gnStructure):

    def __init__(self):
        self.pot = md_pot_data.va_pot.Nb_pbe
        gn_config.gnStructure.__init__(self, self.pot)

    def build_betaTa(self):
        atoms = betaTa(latticeconstant=(10.184, 10.184, 5.371),
                       size=(1, 1, 1), symbol=('Nb'))
        ase.io.write("POSCAR_betaTa", atoms, "vasp")
        self.write_lmp_config_data(atoms)

if __name__ == '__main__':
    drv = cal_betaTa()
    drv.build_betaTa()

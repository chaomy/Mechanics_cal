#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-06-25 14:28:58
# @Last Modified by:   chaomy
# @Last Modified time: 2017-07-06 23:21:22

import ase
import ase.io
import numpy as np
import md_pot_data
import tool_elastic_constants
import stroh_solve
import ase.lattice
import cal_md_dislocation


class cal_dis_dipole(object):

    def __init__(self, pot=None):
        if pot is None:
            pot = md_pot_data.qe_pot.vca_W50Re50
        self.pot = pot
        self.mddis_drv = cal_md_dislocation.md_dislocation(self.pot)
        return

    def set_dipole_box(self, sizen=1):
        n = 7 * sizen
        m = 11 * sizen
        t = 1 * sizen

        print self.pot
        alat = self.pot['lattice']
        atoms = ase.lattice.cubic.BodyCenteredCubic(directions=[[1., 1., -2.],
                                                                [-1., 1., 0],
                                                                [0.5, 0.5, 0.5]],
                                                    latticeconstant=alat,
                                                    size=(n, m, t),
                                                    symbol=self.pot['element'])

        atoms = self.mddis_drv.cut_half_atoms_new(atoms, "cuty")
        supercell = atoms.get_cell()
        strain = np.mat([[1.0, 0.0, 0.0],
                         [0.5, 1.0, 0.5],
                         [0.0, 0.0, 1.0]])
        supercell = strain * supercell
        atoms.set_cell(supercell)
        atoms.wrap(pbc=[1, 1, 1])
        return atoms

    def print_dis_constants(self, stroh):
        print "K tensor", stroh.K_tensor
        print "K (biKijbj)", stroh.K_coeff, "eV/A"
        print "pre-ln alpha = biKijbj/4pi", stroh.preln, "ev/A"
        return

    def bcc_screw_dipole_configs_alongz(self, sizen=1):
        c = tool_elastic_constants.elastic_constants(C11=self.pot['c11'],
                                                     C12=self.pot['c12'],
                                                     C44=self.pot['c44'])
        # c = tool_elastic_constants.elastic_constants(C11=502, C12=173, C44=138)

        axes = np.array([[1, 1, -2],
                         [-1, 1, 0],
                         [1, 1, 1]])

        burgers = self.pot['lattice'] / 2 * np.array([1., 1., 1.])
        stroh = stroh_solve.Stroh(c, burgers, axes=axes)

        atoms = self.set_dipole_box()
        atoms_perf = atoms.copy()
        pos = atoms.get_positions()

        unitx = np.sqrt(6) / 3. * self.pot['lattice']
        unity = np.sqrt(2) / 2. * self.pot['lattice']
        sx = 10.0 * sizen
        sy = 5 * sizen
        ix = 10.5 * sizen

        # c1 = 1. / 3. * np.sum(self.pot['core1'], axis=0)
        # c2 = 1. / 3. * np.sum(self.pot['core2'], axis=0)
        # shiftc1 = np.ones(np.shape(pos)) * np.array([c1[0, 0], c1[0, 1], 0.0])
        # shiftc2 = np.ones(np.shape(pos)) * np.array([c2[0, 0], c2[0, 1], 0.0])

        c1 = [(sx) * unitx, (sy + 1. / 3.) * unity]
        c2 = [(sx + ix) * unitx, (sy + 2. / 3.) * unity]
        shiftc1 = np.ones(np.shape(pos)) * np.array([c1[0], c1[1], 0.0])
        shiftc2 = np.ones(np.shape(pos)) * np.array([c2[0], c2[1], 0.0])

        disp1 = stroh.displacement(pos - shiftc1)
        disp2 = stroh.displacement(pos - shiftc2)

        atoms.set_positions(pos + np.real(disp1) - np.real(disp2))
        # ase.io.write('POSCAR', atoms, format='vasp')
        return (atoms, atoms_perf)

if __name__ == '__main__':
    drv = cal_dis_dipole()
    drv.bcc_screw_dipole_configs_alongz()

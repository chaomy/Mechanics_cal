#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2018-03-28 21:31:56
# @Last Modified by:   chaomy
# @Last Modified time: 2018-06-26 00:56:27

import ase
import ase.io
import gn_config
import ase.lattice.cubic as cubic
import md_pot_data
import glob
import os


class A15Factory(cubic.SimpleCubicFactory):
    bravais_basis = [[0., 0., 0.],
                     [0.5, 0.5, 0.5],
                     [0.25, 0.5, 0.],
                     [0.75, 0.5, 0.],
                     [0., 0.25, 0.5],
                     [0., 0.75, 0.5],
                     [0.5, 0., 0.25],
                     [0.5, 0., 0.75]]


a15 = A15Factory()


class cal_a15(gn_config.gnStructure):

    def __init__(self):
        self.pot = md_pot_data.va_pot.Nb_pbe
        gn_config.gnStructure.__init__(self, self.pot)

    def build_a15(self):
        atoms = a15(latticeconstant=5.292, size=(1, 1, 1), symbol=('Nb'))
        ase.io.write("POSCAR_A15", atoms, "vasp")
        self.write_lmp_config_data(atoms)

    def convert(self):
        atoms = ase.io.read("CONTCAR", format='vasp')
        self.write_lmp_config_data(atoms)

    def convert_lmp_to_POSCAR(self):
        fls = glob.glob("dump.*")
        for i in range(len(fls)):
            atoms = ase.io.read(fls[i], format="lammps-dump")
            ase.io.write("poscar.{:02d}".format(
                i), images=atoms, format="vasp")
            mdir = "dir_{:02d}".format(i)
            os.mkdir(mdir)
            os.system("cp INCAR KPOINTS va.pbs POTCAR {}".format(mdir))
            os.system("cp poscar.{:02d} {}/POSCAR".format(i, mdir))

if __name__ == '__main__':
    drv = cal_a15()
    # drv.build_a15()
    drv.convert_lmp_to_POSCAR()

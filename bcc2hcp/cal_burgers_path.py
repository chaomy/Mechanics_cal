#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2018-03-12 13:55:06
# @Last Modified by:   chaomy
# @Last Modified time: 2018-03-15 02:22:05

import os
import ase
import ase.io
import get_data
import gn_config
import md_pot_data
import numpy as np
import ase.lattice.orthorhombic as otho


class BCCOTHO(otho.SimpleOrthorhombicFactory):
    bravais_basis = [[0.0, 0.0, 0.0],
                     [0.5, 0.5, 0.0],
                     [0.5, 0.0, 0.5],
                     [0.0, 0.5, 0.5]]
BCC2HCP = BCCOTHO()


class burgers_path(get_data.get_data, gn_config.bcc):

    def __init__(self):
        self.pot = md_pot_data.va_pot.Nb_pbe
        gn_config.bcc.__init__(self, self.pot)
        la = self.pot['latbcc']
        self.ubcc = BCC2HCP(latticeconstant=(la * np.sqrt(2), la,
                                             la * np.sqrt(2)),
                            size=(1, 1, 1), symbol=self.pot['element'])

    def add_strain(self):
        atoms = self.set_bcc_convention(in_size=(2, 2, 2))
        ase.io.write("bcc0", images=atoms, format='vasp')

        a = 1.0
        t = 3. / (4. * np.sqrt(2))
        strain = np.mat([[0.5 * a + t, 0.5 * a - t, 0.0],
                         [0.5 * a - t, 0.5 * a + t, 0.0],
                         [0.0, 0.0, np.sqrt(3) / 2.]])
        print(strain)
        cell = atoms.get_cell()
        pos = np.mat(atoms.get_positions())
        pos = pos * strain
        cell = strain * cell
        atoms.set_cell(cell)
        atoms.set_positions(pos)
        ase.io.write("bcc1", images=atoms, format='vasp')
        print(np.sqrt(np.sqrt(2)))

    def othos_burgers_path(self):
        la = self.pot['latbcc']
        atoms = self.ubcc
        ase.io.write("bcc0", images=atoms, format='vasp')
        atoms = self.burgers_path_strain(atoms)
        atoms = self.burgers_path_shuffle(atoms)

    def burgers_path_strain(self, atoms):  # strain
        la = self.pot['latbcc']
        ahcp = self.pot['ahcp']
        chcp = self.pot['chcp']
        elong = ahcp * np.sqrt(3) / (la * np.sqrt(2))
        shorten = ahcp / la
        cdir = chcp / (np.sqrt(2) * la)
        print(elong, shorten, cdir)
        strain = np.mat([[elong, 0, 0],
                         [0, shorten, 0],
                         [0, 0, cdir]])
        cell = atoms.get_cell()
        pos = np.mat(atoms.get_positions())
        pos = pos * strain
        cell = strain * cell
        atoms.set_cell(cell)
        atoms.set_positions(pos)
        return atoms

    def burgers_path_shuffle(self, atoms):  # shuffle
        disp = self.pot['ahcp'] * np.sqrt(3) / 2. * 1. / 3.
        disp_matrix = np.mat([[0.0, 0.0, 0.0], [0.0, 0.0, 0.0],
                              [disp, 0.0, 0.0], [disp, 0.0, 0.0]])
        atoms.translate(disp_matrix)
        ase.io.write("bcc2", images=atoms, format="vasp")
        return atoms

    def burgers_path_strain_mesh(self):
        atoms = self.ubcc
        cnt = 0
        dlt = 0.05
        ahcp = self.pot['ahcp'] - 10 * dlt
        chcp = self.pot['chcp'] - 10 * dlt

        for i in range(20):
            self.pot['ahcp'] = ahcp + i * dlt
            for j in range(20):
                self.pot['chcp'] = chcp + j * dlt
                deform_atoms = self.burgers_path_strain(atoms.copy())
                ase.io.write("POSCAR", images=deform_atoms, format="vasp")
                mdir = 'dir_{:03}'.format(cnt)
                self.mymkdir(mdir)
                self.prepare_vasp_inputs(mdir)
                cnt += 1

    def burgers_path_mesh(self):
        atoms = self.ubcc
        atoms = self.burgers_path_strain(atoms)
        npts = 18
        for i in range(npts):
            mdir = 'dir_{:03}'.format(i)
            self.mymkdir(mdir)

            trans_atoms = atoms.copy()
            disp = self.pot['ahcp'] * np.sqrt(3) / 2. * i / npts
            disp_matrix = np.mat([[0.0, 0.0, 0.0], [0.0, 0.0, 0.0],
                                  [disp, 0.0, 0.0], [disp, 0.0, 0.0]])
            trans_atoms.translate(disp_matrix)
            ase.io.write("POSCAR", images=trans_atoms, format="vasp")
            self.prepare_vasp_inputs(mdir)

    def prepare_vasp_inputs(self, mdir):
        os.system("mv POSCAR {}".format(mdir))
        os.system("cp va.pbs {}".format(mdir))
        os.system("cp KPOINTS {}".format(mdir))
        os.system("cp INCAR {}".format(mdir))
        os.system("cp POTCAR {}".format(mdir))


if __name__ == '__main__':
    drv = burgers_path()
    # drv.burgers_path_mesh()
    drv.burgers_path_strain_mesh()

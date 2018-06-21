#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-07-05 08:12:30
# @Last Modified by:   chaomy
# @Last Modified time: 2018-05-17 21:10:38


import numpy as np
from . import cal_add_strain_otho


class cal_add_strain(cal_add_strain_otho.strain_otho):

    def __init__(self):
        cal_add_strain_otho.strain_otho.__init__(self)
        self.addstrain = {'vol': self.add_volumeric_strain,
                          'otho': self.volume_conserving_ortho_strain_atoms,
                          'mono': self.volume_conserving_mono_strain_atoms}

    def volume_conserving_ortho_strain(self, delta):
        atom_number, supercell_base, comment, atom_position =\
            self.read_vasp_poscar()
        OrthoM = np.mat([[1 + delta, 0, 0],
                         [0, 1 - delta, 0],
                         [0, 0, 1 / (1 - delta**2)]])
        Transformed_Base = OrthoM * supercell_base
        with open("POSCAR", 'w') as fid:
            fid.write("# Screw Bcc\n")
            fid.write("%12.6f\n" % (self.lattice_constant))
            fid.write("%12.6f %12.6f %12.6f\n" %
                      (Transformed_Base[0, 0],
                       Transformed_Base[0, 1],
                       Transformed_Base[0, 2]))
            fid.write("%12.6f %12.6f %12.6f\n" %
                      (Transformed_Base[1, 0],
                       Transformed_Base[1, 1],
                       Transformed_Base[1, 2]))
            fid.write("%12.6f %12.6f %12.6f\n" %
                      (Transformed_Base[2, 0],
                       Transformed_Base[2, 1],
                       Transformed_Base[2, 2]))
            fid.write("%d\n" % (atom_number))
            fid.write(comment)
            for i in range(atom_number):
                fid.write("%12.6f %12.6f %12.6f\n" %
                          (atom_position[0, i],
                           atom_position[1, i],
                           atom_position[2, i]))
            fid.close()

    def add_volumeric_strain(self, delta, atoms):
        strain = np.mat([[1 + delta, 0, 0],
                         [0, 1 + delta, 0],
                         [0, 0, 1 + delta]], "float")
        cell = strain * np.mat(atoms.get_cell())
        positions = np.mat(atoms.get_positions()) * strain
        atoms.set_positions(positions)
        atoms.set_cell(cell)
        return atoms

    def volume_conserving_ortho_strain_atoms(self, delta, atoms):
        strain = np.mat([[1 + delta, 0, 0],
                         [0, 1 - delta, 0],
                         [0, 0, 1 / (1 - delta**2)]])
        cell = strain * atoms.get_cell()
        pos = np.mat(atoms.get_positions()) * strain
        atoms.set_positions(pos)
        atoms.set_cell(cell)
        return atoms

    def volume_conserving_mono_strain_atoms(self, delta, atoms):
        strain = np.mat([[1, delta, 0],
                         [delta, 1, 0],
                         [0, 0, 1 / (1 - delta**2)]])
        cell = atoms.get_cell()
        pos = np.mat(atoms.get_positions())
        cell = strain * cell
        pos = pos * strain
        atoms.set_positions(pos)
        atoms.set_cell(cell)
        return atoms

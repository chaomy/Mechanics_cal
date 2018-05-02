#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-07-05 08:12:30
# @Last Modified by:   chaomy
# @Last Modified time: 2018-05-01 22:37:07


import numpy as np
import math
from ase import Atoms
import cal_add_strain
import cal_cut_cell
import cal_intro_iso_dis
import cal_intro_iso_dis_image
import cal_intro_ani_dis


class cubic_cij:
    c11 = None
    c12 = None
    c44 = None


class vasp_change_box(cal_intro_iso_dis.cal_intro_iso_dis,
                      cal_intro_ani_dis.cal_intro_ani_dis,
                      cal_intro_iso_dis_image.cal_intro_iso_dis_image,
                      cal_add_strain.cal_add_strain,
                      cal_cut_cell.cal_cut_cell):

    def __init__(self):
        self.pi = 3.141592653589793
        if self.pot['structure'] in ['bcc']:
            self.burger = np.sqrt(3.) / 2. * self.pot['lattice']
            self.screw_coeff = self.burger / (2. * self.pi)
            self.Edge_coeff = self.burger / (2. * self.pi)
            self.P = 0.33

    def intro_edge_nuclei(self, atoms, center=None,
                          sign=1, orient=None,
                          zzone=None):
        if orient is not None:
            bdir = orient[0]
            ndir = orient[1]
            tdir = orient[2]
        else:
            tdir = 2  # line direction
            bdir = 0  # burger
            ndir = 1  # glide plane normal

        xc0 = center[bdir]
        yc0 = center[ndir]

        atom_position = atoms.get_positions()

        coeff = self.Edge_coeff
        alpha = 0.1

        for i in range(len(atoms)):
            xc1 = xc0
            yc1 = yc0 * (1 + math.tanh(alpha * (atom_position[i, tdir] - zzone[0]))) - \
                yc0 * (1 + math.tanh(alpha *
                                     (atom_position[i, tdir] - zzone[1])))
            yc1 *= 0.5

            dx1 = atom_position[i, bdir] - xc1
            dy1 = atom_position[i, ndir] - yc1

            A1 = (1 - self.P) * (dx1**2 + dy1**2)

            if A1 != 0:
                ux1 = coeff * (np.arctan2(dy1, dx1) + (dx1 * dy1) / (2. * A1))
                uy1 = (1 - 2 * self.P) / (4 * (1 - self.P)) * \
                    math.log(dx1**2 + dy1**2) + (dx1**2 - dy1**2) / \
                    (4 * (1 - self.P) * (dx1**2 + dy1**2))
                uy1 *= -coeff
            else:
                ux1 = uy1 = 0

            atom_position[i, bdir] = atom_position[i, bdir] + sign * ux1
            atom_position[i, ndir] = atom_position[i, ndir] + sign * uy1

        atoms.set_positions(atom_position)
        return atoms

    ###########################################################################
    # Sun Mar 26 13:01:49 2017
    # intro dipole screw used for LAMMPS
    # x along the [-1 1 0];
    # y is along the dislocation line t
    # z along the [1 1 -2];
    ###########################################################################
    def intro_dipole_screw_atoms_LMP(self, atoms, center):

        lxid = 0  # glide plane
        lzid = 1  # dislocation line
        lyid = 2  # y

        atoms.wrap(pbc=[1, 1, 1])

        supercell_base = atoms.get_cell()
        atom_position = atoms.get_positions()

        xc1 = center[0][0]
        yc1 = center[0][1]

        xc2 = center[1][0]
        yc2 = center[1][1]

        print(("center 1 (%g  %g)\n" % (xc1, yc1)))
        print(("center 2 (%g  %g)\n" % (xc2, yc2)))

        for i in range(len(atom_position)):
            dx1, dx2 = atom_position[i, lxid] - \
                xc1, atom_position[i, lxid] - xc2
            dy1, dy2 = atom_position[i, lyid] - \
                yc1, atom_position[i, lyid] - yc2

            theta1 = np.arctan2(dy1, dx1)
            theta2 = np.arctan2(dy2, dx2)

            dz1 = self.screw_coeff * theta1
            dz2 = self.screw_coeff * theta2

            atom_position[i, lzid] = atom_position[i, lzid] + dz1 - dz2

        atoms.set_positions(atom_position)
        return atoms

    def intro_dipole_screw(self):
        atom_number, supercell_base, comment, atom_position = \
            self.read_vasp_poscar()

        #  add_disp = [0.15 * (supercell_base[0, 0] + supercell_base[1, 0]),
        #  0.15 * (supercell_base[0, 1] + supercell_base[1, 1])]

        add_disp = [0.5 * (supercell_base[0, 0] + supercell_base[1, 0]),
                    0.5 * (supercell_base[0, 1] + supercell_base[1, 1])]

        xc1 = 0 + add_disp[0]
        yc1 = 0 + add_disp[1]

        xc2 = 0.5 * supercell_base[0, 0] + add_disp[0]
        yc2 = 0.0 * supercell_base[1, 1] + add_disp[1]

        for i in range(atom_number):
            dx1, dx2 = atom_position[0, i] - xc1, atom_position[0, i] - xc2
            dy1, dy2 = atom_position[1, i] - yc1, atom_position[1, i] - yc2

            theta1 = np.arctan2(dy1, dx1)
            theta2 = np.arctan2(dy2, dx2)

            dz1 = self.screw_coeff * theta1
            dz2 = self.screw_coeff * theta2

            atom_position[2, i] = atom_position[2, i] + dz1 - dz2

    def lmp_map_atoms_list(self, atoms, atoms_tobe_change):
        position_sample = atoms.get_positions()
        position_change = atoms_tobe_change.get_positions()
        num_atoms = len(atoms)
        map_list = [[], []]
        for i in range(num_atoms):
            for k in range(num_atoms):
                if np.linalg.norm(position_sample[i] - position_change[k]) < 0.001:
                    map_list[0].append(i)
                    map_list[1].append(k)

        print(len(map_list[1]))
        return map_list

    def add_perturbation(self, atoms, disp):
        coeff = disp
        for atom in atoms:
            sign_ran = np.random.rand(3)
            disp_ran = np.random.rand(3)
            for i in range(3):
                if sign_ran[i] >= 0.5:
                    atom.position[i] += coeff * disp_ran[i]
                else:
                    atom.position[i] -= coeff * disp_ran[i]
        return atoms

    def map_atoms_list(self, atoms, atoms_tobe_change):
        position_sample = atoms.get_positions()
        position_change = atoms_tobe_change.get_positions()
        num_atoms = len(atoms)
        map_list = [[], []]

        for i in range(num_atoms):
            for k in range(num_atoms):
                if np.linalg.norm(position_sample[i] - position_change[k]) < 0.001:
                    map_list[0].append(i)
                    map_list[1].append(k)

        print(len(map_list[1]))
        return map_list

    def sort_atoms(self, atoms, map_list):
        new_atoms = Atoms()
        print(map_list[0])
        print(map_list[1])

        for i in map_list[1]:
            #  print i
            print(atoms[i].position)
            new_atoms.append(atoms[i])
        new_atoms.set_cell(atoms.get_cell())
        return new_atoms

    def intro_ani_edge_fcc(self):
        self.fcc_edge()

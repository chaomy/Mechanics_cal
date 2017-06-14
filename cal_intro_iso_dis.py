#!/usr/bin/env python
# encoding: utf-8

###################################################################
#
# File Name : cal_intro_iso_dis.py
#
###################################################################
#
# Purpose : introduce isotrpic dislocations
#
# Creation Date :
# Last Modified :
# Created By    : Chaoming Yang
#
###################################################################


import numpy as np
import math


class cal_intro_iso_dis(object):
    def __init__(self, lattice_constant=3.16):
        self.lattice_constant = lattice_constant
        self._set_intro_coeff()
        return

    def _set_intro_coeff(self):
        self.pi = 3.141592653589793
        self.burger = np.sqrt(3.) / 2. * self.lattice_constant
        self.screw_coeff = self.burger / (2. * self.pi)
        self.Edge_coeff = self.burger / (2. * self.pi)
        print self.screw_coeff
        self.P = 0.33
        return

    def _intro_single_screw(self):
        atom_number, supercell_base, comment, atom_position = self.read_vasp_poscar()
        xc = 0.5 * supercell_base[0, 0]
        yc = 0.5 * supercell_base[1, 1]

        for i in range(atom_number):
            dx, dy = atom_position[0, i] - xc, atom_position[1, i] - yc
            theta = np.arctan2(dy, dx)
            dz = self.screw_coeff * theta
            atom_position[2, i] += dz
        with open("POSCAR", 'w') as fid:
            fid.write("# Screw Bcc\n")
            fid.write("%12.6f\n" % (self.lattice_constant))
            fid.write("%12.6f %12.6f %12.6f\n" %
                      (supercell_base[0, 0], supercell_base[0, 1], supercell_base[0, 2]))
            fid.write("%12.6f %12.6f %12.6f\n" %
                      (supercell_base[1, 0], supercell_base[1, 1], supercell_base[1, 2]))
            fid.write("%12.6f %12.6f %12.6f\n" %
                      (supercell_base[2, 0], supercell_base[2, 1], supercell_base[2, 2]))
            fid.write("%d\n" % (atom_number))
            fid.write("Cartesian\n")
            for i in range(atom_number):
                fid.write("%12.6f %12.6f %12.6f\n" %
                          (atom_position[0, i], atom_position[1, i], atom_position[2, i]))
            fid.close()
        os.system("cp POSCAR POSCAR.vasp")
        return

    def _intro_single_screw_atoms(self, atoms,
                                  center=None,
                                  sign=None,
                                  orient=[0, 1, 2]):

        supercell_base = atoms.get_cell()
        atom_position = atoms.get_positions()

        if center is None:
            add_disp = [0.5 * (supercell_base[0, 0] + supercell_base[1, 0]),
                        0.5 * (supercell_base[0, 1] + supercell_base[1, 1])]

            xc = 0 + add_disp[0]
            yc = 0 + add_disp[1]
        else:
            xc = center[0]
            yc = center[1]

        for i in range(len(atoms)):
            dx, dy = atom_position[i, 0] - xc, atom_position[i, 1] - yc

            ######  important !!!!! add shift  ######
            if atom_position[i, 0] > xc:
                atom_position[i, 2] += 0.30 * self.burger  # 0.3 or 0.2 for e1 = 1 -1 0
                # 0.2 e1 = -1, -1, 2

            theta = np.arctan2(dy, dx)
            dz = self.screw_coeff * theta
            atom_position[i, 2] -= dz

        atoms.set_positions(atom_position)
        return atoms

    def _intro_single_screw_with_image_atoms(self, atoms):
        supercell_base = atoms.get_cell()
        add_disp = [0.5 * (supercell_base[0, 0] + supercell_base[1, 0]),
                    0.5 * (supercell_base[0, 1] + supercell_base[1, 1])]

        xc1 = 0 + add_disp[0]
        yc1 = 0 + add_disp[1]

        positive_positions = [[xc1], [yc1]]

        first_neigh_list = [[1, 0], [0, 1], [1, 1], [1, -1], [-1, 0], [0, -1], [-1, -1], [-1, 1]]
        second_neigh_list = [[2, 0], [2, 1], [2, 2], [1, 2], [0, 2], [-1, 2], [-2, 2], [-2, 1]]

        for i in range(8):
            second_neigh_list.append(
                [-1 * second_neigh_list[i][0], -1 * second_neigh_list[i][1]])

        for i in range(len(first_neigh_list)):
            positive_positions[0].append(xc1 + first_neigh_list[i][0] * supercell_base[0, 0])
            positive_positions[1].append(yc1 + first_neigh_list[i][1] * supercell_base[1, 1])

        for i in range(len(second_neigh_list)):
            positive_positions[0].append(xc1 + second_neigh_list[i][0] * supercell_base[0, 0])
            positive_positions[1].append(yc1 + second_neigh_list[i][1] * supercell_base[1, 1])

        for i in range(len(positive_positions[0])):
            atoms = self.intro_single_screw_atoms(atoms,
                                                  center=[positive_positions[0][i],
                                                          positive_positions[1][i]],
                                                  sign=1)
        return atoms

    def _intro_single_edge_atoms(self,
                                 atoms,
                                 center=None,
                                 sign=1,
                                 orient=None):

        # default
        tdir = 2  # line direction
        bdir = 0  # burger
        ndir = 1  # glide plane normal

        # hcp
        tdir = 0
        bdir = 1
        ndir = 2

        if orient is not None:
            bdir = orient[0]
            ndir = orient[1]
            tdir = orient[2]

        supercell_base = atoms.get_cell()
        atom_position = atoms.get_positions()

        if center is None:
            add_disp = [0.5 * (supercell_base[bdir, bdir] + supercell_base[ndir, bdir]),
                        0.5 * (supercell_base[bdir, ndir] + supercell_base[ndir, ndir])]

            xc1 = 0 + add_disp[0]
            yc1 = 0 + add_disp[1]
        else:
            xc1 = center[bdir]
            yc1 = center[ndir]

        coeff = self.Edge_coeff

        for i in range(len(atoms)):
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

    def _intro_dipole_edge_atoms(self, atoms):
        supercell_base = atoms.get_cell()
        atom_position = atoms.get_positions()

        add_disp = [0.5 * (supercell_base[0, 0] + supercell_base[1, 0]),
                    0.25 * (supercell_base[0, 1] + supercell_base[1, 1])]

        xc1 = 0 + add_disp[0]
        yc1 = 0 + add_disp[1]

        xc2 = 0.0 * supercell_base[0, 0] + add_disp[0]
        yc2 = 0.5 * supercell_base[1, 1] + add_disp[1]

        coeff = self.Edge_coeff
        for i in range(len(atoms)):
            dx1, dx2 = atom_position[i, 0] - xc1, atom_position[i, 0] - xc2
            dy1, dy2 = atom_position[i, 1] - yc1, atom_position[i, 1] - yc2

            A1 = (1 - self.P) * (dx1**2 + dy1**2)
            A2 = (1 - self.P) * (dx2**2 + dy2**2)

            if A1 != 0:
                ux1 = coeff * (np.arctan2(dy1, dx1) + (dx1 * dy1) / (2. * A1))
                uy1 = (1 - 2 * self.P) / (4 * (1 - self.P)) * \
                    math.log(dx1**2 + dy1**2) + (dx1**2 - dy1**2) / \
                    (4 * (1 - self.P) * (dx1**2 + dy1**2))
                uy1 *= -coeff
            else:
                ux1 = uy1 = 0

            if A2 != 0:
                ux2 = coeff * (np.arctan2(dy2, dx2) + (dx2 * dy2) / (2. * A2))
                uy2 = (1 - 2 * self.P) / (4 * (1 - self.P)) * \
                    math.log(dx2**2 + dy2**2) + (dx2**2 - dy2**2) / \
                    (4 * (1 - self.P) * (dx2**2 + dy2**2))
                uy2 *= -coeff
            else:
                ux2 = uy2 = 0

            atom_position[i, 0] = atom_position[i, 0] + ux1 - ux2
            atom_position[i, 1] = atom_position[i, 1] + uy1 - uy2

        atoms.set_positions(atom_position)
        return atoms

    def _intro_dipole_edge_with_image_atoms(self, atoms):
        supercell_base = atoms.get_cell()

        add_disp = [0.5 * (supercell_base[0, 0] + supercell_base[1, 0]),
                    0.0 * (supercell_base[0, 1] + supercell_base[1, 1])]

        xc1 = 0 + add_disp[0]
        yc1 = 0 + add_disp[1]

        xc2 = 0.0 * supercell_base[0, 0] + add_disp[0]
        yc2 = 0.5 * supercell_base[1, 1] + add_disp[1]

        positive_positions = [[xc1], [yc1]]
        negative_positions = [[xc2], [yc2]]

        #  positive_positions[0].append(xc2)
        #  positive_positions[1].append(yc2)

        first_neigh_list = [[1, 0], [0, 1], [1, 1], [1, -1], [-1, 0], [0, -1], [-1, -1], [-1, 1]]

        second_neigh_list = [[2, 0], [2, 1], [2, 2], [1, 2], [0, 2], [-1, 2], [-2, 2], [-2, 1]]
        for i in range(8):
            second_neigh_list.append(
                [-1 * second_neigh_list[i][0], -1 * second_neigh_list[i][1]])

        for i in range(len(first_neigh_list)):
            positive_positions[0].append(xc1 + first_neigh_list[i][0] * supercell_base[0, 0])
            positive_positions[1].append(yc1 + first_neigh_list[i][1] * supercell_base[1, 1])

            negative_positions[0].append(xc2 + first_neigh_list[i][0] * supercell_base[0, 0])
            negative_positions[1].append(yc2 + first_neigh_list[i][1] * supercell_base[1, 1])

        for i in range(len(second_neigh_list)):
            positive_positions[0].append(xc1 + second_neigh_list[i][0] * supercell_base[0, 0])
            positive_positions[1].append(yc1 + second_neigh_list[i][1] * supercell_base[1, 1])

            negative_positions[0].append(xc2 + second_neigh_list[i][0] * supercell_base[0, 0])
            negative_positions[1].append(yc2 + second_neigh_list[i][1] * supercell_base[1, 1])

        for i in range(len(positive_positions[0])):
            atoms = self.intro_single_edge_atoms(atoms,
                                                 center=[positive_positions[0][i],
                                                         positive_positions[1][i]],
                                                 sign=1)
            atoms = self.intro_single_edge_atoms(atoms,
                                                 center=[negative_positions[0][i],
                                                         negative_positions[1][i]],
                                                 sign=-1)
        return atoms

    def _intro_dipole_screw_with_image_atoms(self, atoms):
        supercell_base = atoms.get_cell()

        add_disp = [0.5 * (supercell_base[0, 0] + supercell_base[1, 0]),
                    0.5 * (supercell_base[0, 1] + supercell_base[1, 1])]

        xc1 = 0.5 * supercell_base[0, 0]
        yc1 = 0.0 * supercell_base[1, 1]

        xc2 = 0.5 * supercell_base[0, 0]
        yc2 = 0.0 * supercell_base[1, 1] + add_disp[1]

        positive_positions = [[xc1], [yc1]]
        negative_positions = [[xc2], [yc2]]

        first_neigh_list = [[1, 0], [0, 1], [1, 1], [1, -1], [-1, 0], [0, -1], [-1, -1], [-1, 1]]

        second_neigh_list = [[2, 0], [2, 1], [2, 2], [1, 2], [0, 2], [-1, 2], [-2, 2], [-2, 1]]
        for i in range(8):
            second_neigh_list.append(
                [-1 * second_neigh_list[i][0], -1 * second_neigh_list[i][1]])

        for i in range(len(first_neigh_list)):
            positive_positions[0].append(xc1 + first_neigh_list[i][0] * supercell_base[0, 0])
            positive_positions[1].append(yc1 + first_neigh_list[i][1] * supercell_base[1, 1])

            negative_positions[0].append(xc2 + first_neigh_list[i][0] * supercell_base[0, 0])
            negative_positions[1].append(yc2 + first_neigh_list[i][1] * supercell_base[1, 1])

        for i in range(len(second_neigh_list)):
            positive_positions[0].append(xc1 + second_neigh_list[i][0] * supercell_base[0, 0])
            positive_positions[1].append(yc1 + second_neigh_list[i][1] * supercell_base[1, 1])

            negative_positions[0].append(xc2 + second_neigh_list[i][0] * supercell_base[0, 0])
            negative_positions[1].append(yc2 + second_neigh_list[i][1] * supercell_base[1, 1])

        #  check the correctness of the image
        #  import matplotlib.pylab as plt
        #  print positive_positions[1]
        #  plt.scatter(np.array(positive_positions[0][:]),
            #  np.array(positive_positions[1][:]))
        #  plt.scatter(np.array(negative_positions[0][:]),
            #  np.array(negative_positions[1][:]),
            #  color='r')

        for i in range(len(positive_positions[0])):
            atoms = self.intro_single_screw_atoms(atoms,
                                                  center=[positive_positions[0][i],
                                                          positive_positions[1][i]],
                                                  sign=1)
            atoms = self.intro_single_screw_atoms(atoms,
                                                  center=[negative_positions[0][i],
                                                          negative_positions[1][i]],
                                                  sign=-1)
        return atoms

    def _core_center_by_atom_index(self, atoms):
        positions = atoms.get_positions()

        a_list = [102, 107, 103, 105, 104]
        b_list = [123, 122, 124, 84, 89]

        a_atoms = []
        b_atoms = []
        for i in range(len(a_list)):
            a_atoms.append(positions[a_list[i]])
            b_atoms.append(positions[b_list[i]])

        aeasy_core = (a_atoms[0][0:2] + a_atoms[1][0:2] + a_atoms[2][0:2]) / 3.0;
        ahard_core = (a_atoms[1][0:2] + a_atoms[2][0:2] + a_atoms[3][0:2]) / 3.0;
        aeasy_core2 = (a_atoms[2][0:2] + a_atoms[3][0:2] + a_atoms[4][0:2]) / 3.0;
        asplit = a_atoms[2][0:2];
        aM_point = 0.5 * (ahard_core + asplit)

        beasy_core = (b_atoms[0][0:2] + b_atoms[1][0:2] + b_atoms[2][0:2]) / 3.0;
        bhard_core = (b_atoms[1][0:2] + b_atoms[2][0:2] + b_atoms[3][0:2]) / 3.0;
        beasy_core2 = (b_atoms[2][0:2] + b_atoms[3][0:2] + b_atoms[4][0:2]) / 3.0;
        bsplit = b_atoms[2][0:2];
        bM_point = 0.5 * (bhard_core + bsplit)

        return ([aeasy_core, ahard_core, asplit, aM_point, aeasy_core2],
                [beasy_core, bhard_core, bsplit, bM_point, beasy_core2])

    def _intro_dipole_screw_atoms_vasp(self, atoms,
                                       lattice=3.307,
                                       move_x=0,
                                       in_tag='easy_easy',
                                       input_s=1.0):

        atoms.wrap(pbc=[1, 1, 1])
        atom_position = atoms.get_positions()

        ([aeasy_core, ahard_core, asplit, aM_point, aeasy_core2],
         [beasy_core, bhard_core, bsplit, bM_point, beasy_core2]) = self._core_center_by_atom_index(atoms)

        ########### easy to hard ################
        tag = in_tag
        s = input_s
        if tag == 'easy_Mid':
            xc1 = (1 - s) * aeasy_core[0] + s * aM_point[0]
            yc1 = (1 - s) * aeasy_core[1] + s * aM_point[1]

            xc2 = (1 - s) * beasy_core[0] + s * bM_point[0]
            yc2 = (1 - s) * beasy_core[1] + s * bM_point[1]

        elif tag == 'easy_hard':
            xc1 = (1 - s) * aeasy_core[0] + s * ahard_core[0]
            yc1 = (1 - s) * aeasy_core[1] + s * ahard_core[1]

            xc2 = (1 - s) * beasy_core[0] + s * bhard_core[0]
            yc2 = (1 - s) * beasy_core[1] + s * bhard_core[1]

        elif tag == 'easy_easy':
            xc1 = (1 - s) * aeasy_core[0] + s * aeasy_core2[0]
            yc1 = (1 - s) * aeasy_core[1] + s * aeasy_core2[1]

            xc2 = (1 - s) * beasy_core[0] + s * beasy_core2[0]
            yc2 = (1 - s) * beasy_core[1] + s * beasy_core2[1]

        elif tag == 'hard_split':
            xc1 = (1 - s) * ahard_core[0] + s * asplit[0]
            yc1 = (1 - s) * ahard_core[1] + s * asplit[1]

            xc2 = (1 - s) * bhard_core[0] + s * bsplit[0]
            yc2 = (1 - s) * bhard_core[1] + s * bsplit[1]

        elif tag == 'easy_split':
            xc1 = (1 - s) * aeasy_core[0] + s * asplit[0]
            yc1 = (1 - s) * aeasy_core[1] + s * asplit[1]

            xc2 = (1 - s) * beasy_core[0] + s * bsplit[0]
            yc2 = (1 - s) * beasy_core[1] + s * bsplit[1]

        print "xc1 = %g, yc1 = %g \n" % (xc1, yc1)
        print "xc2 = %g, yc2 = %g \n" % (xc2, yc2)

        xc1 += move_x
        xc2 += move_x

        #  dispx = np.sqrt(2)/np.sqrt(3) * lattice
        #  dispy = np.sqrt(2)/2 * lattice;

        for i in range(len(atom_position)):
            dx1, dx2 = atom_position[i, 0] - xc1, atom_position[i, 0] - xc2
            dy1, dy2 = atom_position[i, 1] - yc1, atom_position[i, 1] - yc2

            theta1 = np.arctan2(dy1, dx1)
            theta2 = np.arctan2(dy2, dx2)

            dz1 = self.screw_coeff * theta1
            dz2 = self.screw_coeff * theta2

            atom_position[i, 2] = atom_position[i, 2] + dz1 - dz2

        atoms.set_positions(atom_position)
        return atoms

    def _intro_split_core(self, atoms, lattice=None,
                          move_x=None,
                          in_tag=None, input_s=None):
        atoms.wrap(pbc=[1, 1, 1])

        atom_position = atoms.get_positions()

        ([aeasy_core, ahard_core, asplit, aM_point, aeasy_core2],
         [beasy_core, bhard_core, bsplit, bM_point, beasy_core2]) = self.core_center_by_atom_index(atoms)

        alpha_list = np.array([107, 105])
        beta_list = np.array([122, 84])

        ## introduce hard core ##
        s = 0
        xc1 = (1 - s) * ahard_core[0] + s * asplit[0]
        yc1 = (1 - s) * ahard_core[1] + s * asplit[1]

        xc2 = (1 - s) * bhard_core[0] + s * bsplit[0]
        yc2 = (1 - s) * bhard_core[1] + s * bsplit[1]

        for i in range(len(atom_position)):
            dx1, dx2 = atom_position[i, 0] - xc1, atom_position[i, 0] - xc2
            dy1, dy2 = atom_position[i, 1] - yc1, atom_position[i, 1] - yc2

            theta1 = np.arctan2(dy1, dx1)
            theta2 = np.arctan2(dy2, dx2)

            dz1 = self.screw_coeff * theta1
            dz2 = self.screw_coeff * theta2

            atom_position[i, 2] = atom_position[i, 2] + dz1 - dz2

            if ((i == alpha_list).any()):
                atom_position[i, 2] = atom_position[i, 2] + 1. / 6. * self.burger

            elif ((i == beta_list).any()):
                atom_position[i, 2] = atom_position[i, 2] - 1. / 6. * self.burger

        atoms.set_positions(atom_position)
        return atoms

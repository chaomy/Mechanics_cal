#!/usr/bin/env python
# encoding: utf-8
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-07-05 08:12:30
# @Last Modified by:   chaomy
# @Last Modified time: 2017-09-21 11:28:44


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

    def __init__(self, pot=None):
        print pot
        self.pot = pot
        self.set_intro_coeff()
        return

    def set_intro_coeff(self):
        self.pi = 3.141592653589793
        self.burger = np.sqrt(3.) / 2. * self.pot['lattice']
        self.screw_coeff = self.burger / (2. * self.pi)
        self.Edge_coeff = self.burger / (2. * self.pi)
        print self.screw_coeff
        self.P = 0.33
        return

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

    def intro_kink_screw_dislocations(self, atoms,
                                      center1, center2,
                                      H=None, tilpnt=1. / 4., sign=1):
        # find how many layer of atoms #
        #  (layerList, distance) = \
        #  ase.utils.geometry.get_layers(atoms, [0, 0, 1], tolerance=0.01);
        boundary = 'ppp'

        lxid = 0  # glide plane
        lyid = 2  # y
        lzid = 1  # dislocation line

        cell = atoms.get_cell()
        nleft = tilpnt * cell[lzid, lzid]
        nrigh = (1 - tilpnt) * cell[lzid, lzid]

        if H is None:
            H = 0.5 * (center2[lxid] - center1[lxid])

        alpha = 0.1
        xc0, yc0 = center1[0], center1[1]

        # cut atoms (do not need for peierodic boundary condition ) ###
        if not boundary == 'ppp':
            # not using the ppp bounary condition ####
            # we have some other treaments ####
            lowx = 5 * (0.5 * (2.56656))
            higx = cell[0, 0] - lowx

            lowy = 5 * 2.22271
            higy = cell[1, 1] - lowy

            ##########################################################
            # the vol I cut here is
            # (lowx * y + (lowy * (x - 2 * lowx))) * 2 * z
            ##########################################################
            vol = np.linalg.det(cell)
            areaxy = (cell[0, 0] - 2 * lowx) * cell[1, 1]
            volcut = ((lowx * cell[1, 1] + lowy *
                       (cell[0, 0] - 2 * lowx)) * 2 * cell[2, 2])

            print("vol[%f] - volcut[%f] = %f " % (vol, volcut, vol - volcut))
            print("area is [%f]" % (areaxy))  # 29094.182232  A^2

            index_list = []
        cut = None  # only cut y direction ####
        #  shift   = True;
        shift = False

        # loop over atoms to add displacement #
        for atom in atoms:
            pos = atom.position

            if not boundary == 'ppp':
                # cut  10 A along x and y  #
                if cut == "xy":
                    if ((pos[lxid] < lowx) or (pos[lxid] > higx) or
                            (pos[lxid] > higy) or (pos[lxid] < lowy)):
                        index_list.append(atom.index)

                elif cut == "y":
                    if ((pos[lyid] > higy) or (pos[lyid] < lowy)):
                        index_list.append(atom.index)

            # calculate the h #####
            h = H * (1 + math.tanh(alpha * (pos[lzid] - nleft))) - \
                H * (1 + math.tanh(alpha * (pos[lzid] - nrigh)))

            xc = xc0 + h
            yc = yc0

            # calculate the displacement #####
            dx = pos[lxid] - xc
            dy = pos[lyid] - yc

            # add shift if peierodic along x   ######
            if shift is True:
                if pos[lxid] > xc:
                    # 0.3 or 0.2 for e1 = 1 -1 0
                    pos[lzid] += 0.25 * self.burger

            theta = np.arctan2(dy, dx)
            dz = sign * self.screw_coeff * theta
            pos[lzid] -= dz

            # assign the position ###
            atom.position = pos

        if not boundary == 'ppp':
            # delete boundary atoms ####
            del atoms[index_list]
        return atoms

    def intro_single_screw_atoms_xdirection(self,
                                            atoms):
        supercell_base = atoms.get_cell()
        atom_position = atoms.get_positions()

        add_disp = [0.5 * (supercell_base[1, 0] + supercell_base[1, 1]),
                    0.5 * (supercell_base[2, 0] + supercell_base[2, 2])]

        yc = 0 + add_disp[0]
        zc = 0 + add_disp[1]

        for i in range(len(atoms)):
            dy, dz = atom_position[i, 1] - yc, atom_position[i, 2] - zc
            theta = np.arctan2(dy, dz)
            ux = self.screw_coeff * theta
            atom_position[i, 0] += ux

        atoms.set_positions(atom_position)
        return atoms

    def core_center(self):
        # work for 7 x 11
        boxMod = np.array([0.0, 0.0])
        n = 7

        a1 = np.array([23.0991, 8.89085]) + boxMod
        a2 = np.array([24.3824, 11.1136]) + boxMod      # easy core only
        a3 = np.array([25.6657, 8.89085]) + boxMod  # split core
        a4 = np.array([26.949, 11.1136]) + boxMod      # hard core only
        a5 = np.array([28.2322, 8.89085]) + boxMod      # another easy core

        aeasy_core = (a1 + a2 + a3) / 3.0
        ahard_core = (a2 + a3 + a4) / 3.0
        aeasy_core2 = (a3 + a4 + a5) / 3.0
        asplit = a3
        aM_point = 0.5 * (ahard_core + asplit)

        dx = a3 - a1

        inter = 0.5 * (3 * n - 1)

        b1 = a2 + inter * dx
        b2 = a3 + inter * dx
        b3 = a4 + inter * dx
        b4 = a5 + inter * dx
        b5 = b3 + dx

        beasy_core = (b1 + b2 + b3) / 3.0
        bhard_core = (b2 + b3 + b4) / 3.0
        beasy_core2 = (b3 + b4 + b5) / 3.0
        bsplit = b3
        bM_point = 0.5 * (bhard_core + bsplit)

        return ([aeasy_core, ahard_core, asplit, aM_point, aeasy_core2],
                [beasy_core, bhard_core, bsplit, bM_point, beasy_core2])

    ###########################################################################
    # Sun Mar 26 13:01:49 2017
    # intro dipole screw used for LAMMPS
    # x along the [-1 1 0];
    # y is along the dislocation line t
    # z along the [1 1 -2];
    ###########################################################################
    def intro_dipole_screw_atoms_LMP(self,
                                     atoms,
                                     center=None,
                                     lattice=3.307,
                                     in_tag='hard_split',
                                     input_s=None):

        lxid = 0  # glide plane
        lzid = 1  # dislocation line
        lyid = 2  # y

        atoms.wrap(pbc=[1, 1, 1])

        supercell_base = atoms.get_cell()
        atom_position = atoms.get_positions()
        if center is None:
            add_disp = [0.25 * (supercell_base[lxid, lxid] + supercell_base[lyid, lxid]),
                        0.50 * (supercell_base[lxid, lyid] + supercell_base[lyid, lyid])]

            xc1 = 0 + add_disp[0]
            yc1 = 0 + add_disp[1]

            xc2 = 0.5 * supercell_base[lxid, lxid] + add_disp[0]
            yc2 = 0.0 * supercell_base[lyid, lyid] + add_disp[1]
        else:
            xc1 = center[0][0]
            yc1 = center[0][1]

            xc2 = center[1][0]
            yc2 = center[1][1]

        print("center 1 (%g  %g)\n" % (xc1, yc1))
        print("center 2 (%g  %g)\n" % (xc2, yc2))

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
        return

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

        print len(map_list[1])
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

        print len(map_list[1])
        return map_list

    def sort_atoms(self, atoms, map_list):
        new_atoms = Atoms()
        print map_list[0]
        print map_list[1]

        for i in map_list[1]:
            #  print i
            print atoms[i].position
            new_atoms.append(atoms[i])
        new_atoms.set_cell(atoms.get_cell())
        return new_atoms

    def intro_ani_edge_fcc(self):
        self.fcc_edge()

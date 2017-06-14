#!/usr/bin/env python
# encoding: utf-8
#
###################################################################
#
# File Name : ./Intro_vasp.py
#
###################################################################
#
# Purpose :   contains the routine that change the configurations of atoms
#
# Creation Date :
# Last Modified : Sun Apr  9 20:23:17 2017
# Created By    : Chaoming Yang
#
###################################################################

import os
import numpy as np
import math
from ase import Atoms
import cal_add_strain
import cal_intro_iso_dis
import cal_intro_ani_dis


class cubic_cij:
    c11 = None
    c12 = None
    c44 = None


class vasp_change_box(object):
    def __init__(self, pot=None):
        print pot 
        self.pot = pot
        self.screw_coeff = None
        self.Edge_coeff = None
        self.lattice_constant = self.pot['lattice']
        self.add_strain_drv = cal_add_strain.cal_add_strain()
        self.add_iso_dis_drv = cal_intro_iso_dis.cal_intro_iso_dis(self.lattice_constant)
        self.add_ani_dis_drv = cal_intro_ani_dis.cal_intro_iso_dis()
        self.set_intro_coeff()
        return

    def set_intro_coeff(self):
        self.pi = 3.141592653589793
        self.burger = np.sqrt(3.) / 2. * self.lattice_constant
        self.screw_coeff = self.burger / (2. * self.pi)
        self.Edge_coeff = self.burger / (2. * self.pi)
        print self.screw_coeff
        self.P = 0.33
        return

    def set_intro_lattice_constnat(self, lattice_constant):
        self.lattice_constant = lattice_constant
        return

    def volume_conserving_ortho_strain(self, delta):
        self.add_strain_drv._volume_conserving_ortho_strain(delta)
        return

    def add_volumeric_strain(self, atoms, delta):
        atoms = self.add_strain_drv._add_volumeric_strain(atoms, delta)
        return atoms

    def volume_conserving_ortho_strain_atoms(self, delta, atoms):
        atoms = self.add_strain_drv._volume_conserving_ortho_strain_atoms(delta, atoms)
        return atoms

    def volume_conserving_mono_strain_atoms(self, delta, atoms):
        atoms = self.add_strain_drv._volume_conserving_mono_strain_atoms(delta, atoms)
        return atoms

    def volume_conserving_mono_strain(self, delta):
        self.add_strain_drv._volume_conserving_mono_strain(delta)
        return

    def intro_single_screw(self):
        self.add_iso_dis_drv._intro_single_screw()
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
                yc0 * (1 + math.tanh(alpha * (atom_position[i, tdir] - zzone[1])))
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
                                      H=None,
                                      tilpnt=1. / 4., sign=1):
        ### find how many layer of atoms ###
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

        #### cut atoms (do not need for peierodic boundary condition ) ###
        if not boundary == 'ppp':
            #### not using the ppp bounary condition ####
            #### we have some other treaments ####
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
            volcut = ((lowx * cell[1, 1] + lowy * (cell[0, 0] - 2 * lowx)) * 2 * cell[2, 2])

            print ("vol[%f] - volcut[%f] = %f " % (vol, volcut, vol - volcut))
            print ("area is [%f]" % (areaxy))  # 29094.182232  A^2

            index_list = []
        cut = None  # only cut y direction ####
        #  shift   = True;
        shift = False

        #### loop over atoms to add displacement ####
        for atom in atoms:
            pos = atom.position

            if not boundary == 'ppp':
                ### cut  10 A along x and y  ###
                if cut == "xy":
                    if ((pos[lxid] < lowx) or (pos[lxid] > higx) or
                            (pos[lxid] > higy) or (pos[lxid] < lowy)):
                        index_list.append(atom.index)

                elif cut == "y":
                    if ((pos[lyid] > higy) or (pos[lyid] < lowy)):
                        index_list.append(atom.index)

            #### calculate the h #####
            h = H * (1 + math.tanh(alpha * (pos[lzid] - nleft))) - \
                H * (1 + math.tanh(alpha * (pos[lzid] - nrigh)))

            xc = xc0 + h
            yc = yc0

            #### calculate the displacement #####
            dx = pos[lxid] - xc
            dy = pos[lyid] - yc

            ###### add shift if peierodic along x   ######
            if shift is True:
                if pos[lxid] > xc:
                    pos[lzid] += 0.25 * self.burger  # 0.3 or 0.2 for e1 = 1 -1 0

            theta = np.arctan2(dy, dx)
            dz = sign * self.screw_coeff * theta
            pos[lzid] -= dz

            ### assign the position ###
            atom.position = pos

        if not boundary == 'ppp':
            #### delete boundary atoms ####
            del atoms[index_list]
        return atoms

    ################################################
    # add single screw dislocation
    # with shiftment used for (peierodic along x direction)
    # Sun Mar 26 13:03:27 2017
    ################################################
    def intro_single_screw_atoms(self,
                                 atoms,
                                 center=None, sign=None):
        atoms = self.add_iso_dis_drv._intro_single_screw_atoms(atoms,
                                                               center,
                                                               sign)
        return atoms

    def intro_single_screw_with_image_atoms(self,
                                            atoms):
        atoms = self.add_iso_dis_drv._intro_single_screw_with_image_atoms(atoms)
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

    def intro_single_edge_atoms(self,
                                atoms,
                                center=None,
                                sign=1,
                                orient=None):
        atoms = self.add_iso_dis_drv._intro_single_edge_atoms(atoms,
                                                              center,
                                                              sign,
                                                              orient)
        return atoms

    def intro_dipole_edge_atoms(self, atoms):
        atoms = self.add_iso_dis_drv._intro_dipole_edge_atoms(atoms)
        return atoms

    def intro_dipole_edge_with_image_atoms(self, atoms):
        atoms = self.add_iso_dis_drv._intro_dipole_edge_with_image_atoms(atoms)
        return atoms

    def intro_dipole_screw_with_image_atoms(self, atoms):
        atoms = self.add_iso_dis_drv._intro_dipole_screw_with_image_atoms(atoms)
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

    def core_center_by_atom_index(self, atoms):
        ([aeasy_core, ahard_core, asplit, aM_point, aeasy_core2],
         [beasy_core, bhard_core, bsplit, bM_point, beasy_core2]) = \
            self.add_iso_dis_drv._core_center_by_atom_index(atoms)

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

        print ("center 1 (%g  %g)\n" % (xc1, yc1))
        print ("center 2 (%g  %g)\n" % (xc2, yc2))

        for i in range(len(atom_position)):
            dx1, dx2 = atom_position[i, lxid] - xc1, atom_position[i, lxid] - xc2
            dy1, dy2 = atom_position[i, lyid] - yc1, atom_position[i, lyid] - yc2

            theta1 = np.arctan2(dy1, dx1)
            theta2 = np.arctan2(dy2, dx2)

            dz1 = self.screw_coeff * theta1
            dz2 = self.screw_coeff * theta2

            atom_position[i, lzid] = atom_position[i, lzid] + dz1 - dz2

        atoms.set_positions(atom_position)
        return atoms

    ###########################################################################
    # normal way of intro_dipole_screw_atoms used for vasp calculation
    ###########################################################################
    def intro_dipole_screw_atoms(self, 
                                 atoms,
                                 lattice=3.307,
                                 move_x=0,
                                 in_tag='easy_easy',
                                 input_s=1.0):

        atoms = self.add_iso_dis_drv._intro_dipole_screw_atoms_vasp(atoms,
                                                                    lattice,
                                                                    move_x,
                                                                    in_tag,
                                                                    input_s)

        return atoms

    def intro_split_core(self, atoms, lattice=None,
                         move_x=None,
                         in_tag=None, input_s=None):
        atoms = self.add_iso_dis_drv._intro_split_core(atoms,
                                                       lattice,
                                                       move_x,
                                                       in_tag,
                                                       input_s)
        return atoms

    def cut_half_atoms_new(self, atoms, tag="cuty"):

        cell = atoms.get_cell()
        positions = atoms.get_positions()
        scale_positions = atoms.get_scaled_positions()

        if (tag == "cuty"):
            ### cut along y direction ###
            cutIndx = 1

        elif (tag == "cutx"):
            ### cut along x direction ###
            cutIndx = 0

        elif (tag == "cutz"):
            ### cut along z direction ###
            cutIndx = 2

        yy_max = 0.5 * np.max(positions[:, cutIndx])
        y_max = 0.499 * np.max(scale_positions[:, cutIndx])

        index_list = []
        # delele half atoms
        print yy_max
        inv_cell_t = np.linalg.inv(cell.transpose())

        for i in range(len(atoms)):
            atom = atoms[i]
            b = (inv_cell_t * np.mat(atom.position).transpose())[cutIndx]

            if b >= y_max:
                index_list.append(atom.index)

        print len(index_list)
        del atoms[index_list]

        # change the supercell
        cell[cutIndx, :] = cell[cutIndx, :] * 0.5

        atoms.set_cell(cell)
        print len(atoms.get_positions())
        return atoms

    def cut_half_atoms(self, atoms):
        cell = atoms.get_cell()
        positions = atoms.get_positions()
        scale_positions = atoms.get_scaled_positions()

        print np.max(scale_positions[:, 0])
        yy_max = 0.5 * np.max(positions[:, 1])
        y_max = 0.5 * np.max(scale_positions[:, 1])

        index_list = []

        # delele half atoms
        print yy_max
        inv_cell_t = np.linalg.inv(cell.transpose())

        for i in range(len(atoms)):
            atom = atoms[i]
            b = (inv_cell_t * np.mat(atom.position).transpose())[1]

            if b > y_max:
                index_list.append(atom.index)

        print len(index_list)

        del atoms[index_list]

        # change the supercell
        cell[1, :] = cell[1, :] * 0.4999
        atoms.set_cell(cell)
        print len(atoms.get_positions())
        return atoms

    def cut_y_normal_atoms(self, atoms, gp_n=1):

        cell = atoms.get_cell()
        positions = atoms.get_positions()

        # cut along glide plane normal
        print cell

        y_above = 7. / 8. * cell[gp_n, gp_n]  # edge
        y_below = 1. / 8. * cell[gp_n, gp_n]  # edge

        print y_above, y_below

        index_list = []

        for i in range(len(atoms)):
            atom = atoms[i]
            if ((atom.position[gp_n] > y_above) or (atom.position[gp_n] < y_below)):
                index_list.append(atom.index)

        del atoms[index_list]
        return atoms

    def cut_z_normal_top_atoms(self, atoms):
        cell = atoms.get_cell()
        # cut along y
        print cell

        y_above = 7. / 8. * cell[2, 2]  # edge
        y_below = 1. / 8. * cell[2, 2]  # edge

        print y_above, y_below

        index_list = []

        for i in range(len(atoms)):
            atom = atoms[i]
            if ((atom.position[2] > y_above) or (atom.position[2] < y_below)):
                index_list.append(atom.index)

        del atoms[index_list]
        return atoms

    def cut_z_normal_atoms(self, atoms, lattice_constant):
        cell = atoms.get_cell()
        positions = atoms.get_positions()

        index_list = []
        # #####################################################
        # be careful, to delete
        # x_crit =  np.sqrt(3.)/2.5 * lattice constant
        # #####################################################
        z_crit = -0.2 * np.sqrt(3.) / 4. * lattice_constant

        for i in range(len(atoms)):
            atom = atoms[i]
            if (atom.position[2] < z_crit):
                index_list.append(atom.index)

        if index_list is not []:
            print "delete %s atoms" % (len(index_list))
            del atoms[index_list]
        return atoms

    def cut_x_normal_atoms(self,
                           atoms,
                           lattice_constant,
                           bdir=0,
                           ratio=np.sqrt(3) / 2.5):

        index_list = []
        # #####################################################
        # be careful, to delete
        # x_crit =  np.sqrt(3.)/2.5 * lattice constant
        # #####################################################
        #  x_crit = np.sqrt(3.) / 2.5 * lattice_constant
        x_crit = ratio * lattice_constant
        for i in range(len(atoms)):
            atom = atoms[i]
            if (atom.position[bdir] < x_crit):
                index_list.append(atom.index)

        if index_list is not []:
            print "delete %s atoms" % (len(index_list))
            del atoms[index_list]
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
        self.add_ani_dis_drv._fcc_edge()


if __name__ == "__main__":
    N = vasp_change_box()

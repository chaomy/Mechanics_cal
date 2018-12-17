#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-07-05 08:12:30
# @Last Modified by:   chaomy
# @Last Modified time: 2018-12-12 14:24:11

import numpy as np
import math


class cal_cut_cell(object):

    def assign_cubic(self, atoms, opt, symbol, *args):
        if args[0] == []:
            cell = atoms.get_cell()
            lob = np.array([0.0, 0.0, 0.0])
            hib = np.array([cell[0, 0], cell[1, 1], cell[2, 2]])
        elif len(args) == 2:
            lob = args[0]
            hib = args[1]

        if opt in ['in']:
            for i in range(len(atoms)):
                if all(((atoms[i].position - lob) * (atoms[i].position - hib)) <= 0.0):
                    atoms[i].symbol = symbol
        else:
            for i in range(len(atoms)):
                if any(((atoms[i].position - lob) * (atoms[i].position - hib)) > 0.0):
                    atoms[i].symbol = symbol
        return atoms

    def assign_ynormal_fixatoms(self, atoms, gp_n=1, dh=20):
        cell = atoms.get_cell()
        y_above = cell[gp_n, gp_n] - dh
        y_below = dh
        nup = 0
        ndn = 0
        for i in range(len(atoms)):
            atom = atoms[i]
            # atom.symbol = 'W'
            if (atom.position[gp_n] > y_above):
                atom.symbol = 'Mo'
                nup += 1
            if (atom.position[gp_n] < y_below):
                atom.symbol = 'Nb'
                ndn += 1
        print(("up = ", nup, "dn = ", ndn))
        return atoms

    def make_cylinder(self, opt, atoms, *args):
        lob = args[0]
        hib = args[1]

        center = 0.5 * (hib + lob)
        radius = hib - center

        index = []
        Invaa = 1. / (radius[0] * radius[0])
        Invbb = 1. / (radius[1] * radius[1])
        Invcc = 1. / (radius[2] * radius[2])

        if opt in ["in"]:
            for i in range(len(atoms)):
                x = atoms[i].position[0] - center[0]
                y = atoms[i].position[1] - center[1]
                z = atoms[i].position[2] - center[2]
                if (abs(y) < radius[1] and (x * x * Invaa + z * z * Invcc < 1)):
                    index.append(atoms[i].index)
        else:
            for i in range(len(atoms)):
                x = atoms[i].position[0] - center[0]
                y = atoms[i].position[1] - center[1]
                z = atoms[i].position[2] - center[2]
                if (abs(y) >= radius[1] or (x * x * Invaa + z * z * Invcc >= 1)):
                    index.append(atoms[i].index)
        del atoms[index]
        return atoms

    def make_cubic(self, opt, atoms, *args):
        if args[0] == []:
            cell = atoms.get_cell()
            lob = np.array([0.0, 0.0, 0.0])
            hib = np.array([cell[0, 0], cell[1, 1], cell[2, 2]])
        elif len(args) == 2:
            lob = args[0]
            hib = args[1]
        elif len(args) == 1:
            arr = args[0]
            lob = arr[:3]
            hib = arr[3:]
        index = []
        if opt in ["in"]:
            for i in range(len(atoms)):
                if all(((atoms[i].position - lob) * (atoms[i].position - hib)) <= 0.0):
                    index.append(atoms[i].index)
        else:
            for i in range(len(atoms)):
                if any(((atoms[i].position - lob) * (atoms[i].position - hib)) > 0.0):
                    index.append(atoms[i].index)
        del atoms[index]
        return atoms

    def make_sphere(self, atoms):
        cell = atoms.get_cell()
        xc, yc, zc = 0.75 * cell[0, 0], 0.5 * cell[1, 1], 0.5 * cell[2, 2]
        center = np.array([xc, yc, zc])
        rc = 25
        fid = open('center.txt', 'w')
        fid.write("{}  {}  {}  {}".format(xc, yc, zc, rc))
        fid.close()
        cnt = 0
        for i in range(len(atoms)):
            if np.linalg.norm(atoms[i].position - center) < rc:
                if cnt % 2 == 0:
                    atoms[i].symbol = 'Y'
                if cnt % 2 == 1:
                    atoms[i].symbol = 'O'
                cnt += 1
        return atoms

    def cut_y_normal_atoms(self, atoms, gp_n=1):
        cell = atoms.get_cell()
        # cut along glide plane normal
        dh = 20.
        y_above = cell[gp_n, gp_n] - dh
        y_below = dh
        # y_above = 7. / 8. * cell[gp_n, gp_n]  # edge
        # y_below = 1. / 8. * cell[gp_n, gp_n]  # edge
        index = []
        for i in range(len(atoms)):
            atom = atoms[i]
            if ((atom.position[gp_n] > y_above) or (atom.position[gp_n] < y_below)):
                index.append(atom.index)
        del atoms[index]
        return atoms

    def cut_z_normal_top_atoms(self, atoms):
        cell = atoms.get_cell()
        # cut along y
        print(cell)
        y_above = 7. / 8. * cell[2, 2]  # edge
        y_below = 1. / 8. * cell[2, 2]  # edge

        print(y_above, y_below)

        index_list = []

        for i in range(len(atoms)):
            atom = atoms[i]
            if ((atom.position[2] > y_above) or (atom.position[2] < y_below)):
                index_list.append(atom.index)

        del atoms[index_list]
        return atoms

    def cut_z_normal_atoms(self, atoms):
        index_list = []
        # z_crit = -0.2 * np.sqrt(3.) / 4.0 * self.pot['lattice']
        z_crit = -0.2 * np.sqrt(3.) / 0.7 * self.pot['lattice']
        for i in range(len(atoms)):
            atom = atoms[i]
            if (atom.position[2] < z_crit):
                index_list.append(atom.index)
        if index_list is not []:
            print("delete %s atoms" % (len(index_list)))
            del atoms[index_list]
        return atoms

    def cut_x_normal_atoms(self, atoms, bdir=0):
        index_list = []
        # #####################################################
        # be careful, to delete
        # x_crit =  np.sqrt(3.)/2.5 * lattice constant
        # #####################################################
        #  x_crit = np.sqrt(3.) / 2.5 * lattice_constant
        # ratio = np.sqrt(3) / 2. * 1. / 3.    # for bcc Fe
        # x_crit = ratio * self.pot['lattice']
        # x_crit = np.sqrt(3) / 6 * 3.2
        numDis = 1
        x_crit = np.sqrt(3) / 6 * 3.2 * numDis
        for i in range(len(atoms)):
            atom = atoms[i]
            if (atom.position[bdir] < x_crit):
                index_list.append(atom.index)
        if index_list is not []:
            print("delete %s atoms" % (len(index_list)))
            del atoms[index_list]
        return atoms

    def cut_half_atoms_new(self, atoms, tag="cuty"):
        cell = atoms.get_cell()
        positions = atoms.get_positions()
        scale_positions = atoms.get_scaled_positions()

        if (tag in ["cuty"]):
            # cut along y direction #
            cutIndx = 1

        elif (tag in ["cutx"]):
            # cut along x direction #
            cutIndx = 0

        elif (tag in ["cutz"]):
            # cut along z direction #
            cutIndx = 2

        yy_max = 0.5 * np.max(positions[:, cutIndx])
        y_max = 0.499 * np.max(scale_positions[:, cutIndx])

        index_list = []
        # delele half atoms
        inv_cell_t = np.linalg.inv(cell.transpose())

        for i in range(len(atoms)):
            atom = atoms[i]
            b = (inv_cell_t * np.mat(atom.position).transpose())[cutIndx]

            if b >= y_max:
                index_list.append(atom.index)

        del atoms[index_list]

        # change the supercell
        cell[cutIndx, :] = cell[cutIndx, :] * 0.5

        atoms.set_cell(cell)
        return atoms

    def cut_half_atoms(self, atoms):
        cell = atoms.get_cell()
        positions = atoms.get_positions()
        scale_positions = atoms.get_scaled_positions()

        print(np.max(scale_positions[:, 0]))
        yy_max = 0.5 * np.max(positions[:, 1])
        y_max = 0.5 * np.max(scale_positions[:, 1])

        index_list = []

        # delele half atoms
        print(yy_max)
        inv_cell_t = np.linalg.inv(cell.transpose())

        for i in range(len(atoms)):
            atom = atoms[i]
            b = (inv_cell_t * np.mat(atom.position).transpose())[1]

            if b > y_max:
                index_list.append(atom.index)

        print(len(index_list))

        del atoms[index_list]

        # change the supercell
        cell[1, :] = cell[1, :] * 0.4999
        atoms.set_cell(cell)
        print(len(atoms.get_positions()))
        return atoms

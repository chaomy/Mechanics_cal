#!/usr/bin/env python
# encoding: utf-8
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-07-05 08:12:30
# @Last Modified by:   chaomy
# @Last Modified time: 2017-09-21 11:27:54


class cal_intro_iso_dis_image(object):

    def intro_single_screw_with_image_atoms(self, atoms):
        ucell = atoms.get_cell()
        add_disp = [0.5 * (ucell[0, 0] + ucell[1, 0]),
                    0.5 * (ucell[0, 1] + ucell[1, 1])]

        xc1 = 0 + add_disp[0]
        yc1 = 0 + add_disp[1]

        positive_positions = [[xc1], [yc1]]

        first_neigh_list = [[1, 0], [0, 1], [1, 1], [
            1, -1], [-1, 0], [0, -1], [-1, -1], [-1, 1]]
        second_neigh_list = [[2, 0], [2, 1], [2, 2], [
            1, 2], [0, 2], [-1, 2], [-2, 2], [-2, 1]]

        for i in range(8):
            second_neigh_list.append(
                [-1 * second_neigh_list[i][0], -1 * second_neigh_list[i][1]])

        for i in range(len(first_neigh_list)):
            positive_positions[0].append(
                xc1 + first_neigh_list[i][0] * ucell[0, 0])
            positive_positions[1].append(
                yc1 + first_neigh_list[i][1] * ucell[1, 1])

        for i in range(len(second_neigh_list)):
            positive_positions[0].append(
                xc1 + second_neigh_list[i][0] * ucell[0, 0])
            positive_positions[1].append(
                yc1 + second_neigh_list[i][1] * ucell[1, 1])

        for i in range(len(positive_positions[0])):
            atoms = self.intro_single_screw_atoms(atoms,
                                                  center=[positive_positions[0][i],
                                                          positive_positions[1][i]],
                                                  sign=1)
        return atoms

    def intro_dipole_edge_with_image_atoms(self, atoms):
        ucell = atoms.get_cell()
        # add_disp = [0.5 * (ucell[0, 0] + ucell[1, 0]),
        #             0.0 * (ucell[0, 1] + ucell[1, 1])]
        add_disp = [0.25 * (ucell[0, 0] + ucell[1, 0]),
                    0.0 * (ucell[0, 1] + ucell[1, 1])]

        xc1 = 0 + add_disp[0]
        yc1 = 0 + add_disp[1]

        xc2 = 0.0 * ucell[0, 0] + add_disp[0]
        yc2 = 0.5 * ucell[1, 1] + add_disp[1]

        positive_positions = [[xc1], [yc1]]
        negative_positions = [[xc2], [yc2]]

        #  positive_positions[0].append(xc2)
        #  positive_positions[1].append(yc2)

        first_neigh_list = [[1, 0], [0, 1], [1, 1], [
            1, -1], [-1, 0], [0, -1], [-1, -1], [-1, 1]]
        second_neigh_list = [[2, 0], [2, 1], [2, 2], [
            1, 2], [0, 2], [-1, 2], [-2, 2], [-2, 1]]
        for i in range(8):
            second_neigh_list.append(
                [-1 * second_neigh_list[i][0], -1 * second_neigh_list[i][1]])

        for i in range(len(first_neigh_list)):
            positive_positions[0].append(
                xc1 + first_neigh_list[i][0] * ucell[0, 0])
            positive_positions[1].append(
                yc1 + first_neigh_list[i][1] * ucell[1, 1])

            negative_positions[0].append(
                xc2 + first_neigh_list[i][0] * ucell[0, 0])
            negative_positions[1].append(
                yc2 + first_neigh_list[i][1] * ucell[1, 1])

        for i in range(len(second_neigh_list)):
            positive_positions[0].append(
                xc1 + second_neigh_list[i][0] * ucell[0, 0])
            positive_positions[1].append(
                yc1 + second_neigh_list[i][1] * ucell[1, 1])

            negative_positions[0].append(
                xc2 + second_neigh_list[i][0] * ucell[0, 0])
            negative_positions[1].append(
                yc2 + second_neigh_list[i][1] * ucell[1, 1])

        for i in range(len(positive_positions[0])):
            atoms = self.intro_single_edge_atoms(
                atoms,
                center=[positive_positions[0][i],
                        positive_positions[1][i], 0.0],
                sign=1)
            atoms = self.intro_single_edge_atoms(
                atoms,
                center=[negative_positions[0][i],
                        negative_positions[1][i], 0.0],
                sign=-1)
        return atoms

    def intro_dipole_screw_with_image_atoms(self, atoms):
        ucell = atoms.get_cell()

        add_disp = [0.5 * (ucell[0, 0] + ucell[1, 0]),
                    0.5 * (ucell[0, 1] + ucell[1, 1])]

        xc1 = 0.5 * ucell[0, 0]
        yc1 = 0.0 * ucell[1, 1]

        xc2 = 0.5 * ucell[0, 0]
        yc2 = 0.0 * ucell[1, 1] + add_disp[1]

        positive_positions = [[xc1], [yc1]]
        negative_positions = [[xc2], [yc2]]

        first_neigh_list = [[1, 0], [0, 1], [1, 1], [
            1, -1], [-1, 0], [0, -1], [-1, -1], [-1, 1]]

        second_neigh_list = [[2, 0], [2, 1], [2, 2], [
            1, 2], [0, 2], [-1, 2], [-2, 2], [-2, 1]]
        for i in range(8):
            second_neigh_list.append(
                [-1 * second_neigh_list[i][0], -1 * second_neigh_list[i][1]])

        for i in range(len(first_neigh_list)):
            positive_positions[0].append(
                xc1 + first_neigh_list[i][0] * ucell[0, 0])
            positive_positions[1].append(
                yc1 + first_neigh_list[i][1] * ucell[1, 1])

            negative_positions[0].append(
                xc2 + first_neigh_list[i][0] * ucell[0, 0])
            negative_positions[1].append(
                yc2 + first_neigh_list[i][1] * ucell[1, 1])

        for i in range(len(second_neigh_list)):
            positive_positions[0].append(
                xc1 + second_neigh_list[i][0] * ucell[0, 0])
            positive_positions[1].append(
                yc1 + second_neigh_list[i][1] * ucell[1, 1])

            negative_positions[0].append(
                xc2 + second_neigh_list[i][0] * ucell[0, 0])
            negative_positions[1].append(
                yc2 + second_neigh_list[i][1] * ucell[1, 1])

        #  check the correctness of the image
        #  import matplotlib.pylab as plt
        #  print positive_positions[1]
        #  plt.scatter(np.array(positive_positions[0][:]),
            #  np.array(positive_positions[1][:]))
        #  plt.scatter(np.array(negative_positions[0][:]),
            #  np.array(negative_positions[1][:]),
            #  color='r')

        for i in range(len(positive_positions[0])):
            atoms = self.intro_single_screw_atoms(
                atoms,
                center=[positive_positions[0][i],
                        positive_positions[1][i]],
                sign=1)
            atoms = self.intro_single_screw_atoms(
                atoms,
                center=[negative_positions[0][i],
                        negative_positions[1][i]],
                sign=-1)
        return atoms

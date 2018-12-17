# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2018-02-20 14:11:07
# @Last Modified by:   chaomy
# @Last Modified time: 2018-12-16 20:24:56

import ase
import ase.lattice.orthorhombic as otho
import ase.lattice.cubic as cubic
import ase.io
from numpy import sqrt
import numpy as np
import os


class othoHCPFractory(otho.SimpleOrthorhombicFactory):
    bravais_basis = [[0.0, 0.0, 0.0],
                     [0.5, 0.0, 0.5],
                     [0.0, 1. / 2., 1. / 3.],
                     [1. / 2., 1. / 2., 5. / 6.]]


class othoHCPFractoryB(otho.SimpleOrthorhombicFactory):
    bravais_basis = [[0.0, 0.0, 0.0],
                     [0.5, 0.0, 0.5],
                     [1. / 3., 1. / 2., 0.0],
                     [5. / 6., 1. / 2., 1. / 2.]]

othoHCP = othoHCPFractory()
othoHCPB = othoHCPFractoryB()
# unit cell


class D03Factory(cubic.SimpleCubicFactory):
    "A factory for creating Mg3Nd (D03) lattices."
    bravais_basis = [[0, 0, 0],
                     [0, 0.5, 0.5],
                     [0.5, 0, 0.5],
                     [0.5, 0.5, 0],
                     [0.25, 0.25, 0.75],
                     [0.25, 0.75, 0.25],
                     [0.75, 0.25, 0.25],
                     [0.25, 0.25, 0.25],
                     [0.75, 0.75, 0.25],
                     [0.75, 0.25, 0.75],
                     [0.25, 0.75, 0.75],
                     [0.75, 0.75, 0.75],
                     [0.5, 0.0, 0.0],
                     [0.0, 0.5, 0.0],
                     [0.0, 0.0, 0.5],
                     [0.5, 0.5, 0.5]]
    element_basis = (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1)


# type B  x: [-1, 1, 0]  y:[0, 0, 1] z: [1, 1, 0]
class D03FactoryP110B(otho.SimpleOrthorhombicFactory):
    bravais_basis = [[0.0, 0.0, 0.0],
                     [0.5, 0.25, 0.0],
                     [0.5, 0.75, 0.0],
                     [0.5, 0.5, 0.5],
                     [0.0, 0.25, 0.5],
                     [0.0, 0.75, 0.5],
                     [0.0, 0.5, 0.0],
                     [0.5, 0.0, 0.5]]
    element_basis = (0, 0, 0, 0, 0, 0,
                     1, 1)

# type C  x: [1, 1, 1], y[-1  1  0], z [-1 -1  2]
dx = 0.25 / 3.
dz = 1. / 6.


class D03FactoryP211(otho.SimpleOrthorhombicFactory):
    bravais_basis = [[0.0, 0.0, 0.0], [0.25, 0.0, 0.0],
                     [0.5, 0.0, 0.0], [0.75, 0.0, 0.0],

                     [0.0 + dx, 0.5, dz], [0.25 + dx, 0.5, dz],
                     [0.5 + dx, 0.5, dz], [0.75 + dx, 0.5, dz],

                     [0.0 + 2 * dx, 0.0, 2 * dz], [0.25 + 2 * dx, 0.0, 2 * dz],
                     [0.5 + 2 * dx, 0.0, 2 * dz], [0.75 + 2 * dx, 0.0, 2 * dz],

                     [0.0, 0.5, 3 * dz], [0.25, 0.5, 3 * dz],
                     [0.5, 0.5, 3 * dz], [0.75, 0.5, 3 * dz],

                     [0.0 + dx, 0.0, 4 * dz], [0.25 + dx, 0.0, 4 * dz],
                     [0.5 + dx, 0.0, 4 * dz], [0.75 + dx, 0.0, 4 * dz],

                     [0.0 + 2 * dx, 0.5, 5 * dz], [0.25 + 2 * dx, 0.5, 5 * dz],
                     [0.5 + 2 * dx, 0.5, 5 * dz], [0.75 + 2 * dx, 0.5, 5 * dz]]

    element_basis = (1, 0, 0, 0,
                     0, 1, 0, 0,
                     0, 0, 1, 0,
                     1, 0, 0, 0,
                     0, 1, 0, 0,
                     0, 0, 1, 0)


class D03FactoryP211B(otho.SimpleOrthorhombicFactory):
    bravais_basis = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.25],
                     [0.0, 0.0, 0.5], [0.0, 0.0, 0.75],

                     [dz, 0.5, 0.0 + dx], [dz, 0.5, 0.25 + dx],
                     [dz, 0.5, 0.5 + dx], [dz, 0.5, 0.75 + dx],

                     [2 * dz, 0.0, 0.0 + 2 * dx], [2 * dz, 0.0, 0.25 + 2 * dx],
                     [2 * dz, 0.0, 0.5 + 2 * dx], [2 * dz, 0.0, 0.75 + 2 * dx],

                     [3 * dz, 0.5, 0.0], [3 * dz, 0.5, 0.25],
                     [3 * dz, 0.5, 0.5], [3 * dz, 0.5, 0.75],

                     [4 * dz, 0.0, 0.0 + dx], [4 * dz, 0.0, 0.25 + dx],
                     [4 * dz, 0.0, 0.5 + dx], [4 * dz, 0.0, 0.75 + dx],

                     [5 * dz, 0.5, 0.0 + 2 * dx], [5 * dz, 0.5, 0.25 + 2 * dx],
                     [5 * dz, 0.5, 0.5 + 2 * dx], [5 * dz, 0.5, 0.75 + 2 * dx]]

    element_basis = (1, 0, 0, 0,
                     0, 1, 0, 0,
                     0, 0, 1, 0,
                     1, 0, 0, 0,
                     0, 1, 0, 0,
                     0, 0, 1, 0)

Mg3NdCubic = D03Factory()
Mg3Nd = D03FactoryP211()
Mg3NdB = D03FactoryP211B()

# vol = (342.169 - 20.276) = 321.893 * 55.582636 * 449.269099


class md_prec(object):

    # maka screw dislocation iteracts with the precipitates
    def make_screw_prec(self):
        ux = self.pot["ahcp"] * sqrt(3)
        uy = self.pot["chcp"]
        uz = self.pot['ahcp']

        # sz = [80, 80, 2]
        sz = [150, 120, 100]
        atoms = othoHCPB(latticeconstant=(ux, uy, uz),
                         size=sz, symbol=self.pot['element'])

        lata, latc = self.pot["ahcp"], self.pot["chcp"]
        cell = atoms.get_cell()

        lob = np.array([ux * (sz[0] - 75), 0.0 + 55,
                        uz * (1. / 2. * sz[2] - 15)])
        hib = np.array([ux * (sz[0] - 50), uy * sz[1] -
                        56, uz * (1. / 2. * sz[2] + 15)])

        atoms = self.make_cubic("in", atoms, lob, hib)

        atoms2 = self.buildd03B()
        atoms2.set_positions(lob + atoms2.get_positions())
        atoms2 = self.make_cubic("out", atoms2, lob - 1.0, hib + 1.0)

        # atoms = self.build_screw_basal_hcp_atoms_dipole(atoms)
        atoms = self.build_screw_basal_hcp_atoms(atoms)

        atoms = self.cut_y_normal_atoms(atoms, 1)
        atoms = self.cut_y_normal_atoms(atoms, 0)

        # index = []
        # cryticalLen = 0.5 * cell[1, 1]

        # for i in range(len(atoms)):
        #     atom = atoms[i]
        #     if (atom.position[1] > cryticalLen):
        #         index.append(atom.index)
        # del atoms[index]

        # cell[1, 1] = cryticalLen + 30.0

        # atoms.translate([0.0, 15, 0.0])
        # atoms.set_cell(cell)

        # atoms = self.cut_z_normal_atoms(atoms)
        # atoms = self.assign_ynormal_fixatoms(atoms)

        atoms.extend(atoms2)
        # print("volume", cell[0, 0] * (cell[1, 1] - 40.0) * cell[2, 2])
        self.write_lmp_config_data(atoms, "lmp_init.txt")

    def make_dipole_prec(self):
        ux, uy, uz = self.pot['ahcp'], self.pot[
            'chcp'], self.pot['ahcp'] * sqrt(3.)

        sz = [150, 90, 60]  # z = 60 gives around 30 nm  100->55. 194->100
        # sz = [150, 90, 3]

        atoms = othoHCP(latticeconstant=(ux, uy, uz),
                        size=sz, symbol=self.pot['element'])
        lata, latc = self.pot["ahcp"], self.pot["chcp"]
        self.burger = self.pot["lattice"]
        cell = atoms.get_cell()

        # + 14 * uy
        # - 14 * uy
        lob = np.array([ux * (sz[0] - 68), 0.0,
                        uz * (1. / 2. * sz[2] - 6)])
        hib = np.array([ux * (sz[0] - 28), uy * sz[1],
                        uz * (1. / 2. * sz[2] + 6)])
        atoms = self.make_cubic("in", atoms, lob, hib)
        # atoms = self.make_cylinder("in", atoms, lob, hib)
        atoms2 = self.buildd03()

        pos = atoms2.get_positions()
        pos += lob
        atoms2.set_positions(pos)

        # atoms2 = self.make_cylinder("out", atoms2, lob - 1.0, hib + 1.0)
        atoms = self.make_cubic("in", atoms, lob, hib)

        atoms = self.build_edge_basal_hcp_atoms(
            atoms, center=[ux * 40, 0.25 * sz[1] * uy], sign=1)

        atoms = self.build_edge_basal_hcp_atoms(
            atoms, center=[ux * 40, 0.75 * sz[1] * uy], sign=-1)

        # atoms = self.cut_y_normal_atoms(atoms)
        atoms = self.cut_x_normal_atoms(atoms)

        atoms.extend(atoms2)
        cml = "volume = {}".format(
            cell[0, 0] * (cell[1, 1] - 40.0) * cell[2, 2])
        self.write_lmp_config_data(atoms, "lmp_init.txt", cml)

    def make_prec_read_screw(self):
        atoms = ase.io.read("dump.422", format="lammps-dump")
        ux, uy, uz = self.pot['ahcp'] * sqrt(3.), self.pot[
            'chcp'], self.pot['ahcp']

        sz = [84, 90, 104]
        cell = atoms.get_cell()
        # print(cell[2, 2])
        # print(sz[2] * uz)

        lob = np.array([ux * (sz[0] - 46), 0.0 + 14 * uy,
                        uz * (1. / 2. * sz[2] - 12)])

        hib = np.array([ux * (sz[0] - 24), cell[1, 1] - 14 * uy,
                        uz * (1. / 2. * sz[2] + 12)])

        atoms = self.make_cubic("in", atoms, lob, hib)

        atoms2 = self.buildd03B()
        atoms2.set_positions(lob + atoms2.get_positions())
        atoms2 = self.make_cubic("out", atoms2, lob - 1.0, hib + 1.0)

        atoms.extend(atoms2)

        # calculate volume
        atoms = self.assign_ynormal_fixatoms(atoms, 1, 20)
        cml = "volume = {}".format(
            cell[0, 0] * (cell[1, 1] - 20 * 2) * cell[2, 2])
        self.write_lmp_config_data(atoms, "lmp_init.txt", cml)

    def make_prec(self):  # FORMALLY USED # [1-210], [0001], [10-10]
        ux, uy, uz = self.pot['ahcp'], self.pot[
            'chcp'], self.pot['ahcp'] * sqrt(3.)

        sz = [150, 90, 60]  # z = 60 gives around 30 nm  100->55. 194->100
        # sz = [200, 100, 5]  # z = 60 gives around 30 nm  100->55. 194->100

        atoms = othoHCP(latticeconstant=(ux, uy, uz),
                        size=sz, symbol=self.pot['element'])

        lata, latc = self.pot["ahcp"], self.pot["chcp"]
        self.burger = self.pot["lattice"]

        lob = np.array([ux * (sz[0] - 68), 0.0 + 14 * uy,
                        uz * (1. / 2. * sz[2] - 6)])
        hib = np.array([ux * (sz[0] - 28), uy * sz[1] - 14 * uy,
                        uz * (1. / 2. * sz[2] + 6)])

        atoms = self.make_cubic("in", atoms, lob, hib)
        atoms2 = self.buildd03()

        pos = atoms2.get_positions()
        pos += lob
        atoms2.set_positions(pos)
        atoms2 = self.make_cubic("out", atoms2, lob - 1.0, hib + 1.0)

        atoms = self.build_edge_basal_hcp_atoms(
            atoms, center=[ux * 40, 0.5 * sz[1] * uy])

        atoms = self.cut_y_normal_atoms(atoms)
        atoms = self.cut_x_normal_atoms(atoms)

        atoms = self.assign_ynormal_fixatoms(atoms)

        atoms.extend(atoms2)

        cell = atoms.get_cell()
        cml = "volume = {}".format(
            cell[0, 0] * (cell[1, 1] - 40.0) * cell[2, 2])
        self.write_lmp_config_data(atoms, "lmp_init.txt", cml)

    def make_double_prec_test(self):
        ux, uy, uz = self.pot['ahcp'], self.pot[
            'chcp'], self.pot['ahcp'] * sqrt(3.)

        sz = [150, 90, 120]  # z = 60 gives around 30 nm  100->55. 194->100
        lata, latc = self.pot["ahcp"], self.pot["chcp"]

        atoms3 = self.buildd03()
        tmpc = atoms3.get_cell()
        self.write_lmp_config_data(atoms3, "atm3_1.txt")

        # mirrow along z direction
        image = tmpc[2, 2]
        for atom in atoms3:
            atom.position[2] = image - atom.position[2]
        self.write_lmp_config_data(atoms3, "atm3_2.txt")

    def make_double_prec(self):
        ux, uy, uz = self.pot['ahcp'], self.pot[
            'chcp'], self.pot['ahcp'] * sqrt(3.)

        sz = [150, 90, 120]  # z = 60 gives around 30 nm  100->55. 194->100
        atoms = othoHCP(latticeconstant=(ux, uy, uz),
                        size=sz, symbol=self.pot['element'])
        lata, latc = self.pot["ahcp"], self.pot["chcp"]
        self.burger = self.pot["lattice"]
        cell = atoms.get_cell()
        lob = np.array([ux * (sz[0] - 68), 0.0 + 14 * uy,
                        uz * (1. / 4. * sz[2] - 6)])
        hib = np.array([ux * (sz[0] - 28), uy * sz[1] - 14 * uy,
                        uz * (1. / 4. * sz[2] + 6)])
        atoms = self.make_cubic("in", atoms, lob, hib)
        atoms2 = self.buildd03()
        atoms2.translate(lob)
        # pos = atoms2.get_positions()
        # pos += lob
        # atoms2.set_positions(pos)
        atoms2 = self.make_cubic("out", atoms2, lob - 1.0, hib + 1.0)

        lob2 = np.array([ux * (sz[0] - 68), 0.0 + 14 * uy,
                         uz * (3. / 4. * sz[2] - 6)])
        hib2 = np.array([ux * (sz[0] - 28), uy * sz[1] - 14 * uy,
                         uz * (3. / 4. * sz[2] + 6)])
        atoms = self.make_cubic("in", atoms, lob2, hib2)
        atoms3 = self.buildd03()
        tmpc = atoms3.get_cell()

        # mirror along z direction
        image = tmpc[2, 2]
        for atom in atoms3:
            atom.position[2] = image - atom.position[2]
        atoms3.translate(lob2)
        # pos = atoms3.get_positions()
        # pos += lob2
        # atoms3.set_positions(pos)
        atoms3 = self.make_cubic("out", atoms3, lob2 - 1.0, hib2 + 1.0)

        atoms = self.build_edge_basal_hcp_atoms(
            atoms, center=[ux * 40, 0.5 * sz[1] * uy])

        atoms = self.cut_y_normal_atoms(atoms)
        atoms = self.cut_x_normal_atoms(atoms)
        atoms = self.assign_ynormal_fixatoms(atoms)

        atoms.extend(atoms2)
        atoms.extend(atoms3)
        cml = "volume = {}".format(
            cell[0, 0] * (cell[1, 1] - 40.0) * cell[2, 2])
        self.write_lmp_config_data(atoms, "lmp_init.txt", cml)

    def make_r60_prec(self):
        ux, uy, uz = self.pot['ahcp'], self.pot[
            'chcp'], self.pot['ahcp'] * sqrt(3.)

        sz = (150, 90, 60)
        atoms = othoHCP(latticeconstant=(ux, uy, uz),
                        size=sz, symbol=self.pot['element'])
        lata, latc = self.pot["ahcp"], self.pot["chcp"]
        self.burger = self.pot["lattice"]

        cell = atoms.get_cell()
        lob = np.array([ux * 112, 0.0 + 14 * uy, -16 * uz])
        hib = np.array([ux * 132, uy * sz[1] - 14 * uy, 10 * uz])

        atoms.rotate(30, 'y')
        atoms = self.make_cubic("in", atoms, lob, hib)
        atoms.rotate(-30, 'y')

        # atoms = self.intro_single_edge_atoms(
        #     atoms, center=[ux * 30, 0.5 * sz[1] * uy, 15 * uz])
        atoms = self.build_edge_basal_hcp_atoms(
            atoms, center=[ux * 40, 0.5 * sz[1] * uy])

        atoms = self.cut_y_normal_atoms(atoms)
        atoms = self.cut_x_normal_atoms(atoms)
        atoms = self.assign_ynormal_fixatoms(atoms)

        atoms.rotate(30, 'y')
        # make precipitates
        atoms2 = self.buildd03()
        atoms2.set_positions(atoms2.get_positions() + lob)
        atoms2 = self.make_cubic("out", atoms2, lob - 1.0, hib + 1.0)
        atoms.extend(atoms2)

        # self.write_lmp_config_data(atoms, "lmp_rotate.txt")
        atoms.rotate(-30, 'y')
        self.write_lmp_config_data(atoms, "lmp_init.txt")

    def make_only_prec(self):
        ux, uy, uz = self.pot['ahcp'], self.pot[
            'chcp'], self.pot['ahcp'] * sqrt(3.)

        sz = (40, 10, 10)
        atoms = othoHCP(latticeconstant=(ux, uy, uz),
                        size=sz, symbol=self.pot['element'])
        lata, latc = self.pot["ahcp"], self.pot["chcp"]
        self.burger = self.pot["lattice"]

        cell = atoms.get_cell()
        lob = np.array([ux * 15, 0.0, -13 * uz])
        hib = np.array([ux * 23, uy * 10, -8 * uz])

        atoms.rotate(60, 'y')
        atoms = self.make_cubic("in", atoms, lob, hib)

        # make precipitates
        atoms2 = self.buildd03()
        atoms2.set_positions(atoms2.get_positions() + lob)
        atoms2 = self.make_cubic("out", atoms2, lob - 1.0, hib + 1.0)

        atoms.extend(atoms2)
        self.write_lmp_config_data(atoms, "lmp_rotate.txt")
        atoms.rotate(-60, 'y')
        self.write_lmp_config_data(atoms, "lmp_init.txt")

    def buildd03(self):
        la = self.pot['latD03']
        SZ = (30, 110, 30)
        # SZ = (30, 5, 30)
        # type A unit cell
        # atoms = Mg3Nd(latticeconstant=la, size=(1, 1, 1),
        #               symbol=('Mg', 'Nd'))

        # type B  x: [-1, 1, 0]  y:[0, 0, 1] z: [1, 1, 0]
        # atoms = Mg3Nd(latticeconstant=(la * sqrt(2) / 2, la,
        #                                la * sqrt(2) / 2.),
        #               size=(1, 1, 18), symbol=('Mg', 'Nd'))

        # type C  x: [1, 1, 1], y[-1  1  0], z [-1 -1  2]
        atoms = Mg3Nd(latticeconstant=(la * sqrt(3), la * sqrt(2) / 2.,
                                       la * sqrt(6) / 2.), size=SZ, symbol=('Mg', 'Nd'))
        # U = np.mat([[-1, 1, 0], [0, 0, 1], [0.5, 0.5, 0]])
        # Uinv = np.linalg.inv(U)
        # pos = atoms.get_scaled_positions()
        # print np.linalg.det(U)
        self.write_lmp_config_data(atoms, "lmp_d03.txt")
        return atoms

    def buildd03B(self):
        la = self.pot['latD03']
        # type C  x: [1, 1, 1], y[-1  1  0], z [-1 -1  2]
        atoms = Mg3NdB(latticeconstant=(la * sqrt(6) / 2.,
                                        la * sqrt(2) / 2.,
                                        la * sqrt(3)), size=(20, 100, 20), symbol=('Mg', 'Nd'))
        self.write_lmp_config_data(atoms, "lmp_d03.txt")
        return atoms

    def buildHCP(self):
        ux, uy, uz = self.pot['ahcp'], self.pot[
            'chcp'], self.pot['ahcp'] * sqrt(3.)
        sz = [1, 1, 16]
        atoms = othoHCP(latticeconstant=(ux, uy, uz),
                        size=sz, symbol=self.pot['element'])
        cell = atoms.get_cell()
        print(cell[2, 2])
        cell[2, 2] = cell[2, 2] + 20
        atoms.translate([0, 0, 10])
        atoms.set_cell(cell)
        cryticalLen = 0.5 * cell[2, 2]

        fname = "lmp_init{:03d}.txt".format(0)
        self.write_lmp_config_data(atoms, fname)
        os.system("cp {} lmp_init.txt".format(fname))
        os.system("lmp_mpi -i in.min_hcp")

        numP = 20
        for i in range(numP):
            for atom in atoms:
                if (atom.position[2] > cryticalLen):
                    atom.position[1] = atom.position[
                        1] + (1. / numP) * cell[1, 1]
                fname = "lmp_init{:03d}.txt".format(i + 1)
                self.write_lmp_config_data(atoms, fname)
            os.system("cp {} lmp_init.txt".format(fname))
            os.system("lmp_mpi -i in.min_hcp")

    def calculateGSFInterface(self):
        atoms = ase.io.read("dump.rst", format='lammps-dump')
        cell = atoms.get_cell()
        numP = 20
        for i in range(20):
            for atom in atoms:
                if (atom.position[2] > 57.65):
                    atom.position[1] = atom.position[
                        1] + (1. / numP) * cell[1, 1]
            fname = "lmp_init{:03d}.txt".format(i)
            self.write_lmp_config_data(atoms, fname)
            os.system("cp {} lmp_init.txt".format(fname))
            os.system("lmp_mpi -i in.min")

    def buildHCP_D03_Interface(self):
        ux, uy, uz = self.pot['ahcp'], self.pot[
            'chcp'], self.pot['ahcp'] * sqrt(3.)

        sz = [4, 1, 18]
        atoms1 = othoHCP(latticeconstant=(ux, uy, uz),
                         size=sz, symbol='Na')
        cell1 = atoms1.get_cell()
        shift = cell1[2, 2] * 0.5

        idx = []
        for atom in atoms1:
            if atom.position[2] > shift:
                idx.append(atom.index)
        del atoms1[idx]

        la = self.pot['d03']
        SZ = (1, 1, 4)
        atoms2 = Mg3Nd(latticeconstant=(la * sqrt(3), la * sqrt(2) / 2.,
                                        la * sqrt(6) / 2.), size=SZ, symbol=('Mg', 'Nd'))
        cell2 = atoms2.get_cell()
        atoms2.translate([0.0, 0.0, shift])
        atoms1.extend(atoms2)

        atoms1.translate([0.0, 0.0, 0.125 * (cell1[2, 2] - cell2[2, 2])])
        self.write_lmp_config_data(atoms1, "lmp_init.txt")

    def buildd03small(self):
        la = latd03 = self.pot["latD03"]
        # atoms = Mg3NdCubic(latticeconstant=la, size=(1, 1, 1),
        #                    symbol=('Mg', 'Nd')) 
        # x [1, 1, 1],  z[1, 1, 2]
        atoms = Mg3Nd(latticeconstant=(la * sqrt(3), la * sqrt(2) / 2.,
                                       la * sqrt(6) / 2.), size=(2, 1, 2), symbol=('Mg', 'Nd'))
        cell = atoms.get_cell() 
        for atm in atoms: 
            if (atm.position[1] < 0.5 * cell[1, 1]): 
                if (atm.symbol == 'Mg'):
                    atm.symbol = 'Al' 
                else:
                    atm.symbol = 'Y'  
        self.write_lmp_config_data(atoms, "lmp_d03.txt")

    def genGSF2D(self, atoms, org_atoms):
        la = latd03 = self.pot["latD03"]

        unitx = la * sqrt(3)
        unity = la * la * sqrt(2) / 2. 
        unitz = la * sqrt(6) / 2. 
        delta = 1. / 24  

        cell = atoms.get_cell() 
        
        for j in range(25):
            for i in range(25):
                jobid = j * 25 + i 
                fname = "lmp_gsf_{}.txt".format(jobid) 
                self.write_lmp_config_data(atoms, "CONFIGS/" + fname)

                for k in range(len(atoms)):  
                    if (atoms[k].position[1] < 0.5 * cell[1, 1]): 
                        atoms[k].position[0] = org_atoms[k].position[0] + i * delta * unitx 
                        atoms[k].position[2] = org_atoms[k].position[2] + j * delta * unitz

                os.system("cp CONFIGS/{} lmp_init.txt".format(fname)) 
                os.system("lmp_mpi -i in.init")  

    def genGSF1D(self):
        la = latd03 = self.pot["latD03"]

        unitx = la * sqrt(3)
        unity = la * la * sqrt(2) / 2. 
        unitz = la * sqrt(6) / 2. 
        delta = 1. / 24  

        # gen1D
        delta = 1. / 12 
        for i in range(13):
            self.write_lmp_config_data(atoms, "lmp_gsf_{}.txt".format(i))
            for atm in atoms: 
                if (atm.position[1] < 0.5 * cell[1, 1]): 
                    atm.position[0] += delta * unitx 

                    # atm.position[0] += np.sqrt(3) / 2. * burg 
                    # atm.position[2] += 1./2. * burg      

    def buildd03small_apb(self):
        la = latd03 = self.pot["latD03"]
        unitx = la * sqrt(3)
        unity = la * la * sqrt(2) / 2. 
        unitz = la * sqrt(6) / 2. 

        atoms = Mg3Nd(latticeconstant=(la * sqrt(3), la * sqrt(2) / 2.,
                                       la * sqrt(6) / 2.), size=(2, 12, 2), symbol=('Mg', 'Nd')) 

        burg = self.pot["ahcp"]

        cell = atoms.get_cell() 
        cell[1, 1] += 20.0  
        atoms.set_cell(cell)
        atoms.translate([0.0, 10.0, 0.0])

        self.genGSF2D(atoms, atoms.copy()) 
        
        # run 1 D  
        # for i in range(13): 
        #     os.system("cp lmp_gsf_{}.txt lmp_init.txt".format(i))  
        #     os.system("lmp_mpi -i in.init") 

        # run 2 D 

    def cal_thermo(self):
        ux, uy, uz = self.pot['ahcp'], self.pot[
            'chcp'], self.pot['ahcp'] * sqrt(3.)
        sz = (10, 10, 10)
        atoms = othoHCP(latticeconstant=(ux, uy, uz),
                        size=sz, symbol=self.pot['element'])
        lata, latc = self.pot["ahcp"], self.pot["chcp"]
        self.burger = self.pot["lattice"]
        self.write_lmp_config_data(atoms, "thermo2.txt")

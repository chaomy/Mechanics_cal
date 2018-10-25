# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2018-02-20 14:11:07
# @Last Modified by:   chaomy
# @Last Modified time: 2018-08-24 16:07:51

import ase.lattice.orthorhombic as otho
import ase.lattice.cubic as cubic
from numpy import sqrt
import numpy as np


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

# Mg3Nd = D03Factory()
# Mg3Nd = D03FactoryP110B()
Mg3Nd = D03FactoryP211()
Mg3NdB = D03FactoryP211B()

# vol = (342.169 - 20.276) = 321.893 * 55.582636 * 449.269099


class md_prec(object):
    # def make_screw_prec(self):

    def make_screw_prec(self):
        ux = self.pot["ahcp"] * sqrt(3)
        uy = self.pot["chcp"]
        uz = self.pot['ahcp']

        # sz = [80, 80, 100]
        sz = [80, 80, 1]
        atoms = othoHCPB(latticeconstant=(ux, uy, uz),
                         size=sz, symbol=self.pot['element'])

        lata, latc = self.pot["ahcp"], self.pot["chcp"]
        cell = atoms.get_cell()

        lob = np.array([ux * (sz[0] - 30), 0.0 + 55,
                        uz * (1. / 2. * sz[2] - 5)])
        hib = np.array([ux * (sz[0] - 20), uy * sz[1] - 56,
                        uz * (1. / 2. * sz[2] + 5)])

        atoms = self.make_cubic("in", atoms, lob, hib)

        # atoms2 = self.buildd03B()
        # atoms2.set_positions(lob + atoms2.get_positions())
        # atoms2 = self.make_cubic("out", atoms2, lob - 1.0, hib + 1.0)

        atoms = self.build_screw_basal_hcp_atoms(atoms)
        atoms = self.cut_y_normal_atoms(atoms)
        # atoms = self.cut_z_normal_atoms(atoms)
        atoms = self.assign_ynormal_fixatoms(atoms)
        # atoms.extend(atoms2)
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
        # atoms = self.make_cubic("in", atoms, lob, hib)
        atoms = self.make_cylinder("in", atoms, lob, hib)
        atoms2 = self.buildd03()

        pos = atoms2.get_positions()
        pos += lob
        atoms2.set_positions(pos)

        atoms2 = self.make_cylinder("out", atoms2, lob - 1.0, hib + 1.0)

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

    def make_prec(self):  # FORMALLY USED # [1-210], [0001], [10-10]
        ux, uy, uz = self.pot['ahcp'], self.pot[
            'chcp'], self.pot['ahcp'] * sqrt(3.)

        # sz = [150, 90, 60]  # z = 60 gives around 30 nm  100->55. 194->100
        sz = [2000, 100, 5]  # z = 60 gives around 30 nm  100->55. 194->100

        atoms = othoHCP(latticeconstant=(ux, uy, uz),
                        size=sz, symbol=self.pot['element'])
        lata, latc = self.pot["ahcp"], self.pot["chcp"]
        self.burger = self.pot["lattice"]
        cell = atoms.get_cell()

        lob = np.array([ux * (sz[0] - 68), 0.0 + 14 * uy,
                        uz * (1. / 2. * sz[2] - 6)])
        hib = np.array([ux * (sz[0] - 28), uy * sz[1] - 14 * uy,
                        uz * (1. / 2. * sz[2] + 6)])
        # atoms = self.make_cubic("in", atoms, lob, hib)
        # atoms2 = self.buildd03()

        # pos = atoms2.get_positions()
        # pos += lob
        # atoms2.set_positions(pos)
        # atoms2 = self.make_cubic("out", atoms2, lob - 1.0, hib + 1.0)

        # atoms = self.build_edge_basal_hcp_atoms(
        #     atoms, center=[ux * 40, 0.5 * sz[1] * uy])

        # atoms = self.cut_y_normal_atoms(atoms)
        # atoms = self.cut_x_normal_atoms(atoms)
        # atoms = self.assign_ynormal_fixatoms(atoms)

        # atoms.extend(atoms2)

        cml = "volume = {}".format(
            cell[0, 0] * (cell[1, 1] - 40.0) * cell[2, 2])
        self.write_lmp_config_data(atoms, "lmp_init.txt", cml)

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

        pos = atoms2.get_positions()
        pos += lob
        atoms2.set_positions(pos)
        atoms2 = self.make_cubic("out", atoms2, lob - 1.0, hib + 1.0)

        lob = np.array([ux * (sz[0] - 68), 0.0 + 14 * uy,
                        uz * (3. / 4. * sz[2] - 6)])
        hib = np.array([ux * (sz[0] - 28), uy * sz[1] - 14 * uy,
                        uz * (3. / 4. * sz[2] + 6)])
        atoms = self.make_cubic("in", atoms, lob, hib)
        atoms3 = self.buildd03()
        tmpc = atoms3.get_cell()
        atoms3.rotate(180, 'z', center=(0.5 * tmpc[0, 0], 0.5 * tmpc[1, 1], 0))
        pos = atoms3.get_positions()
        pos += lob
        atoms3.set_positions(pos)
        atoms3 = self.make_cubic("out", atoms3, lob - 1.0, hib + 1.0)

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
        la = latd03 = 7.46627803307887
        SZ = (120, 110, 30)
        # unit cell
        # atoms = Mg3Nd(latticeconstant=latd03, size=(1, 1, 1),
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
        la = latd03 = 7.46627803307887
        # type C  x: [1, 1, 1], y[-1  1  0], z [-1 -1  2]
        atoms = Mg3NdB(latticeconstant=(la * sqrt(6) / 2.,
                                        la * sqrt(2) / 2.,
                                        la * sqrt(3)), size=(20, 85, 15), symbol=('Mg', 'Nd'))
        self.write_lmp_config_data(atoms, "lmp_d03.txt")
        return atoms

    def buildHCP(self):
        ux, uy, uz = self.pot['ahcp'], self.pot[
            'chcp'], self.pot['ahcp'] * sqrt(3.)
        sz = [2, 2, 2]
        atoms = othoHCP(latticeconstant=(ux, uy, uz),
                        size=sz, symbol=self.pot['element'])
        self.write_lmp_config_data(atoms, "lmp_HCP.txt")

    def buildd03small(self):
        la = latd03 = 7.46627803307887
        atoms = Mg3Nd(latticeconstant=(la * sqrt(3),
                                       la * sqrt(2) / 2.,
                                       la * sqrt(6) / 2.), size=(2, 2, 2), symbol=('Mg', 'Nd'))
        self.write_lmp_config_data(atoms, "lmp_d03.txt")

    def cal_thermo(self):
        ux, uy, uz = self.pot['ahcp'], self.pot[
            'chcp'], self.pot['ahcp'] * sqrt(3.)
        sz = (10, 10, 10)
        atoms = othoHCP(latticeconstant=(ux, uy, uz),
                        size=sz, symbol=self.pot['element'])

        lata, latc = self.pot["ahcp"], self.pot["chcp"]
        self.burger = self.pot["lattice"]
        self.write_lmp_config_data(atoms, "thermo2.txt")

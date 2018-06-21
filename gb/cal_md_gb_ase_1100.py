# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-12-03 11:07:29
# @Last Modified by:   chaomy
# @Last Modified time: 2018-06-20 11:41:17

import ase.lattice.orthorhombic as otho
import ase.io
import os
import numpy as np
import atomman as am
from utils import stroh_solve
from ase import Atoms
from ase.lattice.orthorhombic import SimpleOrthorhombicFactory
from numpy import sqrt, deg2rad, floor, cos, sin


class othoHCPFractory(otho.SimpleOrthorhombicFactory):
    bravais_basis = [[0.0, 0.0, 0.0],
                     [0.5, 0.0, 0.5],
                     [0.0, 0.5, 1. / 3.],
                     [0.5, 0.5, 5. / 6.]]
othoHCP = othoHCPFractory()


class othoHCPFractoryB(otho.SimpleOrthorhombicFactory):
    bravais_basis = [[0.0, 0.0, 0.0],
                     [0.5, 0.0, 0.5],
                     [0.0, 0.5, 2. / 3.],
                     [0.5, 0.5, 1. / 6.]]
othoHCPB = othoHCPFractoryB()


class md_gb_ase_1100(object):

    def auto(self):
        self.build_hcp_ase_1100_small()
        # minimie
        self.make_perf()
        # minimie
        self.intro_edge_dipole()

    def loop_init_1100(self):
        self.find_angles_1100()
        cn = 0
        for e in self.ag[:]:
            mdir = "1100_{:.2f}".format(e[0])
            self.mymkdir(mdir)
            self.write_1100_small(e)
            # self.write_1100_large(e)
            os.system("cp lmp_init.txt pos_{:02d}".format(cn))
            os.system("mv lmp_init.txt {}".format(mdir))
            cn += 1

    def write_1100_small(self, ag):
        ux = self.pot['ahcp']
        uy = self.pot['chcp']
        uz = self.pot['ahcp'] * sqrt(3.)

        VACUMM = 10.0

        # angle, length, i, j
        atoms = othoHCP(latticeconstant=(ux, uy, uz), size=(
            130, 130, 2), symbol=self.pot['element'])

        atoms.rotate(ag[0], 'z')
        cell = atoms.get_cell()
        cell[0, 0], cell[1, 1] = ag[1], 160
        atoms.translate(
            np.array([cell[0, 0] - floor(15 * cos(deg2rad(ag[0]))) * ux,
                      -floor(47 * cos(deg2rad(ag[0]))) * uy, 0]))

        lob = np.array([0.0, 0.0, 0.0])
        hib = np.array([cell[0, 0], 0.5 * cell[1, 1], cell[2, 2]])
        atoms = self.make_cubic('out', atoms, lob, hib)

        # the other grain
        atoms2 = othoHCPB(latticeconstant=(ux, uy, uz), size=(
            140, 140, 2), symbol='Al')     # for 1100 72.877
        # atoms2 = othoHCP(latticeconstant=(ux, uy, uz), size=(
        #     80, 80, 3), symbol='Nb')   # for 1100 58.361

        lob = np.array([0.0, 0.5 * cell[1, 1], 0.0])
        hib = np.array([cell[0, 0], cell[1, 1] - 0.2, cell[2, 2]])
        # hib = np.array([cell[0, 0], cell[1, 1] - VACUMM, cell[2, 2]])

        atoms2.rotate(-ag[0], 'z')
        # atoms2.translate(np.array([-10 * ux, cell[1, 1], 0]))  # for
        # 72.877
        atoms2.translate(np.array(
            [floor(60 * cos(deg2rad(ag[0]))) * -ux,
             cell[1, 1] - floor(12 * cos(deg2rad(ag[0]))) * uy, 0]))  # for 72.877

        atoms2 = self.make_cubic('out', atoms2, lob, hib)
        atoms2.translate(np.array([0.0, 0.2, 0.25 * uz]))
        atoms.extend(atoms2)

        # assign low grain
        assign_gb = 1
        if assign_gb == 1:
            for atom in atoms:
                if atom.position[1] <= 20:
                    atom.symbol = 'W'
                if atom.position[1] >= cell[1, 1] - 20:
                    atom.symbol = 'Mo'

        vacumm = 1
        if vacumm == 1:
            cell[1, 1] += 40.0
        atoms.set_cell(cell)
        if vacumm == 1:
            atoms.translate(np.array([0.0, 20.0, 0.0]))
        self.write_lmp_config_data(atoms, "lmp_init.txt")

    def write_1100_large(self, ag):
        ux = self.pot['ahcp']
        uy = self.pot['chcp']
        uz = self.pot['ahcp'] * sqrt(3.)

        # angle, length, i, j
        atoms = othoHCP(latticeconstant=(ux, uy, uz), size=(
            400, 400, 1), symbol=self.pot['element'])

        atoms.rotate(ag[0], 'z')
        cell = atoms.get_cell()
        cell[0, 0], cell[1, 1] = ag[1], 280
        atoms.translate(
            np.array([cell[0, 0] - floor(15 * cos(deg2rad(ag[0]))) * ux,
                      -floor(47 * cos(deg2rad(ag[0]))) * uy, 0]))

        lob = np.array([0.0, 0.0, 0.0])
        hib = np.array([cell[0, 0], 0.5 * cell[1, 1], cell[2, 2]])
        atoms = self.make_cubic('out', atoms, lob, hib)

        # the other grain
        atoms2 = othoHCPB(latticeconstant=(ux, uy, uz), size=(
            400, 400, 1), symbol='Al')     # for 1100 72.877

        # atoms2 = othoHCP(latticeconstant=(ux, uy, uz), size=(
        #     80, 80, 3), symbol='Nb')   # for 1100 58.361

        lob = np.array([0.0, 0.5 * cell[1, 1], 0.0])
        hib = np.array([cell[0, 0], cell[1, 1] - 0.2, cell[2, 2]])
        # hib = np.array([cell[0, 0], cell[1, 1] - VACUMM, cell[2, 2]])

        atoms2.rotate(-ag[0], 'z')
        # atoms2.translate(np.array([-10 * ux, cell[1, 1], 0]))  # for
        # 72.877
        atoms2.translate(np.array(
            [floor(60 * cos(deg2rad(ag[0]))) * -ux,
             cell[1, 1] - floor(12 * cos(deg2rad(ag[0]))) * uy, 0]))  # for 72.877

        atoms2 = self.make_cubic('out', atoms2, lob, hib)
        atoms2.translate(np.array([0.0, 0.2, 0.25 * uz]))
        atoms.extend(atoms2)

        # lob = np.array([0.0, 40, 0.0])
        # hib = np.array([cell[0, 0], cell[1, 1] - 40, cell[2, 2]])
        # atoms = self.assign_cubic(atoms, 'out', 'W', lob, hib)

        # assign gb region
        assign_gb = 1
        if assign_gb == 1:
            for atom in atoms:
                if atom.position[1] <= 17:
                    atom.symbol = 'W'
                if atom.position[1] >= cell[1, 1] - 17:
                    atom.symbol = 'Mo'

        # to add vacancy optional
        cell[1, 1] += 30
        atoms.translate(np.array([0.0, 15, 0.0]))
        atoms.set_cell(cell)
        self.write_lmp_config_data(atoms, "lmp_init.txt")

    def build_hcp_ase_1100_small(self):  # to examine the GB structures
        self.find_angles_1100(il=[[], [1]], jl=[2])   # 72.877    ABAB --
        self.write_1100_small(self.ag[0])

    def build_hcp_ase_1100(self):
        self.find_angles_1100(il=[[1], [1]], jl=[1])    # 58.361
        self.write_1100_large(self.ag[0])

    def build_hcp_ase_1100_3ABA(self):
        self.find_angles_1100(il=[[1], [1]], jl=[1])    # 58.361
        ux = self.pot['ahcp']
        uy = self.pot['chcp']
        uz = self.pot['ahcp'] * sqrt(3.)

        atoms = othoHCP(latticeconstant=(ux, uy, uz), size=(
            260, 160, 1), symbol=self.pot['element'])
        cell = atoms.get_cell()
        cell[0, 0], cell[1, 1] = 8 * self.ag[0][1], 140 * \
            uy * np.cos(self.ag[0][0] / np.pi)

        atoms.rotate(self.ag[0][0], 'z')
        atoms.translate(np.array([cell[0, 0], -33 * uy, 0]))

        lob = np.array([0.0, 0.75 * cell[1, 1], 0.0])  # double twins
        hib = np.array([cell[0, 0], 0.90 * cell[1, 1], cell[2, 2]])
        atoms = self.make_cubic('in', atoms, lob, hib)

        # the other grain
        atoms2 = othoHCPB(latticeconstant=(ux, uy, uz), size=(
            200, 200, 1), symbol='Ta')

        lob = np.array([0.0, 0.75 * cell[1, 1], 0.0])
        hib = np.array([cell[0, 0], 0.9 * cell[1, 1], cell[2, 2]])

        atoms2.rotate(-self.ag[0][0], 'z')
        atoms2.translate(np.array([-78 * ux, cell[1, 1], 0]))
        atoms2 = self.make_cubic('out', atoms2, lob, hib)
        atoms.extend(atoms2)

        lob = np.array([0.0, 0.0, 0.0])  # double twins
        hib = np.array([cell[0, 0], cell[1, 1], cell[2, 2]])
        atoms = self.make_cubic('out', atoms, lob, hib)
        atoms.set_cell(cell)
        # atoms = atoms.repeat((2, 1, 1))
        self.write_lmp_config_data(atoms, "lmp_init.txt")

    def make_gb_cut(self):
        atoms = ase.io.read("dump/dump.00000", format='lammps-dump')
        cell = atoms.get_cell()
        usxpex = 8.0
        idx = []
        miny = cell[1, 1]
        maxy = 0.0
        ix = np.floor(cell[0, 0] / 10.0)
        for atom in atoms:
            if (atom.position[1] >= 0.5 * cell[1, 1] - 0.5 * usxpex) and (atom.position[1] <= 0.5 * cell[1, 1] + 0.5 * usxpex):
                idx.append(atom.index)
        del atoms[idx]
        print("num of images along x ", ix)
        print("atoms to be predict", len(idx) / ix)
        cnt = 0
        for atom in atoms:
            if (atom.position[1] > 0.5 * cell[1, 1]):
                miny = min(miny, atom.position[1])
                maxy = max(maxy, atom.position[1])
                cnt += 1
        cml = "{:.5f} {:.5f} {:.5f} {:.5f} {} {} {}".format(
            0.5 * cell[1, 1] - 0.5 * usxpex,
            usxpex, miny, maxy, len(atoms) - cnt, cnt, ix)
        self.write_lmp_config_data(atoms, "gb.txt", cml)

    def make_repeat(self):
        atoms = ase.io.read("dump.restart", format='lammps-dump')
        cell = atoms.get_cell()

        # add thermoexpansion
        cell[0, 0] *= self.pot['ahcp200'] / self.pot['ahcp']
        cell[1, 1] *= self.pot['chcp200'] / self.pot['chcp']
        cell[2, 2] *= self.pot['ahcp200'] / self.pot['ahcp']
        atoms.set_cell(cell, scale_atoms=True)

        # width = 300
        rx = int(np.around(180. / cell[0, 0]))
        rz = int(np.around(180. / cell[2, 2]))
        atoms = atoms.repeat([rx, 1, rz])

        vol = np.linalg.det(atoms.get_cell())
        cml = "vol = {}".format(vol * (cell[1, 1] - 30.0) / cell[1, 1])
        self.write_lmp_config_data(atoms, "large.txt", cml)

    def make_perf_HCP(self):
        ux = self.pot['ahcp']
        uy = self.pot['chcp']
        uz = self.pot['ahcp'] * sqrt(3.)
        atoms = othoHCP(latticeconstant=(ux, uy, uz), size=(
            7, 7, 7), symbol=self.pot['element'])
        self.write_lmp_config_data(atoms, "pos.txt")

    def intro_edge_dipole(self):
        self.find_angles_1100(il=[[1], [1]], jl=[1])    # 58.361
        ux = self.pot['ahcp']
        uy = self.pot['chcp']
        uz = self.pot['ahcp'] * sqrt(3.)
        atoms = ase.io.read("dump_twin/dump.00002", format='lammps-dump')

        atoms.rotate(self.ag[0][0], 'z')  # rotate
        atoms.translate(np.array([80 * ux, 0.0, 0]))

        axes = np.array([[1, 0, 0], [0, 0, 1], [0, -1, 0]])
        burgers = self.pot['lattice'] * np.array([1, 0, 0])
        c = am.ElasticConstants()
        c.hexagonal(C11=self.pot['C11'], C12=self.pot['C12'], C33=self.pot[
                    'C33'], C13=self.pot['C13'], C44=self.pot['C44'])
        stroh = stroh_solve.Stroh(c, burgers, axes=axes)
        pos = atoms.get_positions()
        sh1 = np.ones(pos.shape) * \
            np.array([ux + 0.01, 80 * uy + 0.01, 0.0])
        sh2 = np.ones(pos.shape) * \
            np.array([ux + 0.01, 100 * uy + 0.01, 0.0])
        d1 = stroh.displacement(pos - sh1)
        d2 = stroh.displacement(pos - sh2)
        atoms.set_positions(pos + np.real(d1) - np.real(d2))

        atoms.translate(np.array([-80 * ux, 0.0, 0]))
        atoms.rotate(-self.ag[0][0], 'z')  # rotate back

        idx = []
        cy = 699
        for i in range(len(atoms)):
            if atoms[i].position[1] >= cy:
                idx.append(atoms[i].index)
        del atoms[idx]
        self.write_lmp_config_data(atoms, "lmp_init.txt")

    def intro_edge(self):
        self.find_angles_1100(il=[[1], [1]], jl=[1])    # 58.361
        ux = self.pot['ahcp']
        uy = self.pot['chcp']
        uz = self.pot['ahcp'] * sqrt(3.)
        atoms = ase.io.read("dump/dump.00017", format='lammps-dump')

        atoms.rotate(self.ag[0][0], 'z')  # rotate
        atoms.translate(np.array([80 * ux, 0.0, 0]))
        # # introduce dislocation
        axes = np.array([[1, 0, 0], [0, 0, 1], [0, -1, 0]])
        burgers = self.pot['lattice'] * np.array([1, 0, 0])
        c = am.ElasticConstants()
        c.hexagonal(C11=self.pot['C11'], C12=self.pot['C12'],
                    C33=self.pot['C33'], C13=self.pot['C13'], C44=self.pot['C44'])
        stroh = stroh_solve.Stroh(c, burgers, axes=axes)
        pos = atoms.get_positions()
        cx = 60 * ux + 0.01
        cy = 60 * uy + 0.01
        shift = np.ones(pos.shape) * np.array([cx, cy, 0.0])
        disp = stroh.displacement(pos - shift)
        atoms.set_positions(pos + np.real(disp))

        atoms.translate(np.array([-80 * ux, 0.0, 0]))
        atoms.rotate(-self.ag[0][0], 'z')  # rotate back

        # try with non periodict boundary conditions
        cell = atoms.get_cell()
        cell[0, 0] += 30
        cell[1, 1] += 30
        atoms.translate(np.array([15, 15, 0.0]))
        atoms.set_cell(cell)
        # ase.io.write("pos02", images=atoms, format='cfg')
        self.write_lmp_config_data(atoms, "pos02")

    def build_hcp_ase_1100_with_edge(self):
        self.find_angles_1100(il=[[1], [1]], jl=[1])    # 58.361
        ux = self.pot['ahcp']
        uy = self.pot['chcp']
        uz = self.pot['ahcp'] * sqrt(3.)

        print(self.ag)
        atoms = othoHCP(latticeconstant=(ux, uy, uz),
                        size=(80, 80, 2), symbol=self.pot['element'])

        # introduce dislocation
        axes = np.array([[1, 0, 0], [0, 0, 1], [0, -1, 0]])
        burgers = self.pot['lattice'] * np.array([1, 0, 0])
        c = am.ElasticConstants()
        c.hexagonal(C11=self.pot['C11'],
                    C12=self.pot['C12'],
                    C33=self.pot['C33'],
                    C13=self.pot['C13'],
                    C44=self.pot['C44'])

        stroh = stroh_solve.Stroh(c, burgers, axes=axes)
        pos = atoms.get_positions()

        cell = atoms.get_cell()
        cell[0, 0], cell[1, 1] = 8 * self.ag[0][1], 300

        cx = 30 * ux + 0.01
        cy = 30 * uy + 0.01

        shift = np.ones(pos.shape) * np.array([cx, cy, 0.0])
        disp = stroh.displacement(pos - shift)
        atoms.set_positions(pos + np.real(disp))

        atoms.rotate(self.ag[0][0], 'z')
        atoms.translate(np.array([cell[0, 0], -100, 0]))
        lob = np.array([0.0, 10, 0.0])
        hib = np.array([cell[0, 0], 0.5 * cell[1, 1] + 0.3, cell[2, 2]])
        atoms = self.make_cubic('out', atoms, lob, hib)
        atoms2 = othoHCP(latticeconstant=(ux, uy, uz),
                         size=(80, 80, 2), symbol=self.pot['element'])

        lob = np.array([0.0, 0.5 * cell[1, 1], 0.0])
        hib = np.array([cell[0, 0], cell[1, 1] - 10, cell[2, 2]])
        atoms2.rotate(-self.ag[0][0], 'z')
        atoms2.translate(np.array([-90, cell[1, 1], 0]))  # for 72.877
        atoms2 = self.make_cubic('out', atoms2, lob, hib)
        atoms.extend(atoms2)

        # try with non periodict boundary conditions
        atoms.set_cell(cell)
        self.write_lmp_config_data(atoms, "lmp_init.txt")

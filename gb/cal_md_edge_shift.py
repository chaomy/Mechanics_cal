# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2018-05-14 12:59:00
# @Last Modified by:   chaomy
# @Last Modified time: 2018-05-18 01:39:24

import ase.io
import numpy as np
import atomman as am
from utils import stroh_solve
import glob


class md_edge_shift(object):

    def shift_make_edge(self):
        self.find_angles_1100(il=[[1], [1]], jl=[1])    # 58.361
        ux = self.pot['ahcp']
        uy = self.pot['chcp']
        uz = self.pot['ahcp'] * np.sqrt(3.)
        atoms = ase.io.read("dump_init/dump.00549", format='lammps-dump')
        print(len(atoms))

        atoms.rotate(self.ag[0][0], 'z')  # rotate
        atoms.translate(np.array([83.5 * ux, -80 * uy, 0]))

        # assign shift region
        cell = atoms.get_cell()
        lob = np.array([15 * ux, 0.0, 0.0])
        hib = np.array([29 * ux, 0.5 * uy, cell[2, 2]])

        for i in range(len(atoms)):
            if all(((atoms[i].position - lob) * (atoms[i].position - hib)) <= 0.0):
                atoms[i].symbol = 'Mg'

        lob2 = np.array([15 * ux, 0.5 * uy, 0.0])
        hib2 = np.array([29 * ux, 1.0 * uy, cell[2, 2]])
        for i in range(len(atoms)):
            if all(((atoms[i].position - lob2) * (atoms[i].position - hib2)) <= 0.0):
                atoms[i].symbol = 'Cr'

        # exert shift on 'upper layer'
        for i in range(len(atoms)):
            if atoms[i].symbol == 'Cr':
                atoms[i].position += np.array([0.3 * ux, 0.0, 0.0])

        atoms.translate(np.array([-83.5 * ux, 80 * uy, 0]))
        atoms.rotate(-self.ag[0][0], 'z')  # rotate back
        self.write_lmp_config_data(atoms, "lmp_init.txt")

    def shift_make_edge_cnt(self):
        self.find_angles_1100(il=[[1], [1]], jl=[1])    # 58.361
        ux = self.pot['ahcp']
        uy = self.pot['chcp']
        uz = self.pot['ahcp'] * np.sqrt(3.)
        files = glob.glob("dump/dump.*")
        atoms = ase.io.read(files[-1], format='lammps-dump')

        atoms.rotate(self.ag[0][0], 'z')  # rotate
        atoms.translate(np.array([83.5 * ux, -80 * uy, 0]))
        # exert shift on 'upper layer'
        c1, c2 = 0, 0
        allsymbos = atoms.get_chemical_symbols()
        symbols = np.unique(allsymbos)
        print(symbols[1], symbols[0])
        for i in range(len(atoms)):
            if atoms[i].symbol == symbols[4]:
                atoms[i].position += np.array([0.1 * ux, 0.0, 0.0])
                c1 += 1
            elif atoms[i].symbol == symbols[2]:
                c2 += 1
        atoms.translate(np.array([-83.5 * ux, 80 * uy, 0]))
        atoms.rotate(-self.ag[0][0], 'z')  # rotate back
        self.write_lmp_config_data(atoms, "lmp_init02.txt")

    def duplicate_make_edge(self):
        atoms = ase.io.read(glob.glob("dump_init/dump.*")
                            [-1], format='lammps-dump')
        atoms = atoms.repeat((1, 1, 30))
        self.write_lmp_config_data(atoms, "lmp_init01.txt")

    def formular_make_edge_large(self):
        self.find_angles_1100(il=[[1], [1]], jl=[1])    # 58.361
        ux = self.pot['ahcp']
        uy = self.pot['chcp']
        uz = self.pot['ahcp'] * np.sqrt(3.)

        atoms = ase.io.read(glob.glob("dump_init_v2/dump.*")
                            [-1], format='lammps-dump')

        atoms.rotate(180 - self.ag[0][0], 'z')  # rotate
        atoms.translate(np.array([80 * ux, 0, 0]))
        axes = np.array([[1, 0, 0], [0, 0, 1], [0, -1, 0]])
        burgers = self.pot['lattice'] * np.array([1, 0, 0])

        c = am.ElasticConstants()
        c.hexagonal(C11=self.pot['C11'], C12=self.pot['C12'], C33=self.pot[
                    'C33'], C13=self.pot['C13'], C44=self.pot['C44'])

        stroh = stroh_solve.Stroh(c, burgers, axes=axes)
        pos = atoms.get_positions()
        cell = atoms.get_cell()
        sh1 = np.ones(pos.shape) * \
            np.array([-30 * ux + 0.1, 10 * uy, cell[2, 2]])
        sh2 = np.ones(pos.shape) * \
            np.array([-30 * ux + 0.1, 30 * uy, cell[2, 2]])

        d1 = stroh.displacement(pos - sh1)
        d2 = stroh.displacement(pos - sh2)
        atoms.set_positions(pos + np.real(d1) - np.real(d2))

        atoms.translate(np.array([-80 * ux, 0.0, 0]))
        atoms.rotate(self.ag[0][0] - 180, 'z')  # rotate back
        self.write_lmp_config_data(atoms, "lmp_init02.txt")

    def formular_make_edge(self):
        self.find_angles_1100(il=[[1], [1]], jl=[1])    # 58.361
        ux = self.pot['ahcp']
        uy = self.pot['chcp']
        uz = self.pot['ahcp'] * np.sqrt(3.)

        atoms = ase.io.read(glob.glob("dump_init/dump.*")
                            [-1], format='lammps-dump')
        print(len(atoms))
        cell = atoms.get_cell()

        atoms.rotate(self.ag[0][0], 'z')  # rotate
        atoms.translate(np.array([83.5 * ux, -80 * uy, 0]))

        axes = np.array([[1, 0, 0], [0, 0, 1], [0, -1, 0]])
        burgers = self.pot['lattice'] * np.array([1, 0, 0])
        c = am.ElasticConstants()
        c.hexagonal(C11=self.pot['C11'], C12=self.pot['C12'], C33=self.pot[
                    'C33'], C13=self.pot['C13'], C44=self.pot['C44'])
        stroh = stroh_solve.Stroh(c, burgers, axes=axes)
        pos = atoms.get_positions()

        # sh1 = np.ones(pos.shape) * \
        #     np.array([30 * ux + 0.1, 0.5 * uy, cell[2, 2]])

        sh1 = np.ones(pos.shape) * \
            np.array([34.5 * ux + 0.1, 0.5 * uy, cell[2, 2]])

        sh2 = np.ones(pos.shape) * \
            np.array([-20 * ux + 0.1, 0.5 * uy, cell[2, 2]])

        d1 = stroh.displacement(pos - sh1)
        d2 = stroh.displacement(pos - sh2)
        atoms.set_positions(pos + np.real(d1) - np.real(d2))

        # lob = np.array([0.0, 0.0, 0.0])
        # hib = np.array([29.5 * ux, 0.5 * uy, cell[2, 2]])
        # for i in range(len(atoms)):
        #     if all(((atoms[i].position - lob) * (atoms[i].position - hib)) <= 0.0):
        #         atoms[i].symbol = 'Mg'

        # lob2 = np.array([0.0, 0.5 * uy, 0.0])
        # hib2 = np.array([30 * ux, 1.0 * uy, cell[2, 2]])
        # for i in range(len(atoms)):
        #     if all(((atoms[i].position - lob2) * (atoms[i].position - hib2)) <= 0.0):
        #         atoms[i].symbol = 'Cr'

        atoms.translate(np.array([-83.5 * ux, 80 * uy, 0]))
        atoms.rotate(-self.ag[0][0], 'z')  # rotate back
        self.write_lmp_config_data(atoms, "lmp_init01.txt")

        # add shear
        strain = np.mat([[1.0, 0.0, 0.0],
                         [-0.015, 1.0, 0.0],
                         [0.0, 0.0, 1.0]])

        pos = atoms.get_positions() * strain
        cell = strain * atoms.get_cell()
        atoms.set_cell(cell)
        atoms.set_positions(pos)
        self.write_lmp_config_data(atoms, "lmp_init02.txt")
        self.write_lmp_config_data(atoms, "lmp_init.txt")

    def add_shear(self):
        atoms = ase.io.read("dump/dump.00000", format='lammps-dump')

        # add shear
        strain = np.mat([[1.0, -0.005, 0.0],
                         [0.0, 1.0, 0.0],
                         [0.0, 0.0, 1.0]])
        pos = atoms.get_positions() * strain
        cell = strain * atoms.get_cell()
        atoms.set_cell(cell)
        atoms.set_positions(pos)

        self.write_lmp_config_data(atoms, "lmp_init02.txt")

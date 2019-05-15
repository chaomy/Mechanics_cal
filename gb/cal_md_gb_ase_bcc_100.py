# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-12-03 11:07:29
# @Last Modified by:   1mingfei
# @Last Modified time: 2018-12-03 14:12:22

#import ase.lattice.orthorhombic as otho
from ase.lattice.cubic import BodyCenteredCubic
import ase.io
import os
import numpy as np
import atomman as am
#from utils import stroh_solve
from ase import Atoms
#from ase.lattice.orthorhombic import SimpleOrthorhombicFactory
from numpy import sqrt, deg2rad, floor, cos, sin
import md_pot_data
import glob

class md_gb_ase_bcc_100(object):

    def auto(self):
        self.build_bcc_ase_100_small()
        # minimie
        self.make_perf()
        # minimie
        self.intro_edge_dipole()

    def loop_dup_structures(self):
        files = glob.glob("STRUCT.*")
        for i in range(len(files)):
            atoms = ase.io.read("STRUCT.{}".format(i), format="lammps-dump")
            # rep = int(np.ceil(50.0 / atoms.get_cell()[0, 0]))
            # atoms = atoms.repeat((rep, 1, 1))
            self.write_lmp_config_data(atoms, "CAND.{}".format(i))

    def extend_along_y(self):
        atoms = ase.io.read("CAND.12", format="lammps-data")
        print(len(atoms))

    def loop_clc_init_structures(self):
        dirs = glob.glob("1100_*")
        for i in range(len(dirs)):
            mdir = dirs[i]
            files = glob.glob("{}/dump/*".format(mdir))
            os.system("cp {} ../GENS/STRUCT.{}".format(files[1], i))

    def loop_plot_energy(self):
        data_evo = np.loadtxt('data_evo.txt')
        # data_org = np.loadtxt('data_dirct.txt')
        self.set_111plt((13, 10))
        next(self.keysiter)
        self.ax.plot(data_evo[:, 0], data_evo[:, 3], **next(self.keysiter))
        # self.ax.plot(data_org[:, 0], data_org[:, 1], **next(self.keysiter))
        self.fig.savefig("Fig_evo.png", **self.figsave)
        self.closefig()

    def loop_plt_each(self):
        files = glob.glob("*.txt")
        for file in files:
            data = np.loadtxt(file)
            tags = file.split('.')
            self.set_111plt()
            self.ax.plot(data[:, 0], data[:, 1])
            self.ax.plot(data[:, 0], data[:, 2], lw='4')
            self.fig.savefig("fig_{}.png".format(
                tags[1] + tags[2]), **self.figsave)
            self.closefig()

    def loop_collect_energy(self):
        files = glob.glob("*.txt")
        total = np.ndarray([len(files) + 2, 4])
        total[0, :] = 0
        for i in range(len(files)):
            file = files[i]
            tags = file.split('.')
            angle = float(tags[1].split('_')[1] + '.' + tags[2])
            data = np.loadtxt(file)
            print(angle)
            print(data[-1])
            total[i + 1, 0] = angle
            total[i + 1, 1:] = data[-1]
        total[-1, 0] = 90.0
        total = total[total[:, 0].argsort()]
        np.savetxt('data_evo.txt', total, fmt='%1.8f')
        # for file in files:
        # self.set_111plt()
        # self.ax.plot(data[:, 0], data[:, 1])
        # self.ax.plot(data[:, 0], data[:, 2], lw='4')
        # self.fig.savefig("fig_{}.png".format(
        #     tags[0] + tags[1] + tags[2]), **self.figsave)
        # self.closefig()

    def loop_combine(self):
        self.find_angles_100()
        for e in self.ag[:]:
            mdir = "100_{:.2f}".format(e[0])
            os.chdir(mdir)
            os.system("gb_bound.exe -p ../../gb.param")
            os.chdir(os.pardir)

    def print_angles(self):
        self.find_angles_100()
        for e in self.ag:
            print(e[0], e[1], e[2], e[3])

    def loop_init_bcc100(self):
        self.find_angles_100()
        cn = 0
        for e in self.ag:
            mdir = "100_{:.2f}".format(e[0])
            print(mdir)
            self.mymkdir(mdir)
            self.write_100_small(e)
            # self.write_1100_DFT(e)
            # self.write_1100_DFT_Surf(e)
            # self.write_1100_large(e)
            # self.write_1100_long_thin(e)
            # os.system("cp POSCAR pos_{:02d}".format(cn))
            # os.system("cp INPUTS/* {}".format(mdir))
            os.system("mv lmp.init {}".format(mdir))
            # os.system("mv POSCAR {}/".format(mdir))
            cn += 1

    # to generate surfaces
    '''
    def write_100_DFT_Surf(self, ag):
        ux = self.pot['lattice']
        uy = self.pot['lattice']
        uz = self.pot['lattice'] #check this -1mingfei


        # angle, length, i, j
        atoms = BodyCenteredCubic(latticeconstant=self.pot['lattice'], size=(
            140, 140, 1), symbol=self.pot['element'])

        atoms.rotate(ag[0], 'z')
        cell = atoms.get_cell()
        cell[0, 0], cell[1, 1] = ag[1], 70
        atoms.translate(
            np.array([cell[0, 0] - floor(15 * cos(deg2rad(ag[0]))) * ux,
                      -floor(47 * cos(deg2rad(ag[0]))) * uy, 0]))

        lob = np.array([0.0, 0.0, 0.0])
        hib = np.array([cell[0, 0], 70, cell[2, 2]])
        atoms = self.make_cubic('out', atoms, lob, hib)

        # assign low grain
        vacumm = 1
        if vacumm == 1:
            cell[1, 1] += 20.0
        atoms.set_cell(cell)
        if vacumm == 1:
            atoms.translate(np.array([0.0, 10.0, 0.0]))
        ase.io.write("POSCAR", images=atoms, format="vasp")
    def write_100_DFT(self, ag):
        ux = self.pot['lattice']
        uy = self.pot['lattice']
        uz = self.pot['lattice'] #check this -1mingfei

        VACUMM = 10.0

        # angle, length, i, j
        atoms = BodyCenteredCubic(latticeconstant=self.pot['lattice'], size=(
            140, 140, 1), symbol=self.pot['element'])

        atoms.rotate(ag[0], 'z')
        cell = atoms.get_cell()
        cell[0, 0], cell[1, 1] = ag[1], 120

        atoms.translate(
            np.array([cell[0, 0] - floor(15 * cos(deg2rad(ag[0]))) * ux,
                      -floor(47 * cos(deg2rad(ag[0]))) * uy, 0]))

        lob = np.array([0.0, 0.0, 0.0])
        hib = np.array([cell[0, 0], 0.5 * cell[1, 1], cell[2, 2]])
        atoms = self.make_cubic('out', atoms, lob, hib)

        # the other grain
        atoms2 = BodyCenteredCubic(latticeconstant=self.pot['lattice'], size=(
            140, 140, 1), symbol=self.pot['element'])     # for 1100 72.877

        # atoms2 = othoHCP(latticeconstant=(ux, uy, uz), size=(
        #     80, 80, 3), symbol='Nb')   # for 1100 58.361

        lob = np.array([0.0, 0.5 * cell[1, 1], 0.0])
        hib = np.array([cell[0, 0], cell[1, 1] - 0.2, cell[2, 2]])

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
        m = 0.5 * cell[1, 1]
        assign_gb = 1
        if assign_gb == 1:
            for atom in atoms:
                if atom.position[1] <= m - 35:
                    atom.symbol = 'Re'
                if atom.position[1] >= m + 35:
                    atom.symbol = 'Ta'
                if atom.position[1] >= m + 40 or atom.position[1] <= m - 40:
                    atom.symbol = 'Mo'

        vacumm = 0
        if vacumm == 1:
            cell[1, 1] += 40.0
        atoms.set_cell(cell)
        if vacumm == 1:
            atoms.translate(np.array([0.0, 20.0, 0.0]))

        idx = []
        for atom in atoms:
            if atom.symbol in ['Mo']:
                idx.append(atom.index)
        del atoms[idx]

        self.write_lmp_config_data(atoms, "lmp.init")
        '''

    def write_100_small(self, ag):
        ux = self.pot['lattice']
        uy = self.pot['lattice']
        uz = self.pot['lattice'] 

        # angle, length, i, j
        atoms = BodyCenteredCubic(latticeconstant=self.pot['lattice'], size=(
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
        atoms2 = BodyCenteredCubic(latticeconstant=self.pot['lattice'], size=(
            130, 130, 2), symbol=self.pot['element'])     

        lob = np.array([0.0, 0.5 * cell[1, 1], 0.0])
        hib = np.array([cell[0, 0], cell[1, 1] - 0.2, cell[2, 2]])
        # hib = np.array([cell[0, 0], cell[1, 1] - VACUMM, cell[2, 2]])

        atoms2.rotate(-ag[0], 'z')
        # atoms2.translate(np.array([-10 * ux, cell[1, 1], 0]))  # for
        # 72.877
        atoms2.translate(np.array(
            [floor(60 * cos(deg2rad(ag[0]))) * -ux,
             cell[1, 1] - floor(12 * cos(deg2rad(ag[0]))) * uy, 0]))  # for 72.877

        atoms2.translate(np.array([-cell[0,0]*2.0, -cell[1,1]*0.1, 0.0]))
        atoms2 = self.make_cubic('out', atoms2, lob, hib)
        #atoms2.translate(np.array([0.0, 0.2, 0.25 * uz]))
        atoms2.translate(np.array([0.0, 0.2, 0.0]))
        atoms.extend(atoms2)

        # assign low grain
        m = 0.5 * cell[1, 1]
        assign_gb = 1
        if assign_gb == 1:
            for atom in atoms:
                # if atom.position[1] <= m - 10 or atom.position[1] >= m + 10:    # buff
                #     atom.symbol = 'Mo'
                if atom.position[1] <= m - 40:
                    atom.symbol = 'Re'
                if atom.position[1] >= m + 40:
                    atom.symbol = 'Ta'
                if atom.position[1] >= m + 45 or atom.position[1] <= m - 45:
                    atom.symbol = 'Mo'

        vacumm = 0
        if vacumm == 1:
            cell[1, 1] += 40.0
        atoms.set_cell(cell)
        if vacumm == 1:
            atoms.translate(np.array([0.0, 20.0, 0.0]))

        idx = []
        for atom in atoms:
            if atom.symbol in ['Mo']:
                idx.append(atom.index)
        del atoms[idx]

        rep = int(np.ceil(50 / cell[0, 0]))
        if rep > 1:
            atoms = atoms.repeat((rep, 1, 1))
        self.write_lmp_config_data(atoms, "lmp.init")

    def write_100_long_thin(self, ag):
        ux = self.pot['lattice']
        uy = self.pot['lattice']
        uz = self.pot['lattice'] #check this -1mingfei

        # angle, length, i, j
        atoms = BodyCenteredCubic(latticeconstant=self.pot['lattice'], size=(
            500, 500, 2), symbol=self.pot['element'])

        atoms.rotate(ag[0], 'z')
        cell = atoms.get_cell()
        cell[0, 0], cell[1, 1] = 2 * ag[1], 400

        atoms.translate(
            np.array([cell[0, 0] - floor(15 * cos(deg2rad(ag[0]))) * ux,
                      -floor(47 * cos(deg2rad(ag[0]))) * uy, 0]))

        lob = np.array([0.0, 0.0, 0.0])
        hib = np.array([cell[0, 0], 0.5 * cell[1, 1], cell[2, 2]])
        atoms = self.make_cubic('out', atoms, lob, hib)

        # the other grain
        atoms2 = BodyCenteredCubic(latticeconstant=self.pot['lattice'], size=(
            500, 500, 2), symbol=self.pot['element'])     # for 1100 72.877

        lob = np.array([0.0, 0.5 * cell[1, 1], 0.0])
        hib = np.array([cell[0, 0], cell[1, 1] - 0.2, cell[2, 2]])

        atoms2.rotate(-ag[0], 'z')
        atoms2.translate(np.array(
            [floor(60 * cos(deg2rad(ag[0]))) * -ux,
             cell[1, 1] - floor(22 * cos(deg2rad(ag[0]))) * uy, 0]))  # for 72.877

        atoms2 = self.make_cubic('out', atoms2, lob, hib)
        #atoms2.translate(np.array([0.0, 0.2, 0.25 * uz]))
        atoms.extend(atoms2)

        # assign low grain
        m = 0.5 * cell[1, 1]
        assign_gb = 1
        if assign_gb == 1:
            for atom in atoms:
                if atom.position[1] <= m - 40:
                    atom.symbol = 'Re'
                if atom.position[1] >= m + 40:
                    atom.symbol = 'Ta'
                if atom.position[1] >= m + 180 or atom.position[1] <= m - 180:
                    atom.symbol = 'Mo'

        atoms.set_cell(cell)
        idx = []
        for atom in atoms:
            if atom.symbol in ['Mo']:
                idx.append(atom.index)
        del atoms[idx]

        # rep = int(np.ceil(50 / cell[0, 0]))
        # if rep > 1:
        #     atoms = atoms.repeat((rep, 1, 1))
        self.write_lmp_config_data(atoms, "lmp.init")

    # def write_1100_large(self, ag):
    #     ux = self.pot['ahcp']
    #     uy = self.pot['chcp']
    #     uz = self.pot['ahcp'] * sqrt(3.)

    #     # angle, length, i, j
    #     atoms = othoHCP(latticeconstant=(ux, uy, uz), size=(
    #         400, 400, 1), symbol=self.pot['element'])

    #     atoms.rotate(ag[0], 'z')
    #     cell = atoms.get_cell()
    #     cell[0, 0], cell[1, 1] = 1 * ag[1], 280
    #     atoms.translate(
    #         np.array([cell[0, 0] - floor(15 * cos(deg2rad(ag[0]))) * ux,
    #                   -floor(47 * cos(deg2rad(ag[0]))) * uy, 0]))

    #     lob = np.array([0.0, 0.0, 0.0])
    #     hib = np.array([cell[0, 0], 0.5 * cell[1, 1], cell[2, 2]])
    #     atoms = self.make_cubic('out', atoms, lob, hib)

    #     # the other grain
    #     atoms2 = othoHCPB(latticeconstant=(ux, uy, uz), size=(
    #         400, 400, 1), symbol='Al')     # for 1100 72.877

    #     # atoms2 = othoHCP(latticeconstant=(ux, uy, uz), size=(
    #     #     80, 80, 3), symbol='Nb')   # for 1100 58.361

    #     lob = np.array([0.0, 0.5 * cell[1, 1], 0.0])
    #     hib = np.array([cell[0, 0], cell[1, 1] - 0.2, cell[2, 2]])

    #     atoms2.rotate(-ag[0], 'z')
    #     # atoms2.translate(np.array([-10 * ux, cell[1, 1], 0]))  # for
    #     # 72.877
    #     atoms2.translate(np.array(
    #         [floor(60 * cos(deg2rad(ag[0]))) * -ux,
    # cell[1, 1] - floor(12 * cos(deg2rad(ag[0]))) * uy, 0]))  # for 72.877

    #     atoms2 = self.make_cubic('out', atoms2, lob, hib)
    #     atoms2.translate(np.array([0.0, 0.2, 0.25 * uz]))
    #     atoms.extend(atoms2)

    #     # lob = np.array([0.0, 40, 0.0])
    #     # hib = np.array([cell[0, 0], cell[1, 1] - 40, cell[2, 2]])
    #     # atoms = self.assign_cubic(atoms, 'out', 'W', lob, hib)

    #     # assign gb region
    #     assign_gb = 1
    #     if assign_gb == 1:
    #         for atom in atoms:
    #             if atom.position[1] <= 17:
    #                 atom.symbol = 'W'
    #             if atom.position[1] >= cell[1, 1] - 17:
    #                 atom.symbol = 'Mo'

    #     # to add vacancy optional
    #     cell[1, 1] += 30
    #     atoms.translate(np.array([0.0, 15, 0.0]))
    #     atoms.set_cell(cell)
    #     self.write_lmp_config_data(atoms, "lmp_init.txt")
    #     self.make_repeat(atoms)

    def build_bcc_ase_100_small(self):  # to examine the GB structures
        self.find_angles_100(il=[[], [1]], jl=[2])    # 72.877    ABAB --
        self.write_100_small(self.ag[0])

    def build_bcc_ase_100(self):
        self.fnnd_angles_100(il=[[1], [1]], jl=[1])    # 58.361
        self.write_100_large(self.ag[0])


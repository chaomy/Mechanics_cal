#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: yang37
# @Date:   2017-06-12 17:03:43
# @Last Modified by:   chaomy
# @Last Modified time: 2018-05-12 20:22:50


from multiprocessing import Pool
from scipy.interpolate import InterpolatedUnivariateSpline
from optparse import OptionParser
from utils import Intro_vasp
import ase.utils.geometry
import md_pot_data
import os
import ase
import glob
import numpy as np
import matplotlib.pylab as plt
import ase.lattice
import gn_config
import get_data
import gn_lmp_infile
import gn_pbs

#  import atomman as am
#  import atomman.lammps as lmp
#  import atomman.unitconvert as uc


class cal_barrier(object):

    def __init__(self):
        self._burger = self.pot['lattice'] * np.sqrt(3.) / 2.
        self.num_configures = 15
        self.set_keys()
        self.mfigsize = (8.5, 4.3)

    def get_init_final_state(self):  # state initial #
        poscar_atoms = ase.io.read("./state_init/lmp_init.cfg", format='cfg')
        atoms_tobe_changed = ase.io.read("./state_init/W.0.cfg", format='cfg')

        map_list = self.map_atoms_list(poscar_atoms, atoms_tobe_changed)
        new_atoms = self.sort_atoms(atoms_tobe_changed, map_list)

        ase.io.write("peierls_init.cfg", new_atoms, format='cfg')

        # state final #
        poscar_atoms = ase.io.read("./state_final/lmp_init.cfg",  format='cfg')
        atoms_tobe_changed2 = ase.io.read(
            "./state_final/W.0.cfg", format='cfg')
        map_list = self.map_atoms_list(poscar_atoms, atoms_tobe_changed2)
        new_atoms2 = self.sort_atoms(atoms_tobe_changed2, map_list)

        ase.io.write("peierls_final.cfg", new_atoms2, format='cfg')

    def aniso_dipole_peierls_barrier_image(self, sizen=1):
        sx = 10.0 * sizen
        sy = 5 * sizen
        ix = 10.5 * sizen

        unitx = np.sqrt(6) / 3. * self.pot['lattice']
        unity = np.sqrt(2) / 2. * self.pot['lattice']

        # dx, dy = 0.1, -0.133333333
        dx, dy = 0.0, 0.0

        ci1 = [(sx + dx) * unitx, (sy + 1. / 3. + dy) * unity]
        atoms = self.bcc_screw_dipole_alongz_with_image(
            self.set_dipole_box(1), ci1)
        atoms = self.convert_alongz_to_alongy(atoms)
        self.write_lmp_config_data(atoms, "init.txt")
        ase.io.write("init.ref", images=atoms, format="vasp")

        ci2 = [(sx + 1.0) * unitx, (sy + 1. / 3.) * unity]
        atoms = self.bcc_screw_dipole_alongz_with_image(
            self.set_dipole_box(1), ci2)
        atoms = self.convert_alongz_to_alongy(atoms)
        self.write_lmp_config_data(atoms, "final.txt")
        np.savetxt("dis_center.txt", np.array(
            [ci1[0], ci1[1], ci1[0] + unitx, ci1[1]]))

    def aniso_dipole_peierls_barrier(self, sizen=1):
        sx = 10.0 * sizen
        sy = 5 * sizen
        ix = 10.5 * sizen

        unitx = np.sqrt(6) / 3. * self.pot['lattice']
        unity = np.sqrt(2) / 2. * self.pot['lattice']

        ci1 = [(sx) * unitx, (sy + 1. / 3.) * unity]
        ci2 = [(sx + ix) * unitx, (sy + 2. / 3.) * unity]

        atoms = self.set_dipole_box(sizen)
        atoms = self.bcc_screw_dipole_alongz_atoms(atoms, ci1, ci2)
        atoms = self.convert_alongz_to_alongy(atoms)
        self.write_lmp_config_data(atoms, "init.txt")

        movex = 1.0
        cf1 = [(sx + movex) * unitx, (sy + 1. / 3.) * unity]
        cf2 = [(sx + ix + movex) * unitx, (sy + 2. / 3.) * unity]

        atoms = self.set_dipole_box(sizen)
        atoms = self.bcc_screw_dipole_alongz_atoms(atoms, cf1, cf2)
        atoms = self.convert_alongz_to_alongy(atoms)
        self.write_lmp_config_data(atoms, "final.txt")

        np.savetxt("dis_center.txt", np.array(
            [ci1[0], ci1[1], cf1[0], cf1[1]]))
        print((sx) * unitx, (sy) * unity)

    def dipole_peierls_barrier(self, sizen=1):
        atoms = self.set_dipole_box_alongy(sizen)

        sx = 10.0 * sizen
        sy = 5 * sizen
        ix = 10.5 * sizen

        unitx = np.sqrt(6) / 3. * self.pot['lattice']
        unity = np.sqrt(2) / 2. * self.pot['lattice']

        ci1 = [(sx) * unitx, (sy + 1. / 3.) * unity]
        ci2 = [(sx + ix) * unitx, (sy + 2. / 3.) * unity]
        self.write_lmp_config_data(self.intro_dipole_screw_atoms_LMP(
            atoms.copy(), center=[ci1, ci2]), "init.txt")

        movex = 1.0
        cf1 = [(sx + movex) * unitx, (sy + 1. / 3.) * unity]
        cf2 = [(sx + ix + movex) * unitx, (sy + 2. / 3.) * unity]
        self.write_lmp_config_data(self.intro_dipole_screw_atoms_LMP(
            atoms.copy(), center=[cf1, cf2]), "final.txt")

        np.savetxt("dis_center.txt", np.array(
            [ci1[0], ci1[1], cf1[0], cf1[1]]))

        # write in the middle
        ci1[0] += 0.5 * (cf1[0] - ci1[0])
        ci2[0] += 0.5 * (cf2[0] - ci2[0])
        ci1[1] -= 0.1
        ci2[1] += 0.1
        ase.io.write("test.txt", images=self.intro_dipole_screw_atoms_LMP(
            atoms.copy(), center=[ci1, ci2]), format="vasp")

        # for i in range(16):
        #     movex = (i + 1) / 16.
        #     c1 = [(sx + movex) * unitx, (sy + 1. / 3.) * unity]
        #     c2 = [(sx + ix + movex) * unitx, (sy + 2. / 3.) * unity]
        #     self.write_lmp_coords(self.intro_dipole_screw_atoms_LMP(
        #         atoms.copy(), center=[c1, c2], lattice=self.pot['lattice']),
        #         "lmp_{:d}.txt".format(i + 1))

        # cluster method  (large supercell)
    def cluster_peierls_barrier(self):
        alat = self.pot['lattice']
        move1 = [-alat * np.sqrt(6.) / 3., 0, 0]
        self.make_screw_plate(size=[200, 260, 2], rad=[340, 380],
                              move=move1, tag='[211]',
                              filename="final.txt", opt='neb')
        os.system("cp init.txt  init0.txt")
        os.system("cp final.txt init1.txt")

        # Mid is 150 160
        # Mid2   200 210
        # Lar    400 410

    def interp_peierls(self):
        atoms_init = ase.io.read(
            glob.glob("bcc.init.*")[-1], format='lammps-dump')
        atoms_final = ase.io.read(
            glob.glob("bcc.final.*")[-1], format='lammps-dump')

        self.write_lmp_config_data(atoms_init, "init.data")
        self.write_lmp_config_data(atoms_final, "final.data")

        atoms_interp = atoms_init.copy()

        position_init = atoms_init.get_positions()
        position_final = atoms_final.get_positions()
        delta = position_final - position_init

        num_list = [47, 78, 79, 49, 50, 181, 180, 178, 179, 149]
        delta_copy = np.zeros(np.shape(delta))
        for i in range(len(num_list)):
            delta_copy[num_list[i]] = delta[num_list[i]]

        print(delta_copy)
        delta = delta_copy
        inter_num = self.num_configures
        delta = delta / inter_num
        for i in range(inter_num):
            dir_name = 'dir-%03d' % (i)
            self.mymkdir(dir_name)
            os.chdir(dir_name)
            filename_data = "lmp_init_%d.txt" % (i)
            atoms_interp = position_init + i * delta
            atoms_init.set_positions(atoms_interp)
            self.write_lmp_config_data(atoms_init, filename_data)
            os.chdir(os.pardir)

    #  reaction coordinate method #
    def peierls(self):
        lattice = 3.143390
        e1 = np.array([1.,   1.,  -2.])
        e2 = np.array([-1.,  1.,   0])
        e3 = np.array([0.5,  0.5,  0.5])

        # 3, 5;   5, 9;  7, 11
        r, s = 7,  11
        v1 = r * e1
        v2 = 0.5 * (r * e1 + s * e2) + 0.5 * e3
        v3 = e3

        inter_num = self.num_configures
        delta = 1. / inter_num

        for i in range(inter_num):
            dir_name = 'dir-%03d' % (i)
            self.mymkdir(dir_name)
            os.system("cp  ./w_eam4.fs  %s" % (dir_name))
            os.chdir(dir_name)

            # generate config #
            atoms = self.set_bcc_convention(
                [v1, v2, v3], (r, 1, 1))  # z periodic 12
            atoms2 = atoms.copy()
            atoms2 = self.cut_half_atoms(atoms2)

            #  movex = np.sqrt(6.)/3. * lattice * delta * i;
            ase.io.write("lmp_perf.cfg", atoms2, "cfg")
            s = i * delta
            atoms = self.intro_dipole_screw_atoms(
                atoms, lattice, move_x=None, input_s=s)
            atoms = self.cut_half_atoms(atoms)

            ase.io.write("lmp_init.cfg", atoms, "cfg")
            filename_data = "lmp_init_%d.txt" % (i)
            # self.gn_md_minimize_cfg(filename_data, "w_eam4.fs", "W", fix=True)

            self.write_lmp_config_data(atoms, filename_data)
            os.chdir(os.pardir)

    def collect_peierls_energy(self):
        num = self.num_configures
        energy_list = []
        count_list = []
        for i in range(0, num):
            dir_name = 'dir-%03d' % (i)
            config_name = 'config-%03d.cfg' % (i)
            config_png = 'init-%03d.png' % (i)
            os.chdir(dir_name)

            energy_list.append(self.md_get_final_energy())
            count_list.append(i)

            filelist = glob.glob("./cfg/*")
            file = filelist[-1]

            os.system("python ../../DD_map_dis.py")
            os.system("cp  %s  ../%s" % (file, config_name))
            os.system("cp  init.png  ../%s" % (config_png))

            os.chdir(os.pardir)

        energy_list = np.array(energy_list)
        energy_list = energy_list - np.min(energy_list)

        with open("peierls.txt", 'w') as fid:
            for i in range(num):
                fid.write("%d %7.6f\n" % (count_list[i], energy_list[i]))
        print(energy_list)
        self.plot()

    def plot_engy(self, **kwargs):
        if 'engy' not in list(kwargs.keys()):
            (distances, engy) = np.loadtxt("./peierls.txt")
        else:
            distances, engy = kwargs['distances'], kwargs['engy']

        if 'ax' not in list(kwargs.keys()):
            self.set_111plt(self.mfigsize)
            ax = self.ax
        else:
            ax = kwargs['ax']

        if 'pltkeys' in list(kwargs.keys()):
            pltkeys = kwargs['pltkeys']
        else:
            pltkeys = self.pltkwargs

        ax.plot(distances, engy, label=kwargs['label'], **pltkeys)
        ax.legend(**self.legendarg)
        return ax

    def read_peierls_barrier_neb(self):
        file_list = glob.glob("log.lammps.*")
        nlogs = len(file_list)
        neb_energy = []
        for i in range(nlogs):
            mfile = "log.lammps.%d" % (i)
            neb_energy.append(self.md_get_final_energy_e(mfile))
        neb_energy = np.array(neb_energy)
        neb_energy -= np.min(neb_energy)
        # neb_energy *= (1. / 4.)     # dislocation dipole (count how many
        # burgers vector)
        disp = np.linspace(0.0, self._burger, len(neb_energy))
        data = np.array([disp, neb_energy])
        np.savetxt('pengy.txt', data)
        return data

    def set_pltkargs(self):
        keyslist = [{'linestyle': '-', 'color': 'b', 'linewidth': 2,
                     'marker': 'o'},
                    {'linestyle': '--', 'color': 'r', 'linewidth': 2,
                        'marker': '<'},
                    {'linestyle': ':', 'color': 'g', 'linewidth': 2,
                        'marker': '>'}]
        return keyslist

    def plot_multi_peierls_barrier(self):
        configlist = ['pengy.txt.eam', 'pengy.txt.adp']
        labellist = ['eam', 'adp']
        self.set_111plt(self.mfigsize)
        keyslist = self.set_pltkargs()
        cnt = 0
        for config in configlist:
            (distances, engy) = np.loadtxt(config)
            inargs = {'distances': distances, 'engy': engy, 'ax': self.ax,
                      'pltkeys': keyslist[cnt],
                      'label': labellist[cnt]}
            cnt += 1
            self.ax = self.plot_engy(**inargs)
        self.ax.set_xlabel('displacement [A]', {'fontsize': self.mlabelsize})
        self.ax.set_ylabel('energy [eV]', {'fontsize': self.mlabelsize})
        plt.yticks(size=self.mlabelsize)
        plt.xticks(size=self.mlabelsize)
        plt.savefig("mpeierls.png", **self.figsave)

    def cal_peierls_stress(self):
        if os.path.isfile('pengy.txt'):
            data = np.loadtxt('pengy.txt')
        else:
            data = self.read_peierls_barrier_neb()
        spl = InterpolatedUnivariateSpline(data[0], data[1], k=3)
        splder1 = spl.derivative()
        interpnts = np.linspace(0, self._burger, 51)
        stress = splder1(interpnts)
        plt.savefig("tmp.png")
        print(np.max(stress))

    def convert_dump_to_poscar(self):
        initfile = glob.glob("bcc.init.*")[-1]
        finalfile = glob.glob("bcc.final.*")[-1]
        print(initfile, finalfile)
        initatoms = ase.io.read(initfile, format='lammps-dump')
        finalatoms = ase.io.read(finalfile, format='lammps-dump')
        ase.io.write("POSCAR00", images=initatoms, format='vasp')
        ase.io.write("POSCAR01", images=finalatoms, format='vasp')


#     (options, args) = parser.parse_args()
#     drv = cal_barrier()
#     dispatcher = {'dipole': drv.dipole_peierls_barrier,
#                   'cluster': drv.cluster_peierls_barrier,
#                   'vasp': drv.convert_dump_to_poscar,
#                   'stress': drv.cal_peierls_stress,
#                   'cmp': drv.plot_multi_peierls_barrier,
#                   'looprcut': drv.loop_rcut,
#                   'minimize': drv.multi_thread_minimize}

    #  drv.peierls()
    #  job.multi_thread_minimize()

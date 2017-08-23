#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-06-28 00:35:14
# @Last Modified by:   chaomy
# @Last Modified time: 2017-08-22 20:14:45

import os
import glob
import copy
import numpy as np
from multiprocessing import Pool
from optparse import OptionParser
from itertools import cycle
import matplotlib.pylab as plt
import plt_drv

try:
    import gn_config
    import get_data
    import gn_pbs
    import gn_lmp_infile
    import output_data

except ImportError:
    print "error during import"


def unwrap_self_run_lammps(arg, **kwarg):
    return cal_md_gsf.run_lmp_gsf(*arg, **kwarg)


class cal_md_gsf(gn_config.bcc,
                 gn_config.fcc,
                 gn_config.hcp,
                 get_data.get_data,
                 gn_pbs.gn_pbs,
                 gn_lmp_infile.gn_md_infile,
                 output_data.output_data,
                 plt_drv.plt_drv):

    def __init__(self, gsf_surface_type='110'):
        get_data.get_data.__init__(self)
        gn_pbs.gn_pbs.__init__(self)
        plt_drv.plt_drv.__init__(self)
        gn_lmp_infile.gn_md_infile.__init__(self)
        output_data.output_data.__init__(self)
        self._pot = self.load_data('pot.dat')
        #  self._pot = md_pot.Nb_adp

        self._gsf_surface_type = gsf_surface_type
        self._surface_element = self._pot['element']  # default
        self._surface_lattice_constant = self._pot['lattice']
        self._structure = self._pot['structure']
        self._gsf_potential = self._pot['file']
        self._pot_type = self._pot['pair_style']
        if self._structure == 'bcc':
            gn_config.bcc.__init__(self, self._pot)

        self.set_lattce_constant(self._surface_lattice_constant)
        self.set_element(self._surface_element)
        self.set_config_file_format("lmp")
        self.set_relax_type()
        self.sample_gsf_num = 21
        self.disp_delta = 1. / (self.sample_gsf_num - 1)
        self.config_file = "lmp_init.txt"

        self.root_dir = os.getcwd()
        return

    def set_relax_type(self, relaxtag='relaxed'):
        self.relaxtag = relaxtag
        return

    def set_md_gsf_potential(self, in_potential):
        self._gsf_potential = in_potential
        return

    def set_gsf_surface_type(self, gsf_surface_type):
        self._gsf_surface_type = gsf_surface_type
        return

    def gn_gsf_atoms(self):
        if self._gsf_surface_type == '100_100':
            self._surface_direction = [[1, 0, 0],
                                       [0, 1, 0],
                                       [0, 0, 1]]
            atoms = \
                self.set_bcc_convention(in_direction=self._surface_direction,
                                        in_size=(1, 1, 18))
            for i in range(8):
                atoms.pop()
            return atoms

        elif self._gsf_surface_type == '110_110':
            self._surface_direction = [[-1, 1, 0],
                                       [0, 0, 1],
                                       [1, 1, 0]]
            atoms = \
                self.set_bcc_convention(in_direction=self._surface_direction,
                                        in_size=(1, 1, 14))
            for i in range(12):
                atoms.pop()
            return atoms

        elif self._gsf_surface_type == '111_110':
            self._surface_direction = [[-1, 1, 1],
                                       [1, -1, 2],
                                       [1, 1, 0]]
            atoms = \
                self.set_bcc_convention(in_direction=self._surface_direction,
                                        in_size=(1, 2, 14))
            for i in range(48):
                atoms.pop()
            return atoms

        elif self._gsf_surface_type == '111_211':
            self._surface_direction = [[-1, 1, 1],
                                       [0, -1, 1],
                                       [2, 1, 1]]
            atoms = \
                self.set_bcc_convention(in_direction=self._surface_direction,
                                        in_size=(1, 1, 14))
            for i in range(12):
                atoms.pop()
            return atoms

        elif self._gsf_surface_type == '112_111':
            self._surface_direction = [[1, 1, -2],
                                       [-1, 1, 0],
                                       [1, 1, 1]]
            atoms = \
                self.set_bcc_convention(in_direction=self._surface_direction,
                                        in_size=(1, 1, 14))
            for i in range(24):
                atoms.pop()
            return atoms

        elif self._gsf_surface_type == '110_111':
            self._surface_direction = [[1, 1, 0],
                                       [-1, 1, 2],
                                       [1, 1, 1]]
            atoms = \
                self.set_bcc_convention(in_direction=self._surface_direction,
                                        in_size=(1, 1, 14))
            for i in range(24):
                atoms.pop()
            return atoms

        elif self._gsf_surface_type == '110_211':
            self._surface_direction = [[1, 1, 0],
                                       [-1, 1, 2],
                                       [1, 1, 1]]
            atoms = \
                self.set_bcc_convention(in_direction=self._surface_direction,
                                        in_size=(1, 1, 12))
            for i in range(12):
                atoms.pop()
            return atoms

    def gn_displacement(self,
                        atoms,
                        displacement_vector):
        positions = atoms.get_positions()
        atom_num = len(positions)
        displacement = copy.deepcopy(positions)
        cut = 0.5 * np.max(positions[:, 2])
        for i in range(atom_num):
            print positions[i, 2]
            if positions[i, 2] < cut:
                displacement[i] = [0, 0, 0]
            else:
                displacement[i] = displacement_vector
        return displacement

    def prepare_md_inputs(self,
                          config_file=None,
                          in_potential=None):
        if config_file is not None:
            self.config_file = config_file
        if in_potential is not None:
            self._gsf_potential = in_potential
        self.gn_gsf_minimize(config_file=self.config_file,
                             potential_type=self._pot_type,
                             potential_file=self._gsf_potential,
                             element=self._surface_element,
                             tag='relaxed')
        os.system("cp  ../{}  .".format(self._gsf_potential))
        return

    def run_lmp_gsf(self, loc_dir):
        os.chdir(loc_dir)
        os.system("lmp_mpi -in in.md_gsf")
        os.chdir(self.root_dir)
        return

    def multi_thread_gsf(self):
        dir_list = glob.glob("dir-*")
        num_threads = len(dir_list)
        pool = Pool(processes=num_threads)
        pool.map(unwrap_self_run_lammps,
                 zip([self] * num_threads,
                     dir_list))
        return

    def collect_gsf_energy(self):
        disp_list, energylist, area_list = [], [], []
        for i in range(0, self.sample_gsf_num):
            dir_name = 'dir-x-%03d-%s' % (i, self._gsf_surface_type)
            if os.path.isdir(dir_name):
                os.chdir(dir_name)

                disp_list.append(i * self.disp_delta)
                energylist.append(self.md_get_final_energy())
                super_cell = self.md_get_cell()
                area_list.append(self.cal_poscar_xy_area(super_cell))
                os.chdir(self.root_dir)

        energylist = np.array(energylist)
        energylist = energylist / np.average(np.array(area_list))
        energylist = energylist - np.min(energylist)
        disp_list = np.array(disp_list)
        self.output_disp_energy(disp_list,
                                energylist,
                                "gsf_{}_{}.txt".format(self._pot['pair_type'],
                                                       self._gsf_surface_type))
        self.plot_md_gsf(disp_list, energylist)
        return (disp_list, energylist)

    def plot_md_gsf(self,
                    delta=None,
                    energy=None,
                    filename='md_gsf.png'):
        self.set_keys()
        self.set_111plt((8, 4))
        self.ax.plot(delta, energy,
                     label="$displacement-energy$",
                     **self.pltkwargs)
        plt.savefig(filename, **self.figsave)
        return

    def plot_multi_gsf_curv(self,
                            potlist,
                            typelist,
                            fname='gsf_compare.png'):
        self.set_keys()
        self.set_211plt()
        axlist = [self.ax1, self.ax2]
        ylabiter = cycle([
            r"$\gamma$[{}]({}) [eV/$\AA^2$]".format(typelist[0][:3],
                                                    typelist[0][-3:]),
            r"$\gamma$[{}]({}) [eV/$\AA^2$]".format(typelist[1][:3],
                                                    typelist[1][-3:])
        ])
        for pottype in potlist:
            filename = "gsf_{}_{}.txt".format(pottype,
                                              typelist[0])
            pltlabel = "{}".format(pottype)
            data = np.loadtxt(filename)
            self.ax1.plot(data[:, 0], data[:, 1],
                          label=pltlabel, **next(self.keysiter))

            filename = "gsf_{}_{}.txt".format(pottype,
                                              typelist[1])
            pltlabel = "{}".format(pottype)
            data = np.loadtxt(filename)
            self.ax2.plot(data[:, 0], data[:, 1],
                          label=pltlabel, **next(self.keysiter))
            self.ax2.set_xlabel("normalized displacement along $[{}]$".format(self._gsf_surface_type[:3]),
                                {'fontsize': (self.myfontsize - 3)})   # (110): -110  (11-2) -110
        self.set_tick_size(*axlist)
        self.add_y_labels(ylabiter, *axlist)
        self.add_legends(*axlist)
        plt.savefig(fname, **self.figsave)
        return

    def plot_multi_type_gsf_curv(self, typelist, fname='gsf_compare.png'):
        pltdrv = plt_drv.plt_drv()
        self.set_keys()
        self.set_111plt((8, 4))
        plt.rc('xtick', labelsize='large')
        plt.rc('ytick', labelsize='large')
        cnt = 0
        for gsftype in typelist:
            filename = "gsf_{}_{}.txt".format(self._pot['pair_type'],
                                              gsftype)
            pltlabel = "{}_{}".format(self._pot['pair_type'], gsftype)
            data = np.loadtxt(filename)
            self.ax.plot(data[:, 0], data[:, 1],
                         label=pltlabel,
                         **next(self.keysiter))
            cnt += 1
            self.ax.legend(**self.legendarg)
        plt.xlabel("normalized displacement along $[{}]$".format(gsftype[:3]),
                   {'fontsize': (self.myfontsize)})   # (110): -110  (11-2) -110
        plt.ylabel("stacking fault energy $[eV/\AA^{2}]$",
                   {'fontsize': (self.myfontsize)})   # (110): -110  (11-2) -110
        plt.savefig(fname, **self.figsave)
        return

    # trans dft to md
    def trans_data_format(self):
        data = np.loadtxt("DATA")
        disp = data[:, 1]
        area = data[:, 2]
        energy = data[:, -1]
        energy = energy / area[0]
        energy = energy - np.min(energy)
        fid = open("data_out.txt", 'w')
        for i in range(len(disp)):
            fid.write("{:04f} {:04f} \n".format(disp[i], energy[i]))
        fid.close()
        return

    def md_single_dir_gsf(self):
        atoms = self.gn_gsf_atoms()
        perf_cells = copy.deepcopy(atoms.get_cell())

        for i in range(0, self.sample_gsf_num):
            dir_name = 'dir-x-%03d-%s' % (i, self._gsf_surface_type)
            self.mymkdir(dir_name)
            os.chdir(dir_name)

            disp_vector = [i * self.disp_delta, 0, 0]
            disp_matrix_direct = self.gn_displacement(atoms.copy(),
                                                      disp_vector)

            disp_matrix = copy.deepcopy(disp_matrix_direct)

            cell_length_x = perf_cells[0, 0]
            disp_matrix[:, 0] = disp_matrix_direct[:, 0] * cell_length_x

            local_atoms = atoms.copy()
            local_atoms.translate(disp_matrix)

            self.write_lmp_config_data(local_atoms)

            self.prepare_md_inputs()
            os.system("cp lmp_init.txt ../lmp_init_{0:03d}.txt".format(i))
            os.chdir(self.root_dir)
        return


if __name__ == '__main__':
    usage = "usage:%prog [options] arg1 [options] arg2"
    parser = OptionParser(usage=usage)
    parser.add_option("-t", "--mtype",
                      action="store",
                      type="string",
                      dest="mtype",
                      default="prp_r")
    (options, args) = parser.parse_args()
    drv = cal_md_gsf(gsf_surface_type='111_211')
    if options.mtype.lower() == 'prep':
        drv.md_single_dir_gsf()

    elif options.mtype.lower() == 'run':
        drv.multi_thread_gsf()

    elif options.mtype.lower() == 'collect':
        drv.collect_gsf_energy()

    elif options.mtype.lower() == 'trans':
        drv.trans_data_format()

    elif options.mtype.lower() in ['compare', 'cmp']:
        potlist = ['adp', 'pbe']
        typelist = ['111_211', '111_110']
        drv.plot_multi_gsf_curv(potlist, typelist)

    elif options.mtype.lower() == 'twopath':
        typelist = ['111_211', '111_110']
        drv.plot_multi_type_gsf_curv(typelist)

    elif options.mtype.lower() == 'relaxed':
        potlist = ['adp', 'pbe']
        typelist = ['111_211', '111_110']
        drv.set_relax_type('relaxed')
        for gsftype in typelist:
            drv.set_gsf_surface_type(gsftype)
            drv.md_single_dir_gsf()
            drv.multi_thread_gsf()
            drv.collect_gsf_energy()
        drv.plot_multi_gsf_curv(potlist, typelist, 'gsf_relaxed.png')

    elif options.mtype.lower() == 'unrelaxed':
        potlist = ['adp', 'pbe']
        typelist = ['111_211', '111_110']
        drv.set_relax_type('unrelaxed')
        for gsftype in typelist:
            drv.set_gsf_surface_type(gsftype)
            drv.md_single_dir_gsf()
            drv.multi_thread_gsf()
            drv.collect_gsf_energy()
        #  drv.plot_multi_type_gsf_curv(typelist)
        drv.plot_multi_gsf_curv(potlist, typelist, 'gsf_unrelaxed.png')

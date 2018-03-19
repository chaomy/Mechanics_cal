#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-06-28 00:35:14
# @Last Modified by:   chaomy
# @Last Modified time: 2018-03-17 07:47:12


import os
import glob
import copy
import numpy as np
from multiprocessing import Pool
from optparse import OptionParser
from itertools import cycle
import output_data


def unwrap_self_run_lammps(arg, **kwarg):
    return cal_md_gsf.run_lmp_gsf(*arg, **kwarg)


class cal_md_gsf(output_data.output_data):

    def __init__(self):
        output_data.output_data.__init__(self)
        self.set_relax_type()
        self.config_file = "lmp_init.txt"
        self.infile = "in.md_gsf"

    def set_relax_type(self, relaxtag='relaxed'):
        self.relaxtag = relaxtag

    def set_mgsf(self, mgsf):
        self.mgsf = mgsf

    def prepare_md_inputs(self, config_file=None,
                          in_potential=None):
        if config_file is not None:
            self.config_file = config_file
        self.gn_gsf_minimize(config_file=self.config_file,
                             potential_type=self.pot['pair_style'],
                             potential_file=self.pot['file'],
                             element=self.pot['element'], tag='relaxed')
        os.system("cp  ../{}  .".format(self.pot['file']))

    def run_lmp_gsf(self, loc_dir):
        os.chdir(loc_dir)
        os.system("lmp_mpi -in {}".format(self.infile))
        os.chdir(os.pardir)

    def loop_md_gsf(self):
        for i in range(0, self.sample_gsf_num):
            dir_name = 'dir-x-%03d-%s' % (i, self.mgsf)
            self.run_lmp_gsf(dir_name)

    def multi_thread_gsf(self):
        dir_list = glob.glob("dir-*")
        num_threads = len(dir_list)
        pool = Pool(processes=num_threads)
        pool.map(unwrap_self_run_lammps,
                 zip([self] * num_threads, dir_list))

    def collect_gsf_energy(self):
        disp_list, energylist, area_list = [], [], []
        data = np.ndarray([self.sample_gsf_num, 5])
        for i in range(0, self.sample_gsf_num):
            dir_name = 'dir-x-%03d-%s' % (i, self.mgsf)
            if os.path.isdir(dir_name):
                os.chdir(dir_name)
                super_cell = self.md_get_cell()
                data[i, 0] = i
                data[i, 1] = i * self.disp_delta
                data[i, 2] = self.cal_poscar_xy_area(super_cell)
                data[i, 3] = np.loadtxt('out.txt')
                os.chdir(os.pardir)
        data[:, 4] = (data[:, 3] - data[0, 3]) / data[:, 2]
        np.savetxt('gsf.{}.dat'.format(self.mgsf), data, fmt="%.6f")

    def plot_md_gsf(self, delta=None, energy=None, filename='md_gsf.png'):
        self.set_keys()
        self.set_111plt((8, 4))
        self.ax.plot(delta, energy,
                     label="$displacement-energy$", **self.pltkwargs)
        self.fig.savefig(filename, **self.figsave)

    def plot_multi_gsf_curv(self, potlist, typelist, fname='gsf_compare.png'):
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
            self.ax2.set_xlabel("normalized displacement along $[{}]$".format(self.mgsf[:3]),
                                {'fontsize': (self.myfontsize - 3)})   # (110): -110  (11-2) -110
        self.set_tick_size(*axlist)
        self.add_y_labels(ylabiter, *axlist)
        self.add_legends(*axlist)
        self.fig.savefig(fname, **self.figsave)

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

    def md_single_dir_gsf(self):
        atoms = self.gn_gsf_atoms()
        perf_cells = copy.deepcopy(atoms.get_cell())
        for i in range(0, self.sample_gsf_num):
            dir_name = 'dir-x-%03d-%s' % (i, self.mgsf)
            self.mymkdir(dir_name)
            os.chdir(dir_name)

            disp_vector = [i * self.disp_delta, 0, 0]
            disp_matrix_direct = self.gn_displacement(
                atoms.copy(), disp_vector)
            disp_matrix = copy.deepcopy(disp_matrix_direct)

            cell_length_x = perf_cells[0, 0]
            disp_matrix[:, 0] = disp_matrix_direct[:, 0] * cell_length_x

            local_atoms = atoms.copy()
            local_atoms.translate(disp_matrix)

            self.write_lmp_config_data(local_atoms)
            os.system("cp ../{} .".format(self.infile))
            # os.system("cp lmp_init.txt ../lmp_init_{0:03d}.txt".format(i))
            os.chdir(os.pardir)

    def drv_cmp(self):
        potlist = ['adp', 'pbe']
        typelist = ['111_211', '111_110']
        drv.plot_multi_gsf_curv(potlist, typelist)

    def drv_twopath(self):
        typelist = ['111_211', '111_110']
        drv.plot_multi_type_gsf_curv(typelist)

    def drv_relaxed(self):
        potlist = ['adp', 'pbe']
        typelist = ['111_211', '111_110']
        drv.set_relax_type('relaxed')
        for gsftype in typelist:
            drv.set_mgsf(gsftype)
            drv.md_single_dir_gsf()
            drv.multi_thread_gsf()
            drv.collect_gsf_energy()
        drv.plot_multi_gsf_curv(potlist, typelist, 'gsf_relaxed.png')

    def drv_unrelaxed(self):
        potlist = ['adp', 'pbe']
        typelist = ['111_211', '111_110']
        drv.set_relax_type('unrelaxed')
        for gsftype in typelist:
            drv.set_mgsf(gsftype)
            drv.md_single_dir_gsf()
            drv.multi_thread_gsf()
            drv.collect_gsf_energy()
        drv.plot_multi_gsf_curv(potlist, typelist, 'gsf_unrelaxed.png')


if __name__ == '__main__':
    usage = "usage:%prog [options] arg1 [options] arg2"
    parser = OptionParser(usage=usage)
    parser.add_option("-t", "--mtype", action="store",
                      type="string", dest="mtype")
    parser.add_option('-p', "--param", action="store",
                      type='string', dest="fargs")
    (options, args) = parser.parse_args()
    drv = cal_md_gsf()
    dispatcher = {'prep': drv.md_single_dir_gsf,
                  'run': drv.multi_thread_gsf,
                  'clc': drv.collect_gsf_energy,
                  'trans': drv.trans_data_format,
                  'cmp': drv.drv_cmp,
                  'relaxed': drv.drv_relaxed,
                  'unrelaxed': drv.drv_unrelaxed}

    if options.fargs is not None:
        dispatcher[options.mtype.lower()](options.fargs)
    else:
        dispatcher[options.mtype.lower()]()

#!/usr/bin/env python
# encoding: utf-8
# -*- coding: utf-8 -*-
# @Author: yang37
# @Date:   2017-06-12 17:03:43
# @Last Modified by:   chaomy
# @Last Modified time: 2018-03-06 03:11:33


import os
import glob
import shutil
import numpy as np
import matplotlib.pylab as plt
import ase
from optparse import OptionParser
from copy import deepcopy
import gsf_data


class cal_va_gsf(object):

    def prepare_vasp_inputs(self, mdir):
        self.set_pbs_type('va')
        self.set_wall_time(30)
        self.set_job_title(mdir[:3])
        self.set_nnodes(1)
        self.set_main_job("mpirun vasp")
        self.write_pbs(od=True)
        os.system("mv va.pbs {}".format(mdir))
        os.system("cp INCAR {}".format(mdir))
        os.system("cp POTCAR {}".format(mdir))
        os.system("cp KPOINTS {}".format(mdir))
        os.system("cp POSCAR {}".format(mdir))

    def gn_va_single_dir_gsf(self, mtype='relax'):
        # if self.mgsf in ['x111z110']:
        #     atoms = self.gn_bcc110()
        # elif self.mgsf in ['x111z112']:

        atoms = self.gn_gsf_atoms()
        atoms.wrap()
        perf_cells = deepcopy(atoms.get_cell())
        ase.io.write('perf_poscar', images=atoms, format='vasp')

        disps = np.linspace(0, 1.0, 21)
        disps = np.append(disps, 0.0)
        npts = len(disps)

        for i, disp in zip(range(npts), disps):
            dirname = 'dir-{}-{:4.3f}'.format(self.mgsf, disp)
            self.mymkdir(dirname)

            if self.mgsf in ['x111z110']:
                disp_vector = [disp, disp, 0]
            else:
                disp_vector = [disp, 0.0, 0.0]

            disp_matrix_direct = self.gn_displacement(
                atoms.copy(), disp_vector)
            disp_matrix = deepcopy(disp_matrix_direct)
            disp_matrix[:, 0] = disp_matrix_direct[:, 0] * perf_cells[0, 0]

            local_atoms = atoms.copy()
            local_atoms.translate(disp_matrix)
            ase.io.write('POSCAR', images=local_atoms, format='vasp')
            self.prepare_vasp_inputs(dirname)
            os.system("cp POSCAR poscar.{:03d}".format(i))

    def vasp_single_dir_gsf(self):
        atoms = self.gn_gsf_atoms()
        perf_cells = copy.deepcopy(atoms.get_cell())
        npts = 21
        delta = 1. / (npts - 1)
        for i in range(0, npts):
            mdir = 'dir-x-%03d-%s' % (i, self.mgsf)
            self.mymkdir(mdir)
            os.chdir(mdir)

            disp_vector = [i * delta, 0, 0]
            disp_matrix_direct = self.gn_displacement(atoms.copy(),
                                                      disp_vector)

            disp_matrix = copy.deepcopy(disp_matrix_direct)

            cell_length_x = perf_cells[0, 0]
            disp_matrix[:, 0] = disp_matrix_direct[:, 0] * cell_length_x

            local_atoms = atoms.copy()
            local_atoms.translate(disp_matrix)

            self.write_poscar(local_atoms)
            self.add_selective_dynamics(len(disp_matrix))

            self.prepare_vasp_inputs(mdir)
            os.system("cp POSCAR.vasp ../POSCAR%03d.vasp" % (i))

            os.chdir(os.pardir)

    def collect_vasp_gsf_energy(self):  # to be done
        disp_list, energy_list, area_list = [], [], []
        for i in range(0, self.sample_gsf_num):
            mdir = 'dir-x-%03d-%s' \
                % (i, self.mgsf)

            print "dir is", mdir
            disp_list.append(i * self.disp_delta)

            if os.path.isdir(mdir):
                os.chdir(mdir)
                energy_list.append(self.vasp_energy_stress_vol_quick()[0])
                area_list.append(self.cal_xy_area_read_poscar())
                os.chdir(os.pardir)
            else:
                energy_list.append(0.0)
                area_list.append(0.0)

        with open("DATA", 'w') as fid:
            for i in range(len(disp_list)):
                fid.write("%d  %f  %f  %f \n" %
                          (i, disp_list[i],  area_list[i], energy_list[i]))

    def loop_vasp_gsf_surface(self):
        for j in range(0, 1):
            for i in range(0, 11):
                mdir = 'dir-x-%03d-y-%03d' % (i, j)
                self.mymkdir(mdir)
                shutil.copy('POTCAR', mdir)

                os.chdir(mdir)
                disp_vector = [i * 0.05, j * 0, 0]
                print disp_vector
                disp_matrix_direct = self.gn_displacement(atoms.copy(),
                                                          disp_vector)
                disp_matrix = copy.deepcopy(disp_matrix_direct)
                x1 = self._atoms.get_cell()[0, 0]
                x2 = self._atoms.get_cell()[1, 1]

                disp_matrix[:, 0] = disp_matrix_direct[:, 0] * x1
                disp_matrix[:, 1] = disp_matrix_direct[:, 1] * x2

                print "disp_matrix is", disp_matrix
                print "Length is", len(disp_matrix)

                Localatoms = self._atoms.copy()
                Localatoms.translate(disp_matrix)

                self.add_selective_dynamics(len(disp_matrix))
                self.add_alloy(len(disp_matrix))

                os.system("cp POSCAR POSCAR-{0:03d}.vasp".format(i))
                os.chdir(os.pardir)

    def plot_gsf_data(self, tag="full"):
        data = np.loadtxt("DATA")
        fig = plt.figure()
        ax = fig.add_subplot(111)

        ax.get_yaxis().get_major_formatter().set_powerlimits((-1, 2))
        ax.get_xaxis().get_major_formatter().set_powerlimits((-1, 2))

        gsf_bcc_typelist = ['gsf[111](112)', 'gsf[111](110)', 'gsf[111](123)']
        #  area = 2.8772862656350346 * 4.6985887964826203
        #  data :   i    disp    area   energy

        disp = data[:, 1]
        area = data[:, 2]
        energy = data[:, -1]
        energy = energy / area[0]
        energy = energy - np.min(energy)

        print "disp, energy", disp, energy
        print "unstable stacking fault {}: {}".format(self.mgsf, max(energy))

        if tag.lower() == "half":
            # get the data of half then use symmetry to plt full  #
            npts = len(energy)
            fullnpts = 2 * npts - 1
            delta = 1. / (fullnpts - 1)

            energy_full = np.zeros(fullnpts)
            disp_full = np.zeros(fullnpts)

            for i in range(npts):
                energy_full[i] = energy[i]
                energy_full[fullnpts - 1 - i] = energy[i]
                disp_full[i] = delta * i
                disp_full[fullnpts - 1 - i] = 1 - delta * i

            energy = energy_full
            disp = disp_full

        else:
            ax.plot(disp, energy,
                    linestyle='--',
                    linewidth=2.0,
                    markersize=5.5,
                    marker='o',
                    label=gsf_bcc_typelist[1])

        plt.legend(loc="upper left", bbox_to_anchor=[0.055, 0.95],
                   ncol=1, shadow=True)
        plt.xlabel('frac displacement ', {'fontsize': self.fontsize})
        plt.ylabel('energy eV / A^2 ',  {'fontsize': self.fontsize})

        fig.savefig("gsf.png",
                    bbox_inches='tight', pad_inches=0.01)


if __name__ == '__main__':
    usage = "usage:%prog [options] arg1 [options] arg2"
    parser = OptionParser(usage=usage)
    parser.add_option('-t', "--mtype", action="store",
                      type="string", dest="mtype")
    parser.add_option('-p', "--param", action="store",
                      type='string', dest="fargs")
    (options, args) = parser.parse_args()

    # lattice W    3.1711
    #         Nb   3.3224040
    gsf111_211 = {'type': '111_211', 'kpoints': [23, 21, 1]}
    gsf111_110 = {'type': '111_110', 'kpoints': [23, 11, 1]}

    ingsf = gsf111_211
    Job = cal_va_gsf()
    # cal total lists of given shift direction #
    if options.mtype.lower() == "single":
        Job.vasp_single_dir_gsf()

    if options.mtype == "clc":
        Job.collect_vasp_gsf_energy()

    if options.mtype == "plt":
        Job.plot_gsf_data()

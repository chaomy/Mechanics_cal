#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-06-28 00:35:14
# @Last Modified by:   chaomy
# @Last Modified time: 2018-05-25 20:35:36


import os
import copy
import numpy as np
import matplotlib.pylab as plt
from optparse import OptionParser
import md_pot_data

try:
    import gn_config
    import get_data
    import gn_kpoints
    import gn_incar
    import gn_pbs
    import Intro_vasp

except ImportError:
    print("error during import")


class cal_gsf(gn_config.bcc,
              get_data.get_data,
              gn_kpoints.gn_kpoints,
              gn_incar.gn_incar,
              gn_pbs.gn_pbs,
              Intro_vasp.vasp_change_box):

    def __init__(self,
                 in_structure='bcc',
                 gsf_surface_type='100',
                 lattice_constant=3.17,
                 in_element='W',
                 in_kpoints=[18, 18, 1]):

        gn_kpoints.gn_kpoints.__init__(self)
        get_data.get_data.__init__(self)
        gn_incar.gn_incar.__init__(self)
        gn_pbs.gn_pbs.__init__(self)
        Intro_vasp.vasp_change_box.__init__(self)

        self.exe = "mpirun vasp"
        self._pot = md_pot_data.dft_data.Nb_pbe
        self._gsf_surface_type = gsf_surface_type
        self.figsize = (8, 5)
        self.fontsize = 16

        self._surface_element = self._pot['element']
        self._surface_lattice_constant = self._pot['lattice']
        self._surface_kpoints = in_kpoints
        self._structure = self._pot['structure']

        if self.structure == 'bcc':
            gn_config.bcc.__init__(self, self._pot)

        self.root_dir = os.getcwd()
        self.sample_gsf_num = 21
        self.disp_delta = 1. / (self.sample_gsf_num - 1)

    def set_gsf_surface_type(self, gsf_surface_type):
        self._gsf_surface_type = gsf_surface_type

    def gn_gsf_atoms(self):
        if self._gsf_surface_type == '100_100':
            self._surface_direction = [[1, 0, 0],
                                       [0, 1, 0],
                                       [0, 0, 1]]
            atoms = self.set_bcc_convention(in_direction=self._surface_direction,
                                            in_size=(1, 1, 18))
            for i in range(8):
                atoms.pop()
            return atoms

        elif self._gsf_surface_type == '110_110':
            self._surface_direction = [[-1, 1, 0],
                                       [0, 0, 1],
                                       [1, 1, 0]]
            atoms = self.set_bcc_convention(in_direction=self._surface_direction,
                                            in_size=(1, 1, 14))
            for i in range(12):
                atoms.pop()
            return atoms

        elif self._gsf_surface_type == '112_111':
            self._surface_direction = [[1,  1, -2],
                                       [-1, 1,  0],
                                       [1,  1,  1]]
            atoms = self.set_bcc_convention(in_direction=self._surface_direction,
                                            in_size=(1, 1, 14))
            for i in range(24):
                atoms.pop()
            return atoms

        elif self._gsf_surface_type == '110_111':
            self._surface_direction = [[1,  1,  0],
                                       [-1, 1,  2],
                                       [1,  1,  1]]
            atoms = self.set_bcc_convention(in_direction=self._surface_direction,
                                            in_size=(1, 1, 14))
            for i in range(24):
                atoms.pop()
            return atoms

        elif self._gsf_surface_type == '110_211':
            self._surface_direction = [[1,  1,  0],
                                       [-1, 1,  2],
                                       [1,  1,  1]]
            atoms = self.set_bcc_convention(in_direction=self._surface_direction,
                                            in_size=(1, 1, 12))
            for i in range(12):
                atoms.pop()
            return atoms

        ########################################################
        #   shift along  [111] direction  common in bcc
        ########################################################
        elif self._gsf_surface_type == '111_110':  # kpionts   17  8   1
            self._surface_direction = [[-1, 1,  1],
                                       [1, -1,  2],
                                       [1,  1,  0]]
            ### default is (1, 1, 12)
            atoms = self.set_bcc_convention(in_direction=self._surface_direction,
                                            in_size=(1, 1, 14))
            for i in range(48):
                atoms.pop()
            return atoms

        elif self._gsf_surface_type == '111_211':  # kpoints   17  4  1
            self._surface_direction = [[-1,  1,  1],
                                       [0, -1,  1],
                                       [2,  1,  1]]
            atoms = self.set_bcc_convention(in_direction=self._surface_direction,
                                            in_size=(1, 1, 14))
            # default   is 1, 1, 12 , and pop 12
            for i in range(12):
                atoms.pop()
            return atoms

    def gn_displacement(self, atoms, displacement_vector):
        positions = atoms.get_positions()
        atom_num = len(positions)
        displacement = copy.deepcopy(positions)

        cut = 0.5 * np.max(positions[:, 2])
        for i in range(atom_num):
            print(positions[i, 2])
            if positions[i, 2] < cut:
                displacement[i] = [0, 0, 0]
            else:
                displacement[i] = displacement_vector
        return displacement

    def cal_gsf_given_dis(self, given_disp):
        atoms = self.gn_gsf_atoms()
        perf_cells = copy.deepcopy(atoms.get_cell())
        dir_name = 'dir-x-%.3f-%s' % (given_disp, self._gsf_surface_type)

        if not os.path.isdir(dir_name):
            os.mkdir(dir_name)
        os.chdir(dir_name)

        disp_vector = [given_disp, 0, 0]
        disp_matrix_direct = self.gn_displacement(atoms.copy(),
                                                  disp_vector)
        print(disp_matrix_direct)

        disp_matrix = copy.deepcopy(disp_matrix_direct)

        cell_length_x = perf_cells[0, 0]
        disp_matrix[:, 0] = disp_matrix_direct[:, 0] * cell_length_x

        local_atoms = atoms.copy()
        local_atoms.translate(disp_matrix)
        # add small perturbation #
        local_atoms = self.add_perturbation(local_atoms, 1.0)
        self.write_poscar(local_atoms)

        os.chdir(self.root_dir)
        return

    def vasp_single_dir_gsf(self):
        atoms = self.gn_gsf_atoms()
        perf_cells = copy.deepcopy(atoms.get_cell())
        npts = 21
        delta = 1. / (npts - 1)
        for i in range(0, npts):
            dir_name = 'dir-x-%03d-%s' % (i, self._gsf_surface_type)
            if not os.path.isdir(dir_name):
                os.mkdir(dir_name)
            os.chdir(dir_name)

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

            self.prepare_vasp_inputs(dir_name)
            os.system("cp POSCAR.vasp ../POSCAR%03d.vasp" % (i))
            os.chdir(self.root_dir)

    def collect_vasp_gsf_energy(self):  # to be done
        disp_list, energy_list, area_list = [], [], []
        for i in range(0, self.sample_gsf_num):
            dir_name = 'dir-x-%03d-%s' \
                % (i, self._gsf_surface_type)

            print("dir is", dir_name)
            disp_list.append(i * self.disp_delta)

            if os.path.isdir(dir_name):
                os.chdir(dir_name)
                energy_list.append(self.vasp_energy_stress_vol_quick()[0])
                area_list.append(self.cal_xy_area_read_poscar())
                os.chdir(self.root_dir)
            else:
                energy_list.append(0.0)
                area_list.append(0.0)

        with open("DATA", 'w') as fid:
            for i in range(len(disp_list)):
                fid.write("%d  %f  %f  %f \n"
                          % (i, disp_list[i],  area_list[i], energy_list[i]))
        return

    def plot_gsf_data(self, tag="full"):
        data = np.loadtxt("DATA")
        fig = plt.figure(figsize=self.figsize)
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

        print("disp, energy", disp, energy)
        print("unstable stacking fault {}: {}".format(self._gsf_surface_type, max(energy)))

        if tag.lower() == "half":
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
        return


usage = "usage:%prog [options] arg1 [options] arg2"
parser = OptionParser(usage=usage)
parser.add_option("-t",
                  "--mtype",
                  action="store",
                  type="string",
                  dest="mtype",
                  help="",
                  default="prp_r")
(options, args) = parser.parse_args()

if __name__ == '__main__':
    # lattice W    3.1711
    #         Nb   3.3224040
    gsf111_211 = {'type': '111_211', 'kpoints': [23, 21, 1]}
    gsf111_110 = {'type': '111_110', 'kpoints': [23, 11, 1]}
    ingsf = gsf111_211

    Job = cal_gsf(in_structure='bcc',
                  gsf_surface_type=ingsf['type'],
                  lattice_constant=3.3224040,
                  in_element='Nb',
                  in_kpoints=ingsf['kpoints'])

    #cal total lists of given shift direction # 
    if options.mtype.lower() == "single":
        Job.vasp_single_dir_gsf()

    # cal given strain at given direction #
    if options.mtype.lower() == "given":
        Job.cal_gsf_given_dis(0.25)

    if options.mtype == "clc":
        Job.collect_vasp_gsf_energy()

    if options.mtype == "plt":
        Job.plot_gsf_data()

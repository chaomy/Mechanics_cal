#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-06-28 00:35:14
# @Last Modified by:   chaomy
# @Last Modified time: 2017-07-09 18:27:31


from optparse import OptionParser
from copy import deepcopy
import os
import numpy as np
import md_pot_data
import gn_qe_inputs
import ase.io
import gsf_data
import gn_config
import get_data
import gn_kpoints
import gn_incar
import gn_pbs
import Intro_vasp


class cal_gsf(gn_config.bcc,
              get_data.get_data,
              gn_kpoints.gn_kpoints,
              gn_incar.gn_incar,
              gn_pbs.gn_pbs,
              gn_qe_inputs.gn_qe_infile,
              Intro_vasp.vasp_change_box):

    def __init__(self, mgsf='x111z112'):
        self.pot = md_pot_data.qe_pot.vca_W75Re25
        self.mgsf = mgsf

        gn_kpoints.gn_kpoints.__init__(self)
        get_data.get_data.__init__(self)
        gn_incar.gn_incar.__init__(self)
        gn_pbs.gn_pbs.__init__(self)
        Intro_vasp.vasp_change_box.__init__(self, self.pot)
        gn_config.bcc.__init__(self, self.pot)
        gn_qe_inputs.gn_qe_infile.__init__(self, self.pot)

        self.rootdir = os.getcwd()
        self.sample_gsf_num = 21
        self.disp_delta = 1. / (self.sample_gsf_num - 1)
        return

    def gn_gsf_atoms(self):
        mgsf = self.mgsf
        atoms = self.set_bcc_convention(in_direction=gsf_data.gsfbase[mgsf],
                                        in_size=gsf_data.gsfsize[mgsf])
        for i in range(gsf_data.gsfpopn[mgsf]):
            atoms.pop()
        return atoms

    def gn_displacement(self, atoms,
                        displacement_vector):

        positions = atoms.get_positions()
        atom_num = len(positions)
        displacement = deepcopy(positions)
        cut = 0.5 * np.max(positions[:, 2])
        for i in range(atom_num):
            if positions[i, 2] < cut:
                displacement[i] = [0, 0, 0]
            else:
                displacement[i] = displacement_vector
        return displacement

    def setup_qe_scf(self):
        self.set_cal_type('scf')
        self.set_ecut('43')
        self.set_degauss('0.03D0')
        self.set_thr('1.0D-5')
        self.set_kpnts((17, 9, 1))
        self.set_maxseconds(3600 * 70)
        return

    def gn_qe_single_dir_gsf(self):
        atoms = self.gn_gsf_atoms()
        perf_cells = deepcopy(atoms.get_cell())
        npts = 8
        disprng = [0.35, 0.65]
        delta = (disprng[1] - disprng[0]) / npts
        self.setup_qe_scf()
        for i in range(npts):
            disp = disprng[0] + i * delta
            dirname = 'dir-{}-{:4.3f}'.format(self.mgsf, disp)
            self.mymkdir(dirname)
            os.chdir(dirname)
            disp_vector = [disp, 0, 0]
            disp_matrix_direct = self.gn_displacement(atoms.copy(),
                                                      disp_vector)
            disp_matrix = deepcopy(disp_matrix_direct)
            disp_matrix[:, 0] = disp_matrix_direct[:, 0] * perf_cells[0, 0]
            local_atoms = atoms.copy()
            local_atoms.translate(disp_matrix)

            self.gn_qe_scf(local_atoms)
            self.set_pbs(dirname)
            os.system("cp $POTDIR/{} . ".format(self.pot['file']))
            ase.io.write('poscar', images=local_atoms, format='vasp')
            os.system("mv poscar ../poscar.{:03d}".format(i))
            os.chdir(self.rootdir)
        return

    def set_pbs(self, dirname):
        self.set_nnodes(2)
        self.set_ppn(12)
        self.set_job_title("{}".format(dirname))
        self.set_wall_time(80)
        self.set_main_job("""mpirun  qe.x < qe.in > qe.out""")
        self.write_pbs(od=True)
        return

    def gn_infile_gsf_atoms(self, atoms=None, fname='qe.in'):
        self.set_cal_type('relax')
        self.set_ecut('43')
        self.set_degauss('0.03D0')
        self.set_thr('1.0D-5')
        self.set_maxseconds(3600 * 80)
        with open(fname, 'w') as fid:
            fid = self.qe_write_control(fid, atoms)
            fid = self.qe_write_system(fid, atoms)
            fid = self.qe_write_electrons(fid)
            fid = self.qe_write_cell(fid, atoms.get_cell())
            fid = self.qe_write_species(fid, atoms, self.pot)
            fid = self.qe_write_pos(fid, atoms)
            fid = self.qe_write_kpts(fid, (5, 5, 1))
            fid.close()
        return

    def gn_qe_gsf_given_dis(self, disp):
        atoms = self.gn_gsf_atoms()
        perf_cells = copy.deepcopy(atoms.get_cell())
        dirname = 'dir-{}-{4.3f}' % (self.mgsf, disp)
        self.mymkdir(dirname)
        disp_vector = [disp, 0, 0]
        disp_matrix_direct = self.gn_displacement(atoms.copy(),
                                                  disp_vector)

        disp_matrix = copy.deepcopy(disp_matrix_direct)
        disp_matrix[:, 0] = disp_matrix_direct[:, 0] * perf_cells[0, 0]

        local_atoms = atoms.copy()
        local_atoms.translate(disp_matrix)

        self.gn_qe_scf(local_atoms)
        os.system("cp $POTDIR/{} . ".format(self.pot['file']))

        ase.io.write('poscar', images=local_atoms, format='vasp')
        os.system("mv poscar ../poscar.{:03d}".format(disp))
        os.chdir(self.rootdir)
        return

    def collect_qe_gsf_energy(self):
        disp_list, energy_list, area_list = [], [], []
        for i in range(0, self.sample_gsf_num):
            dir_name = 'dir-x-%03d-%s' \
                % (i, self._gsf_surface_type)

            print "dir is", dir_name
            disp_list.append(i * self.disp_delta)

            if os.path.isdir(dir_name):
                os.chdir(dir_name)
                energy_list.append(self.vasp_energy_stress_vol_quick()[0])
                area_list.append(self.cal_xy_area_read_poscar())
                os.chdir(self.rootdir)
            else:
                energy_list.append(0.0)
                area_list.append(0.0)

        with open("DATA", 'w') as fid:
            for i in range(len(disp_list)):
                fid.write("%d  %f  %f  %f \n"
                          % (i, disp_list[i], area_list[i], energy_list[i]))
        return


if __name__ == '__main__':
    usage = "usage:%prog [options] arg1 [options] arg2"
    parser = OptionParser(usage=usage)
    parser.add_option("-t", "--mtype",
                      action="store",
                      type="string",
                      dest="mtype")
    (options, args) = parser.parse_args()
    gsf111_211 = {'type': 'x111z112', 'kpoints': [23, 21, 1]}
    gsf111_110 = {'type': 'x111z110', 'kpoints': [23, 11, 1]}
    ingsf = gsf111_211

    drv = cal_gsf(mgsf=ingsf['type'])
    dispatcher = {'prep': drv.gn_qe_single_dir_gsf}
    dispatcher[options.mtype.lower()]()

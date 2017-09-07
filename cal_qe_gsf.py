#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-06-28 00:35:14
# @Last Modified by:   chaomy
# @Last Modified time: 2017-09-06 14:05:28


from optparse import OptionParser
from copy import deepcopy
import os
import numpy as np
import plt_drv
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
import cal_sub
import cal_qe_gsf_pos
import cal_qe_gsf_pre


class cal_gsf(gn_config.bcc,
              get_data.get_data,
              gn_kpoints.gn_kpoints,
              gn_incar.gn_incar,
              gn_pbs.gn_pbs,
              cal_sub.subjobs,
              plt_drv.plt_drv,
              gn_qe_inputs.gn_qe_infile,
              cal_qe_gsf_pos.cal_qe_gsf_pos,
              Intro_vasp.vasp_change_box,
              cal_qe_gsf_pre.cal_qe_gsf_pre):

    def __init__(self,
                 pot=md_pot_data.qe_pot.pbe_w,
                 mgsf='x111z110'):
        self.pot = pot
        self.mgsf = mgsf
        cal_sub.subjobs.__init__(self)
        gn_kpoints.gn_kpoints.__init__(self)
        get_data.get_data.__init__(self)
        gn_incar.gn_incar.__init__(self)
        gn_pbs.gn_pbs.__init__(self)
        plt_drv.plt_drv.__init__(self)
        Intro_vasp.vasp_change_box.__init__(self, self.pot)
        gn_config.bcc.__init__(self, self.pot)
        gn_qe_inputs.gn_qe_infile.__init__(self, self.pot)
        cal_qe_gsf_pos.cal_qe_gsf_pos.__init__(self)
        cal_qe_gsf_pre.cal_qe_gsf_pre.__init__(self)

        self.rootdir = os.getcwd()
        self.sample_gsf_num = 21
        self.disp_delta = 1. / (self.sample_gsf_num - 1)
        self.vcapots = {
            'WTa50': md_pot_data.qe_pot.vca_W50Ta50,
            'WRe00': md_pot_data.qe_pot.pbe_w,
            'WRe05': md_pot_data.qe_pot.vca_W95Re05,
            'WRe10': md_pot_data.qe_pot.vca_W90Re10,
            'WRe15': md_pot_data.qe_pot.vca_W85Re15,
            'WRe20': md_pot_data.qe_pot.vca_W80Re20,
            'WRe25': md_pot_data.qe_pot.vca_W75Re25,
            'WRe50': md_pot_data.qe_pot.vca_W50Re50}
        return

    def gn_gsf_atoms(self):
        mgsf = self.mgsf
        atoms = self.set_bcc_convention(
            in_direction=gsf_data.gsfbase[mgsf],
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
        self.set_ecut('41')
        self.set_degauss('0.03D0')
        self.set_thr('1.0D-6')
        self.set_kpnts(gsf_data.gsfkpts[self.mgsf])
        self.set_maxseconds(3600 * 24)
        return

    def setup_qe_relax(self):
        self.set_cal_type('relax')
        self.set_ecut('40')
        self.set_degauss('0.03D0')
        self.set_thr('1.0D-6')
        self.set_kpnts(gsf_data.gsfkpts[self.mgsf])
        self.set_maxseconds(3600 * 100)
        return

    def gn_qe_single_dir_gsf(self, mtype='relax'):
        atoms = self.gn_gsf_atoms()
        atoms.wrap()
        perf_cells = deepcopy(atoms.get_cell())
        ase.io.write('perf_poscar', images=atoms, format='vasp')

        # original  npts = 5
        # disps = np.linspace(0.42, 0.58, npts)
        # disps = np.append(disps, 0.0)

        disps = np.arange(0.02, 1.0, 0.04)
        disps = np.append(disps, 0.0)
        npts = len(disps)

        if mtype in ['scf']:
            self.setup_qe_scf()
        elif mtype in ['relax']:
            self.setup_qe_relax()

        for i, disp in zip(range(npts), disps):
            dirname = 'dir-{}-{:4.3f}'.format(self.mgsf, disp)
            self.mymkdir(dirname)
            os.chdir(dirname)
            disp_vector = [disp, 0, 0]
            disp_matrix_direct = self.gn_displacement(atoms.copy(),
                                                      disp_vector)
            disp_matrix = deepcopy(disp_matrix_direct)
            disp_matrix[:, 0] = disp_matrix_direct[:, 0] * perf_cells[0, 0]
            print disp_matrix
            local_atoms = atoms.copy()
            local_atoms.translate(disp_matrix)

            if mtype in ['scf']:
                self.gn_qe_scf_tf(local_atoms)
            elif mtype in ['relax']:
                self.gn_qe_relax_tf(local_atoms, 'xy')

            self.set_pbs(dirname)
            os.system("cp $POTDIR/{} . ".format(self.pot['file']))
            ase.io.write('poscar', images=local_atoms, format='vasp')
            self.write_lmp_config_data(local_atoms)
            os.system("mv poscar ../poscar.{:03d}".format(i))
            os.chdir(self.rootdir)
        return

    def check_gsf(self):
        for key in gsf_data.gsfsize:
            atomss = self.set_bcc_convention(
                in_direction=gsf_data.gsfbase[key],
                in_size=gsf_data.gsfsize[key])

            atomsb = self.set_bcc_convention(
                in_direction=gsf_data.gsfbase[key],
                in_size=gsf_data.bulksize[key])

            print(key, len(atomss) - len(atomsb),
                  gsf_data.gsfpopn[key])
        return

    def set_pbs(self, dirname):
        self.set_nnodes(1)
        self.set_ppn(12)
        self.set_job_title("{}".format(dirname))
        self.set_wall_time(110)
        self.set_main_job("""mpirun pw.x < qe.in > qe.out""")
        self.write_pbs(od=False)
        return

    def cal_usf(self):
        data = np.loadtxt('gsf.dat')
        gsf = (data[:, 3] - data[-1, 3]) / (2 * data[:, 2])
        usf = np.max(gsf)
        print("usf = {} eV/A^2".format(usf))
        return

    def plt_gsf(self):
        data = np.loadtxt('gsf.dat')
        print data
        gsf = (data[:, 3] - np.min(data[:, 3])) / (2 * data[:, 2])
        print gsf
        self.set_111plt()
        self.ax.plot(data[:, 1], gsf,
                     label='gsf', **next(self.keysiter))
        self.fig.savefig('fig_gsf.png', **self.figsave)
        return


if __name__ == '__main__':
    usage = "usage:%prog [options] arg1 [options] arg2"
    parser = OptionParser(usage=usage)
    parser.add_option('-t', "--mtype", action="store",
                      type="string", dest="mtype")
    parser.add_option('-p', "--param", action="store",
                      type='string', dest="fargs")
    (options, args) = parser.parse_args()

    drv = cal_gsf()
    dispatcher = {'prep': drv.gn_qe_single_dir_gsf,
                  'loopprep': drv.loop_pot_gsf,
                  'clcengy': drv.clc_qe_gsf_engy,
                  'loopclc': drv.loop_clcenergy,
                  'usf': drv.cal_usf,
                  'setpbs': drv.loop_set_pbs,
                  'sub': drv.loop_sub,
                  'plt': drv.plt_gsf,
                  'plttol': drv.plt_tol,
                  'loopplt': drv.loop_plt_gsf,
                  'trans': drv.transdata,
                  'movedir': drv.move_dirs}

    if options.fargs is not None:
        dispatcher[options.mtype.lower()](options.fargs)
    else:
        dispatcher[options.mtype.lower()]()

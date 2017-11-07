#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-06-28 00:35:14
# @Last Modified by:   chaomy
# @Last Modified time: 2017-11-02 18:00:28


import gn_lmp_infile
import numpy as np
import plt_drv
import md_pot_data
import gn_qe_inputs
import gsf_data
import gn_config
import get_data
import gn_kpoints
import gn_incar
import gn_pbs
import Intro_vasp
import cal_sub
import cal_qe_gsf
import cal_qe_gsf_pos
import cal_qe_gsf_pre
import cal_md_gsf
import cal_va_gsf
from copy import deepcopy
from math import sqrt
from optparse import OptionParser
import ase.io
import ase.lattice.orthorhombic as otho
import ase.lattice.cubic as cubic

# x [sqrt(2) a ; y = a; z = sqrt(2) a]
class othoBccp110Factory(otho.SimpleOrthorhombicFactory):
    bravais_basis = [[0.0, 0.0, 0.0],
                     [0.5, 0.5, 0.0],
                     [0.5, 0.0, 0.5],
                     [0.0, 0.5, 0.5]]
othoBccp110 = othoBccp110Factory()


class cal_gsf(gn_config.bcc,
              get_data.get_data,
              gn_kpoints.gn_kpoints,
              gn_incar.gn_incar,
              gn_pbs.gn_pbs,
              plt_drv.plt_drv,
              gn_qe_inputs.gn_qe_infile,
              cal_sub.subjobs,
              cal_qe_gsf_pos.cal_qe_gsf_pos,
              cal_qe_gsf_pre.cal_qe_gsf_pre,
              cal_qe_gsf.cal_qe_gsf,
              Intro_vasp.vasp_change_box,
              gn_lmp_infile.gn_md_infile,
              cal_md_gsf.cal_md_gsf,
              cal_va_gsf.cal_va_gsf):

    def __init__(self,
                 pot=md_pot_data.qe_pot.vca_W95Ta05, 
                 mgsf='x111z112'):
        self.pot = pot
        self.mgsf = mgsf
        self.sample_gsf_num = 21
        self.disp_delta = 1. / (self.sample_gsf_num - 1)

        cal_sub.subjobs.__init__(self)
        gn_kpoints.gn_kpoints.__init__(self)
        get_data.get_data.__init__(self)
        gn_incar.gn_incar.__init__(self)
        gn_pbs.gn_pbs.__init__(self)
        plt_drv.plt_drv.__init__(self)
        # config
        gn_config.bcc.__init__(self, self.pot)
        Intro_vasp.vasp_change_box.__init__(self, self.pot)
        # lmp
        gn_lmp_infile.gn_md_infile.__init__(self, self.pot)
        cal_md_gsf.cal_md_gsf.__init__(self)
        # qe
        gn_qe_inputs.gn_qe_infile.__init__(self, self.pot)
        cal_qe_gsf.cal_qe_gsf.__init__(self)
        cal_qe_gsf_pos.cal_qe_gsf_pos.__init__(self)
        cal_qe_gsf_pre.cal_qe_gsf_pre.__init__(self)
        cal_va_gsf.cal_va_gsf.__init__(self)
        return

    def gn_bcc110(self):
        atoms = othoBccp110(latticeconstant=(self.pot['latbcc'] * sqrt(2),
                                             self.pot['latbcc'],
                                             self.pot['latbcc'] * sqrt(2)),
                            size=(1, 1, 15),
                            symbol=self.pot['element'])
        print atoms
        for i in range(12):
          atoms.pop() 
        ase.io.write('poscar', images=atoms, format='vasp')
        return atoms

    def gn_gsf_atoms(self):
        mgsf = self.mgsf
        atoms = self.set_bcc_convention(
            in_direction=gsf_data.gsfbase[mgsf],
            in_size=gsf_data.gsfsize[mgsf])
        for i in range(gsf_data.gsfpopn[mgsf]):
            atoms.pop()
        return atoms

    def set_pbs(self, dirname, opt='qe'):
        self.set_nnodes(2)
        self.set_ppn(12)
        self.set_job_title("{}".format(dirname))
        self.set_wall_time(110)
        if opt in ['qe']:
            self.set_main_job("""mpirun pw.x < qe.in > qe.out""")
        if opt in ['va']:
            self.set_main_job("""mpirun vasp > vasp.log""")
        if self.pot in [md_pot_data.qe_pot.vca_W80Ta20, 
                        md_pot_data.qe_pot.vca_W85Ta15]: 
          self.write_pbs(od=True)
        else:
          self.write_pbs(od=False)
        return

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

    def cal_usf(self):
        data = np.loadtxt('gsf.dat')
        gsf = (data[:, 3] - data[-1, 3]) / (2 * data[:, 2])
        usf = np.max(gsf)
        print("usf = {} eV/A^2".format(usf))
        return

    def plt_gsf(self):
        data = np.loadtxt('gsf.dat')
        print data
        gsf = (data[:, 3] - np.min(data[:, 3])) / (data[:, 2])
        coeff = 16.021766208
        print gsf * coeff
        self.set_111plt()
        self.ax.plot(data[:, 1], gsf, label='gsf', **next(self.keysiter))
        self.fig.savefig('fig_gsf.png', **self.figsave)
        return

    def wrap_prep(self, opt='qe'):
        if opt in ['qe']:
            self.gn_qe_single_dir_gsf()
        if opt in ['md']:
            self.md_single_dir_gsf()
        return

    def wrap_loop(self, opt='md'):
        if opt in ['md']:
            self.loop_md_gsf()
        if opt in ['plt']:
            self.loop_plt_gsf()
        if opt in ['clc']:
            self.loop_clcenergy()
        if opt in ['prep']:
            self.loop_pot_gsf()
        if opt in ['sub']:
            self.loop_sub()


if __name__ == '__main__':
    usage = "usage:%prog [options] arg1 [options] arg2"
    parser = OptionParser(usage=usage)
    parser.add_option('-t', "--mtype", action="store",
                      type="string", dest="mtype")
    parser.add_option('-p', "--param", action="store",
                      type='string', dest="fargs")
    (options, args) = parser.parse_args()
    drv = cal_gsf()
    dispatcher = {'prep': drv.wrap_prep,
                  'loop': drv.wrap_loop,
                  'clc': drv.clc_qe_gsf_engy,
                  'usf': drv.cal_usf,
                  'setpbs': drv.loop_set_pbs,
                  'plt': drv.plt_gsf,
                  'plttol': drv.plt_tol,
                  'trans': drv.transdata,
                  'movedir': drv.move_dirs,
                  'bcc110': drv.gn_bcc110}

    if options.fargs is not None:
        dispatcher[options.mtype.lower()](options.fargs)
    else:
        dispatcher[options.mtype.lower()]()

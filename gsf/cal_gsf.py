#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-06-28 00:35:14
# @Last Modified by:   chaomy
# @Last Modified time: 2018-03-23 14:45:32


import gn_lmp_infile
import plt_drv
import md_pot_data
import gn_qe_inputs
import gsf_data
import gn_config
import get_data
import gn_pbs
from utils import Intro_vasp
from gsf import cal_qe_gsf
from gsf import cal_qe_gsf_pos
from gsf import cal_qe_gsf_pre
from gsf import cal_md_gsf
from gsf import cal_va_gsf
import numpy as np
from copy import deepcopy
from math import sqrt
from optparse import OptionParser
from itertools import cycle
import ase.io
import ase.lattice.orthorhombic as otho
import ase.lattice.cubic as cubic
import os

# x [sqrt(2) a ; y = a; z = sqrt(2) a]


class othoBccp110Factory(otho.SimpleOrthorhombicFactory):
    bravais_basis = [[0.0, 0.0, 0.0],
                     [0.5, 0.5, 0.0],
                     [0.5, 0.0, 0.5],
                     [0.0, 0.5, 0.5]]
othoBccp110 = othoBccp110Factory()


class cal_gsf(gn_config.gnStructure,
              get_data.get_data,
              gn_pbs.gn_pbs,
              plt_drv.plt_drv,
              gn_qe_inputs.gn_qe_infile,
              cal_qe_gsf_pos.cal_qe_gsf_pos,
              cal_qe_gsf_pre.cal_qe_gsf_pre,
              cal_qe_gsf.cal_qe_gsf,
              Intro_vasp.vasp_change_box,
              gn_lmp_infile.gn_md_infile,
              cal_md_gsf.cal_md_gsf,
              cal_va_gsf.cal_va_gsf):

    def __init__(self, pot=md_pot_data.va_pot.Nb_pbe, mgsf='x111z112'):
        self.pot = self.load_data("../BASICS/pot.dat")
        self.mgsf = mgsf
        self.sample_gsf_num = 21
        self.disp_delta = 1. / (self.sample_gsf_num - 1)
        get_data.get_data.__init__(self)
        gn_pbs.gn_pbs.__init__(self)
        plt_drv.plt_drv.__init__(self)
        # config
        gn_config.gnStructure.__init__(self, self.pot)
        Intro_vasp.vasp_change_box.__init__(self)
        # lmp
        gn_lmp_infile.gn_md_infile.__init__(self)
        cal_md_gsf.cal_md_gsf.__init__(self)
        # qe
        gn_qe_inputs.gn_qe_infile.__init__(self, self.pot)
        cal_qe_gsf.cal_qe_gsf.__init__(self)
        cal_qe_gsf_pos.cal_qe_gsf_pos.__init__(self)
        cal_qe_gsf_pre.cal_qe_gsf_pre.__init__(self)
        cal_va_gsf.cal_va_gsf.__init__(self)

    def gn_bcc110(self):
        atoms = othoBccp110(latticeconstant=(self.pot['latbcc'] * sqrt(2),
                                             self.pot['latbcc'],
                                             self.pot['latbcc'] * sqrt(2)),
                            size=(1, 1, 15),
                            symbol=self.pot['element'])
        for i in range(12):
            atoms.pop()
        ase.io.write('poscar', images=atoms, format='vasp')
        return atoms

    def gn_gsf_atoms(self):
        mgsf = self.mgsf
        atoms = self.set_bcc_convention(
            gsf_data.gsfbase[mgsf], gsf_data.gsfsize[mgsf])
        for i in range(gsf_data.gsfpopn[mgsf]):
            atoms.pop()
        return atoms

    def gn_gsf_one_layer(self, opt="211"):
        if opt in ["211"]:  # 211 plane
            atoms = self.set_bcc_convention(
                gsf_data.gsfbase['x111z112'],
                (1, 1, 1))
        if opt in ["110"]:  # 110 plane
            atoms = othoBccp110(latticeconstant=(self.pot['latbcc'] * sqrt(2),
                                                 self.pot['latbcc'],
                                                 self.pot['latbcc'] * sqrt(2)),
                                size=(1, 1, 1), symbol=self.pot['element'])
        if opt in ["100"]:  # 100 plane
            atoms = self.set_bcc_convention()

        if opt in ["111"]:  # 111 plane
            atoms = othoBccp110(latticeconstant=(self.pot['latbcc'] * sqrt(2),
                                                 self.pot['latbcc'],
                                                 self.pot['latbcc'] * sqrt(2)),
                                size=(1, 1, 1), symbol=self.pot['element'])
        npi = 20
        npj = 15
        cn = 0
        cell = atoms.get_cell()
        for i in range(npi):
            for j in range(npj):
                mdir = 'dir_{:03}'.format(cn)
                self.mymkdir(mdir)
                trans_atoms = atoms.copy()
                dx = cell[0, 0] * i / npi
                dy = cell[1, 1] * j / npj
                # disp_matrix = np.mat([[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0],
                #                       [dx, dy, 0.0], [dx, dy, 0.0], [dx, dy, 0.0]])
                # disp_matrix = np.mat([[0.0, 0.0, 0.0], [0.0, 0.0, 0.0],
                #                       [dx, dy, 0.0], [dx, dy, 0.0]])
                disp_matrix = np.mat([[0.0, 0.0, 0.0], [dx, dy, 0.0]])
                trans_atoms.translate(disp_matrix)

                ase.io.write("POSCAR", images=trans_atoms, format="vasp")
                os.system("cp POSCAR pos{:03}".format(cn))
                os.system("mv POSCAR {}".format(mdir))
                os.system("cp va.pbs KPOINTS INCAR POTCAR {}".format(mdir))
                cn += 1

    def set_pbs(self, dirname, opt='qe'):
        self.set_nnodes(2)
        self.set_ppn(12)
        self.set_job_title("{}".format(dirname))
        self.set_wall_time(110)
        if opt in ['qe']:
            self.set_main_job("""mpirun pw.x < qe.in > qe.out""")
        if opt in ['va']:
            self.set_main_job("""mpirun vasp > vasp.log""")
        self.write_pbs(od=False)

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
        print(("usf = {} eV/A^2".format(usf)))

    def plt_gsf(self):
        dftgsfNb111z110 = [0.0, 0.0055, 0.0183, 0.0292, 0.0389, 0.0442, 0.0389,
                           0.0292, 0.0183, 0.0055, 0.0]
        dftgsfNb111z211 = [0.0, 0.0065, 0.0221, 0.0368, 0.0474, 0.0515, 0.0448,
                           0.0334, 0.0214, 0.0066, 0.0000]
        dftdisp = np.linspace(0, 1, len(dftgsfNb111z110))

        self.set_111plt()
        coeff = 16.021766208

        # for kk in ["gsf_112", "gsf_112u"]:
        kk = 'gsf.{}'.format(self.mgsf)
        data = np.loadtxt(kk + ".dat")
        gsf = (data[:, 3] - np.min(data[:, 3])) / (data[:, 2])
        self.ax.plot(data[:, 1], gsf, label=kk, **next(self.keysiter))

        self.add_x_labels(cycle(["Normalized Dispament"]), *self.axls)
        if self.mgsf in ["x111z110"]:
            self.ax.plot(dftdisp, dftgsfNb111z110,
                         label='dft', **next(self.keysiter))
            self.add_y_labels(
                cycle([r"gsf [111](110) [eV/A$^2$]"]), *self.axls)
        elif self.mgsf in ["x111z112"]:
            self.ax.plot(dftdisp, dftgsfNb111z211,
                         label='dft', **next(self.keysiter))
            self.add_y_labels(
                cycle([r"gsf [111](211) [eV/A$^2$]"]), *self.axls)
        self.add_legends(*self.axls)
        self.set_tick_size(*self.axls)
        self.fig.savefig('fig_gsf.{}.png'.format(self.mgsf), **self.figsave)

    def wrap_prep(self, opt='qe'):
        if opt in ['qe']:
            self.gn_qe_single_dir_gsf()
        if opt in ['md']:
            self.md_single_dir_gsf()
        if opt in ['va']:
            self.gn_va_single_dir_gsf()

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

    def wrap_clc(self, opt='qe'):
        if opt in ['md']:
            self.collect_gsf_energy()
        elif opt in ['qe']:
            self.clc_qe_gsf_engy()

    def wrap_auto(self, opt='md'):
        if opt in ['md']:
            for gsf in ['x111z112', 'x111z110']:
                self.mgsf = gsf
                self.md_single_dir_gsf()
                self.loop_md_gsf()
                self.collect_gsf_energy()
                self.plt_gsf()
                self.clean_gsf()

    def clean_gsf(self):
        os.system("rm -rf dir-x-*")

if __name__ == '__main__':
    usage = "usage:%prog [options] arg1 [options] arg2"
    parser = OptionParser(usage=usage)
    parser.add_option('-t', "--mtype", action="store",
                      type="string", dest="mtype")
    parser.add_option('-p', "--param", action="store",
                      type='string', dest="fargs")
    (options, args) = parser.parse_args()
    drv = cal_gsf()
    dispatcher = {'auto': drv.wrap_auto,
                  'prep': drv.wrap_prep,
                  'loop': drv.wrap_loop,
                  'clc': drv.wrap_clc,
                  'usf': drv.cal_usf,
                  'plt': drv.plt_gsf,
                  'plttol': drv.plt_tol,
                  'trans': drv.transdata,
                  'bcc110': drv.gn_bcc110,
                  'one': drv.gn_gsf_one_layer}

    if options.fargs is not None:
        dispatcher[options.mtype.lower()](options.fargs)
    else:
        dispatcher[options.mtype.lower()]()

#!/usr/bin/env python
# encoding: utf-8

###################################################################
#
# File Name : ./cal_bcc_ideal_shear.py
#
###################################################################
#
# Purpose :
#
# Creation Date :
# Last Modified : Sat Apr 22 21:22:17 2017
# Created By    : Chaoming Yang
#
###################################################################

import os
import ase
import ase.io
import ase.lattice
import numpy as np


class cal_bcc_ideal_shear_pre(object):

    def __init__(self):
        return

    def setup_qe_params(self):
        # set qe simulation setup
        self.set_degauss('0.03D0')
        self.set_ecut('45')
        self.set_kpnts((33, 33, 33))
        return

    def gn_primitive_lmps(self,
                          strain=np.mat(np.identity(3)),
                          tag='lmp'):

        alat = self.alat
        bas = np.mat([[-0.5, 0.5, 0.5],
                      [0.5, -0.5, 0.5],
                      [0.5, 0.5, -0.5]])

        if tag == 'vasp':
            va_bas = bas * strain
            cell = alat * va_bas
            atoms = ase.Atoms(self.pot['element'],
                              positions=[[0, 0, 0]],
                              cell=cell,
                              pbc=[1, 1, 1])
            ase.io.write("POSCAR", images=atoms, format='vasp')

        if tag == 'qe':
            qe_bas = bas * strain
            cell = alat * qe_bas
            atoms = ase.Atoms(self.pot['element'],
                              positions=[[0, 0, 0]],
                              cell=cell,
                              pbc=[1, 1, 1])
            self.gn_qe_scf(atoms)

        if tag == 'lmp':
            lmp_bas = bas * strain
            lmp_bas = self.lmp_change_box(lmp_bas)
            cell = alat * lmp_bas
            atoms = ase.Atoms(self.pot['element'],
                              positions=[[0., 0., 0.]],
                              cell=cell,
                              pbc=[1, 1, 1])
            self.write_lmp_config_data(atoms, 'init.txt')
        return

    def gn_shear_twin_path(self):
        e1 = 0.5 * np.array([-1, 1, 1])
        e2 = 0.5 * np.array([1, -1, 1])
        e3 = 0.5 * np.array([1, 1, -1])

        et = np.array([-1, -1, 1])
        npts = 10
        delta = 1. / npts
        oneo6 = 1. / 6.
        for i in range(npts + 1):
            xtos = i * delta
            a1 = e1 + oneo6 * xtos * et
            a2 = e2 + oneo6 * xtos * et
            a3 = e3
            cell = np.mat([a1, a2, a3])
            print cell
            atoms = ase.Atoms('Nb',
                              positions=[[0, 0, 0]],
                              cell=cell,
                              pbc=[1, 1, 1])
            ase.io.write("POSCAR_{:03d}".format(i),
                         images=atoms,
                         format='vasp')
        return

    def set_pbs(self, dirname, opt='vasp'):
        self.set_nnodes(2)
        self.set_ppn(12)
        self.set_job_title("%s" % (dirname))
        self.set_wall_time(70)
        self.set_main_job("""../cal_md_ideal_shear.py  -t  i{}
                          """.format(opt))
        self.write_pbs(od=True)
        return

    def loop_prep_vasp(self):
        npts = self.npts
        for i in range(npts):
            delta = self.delta * i
            dirname = "dir-{:03d}".format(i)
            self.mymkdir(dirname)
            os.system("echo {} > strain.txt".format(delta))
            self.copy_inputs(dirname, 'KPOINTS',
                             'INCAR', 'POTCAR', 'strain.txt')
            self.set_pbs(dirname, delta)
        return

    def loop_prep_restart(self, opt='va'):
        raw = np.mat(np.loadtxt("ishear.txt"))
        for i in range(len(raw)):
            dirname = "dir-{:03d}".format(i)
            self.mymkdir(dirname)
            np.savetxt("restart.txt", raw[i])
            if opt in ['va', 'vasp']:
                self.copy_inputs(dirname, 'KPOINTS',
                                 'INCAR', 'POTCAR', 'restart.txt')
            elif opt in ['qe']:
                os.system("mv restart.txt {}".format(dirname))
                os.system('cp $POTDIR/{}  {}'.format(self.pot['file'],
                                                     dirname))
            self.set_pbs(dirname, raw[i][0], 'qe')
        return

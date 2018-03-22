#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-06-28 00:35:14
# @Last Modified by:   chaomy
# @Last Modified time: 2018-03-20 14:36:08


import os
import ase
import ase.io
import ase.lattice
import numpy as np
from glob import glob


class cal_bcc_ideal_shear_pre(object):

    def setup_qe_params(self):
        self.set_degauss('0.03D0')
        self.set_ecut('45')
        self.set_kpnts((33, 33, 33))

    def gn_primitive_lmps(self,
                          strain=np.mat(np.identity(3)), tag='lmp'):
        alat = self.pot['lattice']
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
            print(cell)
            atoms = ase.Atoms('Nb',
                              positions=[[0, 0, 0]],
                              cell=cell, pbc=[1, 1, 1])
            ase.io.write("POSCAR_{:03d}".format(i),
                         images=atoms, format='vasp')

    def loop_set_pbs(self):
        dirlist = glob('dir-*')
        for mdir in dirlist:
            os.chdir(mdir)
            self.set_pbs(mdir)
            os.chdir(os.pardir)

    def set_pbs(self, dirname, opt='qe'):
        self.set_nnodes(1)
        self.set_ppn(12)
        self.set_job_title("%s" % (dirname))
        self.set_wall_time(50)
        self.set_main_job("""cal_md_ideal_shear.py  -t iva""".format(opt))
        self.write_pbs(od=True)

    def loop_prep_vasp(self):
        npts = self.npts
        for i in range(npts):
            delta = self.delta * i
            dirname = "dir-{:03d}".format(i)
            self.mymkdir(dirname)
            os.system("echo {} > strain.txt".format(delta))
            self.set_pbs(dirname, delta)
            self.copy_inputs(dirname, 'KPOINTS', 'INCAR',
                             'POTCAR', 'strain.txt', 'va.pbs')

    def loop_prep_restart(self, opt='qe'):
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
            self.set_pbs(dirname, 'qe')
            os.system('mv va.pbs {}'.format(dirname))

#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-06-28 00:35:14
# @Last Modified by:   chaomy
# @Last Modified time: 2017-08-31 22:48:49


import numpy as np
import md_pot_data
from itertools import cycle
from numpy import arange
from numpy import loadtxt, savetxt, ndarray
from md_pot_data import fluxdirs
import os
from glob import glob


class cal_qe_gsf_pre(object):

    def loop_set_pbs(self):
        dirlist = glob('dir-*')
        for mdir in dirlist:
            os.chdir(mdir)
            self.set_pbs(mdir)
            os.chdir(os.pardir)
        return

    def gn_infile_gsf_atoms(self, atoms=None, fname='qe.in'):
        self.set_cal_type('relax')
        self.set_ecut('43')
        self.set_degauss('0.03D0')
        self.set_thr('1.0D-5')
        self.set_maxseconds(3600 * 60)
        with open(fname, 'w') as fid:
            fid = self.qe_write_control(fid, atoms)
            fid = self.qe_write_system(fid, atoms)
            fid = self.qe_write_electrons_tf(fid)
            fid = self.qe_write_cell(fid, atoms.get_cell())
            fid = self.qe_write_species(fid, atoms, self.pot)
            fid = self.qe_write_pos(fid, atoms)
            fid = self.qe_write_kpts(fid, (5, 5, 1))
            fid.close()
        return

    def clc_qe_gsf_engy(self, fname='gsf'):
        disps = 0.0
        disps = np.append(disps, np.arange(0.30, 0.70, 0.04))
        disps = np.append(disps, 1.0)
        npts = len(disps)
        data = np.ndarray([npts, 4])
        for i, disp in zip(range(npts), disps):
            if disp == 1:
                dirname = 'dir-{}-{:4.3f}'.format(self.mgsf, 0.0)
            else:
                dirname = 'dir-{}-{:4.3f}'.format(self.mgsf, disp)
            os.chdir(dirname)
            print dirname
            # print(self.qe_get_cell())
            data[i, 0] = i
            data[i, 1] = disp
            data[i, 2] = self.cal_xy_area()
            data[i, 3] = self.qe_get_energy_stress()[0]
            os.chdir(os.pardir)
        np.savetxt('gsf.dat', data)
        return

    def loop_pot_gsf(self, tag='prep'):
        vcapots = {
            'WTa50': md_pot_data.qe_pot.vca_W50Ta50,
            'WRe00': md_pot_data.qe_pot.pbe_w,
            'WRe05': md_pot_data.qe_pot.vca_W95Re05,
            'WRe10': md_pot_data.qe_pot.vca_W90Re10,
            'WRe15': md_pot_data.qe_pot.vca_W85Re15,
            'WRe20': md_pot_data.qe_pot.vca_W80Re20,
            'WRe25': md_pot_data.qe_pot.vca_W75Re25,
            'WRe50': md_pot_data.qe_pot.vca_W50Re50}
        gsfs = ['x111z112', 'x111z110']
        for key in vcapots:
            for gsf in gsfs:
                mdir = 'Bcc_{}_gsf{}'.format(key, gsf)
                self.mymkdir(mdir)
                os.chdir(mdir)
                self.__init__(vcapots[key], gsf)
                self.gn_qe_single_dir_gsf()
                os.chdir(os.pardir)
        return

    def loop_sub(self):
        vcapots = {
            'WTa50': md_pot_data.qe_pot.vca_W50Ta50,
            'WRe00': md_pot_data.qe_pot.pbe_w,
            'WRe05': md_pot_data.qe_pot.vca_W95Re05,
            'WRe10': md_pot_data.qe_pot.vca_W90Re10,
            'WRe15': md_pot_data.qe_pot.vca_W85Re15,
            'WRe20': md_pot_data.qe_pot.vca_W80Re20,
            'WRe25': md_pot_data.qe_pot.vca_W75Re25,
            'WRe50': md_pot_data.qe_pot.vca_W50Re50}
        gsfs = ['x111z112', 'x111z110']
        for key in vcapots.keys():
            for gsf in gsfs:
                mdir = 'Bcc_{}_gsf{}'.format(key, gsf)
                if os.path.isdir(mdir):
                    os.chdir(mdir)
                    os.system('cal_qe_gsf.py -t setpbs')
                    os.system('cal_sub.py -t sub')
                    os.chdir(os.pardir)
        return

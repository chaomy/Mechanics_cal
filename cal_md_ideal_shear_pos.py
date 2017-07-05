#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-07-04 20:53:50
# @Last Modified by:   chaomy
# @Last Modified time: 2017-07-05 08:20:46


from md_pot_data import unitconv
from scipy.interpolate import InterpolatedUnivariateSpline
import os
import ase
import ase.io
import glob
import ase.lattice
import numpy as np


class cal_bcc_ideal_shear_pos(object):

    def __init__(self):
        return

    def load_input_params(self):
        if os.path.isfile('restart.txt'):
            data = np.loadtxt("restart.txt")
            delta = data[0]
            x0 = data[-5:]
            print delta
            print x0
        else:
            data = np.loadtxt("strain.txt")
            delta = data
            x0 = np.array([1., 1., 1., 0.0, 0.0])
        return (delta, x0)

    def trans_coords_to_cartisian(self, stress):
        basis = self.basis
        stress = basis * stress * basis.transpose()
        return stress

    def convert_mtx_to_vec(self, mtx):
        vect = np.zeros(6)
        vect[0], vect[1], vect[2] = mtx[0, 0], mtx[1, 1], mtx[2, 2]
        vect[3], vect[4], vect[5] = mtx[0, 1], mtx[0, 2], mtx[1, 2]
        return vect

    def qe_loop_stress(self, opt='clc'):
        npts = self.npts
        convunit = unitconv.ustress['evA3toGpa']
        data = np.ndarray([npts, 2 + 5 + 6])
        if opt == 'clc':
            for i in range(npts):
                if i != 4 and i != 18:
                    dirname = "dir-{:03d}".format(i)
                    print dirname
                    if os.path.isdir(dirname):
                        os.chdir(dirname)
                        (engy, vol, stress) = self.qe_get_energy_stress('qe.out')
                        stress = self.trans_coords_to_cartisian(np.mat(stress))
                        stress = self.convert_mtx_to_vec(stress)
                        raw = np.loadtxt("ishear.txt")
                        os.chdir(self.root)
                        # vol = vol * (unitconv.ulength['BohrtoA']**3)
                        data[i, :7] = raw
                        data[i, 7:] = stress
            np.savetxt('stress.txt', data)
        elif opt is 'stress':
            data = np.loadtxt('stress.txt')
            spl = InterpolatedUnivariateSpline(data[:, 0], data[:, 1])
            spl.set_smoothing_factor(1.5)
            splder1 = spl.derivative()
            for i in range(len(data)):
                data[i, -1] = splder1(data[i, 0]) * convunit / data[i, 2]
            print data
        return

    def va_loop_stress(self):
        npts = self.npts
        data = np.ndarray([npts, 4])
        convunit = unitconv.ustress['evA3toGpa']
        for i in range(npts):
            dirname = "dir-{:03d}".format(i)
            os.chdir(dirname)
            atoms = ase.io.read("CONTCAR", format='vasp')
            raw = np.loadtxt("ishear.txt")
            os.chdir(self.root)
            data[i, 0] = raw[0]
            data[i, 1] = raw[1]
            data[i, 2] = np.linalg.det(atoms.get_cell())
        spl = InterpolatedUnivariateSpline(data[:, 0], data[:, 1])
        spl.set_smoothing_factor(1.5)
        splder1 = spl.derivative()
        for i in range(len(data)):
            # append the stress to the last column
            data[i, -1] = splder1(data[i, 0]) * convunit / data[i, 2]
        print data
        np.savetxt("stress.txt", data)
        return

    ##########################################################
    # calculate derivative of energy to get stress
    ##########################################################
    def convert_stress(self):
        raw = np.loadtxt("ishear.txt")
        data = np.zeros((len(raw),
                         len(raw[0]) + 1))
        data[:, :-1] = raw
        convunit = unitconv.ustress['evA3toGpa']
        vol = np.zeros(len(raw))
        vperf = 0.5 * self.alat**3
        strmat = np.zeros([3, 3])
        for i in range(len(raw)):
            strmat[0, 0], strmat[1, 1], strmat[2, 2] = \
                raw[i, 2], raw[i, 3], raw[i, 4]
            strmat[1, 0], strmat[2, 0], strmat[2, 1] = \
                raw[i, 0], raw[i, 5], raw[i, 6]
            strmat = np.mat(strmat)
            vol[i] = vperf * np.linalg.det(strmat)
        tag = 'interp'
        if tag == 'interp':
            # interpolate
            spl = InterpolatedUnivariateSpline(raw[:, 0], raw[:, 1])
            spl.set_smoothing_factor(1.5)
            splder1 = spl.derivative()
            for i in range(len(raw)):
                # append the stress to the last column
                print "coeff", convunit / vol[i]
                data[i, -1] = splder1(raw[i, 0]) * convunit / vol[i]
        print data
        np.savetxt("stress.txt", data)
        return

    def prep_restart_from_log(self):
        flist = glob.glob("s*.txt")
        print flist[0]
        data = np.loadtxt(flist[0])
        data_init = np.loadtxt('restart.txt')
        data_init[1] = data[-1][-1]
        data_init[2:] = data[-1][:-1]
        np.savetxt('restart.txt', data_init)
        dirname = os.getcwd().split('/')[-1]
        self.set_pbs(dirname, 'qe')
        return data_init

    def loop_prep_restart_from_log(self):
        npts = self.npts
        data = np.ndarray([npts, 7])
        for i in range(npts):
            dirname = "dir-{:03d}".format(i)
            if os.path.isdir(dirname):
                os.chdir(dirname)
                raw = self.prep_restart_from_log()
                os.chdir(os.pardir)
            data[i, :] = raw
        np.savetxt('ishear.txt', data)
        return

    def get_qe_stress(self):
        (engy, vol, stress) = self.qe_get_energy_stress()
        basis = self.basis
        stress = np.mat(stress)
        stress = basis * stress * basis.transpose()
        print stress
        return

    # for unfinished runs (temporary)
    def get_engy(self, file):
        fid = open(file, 'r')
        raw = fid.readlines()
        fid.close()
        for line in raw[:-1]:
            if len(line.split()) == 1:
                dat = line.split()[0]
        if dat[0] == '[':
            dat = dat.split('\'')[1]
        return dat

    # for unfinished runs
    def read_ofiles(self, opt='makeup'):
        flist = glob.glob('dir-*')
        data = np.ndarray([2, len(flist)])
        if opt == 'clcengy':
            for i in range(len(flist)):
                file = flist[i]
                cnt = int(file[4:7])
                fid = open(file, 'r')
                raw = fid.readlines()
                fid.close()
                for line in raw[:-1]:
                    if len(line.split()) == 1:
                        dat = line.split()[0]
                if dat[0] == '[':
                    dat = dat.split('\'')[1]
                data[0, i] = 0.02 * cnt
                data[1, i] = dat
            print data
            np.savetxt('ishear.txt', data)

        elif opt == 'clccell':
            data = np.ndarray([3, len(flist)])
            for i in range(len(flist)):
                mdir = flist[i]
                cell = self.qe_get_cell('{}/qe.in'.format(mdir))
                data[2, i] = np.linalg.det(cell)
                os.chdir(mdir)
                data[0, i] = np.loadtxt('restart.txt')[0]
                data[1, i] = self.get_engy(glob.glob('dir-*')[0])
                os.chdir(os.pardir)
            np.savetxt('vol.txt', data)

        elif opt == 'clctmp':
            data = np.ndarray([len(flist), 3])
            for i in range(len(flist)):
                mdir = flist[i]
                print mdir
                cell = self.qe_get_cell('{}/qe.in'.format(mdir))
                sfile = glob.glob('{}/s*'.format(mdir))[0]
                data[i, 0] = np.loadtxt('{}/restart.txt'.format(mdir))[0]
                data[i, 1] = np.loadtxt('{}'.format(sfile))[-1][5]
                data[i, 2] = np.linalg.det(cell)
            np.savetxt("itmp.txt", data)

        elif opt == 'convert':
            raw = np.loadtxt('itmp.txt')
            raw = raw[raw[:, 0].argsort()]
            (nrow, ncol) = np.shape(raw)
            data = np.ndarray([nrow, ncol + 1])
            data[:, :-1] = raw
            data[:, 1] = data[:, 1] * unitconv.uengy['rytoeV']
            convunit = unitconv.ustress['evA3toGpa']
            spl = InterpolatedUnivariateSpline(data[:, 0], data[:, 1])
            spl.set_smoothing_factor(1.3)
            splder1 = spl.derivative()
            for i in range(len(data[:, 0])):
                # append the stress to the last column
                data[i, -1] = splder1(data[i, 0]) * convunit / data[i, 2]
                print data[i, -1]
            print data
            np.savetxt("stress.txt", data)
        return

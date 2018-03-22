#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-07-04 20:53:50
# @Last Modified by:   chaomy
# @Last Modified time: 2018-03-20 14:41:33


from md_pot_data import unitconv
from scipy.interpolate import InterpolatedUnivariateSpline
import ase.io
import os
import glob
import numpy as np
from md_pot_data import fluxdirs

dirtree = {'110': {
    '05': 'Bcc_QE_VCA_WRe05_ishear110',
    '10': 'Bcc_QE_VCA_WRe10_ishear110',
    '15': 'Bcc_QE_VCA_WRe15_ishear110',
    '20': 'Bcc_QE_VCA_WRe20_ishear110',
    '25': 'Bcc_QE_VCA_WRe25_ishear110',
    '50': 'Bcc_QE_VCA_WRe50_ishear110'
}, '211': {
    '05': 'Bcc_QE_VCA_WRe05_ishear211',
    '10': 'Bcc_QE_VCA_WRe10_ishear211',
    '15': 'Bcc_QE_VCA_WRe15_ishear211',
    '20': 'Bcc_QE_VCA_WRe20_ishear211',
    '25': 'Bcc_QE_VCA_WRe25_ishear211',
    '50': 'Bcc_QE_VCA_WRe50_ishear211'
}
}


class cal_bcc_ideal_shear_pos(object):

    def load_input_params(self, shtype='110'):
        if os.path.isfile('restart.txt'):
            data = np.loadtxt("restart.txt")
            delta = data[0]
            if shtype in ['110']:
                x0 = np.array([1.004333888227457832e+00,
                               1.024562147487662877e+00,
                               0.993181305132341009e+00,
                               1.079656535712334997e-04,
                               3.744713183531088452e-04])
            print(delta)
            print(x0)
        else:
            data = np.loadtxt("strain.txt")
            delta = data
            x0 = np.array([1., 1., 1., 0.0, 0.0])
        # np.savetxt('init.save', data)
        return (delta, x0)

    def transdata(self, ptype='format'):
        for i in range(20):
            mdir = 'dir-{:03}'.format(i)
            self.mymkdir(mdir)
            if ptype in ['scp']:
                fdir = fluxdirs['QE'] + \
                    'VC_WRe/{}/'.format(dirtree['110']['20'])
                os.system('scp {}/{}/qe.out {}'.format(fdir, mdir, mdir))
                os.system('scp {}/{}/qe.in {}'.format(fdir, mdir, mdir))
                os.system('scp {}/{}/*.txt {}'.format(fdir, mdir, mdir))
                print(fdir)
            elif ptype in ['format']:
                os.chdir(mdir)
                atoms = ase.io.read('qe.out', format='espresso-out')
                ase.io.write(filename='poscar', images=atoms, format='vasp')
                os.system('mv poscar ../poscar_{:03}'.format(i))
                os.chdir(os.pardir)

    def plt_check(self):
        file = glob.glob('s0.*.txt')[0]
        data = np.loadtxt(file)
        self.set_111plt()
        xx = np.linspace(0, 1, len(data))
        self.ax.plot(xx, data[:, -1])
        self.fig.savefig('fig-engy.png')

    def trans_coords_to_cartisian(self, stress):
        basis = self.basis
        stress = basis * stress * basis.transpose()
        return stress

    def convert_mtx_to_vec(self, mtx):
        vect = np.zeros(6)
        vect[0], vect[1], vect[2] = mtx[0, 0], mtx[1, 1], mtx[2, 2]
        vect[3], vect[4], vect[5] = mtx[0, 1], mtx[0, 2], mtx[1, 2]
        return vect

    def qe_loop_stress(self, opt='clc', mtype='out'):
        npts = self.npts
        if mtype in ['out']:
            data = np.ndarray([npts, 2 + 5 + 6])
        elif mtype in ['cal']:
            data = np.ndarray([npts, 7])
        if opt == 'clc':
            for i in range(npts):
                dirname = "dir-{:03d}".format(i)
                print(dirname)
                if os.path.isdir(dirname):
                    os.chdir(dirname)
                    raw = self.load_ishear_txt()
                    data[i, :7] = raw
                    # stress
                    if mtype in ['out']:
                        (engy, vol, stress) = \
                            self.qe_get_energy_stress('qe.out')
                        stress = self.trans_coords_to_cartisian(np.mat(stress))
                        data[i, 7:] = stress = self.convert_mtx_to_vec(stress)
                    os.chdir(self.root)
            if mtype in ['cal']:
                np.savetxt('ishear.txt', data)
            elif mtype in ['out']:
                np.savetxt('stress.txt', data)

    def va_loop_stress(self):
        npts = self.npts
        data = np.ndarray([npts, 2 + 5 + 6])
        for i in range(npts):
            dirname = "dir-{:03d}".format(i)
            os.chdir(dirname)
            raw = self.load_ishear_txt()
            self.get_va_stress()
            print("raw is", raw)
            data[i, :7] = raw
            data[i, 7:] = self.get_va_stress()
            os.chdir(self.root)
        np.savetxt("stress.txt", data)

    def load_sfile_txt(self):
        flist = glob.glob("s0.*.txt")
        if len(flist) >= 1:
            data = np.loadtxt(flist[0])
        indx = np.argmin(data[:, -1])
        print('min engy index', indx)
        tmp = data[indx]
        print(tmp)
        data_init = np.zeros(len(tmp) + 1)
        data_init[1] = tmp[-1]
        data_init[2:] = tmp[:-1]
        return data_init

    def load_ishear_txt(self):
        if os.path.isfile('ishear.txt'):
            raw = np.loadtxt("ishear.txt")
            return raw
        else:
            raw = self.prep_restart_from_log()
            return raw

    # used for transform all
    def convert_stress_vasp(self):
        raw = np.loadtxt("ishear.txt")
        data = np.zeros((len(raw), len(raw[0]) + 1))
        data[:, :-1] = raw
        convunit = unitconv.ustress['evA3toGpa']
        vol = np.zeros(len(raw))
        vperf = 0.5 * self.pot['lattice']**3
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
            # spl.set_smoothing_factor(0.5)
            splder1 = spl.derivative()
            for i in range(len(raw)):
                # append the stress to the last column
                print("coeff", convunit / vol[i])
                data[i, -1] = splder1(raw[i, 0]) * convunit / vol[i]
        np.savetxt("stress.txt", data)

    # used for lammps
    def convert_stress(self, mtype='qe'):
        raw = np.loadtxt("ishear.txt")
        data = np.zeros((len(raw), len(raw[0]) + 1))
        data[:, :-1] = raw
        if mtype in ['qe']:
            data[:, 1] /= unitconv.uengy['rytoeV']
        convunit = unitconv.ustress['evA3toGpa']
        vol = np.zeros(len(raw))
        vperf = 0.5 * self.pot['lattice']**3
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
            # spl.set_smoothing_factor(0.5)
            splder1 = spl.derivative()
            for i in range(len(raw)):
                # append the stress to the last column
                print("coeff", convunit / vol[i])
                data[i, -1] = splder1(raw[i, 0]) * convunit / vol[i]
        print(data)
        np.savetxt("stress.txt", data)

    def prep_restart_from_log(self):
        flist = glob.glob("s0.*.txt")
        if len(flist) >= 1:
            data = np.loadtxt(flist[0])
            data_init = np.loadtxt('restart.txt')
            data_init[1] = data[-1][-1]
            data_init[2:] = data[-1][:-1]
            print(data_init)
        else:
            data_init = np.loadtxt('restart.txt')
        np.savetxt('restart.txt', data_init)
        dirname = os.getcwd().split('/')[-1]
        self.set_pbs(dirname, 'qe')
        return data_init

    def loop_prep_restart_from_log(self):
        npts = self.npts
        data = np.ndarray([npts, 7])
        for i in range(0, npts):
            dirname = "dir-{:03d}".format(i)
            if os.path.isdir(dirname):
                os.chdir(dirname)
                print(dirname)
                raw = self.prep_restart_from_log()
                os.chdir(os.pardir)
                data[i, :] = raw
        np.savetxt('ishear.txt', data)

    def get_va_stress(self):
        basis = self.basis
        (engy, stsvec, vol) = self.vasp_energy_stress_vol()
        vaspmtx = np.mat(np.zeros([3, 3]))

        vaspmtx[0, 0] = stsvec[0]
        vaspmtx[1, 1] = stsvec[1]
        vaspmtx[2, 2] = stsvec[2]

        vaspmtx[1, 0] = stsvec[3]
        vaspmtx[2, 0] = stsvec[4]
        vaspmtx[2, 1] = stsvec[5]

        vaspmtx[0, 1] = stsvec[3]
        vaspmtx[0, 2] = stsvec[4]
        vaspmtx[1, 2] = stsvec[5]

        vaspmtx = basis * vaspmtx * basis.transpose()
        vector = self.convert_mtx_to_vec(vaspmtx)
        return vector

    def get_qe_stress(self):
        (engy, vol, stress) = self.qe_get_energy_stress()
        basis = self.basis
        stress = np.mat(stress)
        stress = basis * stress * basis.transpose()
        print(stress)

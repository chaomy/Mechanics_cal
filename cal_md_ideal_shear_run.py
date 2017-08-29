#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-06-28 00:35:14
# @Last Modified by:   chaomy
# @Last Modified time: 2017-08-28 21:24:15


from scipy.optimize import minimize
import os
import numpy as np


class cal_bcc_ideal_shear_run(object):

    def qe_relax(self):
        (delta, x0) = self.load_input_params()
        data = np.zeros(7)
        res = minimize(self.runqe, x0, delta, tol=1e-4)
        print res
        data[0] = delta
        data[1] = res.fun
        data[2:] = res.x
        np.savetxt("ishear.txt", data)
        return

    def vasp_relax(self):
        (delta, x0) = self.load_input_params()
        data = np.zeros(7)
        res = minimize(self.runvasp, x0, delta, tol=1e-4)
        print res
        data[0] = delta
        data[1] = res.fun
        data[2:] = res.x
        np.savetxt("ishear.txt", data)
        return

    def loop_shear_lmp(self):
        x0 = np.array([1.2, 1.1, 0.9, 0., 0.])
        npts = self.npts
        data = np.ndarray([npts, 7])
        for i in range(npts):
            delta = self.delta * i
            res = minimize(self.runlmp, x0, delta, options={'fatol': 1e-4})
            x0 = res.x
            print res
            data[i][0] = (delta)
            data[i][1] = res.fun
            data[i][2:] = res.x
        np.savetxt("ishear.txt", data)
        return

    def recordstrain(self, delta, x, fval):
        fid = open("s{:4.3f}.txt".format(delta), "a")
        fid.write('{} {} {} {} {} {}\n'.format(x[0], x[1], x[2],
                                               x[3], x[4], fval))
        fid.close()
        return

    def runlmp(self, x, delta):
        basis = self.basis
        # y shear toward x direction
        strain = np.mat([[x[0], 0.0, 0.0],
                         [-delta, x[1], 0.0],
                         [x[3], x[4], x[2]]])
        new_strain = basis.transpose() * strain * basis
        self.gn_primitive_lmps(new_strain, 'lmp')
        os.system("lmp_mpi -i in.init -screen  no")
        raw = np.loadtxt("out.txt")
        engy = raw[0]
        self.recordstrain(delta, x, engy)
        return engy

    def runqe(self, x, delta):
        basis = self.basis
        strain = np.mat([[x[0], 0.0, 0.0],
                         [-delta, x[1], 0.0],
                         [x[3], x[4], x[2]]])
        new_strain = basis.transpose() * strain * basis
        self.gn_primitive_lmps(new_strain, 'qe')
        os.system("mpirun pw.x < qe.in > qe.out")
        (engy, vol, stress) = self.qe_get_energy_stress('qe.out')
        self.recordstrain(delta, x, engy)
        return engy

    def runvasp(self, x, delta):
        basis = self.basis
        strain = np.mat([[x[0], 0.0, 0.0],
                         [-delta, x[1], 0.0],
                         [x[3], x[4], x[2]]])
        new_strain = basis.transpose() * strain * basis
        self.gn_primitive_lmps(new_strain, 'vasp')
        os.system("mpirun vasp > vasp.log")
        (engy, stress, vol) = self.vasp_energy_stress_vol()
        self.recordstrain(delta, x, engy)
        return engy

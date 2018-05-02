#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-06-28 00:35:14
# @Last Modified by:   chaomy
# @Last Modified time: 2018-04-29 22:11:56


from scipy.optimize import minimize
import os
import numpy as np


class cal_bcc_ideal_shear_run(object):

    def qe_relax(self):
        (delta, x0) = self.load_input_params()
        data = np.zeros(7)
        res = minimize(self.runqe, x0, delta, tol=1e-3, method='Nelder-Mead')
        print(res)
        data[0], data[1], data[2:] = delta, res.fun, res.x
        np.savetxt("ishear.txt", data)

    def vasp_relax(self):
        self.cnt = 0
        (delta, x0) = self.load_input_params()
        data = np.zeros(7)
        res = minimize(self.runvasp, x0, delta, tol=1e-3, method='Nelder-Mead')
        print(res)
        data[0], data[1], data[2:] = delta, res.fun, res.x
        np.savetxt("ishear.txt", data)

    def loop_shear_lmp(self):
        # This is a typical starting point coming from running of Nb
        # for 110
        x0 = np.array([9.995455382221108964e-01,
                       1.075933203165313046e+00,
                       1.085678403318995455e+00,
                       2.548558462285926973e-01,
                       3.139171320878977878e-02])

        # for 211
        # x0 = np.array([9.529678394166982702e-01,
        #                1.016847096128853600e+00,
        #                1.048118345240103277e+00,
        #                1.389246099431499289e-04,
        #                -2.142860303806213597e-04])

        npts = self.npts
        # 1 strain  + energy + 5 strains + 6 stresses
        data = np.ndarray([npts, 7 + 6])
        for i in range(npts):
            delta = self.delta * i
            res = minimize(self.runlmp, x0, delta,
                           tol=1e-4, method='Nelder-Mead')
            print(res)
            info = np.loadtxt("out.txt")
            data[i][0], data[i][1], data[i][2:7], data[
                i][7:] = (delta), res.fun, res.x, info[1:]
        np.savetxt("ishear.txt", data)

    def recordstrain(self, delta, x, fval):
        fid = open("s.txt", "a")
        fid.write('{} {} {} {} {} {}\n'.format(x[0], x[1], x[2],
                                               x[3], x[4], fval))
        fid.close()

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
        strain = np.mat(
            [[x[0], 0.0, 0.0], [-delta, x[1], 0.0], [x[3], x[4], x[2]]])
        new_strain = basis.transpose() * strain * basis
        self.gn_primitive_lmps(new_strain, 'vasp')
        os.system("mpirun vasp > vasp.log")
        print("cnt = ", self.cnt)
        (engy, stress, vol) = self.vasp_energy_stress_vol()
        self.recordstrain(delta, x, engy)
        os.system("mv OUTCAR outcar-{:03d}".format(self.cnt))
        self.cnt += 1
        return engy

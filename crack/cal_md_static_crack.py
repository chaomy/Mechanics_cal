#!/usr/bin/env python
# encoding: utf-8
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2018-03-27 16:28:18
# @Last Modified by:   chaomy
# @Last Modified time: 2018-03-27 16:33:12

import Intro as Intr
import os
import numpy as np
import md_pot_data
import cal_compliance_constants as CAL
import lmp_change_configs


# Nb surface_energy_list = [2.04, 2.36, 2.47]

class md_crack(lmp_change_configs.lmp_change_configs):

    def __init__(self):
        lmp_change_configs.lmp_change_configs.__init__(self)
        self.pot = md_pot_data.md_pot.Nb_eam
        self.element = self.pot['element']
        self.surfaceE = self.pot['surf110']
        self.c11, self.c12, self.c44 = \
            self.pot['c11'], self.pot['c12'], self.pot['c44']
        # self.crackcoeff = crack_coeff

    def get_base_mtx(self):
        x = np.array([1, 0, 0], "float")
        y = np.array([0, 1, 0], "float")
        z = x.copy()
        z[0] = x[1] * y[2] - y[1] * x[2]
        z[1] = x[2] * y[0] - y[2] * x[0]
        z[2] = x[0] * y[1] - y[0] * x[1]
        xx, yy, zz = 0., 0., 0.
        for i in range(3):
            xx += x[i] * x[i]
            yy += y[i] * y[i]
            zz += z[i] * z[i]
        x /= np.sqrt(xx)
        y /= np.sqrt(yy)
        z /= np.sqrt(zz)
        lij = np.mat([x, y, z], "float")
        return lij

    def change_sij(self, L, s66):
        s99 = np.zeros([9, 9], "float")
        s66new = np.zeros([6, 6], "float")
        s99new = np.zeros([9, 9], "float")
        for i in range(6):
            for j in range(6):
                s99[i, j] = s66
        for p in range(3):
            for q in range(3):
                for r in range(3):
                    for s in range(3):
                        for i in range(3):
                            for j in range(3):
                                for k in range(3):
                                    for l in range(3):
                                        if i == 0 and j == 0:
                                            m = 0
                                        if k == 0 and l == 0:
                                            n = 0
                                        if p == 0 and q == 0:
                                            u = 0
                                        if r == 0 and s == 0:
                                            v = 0

                                        if i == 1 and j == 1:
                                            m = 1
                                        if k == 1 and l == 1:
                                            n = 1
                                        if p == 1 and q == 1:
                                            u = 1
                                        if r == 1 and s == 1:
                                            v = 1

                                        if i == 2 and j == 2:
                                            m = 2
                                        if k == 2 and l == 2:
                                            n = 2
                                        if p == 2 and q == 2:
                                            u = 2
                                        if r == 2 and s == 2:
                                            v = 2

                                        if i == 1 and j == 2:
                                            m = 3
                                        if k == 1 and l == 2:
                                            n = 3
                                        if p == 1 and q == 2:
                                            u = 3
                                        if r == 1 and s == 2:
                                            v = 3

                                        if i == 2 and j == 0:
                                            m = 4
                                        if k == 2 and l == 0:
                                            n = 4
                                        if p == 2 and q == 0:
                                            u = 4
                                        if r == 2 and s == 0:
                                            v = 4

                                        if i == 0 and j == 1:
                                            m = 5
                                        if k == 0 and l == 1:
                                            n = 5
                                        if p == 0 and q == 1:
                                            u = 5
                                        if r == 0 and s == 1:
                                            v = 5

                                        if i == 2 and j == 1:
                                            m = 6
                                        if k == 2 and l == 1:
                                            n = 6
                                        if p == 2 and q == 1:
                                            u = 6
                                        if r == 2 and s == 1:
                                            v = 6

                                        if i == 0 and j == 2:
                                            m = 7
                                        if k == 0 and l == 2:
                                            n = 7
                                        if p == 0 and q == 2:
                                            u = 7
                                        if r == 0 and s == 2:
                                            v = 7

                                        if i == 1 and j == 0:
                                            m = 8
                                        if k == 1 and l == 0:
                                            n = 8
                                        if p == 1 and q == 0:
                                            u = 8
                                        if r == 1 and s == 0:
                                            v = 8
                                        s99new[u, v] += L[p, i] * L[q, j] * \
                                            L[r, k] * L[s, l] * s99[m, n]
        for i in range(6):
            for j in range(6):
                s66new[i, j] = s99new[i, j]
            print(s99new)
            print(s66new)
        return s66new

    def cal_sij(self):
        M = CAL.cal_para()
        s11, s12, s44 = M.cal_complience(self.c11, self.c12, self.c44)
        sij = np.empty([6, 6])
        sij[0, 0], sij[1, 1], sij[2, 2] = s11, s11, s11
        sij[0, 1], sij[0, 2], sij[1, 2], sij[1, 0], sij[
            2, 0], sij[2, 1] = s12, s12, s12, s12, s12, s12
        sij[3, 3], sij[4, 4], sij[5, 5] = s44, s44, s44
        print(sij)
        # use for coordination trnasformation
        # sij = self.change_sij(L,sij)
        return sij

    def cal_scalarB(self):
        sij = self.cal_sij()
        b11 = (sij[0, 0] * sij[2, 2] - sij[0, 2] * sij[0, 2]) / sij[2, 2]
        b22 = (sij[1, 1] * sij[2, 2] - sij[1, 2] * sij[1, 2]) / sij[2, 2]
        b12 = (sij[0, 1] * sij[2, 2] - sij[0, 2] * sij[1, 2]) / sij[2, 2]
        b16 = (sij[0, 5] * sij[2, 2] - sij[0, 2] * sij[2, 5]) / sij[2, 2]
        b26 = (sij[1, 5] * sij[2, 2] - sij[1, 2] * sij[2, 5]) / sij[2, 2]
        b66 = (sij[5, 5] * sij[2, 2] - sij[2, 5] * sij[2, 5]) / sij[2, 2]
        BB = np.sqrt(0.5 * b11 * b22 * (np.sqrt(b22 / b11) +
                                        (2 * b12 + b66) / (2 * b11)))
        Gg = 2 * self.surfaceE
        ToMpaSqrtM = np.power(10, -1.5)
        self.crackcoeff.Kg = ToMpaSqrtM * np.sqrt(Gg / BB)
        self.crackcoeff.KK = self.crackcoeff.Kg

    def get_coeffs(self):
        sij = self.cal_sij()
        b11 = (sij[0, 0] * sij[2, 2] - sij[0, 2] * sij[0, 2]) / sij[2, 2]
        b22 = (sij[1, 1] * sij[2, 2] - sij[1, 2] * sij[1, 2]) / sij[2, 2]
        b66 = (sij[5, 5] * sij[2, 2] - sij[2, 5] * sij[2, 5]) / sij[2, 2]
        b12 = (sij[0, 1] * sij[2, 2] - sij[0, 2] * sij[2, 1]) / sij[2, 2]
        b16 = (sij[0, 5] * sij[2, 2] - sij[0, 2] * sij[2, 5]) / sij[2, 2]
        b26 = (sij[1, 5] * sij[2, 2] - sij[1, 2] * sij[2, 5]) / sij[2, 2]
        Eq = np.poly1d([b11, -2 * b16, 2 * b12 + b66, -2 * b26, b22])
        print(Eq)
        u1, u3, u2, u4 = Eq.r
        print(u1, u2, u3, u4)
        p1 = b11 * u1**2 - b16 * u1 + b12
        p2 = b11 * u2**2 - b16 * u2 + b12
        q1 = b12 * u1 - b26 + b22 / u1
        q2 = b12 * u2 - b26 + b22 / u2
        self.crackcoeff.q1, self.crackcoeff.q2 = q1, q2
        self.crackcoeff.p1, self.crackcoeff.p2 = p1, p2
        self.crackcoeff.u1, self.crackcoeff.u2 = u1, u2

    def stat_crack(self):
        # execuable = "mpirun lmp_linux -in"
        os.system("python gnStructure.py")
        self.cal_scalarB()
        self.get_coeffs()
        for i in range(2, 3):
            self.crackcoeff.KK = self.crackcoeff.Kg + 0.005 * i
            self.Intro_Crack('cfg', self.crackcoeff)
            # os.system("%s in.read" % (execuable))
            # self.rename_cfg((K))
            # os.system("rm crackxyz/*")

    def md_crack(self):
        execuable = "mpirun lmp_mpi -in"
        os.system("python gnStructure.py")
        self.cal_scalarB()
        self.get_coeffs()
        self.Intro_Crack('cfg', self.crackcoeff)

    def stat_crack_continue(self):
        execuable = "mpirun lmp_linux -in"
        M = Intr.MD_ChangeBox()
        p1, p2, q1, q2, u1, u2 = self.get_coeffs()
        Krestart = 2.4589812
        for i in range(1, 120):
            K = Krestart + 0.007 * i
            M.Intro_Crack_xyz(K, p1, p2, q1, q2, u1, u2)
            os.system("rm crackxyz/*")
            os.system("%s in.read" % (execuable))
            self.rename_cfg(K)

    def Test_init(self):
        M = Intr.MD_ChangeBox()
        Kg = self.cal_scalarB()
        self.get_coeffs()
        for i in range(50):
            K = Kg + 0.01 * i
            M.Intro_Crack_xyz(self.crackcoeff)
            os.system("cp ./Crack.txt Init/Crack_%g" % (0.01 * K))


if __name__ == "__main__":
    N = MD_crack()
    N.stat_crack()

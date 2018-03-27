#!/usr/bin/env python
# encoding: utf-8
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-07-05 08:11:49
# @Last Modified by:   chaomy
# @Last Modified time: 2018-03-27 16:02:40


import numpy as np
from numpy.linalg import inv


class crack_coeff:
    Gg = None
    BB = None
    K1 = None
    Kg = None
    p1 = None
    p2 = None
    q1 = None
    q2 = None
    u1 = None
    u2 = None

    lij = None
    sij = None
    bij_pstrain = {'b11': None, 'b12': None,
                   'b22': None, 'b16': None,
                   'b26': None, 'b66': None}


class md_crack_uti(object):

    def get_crack_coeffs(self):
        return (self.ckcoeff.u1, self.ckcoeff.u2, self.ckcoeff.p1,
                self.ckcoeff.p2, self.ckcoeff.q1, self.ckcoeff.q2)

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
            # print s99new
            # print s66new
        return s66new

    def cubic_cal_complience(self, c11, c12, c44):
        M = np.zeros([6, 6], "float")
        M[0, 0], M[1, 1], M[2, 2] = c11, c11, c11
        M[0, 1], M[0, 2], M[1, 0] = c12, c12, c12
        M[1, 2], M[2, 0], M[2, 1] = c12, c12, c12
        M[3, 3], M[4, 4], M[5, 5] = c44, c44, c44
        N = inv(M)
        s11 = N[0, 0]  # divide(1,float(Youngs_Modulus[i]))
        s12 = N[0, 1]  # divide((-Possion_ratio[i]),Youngs_Modulus[i])
        s44 = N[4, 4]  # divide(1.,float(Shear_Modulus[i]))
        return s11, s12, s44

    def cal_sij(self):
        s11, s12, s44 = self.cubic_cal_complience(self.c11, self.c12, self.c44)
        sij = np.empty([6, 6])
        sij[0, 0], sij[1, 1], sij[2, 2] = s11, s11, s11
        sij[0, 1], sij[0, 2], sij[1, 2], sij[1, 0], sij[
            2, 0], sij[2, 1] = s12, s12, s12, s12, s12, s12
        sij[3, 3], sij[4, 4], sij[5, 5] = s44, s44, s44
        # print sij
        # use for coordination trnasformation
        # sij = self.change_sij(L,sij)
        return sij

    def set_plane_strain_bij(self):
        sij = self.elastic.Sij
        self.ckcoeff.bij_pstrain['b11'] = (
            sij[0, 0] * sij[2, 2] - sij[0, 2] * sij[0, 2]) / sij[2, 2]
        self.ckcoeff.bij_pstrain['b22'] = (
            sij[1, 1] * sij[2, 2] - sij[1, 2] * sij[1, 2]) / sij[2, 2]
        self.ckcoeff.bij_pstrain['b12'] = (
            sij[0, 1] * sij[2, 2] - sij[0, 2] * sij[1, 2]) / sij[2, 2]
        self.ckcoeff.bij_pstrain['b16'] = (
            sij[0, 5] * sij[2, 2] - sij[0, 2] * sij[2, 5]) / sij[2, 2]
        self.ckcoeff.bij_pstrain['b26'] = (
            sij[1, 5] * sij[2, 2] - sij[1, 2] * sij[2, 5]) / sij[2, 2]
        self.ckcoeff.bij_pstrain['b66'] = (
            sij[5, 5] * sij[2, 2] - sij[2, 5] * sij[2, 5]) / sij[2, 2]

    def get_plane_strain_bij(self):
        return (self.ckcoeff.bij_pstrain['b11'],
                self.ckcoeff.bij_pstrain['b22'],
                self.ckcoeff.bij_pstrain['b12'],
                self.ckcoeff.bij_pstrain['b16'],
                self.ckcoeff.bij_pstrain['b26'],
                self.ckcoeff.bij_pstrain['b66'])

    def get_scalarB(self, surfE):
        b11, b22, b12, b16, b26, b66 = self.get_plane_strain_bij()
        tmp1 = np.sqrt(0.5 * b11 * b22)
        tmp2 = np.sqrt(np.sqrt(b22 / b11) + (2 * b12 + b66) / (2 * b11))
        # BB = tmp1 * tmp2
        BB = np.sqrt(0.5 * b11 * b22 * (np.sqrt(b22 / b11) +
                                        (2 * b12 + b66) / (2 * b11)))
        if surfE is None:
            Gg = 2 * self.surfe110
        else:
            Gg = 2 * surfE
        self.ckcoeff.Kg = np.sqrt(Gg / BB * 1e-3)
        self.ckcoeff.K1 = self.ckcoeff.Kg

    def get_coeffs(self):
        b11, b22, b12, b16, b26, b66 = self.get_plane_strain_bij()
        Eq = np.poly1d([b11, -2 * b16, 2 * b12 + b66, -2 * b26, b22])
        u1, u3, u2, u4 = Eq.r
        # print u1, u2, u3, u4
        p1 = b11 * u1**2 - b16 * u1 + b12
        p2 = b11 * u2**2 - b16 * u2 + b12
        q1 = b12 * u1 - b26 + b22 / u1
        q2 = b12 * u2 - b26 + b22 / u2
        self.ckcoeff.q1, self.ckcoeff.q2 = q1, q2
        self.ckcoeff.p1, self.ckcoeff.p2 = p1, p2
        self.ckcoeff.u1, self.ckcoeff.u2 = u1, u2

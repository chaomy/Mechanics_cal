#!/usr/bin/env python
# encoding: utf-8
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-07-05 08:11:49
# @Last Modified by:   chaomy
# @Last Modified time: 2018-03-27 23:19:58


from numpy import sqrt, sin, cos, abs
import axes_check
import numpy as np
import tool_elastic_constants
from utils import stroh_solve
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

    def loop_211(self):
        self.aniso_211()
        K1 = self.ckcoeff.K1
        for i in range(20):
            self.ckcoeff.K1 = K1 + 0.05 * i
            print(self.ckcoeff.K1)
            atoms = self.intro_crack_k1(atoms=self.gn_perf_plate())
            self.write_lmp_config_data(atoms, 'crack_{:03d}.txt'.format(i))

    def aniso_211(self):
        c = tool_elastic_constants.elastic_constants(
            C11=self.pot['c11'],
            C12=self.pot['c12'],
            C44=self.pot['c44'])
        e1 = [1, 0, 0]
        e2 = [0, 1, -1]
        e3 = [0, 1, 1]
        axes = np.array([e1, e2, e3])

        burgers = self.pot["lattice"] / 2. * np.array([1., 1., -1.])
        burgers = axes_check.axes_check(axes).dot(burgers)
        stroh = stroh_solve.Stroh(c, burgers, axes=axes)

        A = np.mat(np.zeros([3, 3]), dtype='complex')
        A[:, 0] = np.mat(stroh.A[0]).transpose()
        A[:, 1] = np.mat(stroh.A[2]).transpose()
        A[:, 2] = np.mat(stroh.A[4]).transpose()

        B = np.mat(np.zeros([3, 3]), dtype='complex')
        B[:, 0] = np.mat(stroh.L[0]).transpose()
        B[:, 1] = np.mat(stroh.L[2]).transpose()
        B[:, 2] = np.mat(stroh.L[4]).transpose()

        Linv = np.real(np.complex(0, 1) * A * np.linalg.inv(B))
        Gamma = 0.5 * Linv

        self.ckcoeff.K1 = sqrt(
            2 * self.pot["surf110"] / abs(Gamma[1, 1] * 1e3))
        # theta = np.deg2rad(54.735610317245346)
        theta = np.deg2rad(35.2643896828)
        omega = np.mat([[cos(theta), sin(theta), 0.0],
                        [-sin(theta), cos(theta), 0.0],
                        [0.0, 0.0, 1.0]])
        Gamma = omega * Gamma * omega.transpose()
        Gamma = np.abs(np.linalg.inv(Gamma))
        Gamma = Gamma * 1e9  # Pa
        print(self.ckcoeff.K1)

        if axes is not None:
            burgers = axes_check.axes_check(axes).dot(burgers)
            c = c.transform(axes)

        G00 = Gamma[0, 0]
        usf = 0.67405172613
        k1e = np.sqrt(G00 * usf) * 1e-6
        print(self.ckcoeff.K1, k1e)

        self.set_plane_strain_bij(c.Sij)
        # self.get_scalarB(self.pot["surf110"]) get the same k1c
        self.get_coeffs()
        # atoms = self.intro_crack_k1(atoms=self.gn_perf_plate())
        # self.write_lmp_config_data(atoms, 'crack.txt')

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
        s11, s12, s44 = self.cubic_cal_complience(self.pot['c11'],
                                                  self.pot['c12'],
                                                  self.pot['c44'])
        sij = np.empty([6, 6])
        sij[0, 0], sij[1, 1], sij[2, 2] = s11, s11, s11
        sij[0, 1], sij[0, 2], sij[1, 2], sij[1, 0], sij[
            2, 0], sij[2, 1] = s12, s12, s12, s12, s12, s12
        sij[3, 3], sij[4, 4], sij[5, 5] = s44, s44, s44
        # print sij
        # use for coordination trnasformation
        # sij = self.change_sij(L,sij)
        return sij

    def set_plane_strain_bij(self, sij):
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

    # Gets the critical K1c
    def get_scalarB(self, surfE):
        b11, b22, b12, b16, b26, b66 = self.get_plane_strain_bij()
        # tmp1 = np.sqrt(0.5 * b11 * b22)
        # tmp2 = np.sqrt(np.sqrt(b22 / b11) + (2 * b12 + b66) / (2 * b11))
        # BB = tmp1 * tmp2
        BB = np.sqrt(0.5 * b11 * b22 * (np.sqrt(b22 / b11) +
                                        (2 * b12 + b66) / (2 * b11)))
        Gg = 2 * surfE
        self.ckcoeff.Kg = np.sqrt(Gg / BB * 1e-3)
        print("Kg is ", self.ckcoeff.Kg)

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

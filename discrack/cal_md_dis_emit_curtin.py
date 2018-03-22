#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-06-25 14:28:58
# @Last Modified by:   chaomy
# @Last Modified time: 2017-11-09 02:01:02


import axes_check
from numpy import cos, sin, sqrt
from collections import OrderedDict
import numpy as np
import tool_elastic_constants
from utils import stroh_solve


matconsts = OrderedDict([('Al1', {'lat': 4.05,
                                  'ugsf': 0.167,
                                  'c11': 113.4,
                                  'c12': 61.5,
                                  'c44': 31.6,
                                  'surf': 0.871}),
                         ('Al2', {'lat': 4.03,
                                  'ugsf': 0.119,
                                  'c11': 118.0,
                                  'c12': 62.2,
                                  'c44': 36.7,
                                  'surf': 0.871}),
                         ('Gold', {'lat': 4.08,
                                   'ugsf': 0.097,
                                   'c11': 183.2,
                                   'c12': 158.7,
                                   'c44': 45.3,
                                   'surf': 0.796}),
                         ('Silver', {'lat': 4.09,
                                     'ugsf': 0.119,
                                     'c11': 129.1,
                                     'c12': 91.7,
                                     'c44': 56.7,
                                     'surf': 0.619})])


class cal_dis_emit_curtin(object):

    def loop_curtin(self):
        for key in list(matconsts.keys())[:]:
            # self.get_cutin_result_mat_k2e(matconsts[key])
            self.get_cutin_result_mat_k1e(matconsts[key])
        return

    def get_cutin_result_mat_k1e(self, param):
        c = tool_elastic_constants.elastic_constants(
            C11=param['c11'],
            C12=param['c12'],
            C44=param['c44'])

        axes = np.array([[-1, -1, 2],
                         [1, 1, 1],
                         [-1, 1, 0]])
        burgers = param['lat'] / np.sqrt(6.) * np.array([-1., -1., 2.])

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

        v1 = np.array([-1, -1, 2])
        v2 = np.array([1, 1, 2])
        theta = (np.arccos(np.dot(v1, v2) /
                           (np.linalg.norm(v1) * np.linalg.norm(v2))))
        omega = np.mat([[cos(theta), sin(theta), 0.0],
                        [-sin(theta), cos(theta), 0.0],
                        [0.0, 0.0, 1.0]])
        Gamma = omega * Gamma * omega.transpose()
        # phi = 0
        # svect = np.mat(np.array([cos(phi), 0.0, sin(phi)]))
        Gamma = np.linalg.inv(Gamma)
        Gamma = np.abs(Gamma)

        Gamma00 = Gamma[0, 0]  # in GPa
        Gamma11 = Gamma[1, 1]

        Gamma00 *= 1e9
        Gamma11 *= 1e9

        if axes is not None:
            T = axes_check.axes_check(axes)
            burgers = T.dot(burgers)
            c = c.transform(axes)

        (coeff, Kg) = self.cal_crack(param, c)
        print((cos(0.5 * theta))**2 * sin(0.5 * theta))
        print(coeff)

        k1e = sqrt(Gamma00 * param['ugsf']) * 1e-6   # Mpa
        print('k2e', k1e)
        print('k2e', k1e / coeff)
        # / Fmat2[1], k2e / cos(theta)
        return


    # def get_cutin_result_mat_k2e(self, param):
    #     c = tool_elastic_constants.elastic_constants(
    #         C11=param['c11'],
    #         C12=param['c12'],
    #         C44=param['c44'])

    #     # A
    #     axes = np.array([[1, 1, 2],
    #                      [1, 1, -1],
    #                      [-1, 1, 0]])

    #     burgers = param['lat'] / np.sqrt(6.) * np.array([1., 1., 2.])

    #     stroh = stroh_solve.Stroh(c, burgers, axes=axes)
    #     A = np.mat(np.zeros([3, 3]), dtype='complex')
    #     A[:, 0] = np.mat(stroh.A[0]).transpose()
    #     A[:, 1] = np.mat(stroh.A[2]).transpose()
    #     A[:, 2] = np.mat(stroh.A[4]).transpose()

    #     B = np.mat(np.zeros([3, 3]), dtype='complex')
    #     B[:, 0] = np.mat(stroh.L[0]).transpose()
    #     B[:, 1] = np.mat(stroh.L[2]).transpose()
    #     B[:, 2] = np.mat(stroh.L[4]).transpose()

    #     # print A * B.transpose() + B * A.transpose()
    #     Linv = np.real(np.complex(0, 1) * A * np.linalg.inv(B))
    #     Gamma = 0.5 * Linv

    #     theta = np.deg2rad(0.0)
    #     omega = np.mat([[cos(theta), sin(theta), 0.0],
    #                     [-sin(theta), cos(theta), 0.0],
    #                     [0.0, 0.0, 1.0]])
    #     Gamma = np.abs(omega * Gamma * omega)

    #     # phi = 0
    #     # svect = np.mat(np.array([cos(phi), 0.0, sin(phi)]))
    #     Gamma = np.linalg.inv(Gamma)
    #     Gamma00 = Gamma[0, 0]  # in GPa
    #     Gamma11 = Gamma[1, 1]
    #     print 'Gamma', Gamma00, Gamma11   # GPa
    #     Gamma00 *= 1e9
    #     Gamma11 *= 1e9

    #     # cal F(theta)
    #     soltrace = np.mat([[stroh.p[0], 0.0, 0.0],
    #                        [0.0, stroh.p[2], 0.0],
    #                        [0.0, 0.0, stroh.p[4]]], dtype='complex')

    #     Fmat1 = soltrace / np.sqrt(np.cos(theta) + soltrace * np.sin(theta))
    #     Fmat1 = np.real(B * Fmat1 * np.linalg.inv(B))

    #     Fmat2 = 1.0 / np.sqrt(np.cos(theta) + soltrace * np.sin(theta))
    #     Fmat2 = np.real(B * Fmat2 * np.linalg.inv(B))

    #     # print Gamma
    #     k2e = np.sqrt(Gamma00 * param['ugsf']) * 1e-6   # Mpa
    #     print 'k2e', k2e
    #     return
   #       np.rad2deg(np.arccos(np.dot(v1, v2)/( np.linalg.norm(v1) * np.linalg.norm(v2) )))

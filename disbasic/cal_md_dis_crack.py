#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2018-05-01 22:40:40
# @Last Modified by:   chaomy
# @Last Modified time: 2018-05-02 16:36:22

from crack import cal_md_crack_ini
from collections import OrderedDict

matconsts = OrderedDict([('Al1', {'lat': 4.05,
                                  'ugsf': 0.167,
                                  'c11': 113.4,
                                  'c12': 61.5,
                                  'c44': 31.6}),
                         ('Al2', {'lat': 4.03,
                                  'ugsf': 0.119,
                                  'c11': 118.0,
                                  'c12': 62.2,
                                  'c44': 36.7}),
                         ('Gold', {'lat': 4.08,
                                   'ugsf': 0.097,
                                   'c11': 183.2,
                                   'c12': 158.7,
                                   'c44': 45.3}),
                         ('Silver', {'lat': 4.09,
                                     'ugsf': 0.119,
                                     'c11': 129.1,
                                     'c12': 91.7,
                                     'c44': 56.7})])


class dis_init_crack(object):

    def cal_crack(self):
        drv = cal_md_crack_ini.md_crack_ini()
        drv.cal_crack_anglecoeff()

    def loop_table(self):  # repeat curtain results
        for key in list(matconsts.keys()):
            print(key)
            self.get_cutin_result(matconsts[key])

    def get_cutin_result(self, param):
        c = tool_elastic_constants.elastic_constants(
            C11=param['c11'],
            C12=param['c12'],
            C44=param['c44'])
        # A
        axes = np.array([[-1, -1, 2],
                         [1, 1, 1],
                         [-1, 1, 0]])
        burgers = param['lat'] * sqrt(2.) / 2. * np.array([-1., 1., 0])
        stroh = stroh_solve.Stroh(c, burgers, axes=axes)

        A = mat(np.zeros([3, 3]), dtype='complex')
        A[:, 0] = mat(stroh.A[0]).transpose()
        A[:, 1] = mat(stroh.A[2]).transpose()
        A[:, 2] = mat(stroh.A[4]).transpose()

        B = mat(np.zeros([3, 3]), dtype='complex')
        B[:, 0] = mat(stroh.L[0]).transpose()
        B[:, 1] = mat(stroh.L[2]).transpose()
        B[:, 2] = mat(stroh.L[4]).transpose()

        Gamma = 0.5 * np.real(np.complex(0, 1) * A * np.linalg.inv(B))
        theta = np.deg2rad(70.0)
        omega = mat([[cos(theta), sin(theta), 0.0],
                     [-sin(theta), cos(theta), 0.0],
                     [0.0, 0.0, 1.0]])
        phi = 0
        svect = mat(np.array([cos(phi), 0.0, sin(phi)]))
        usf = param['ugsf']  # J/m^2
        Gamma = np.abs(np.linalg.inv(Gamma))
        # Gamma = omega * Gamma * omega

        Gamma = (svect * Gamma * svect.transpose())[0, 0]  # in GPa
        print(Gamma)

        Gamma = Gamma * 1e9  # Pa
        ke1 = sqrt(Gamma * usf)
        print(ke1 * 1e-6)

        # A = mat(np.zeros([3, 3]), dtype='complex')
        # A[:, 0] = mat(stroh.A[1]).transpose()
        # A[:, 1] = mat(stroh.A[3]).transpose()
        # A[:, 2] = mat(stroh.A[5]).transpose()

        # B = mat(np.zeros([3, 3]), dtype='complex')
        # B[:, 0] = mat(stroh.L[1]).transpose()
        # B[:, 1] = mat(stroh.L[3]).transpose()
        # B[:, 2] = mat(stroh.L[5]).transpose()
        # Gamma = np.real(np.complex(0, 1) * A * np.linalg.inv(B))

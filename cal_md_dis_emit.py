#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-06-25 14:28:58
# @Last Modified by:   chaomy
# @Last Modified time: 2017-08-31 15:31:41


import axes_check
from numpy import cos, sin, sqrt, abs
from collections import OrderedDict
import numpy as np
import md_pot_data
import tool_elastic_constants
import stroh_solve
import cal_md_dislocation
import cal_md_crack_ini

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

vcaw = OrderedDict([('WRe00', {'lat': 3.17093,
                               'ugsf1': 1.2289,
                               'ugsf2': 1.2616,
                               'c11': 534.2,
                               'c12': 193.8,
                               'c44': 177.1,
                               'surf': 1.5205}),
                    ('WRe05', {'lat': 3.1684646,
                               'ugsf1': 1.1889,
                               'ugsf2': 1.2067,
                               'c11': 537.460016,
                               'c12': 198.714877,
                               'c44': 181.790297,
                               'surf': 1.4964}),
                    ('WRe10', {'lat': 3.16596693,
                               'ugsf1': 1.1436,
                               'ugsf2': 1.1496,
                               'c11': 536.324450,
                               'c12': 203.742842,
                               'c44': 184.820205,
                               'surf': 1.4708}),
                    ('WRe15', {'lat': 3.16346922,
                               'ugsf1': 1.0965,
                               'ugsf2': 1.0926,
                               'c11': 531.897067,
                               'c12': 208.441515,
                               'c44': 187.130407,
                               'surf': 1.4468}),
                    ('WRe20', {'lat': 3.160970015,
                               'ugsf1': 1.0488,
                               'ugsf2': 1.0371,
                               'c11': 527.776820,
                               'c12': 211.722869,
                               'c44': 189.347735,
                               'surf': 1.4227}),
                    ('WRe25', {'lat': 3.158786,
                               'ugsf1': 1.0005,
                               'ugsf2': 0.9841,
                               'c11': 528.780194,
                               'c12': 217.275650,
                               'c44': 191.632656,
                               'surf': 1.4019})])


class cal_dis_emit(object):

    def __init__(self, pot=None):
        if pot is None:
            pot = md_pot_data.qe_pot.vca_W75Re25
        self.pot = pot
        self.mddis_drv = cal_md_dislocation.md_dislocation(self.pot)

        # self.theta = np.deg2rad(70.0)
        # x [1, 0, 0], y[0, 1, -1], z[0, 1, 1]
        # angle between [0, 1, -1] and [2, 1, -1]
        self.theta = 54.735610317245346
        # self.theta = 90 - self.theta
        return

    def loop_curtin(self):
        for key in matconsts.keys()[:]:
            # self.get_cutin_result_mat_k2e(matconsts[key])
            self.get_cutin_result_mat_k1e(matconsts[key])
        return

    def loop_table(self):
        data = np.ndarray([6, 4])
        for key, i in zip(vcaw.keys(), range(6)):
            print key
            data[i, 0] = 0.05 * i
            data[i, 1:] = self.get_bcc_w_result(vcaw[key])
        np.savetxt('vcaw_112_Ke.txt', data)
        return

    def cal_crack(self, param, c):
        drv = cal_md_crack_ini.md_crack_ini()

        drv.set_params(c)
        drv.set_plane_strain_bij()
        drv.get_scalarB(param['surf'])
        drv.get_coeffs()

        s1, s2 = drv.ckcoeff.u1, drv.ckcoeff.u2
        theta = self.theta
        tmp = 1. / (s1 - s2)
        costhe = cos(theta)
        sinthe = sin(theta)
        coeff = np.real(tmp * (s1 / sqrt(costhe + s1 * sinthe) -
                               s2 / sqrt(costhe + s2 * sinthe)))
        return (coeff, drv.ckcoeff.Kg)

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

        # Gamma = omega * Gamma * omega
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

        print(cos(0.5 * theta))**2 * sin(0.5 * theta)
        print coeff

        k1e = sqrt(Gamma00 * param['ugsf']) * 1e-6   # Mpa
        print 'k2e', k1e
        print 'k2e', k1e / coeff
        # / Fmat2[1], k2e / cos(theta)
        return

    def get_bcc_w_result(self, param):
        c = tool_elastic_constants.elastic_constants(
            C11=param['c11'],
            C12=param['c12'],
            C44=param['c44'])

        e1 = [1, 0, 0]
        e2 = [0, 1, -1]
        e3 = [0, 1, 1]

        # glide plane [2, 1, -1]
        axes = np.array([e1, e2, e3])

        # x [1, 0, 0], y[0, 1, -1], z[0, 1, 1]
        burgers = param['lat'] / 2. * np.array([-1., 1., -1.])
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

        theta = np.deg2rad(54.735610317245346)

        Gamma = np.abs(np.linalg.inv(Gamma))
        Gamma = Gamma * 1e9  # Pa

        # K1c
        G11 = Gamma[1, 1]
        surf = param['surf']
        k1c = sqrt(2 * surf / G11) * 1e6

        # Ke
        if axes is not None:
            T = axes_check.axes_check(axes)
            burgers = T.dot(burgers)
            c = c.transform(axes)

        (coeff, Kg) = self.cal_crack(param, c)
        usf = param['ugsf1']  # J/m^2

        # Gamma = (svect * Gamma * svect.transpose())[0, 0]  # in GPa
        G00 = Gamma[0, 0]
        k2e = np.sqrt(G00 * usf) * 1e-6
        k1e = k2e / coeff
        print k1e, k2e, k1c
        return (k1e, k2e, k1c)

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


if __name__ == '__main__':
    drv = cal_dis_emit()
    # drv.get_cutin_result()
    drv.loop_table()
    # drv.loop_curtin()
    # drv.cal_crack()

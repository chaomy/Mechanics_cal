#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-06-25 14:28:58
# @Last Modified by:   chaomy
# @Last Modified time: 2017-11-08 15:48:58


from optparse import OptionParser
from numpy import cos, sin, sqrt
from collections import OrderedDict
import axes_check
import plt_drv
import numpy as np
import md_pot_data
import tool_elastic_constants
import stroh_solve
import cal_md_dislocation
import cal_md_crack_ini
import cal_md_dis_emit_curtin
import cal_md_dis_emit_plt


vcaw = OrderedDict([('WTa05', {'ugsf1': 2.2446,   # 
                               'ugsf2': 2.0327}), # 3.1945 
                    ('WTa10', {'ugsf1': 2.2584,   # 2.239 
                               'ugsf2': 1.9952}), # 3.157
                    ('WTa15', {'ugsf1': 2.241,    # 2.241
                               'ugsf2': 1.9472}), # 3.109 
                    ('WTa20', {'ugsf1': 2.191,    # 2.191
                               'ugsf2': 1.8815}), # 3.0433
                    ('WTa25', {'ugsf1': 2.1815,   # 
                               'ugsf2': 1.7912}), # 2.953 
                    ('WTa50', {'lat':   3.2502,
                               'ugsf1': 1.5988,
                               'ugsf2': 1.2763,   # 2.3793
                               'c11': 343.7798,
                               'c12': 187.0044,
                               'c44': 91.8360,
                               'surf': 3.941229200}),
                    ('WRe00', {'lat': 3.17093,
                               'ugsf1': 2.1203,  # 2.1147 
                               'ugsf2': 2.0213,  # 110   # 3.1831   => 1.5747786078266461 (shift = 1.1618) 
                               'c11': 534.2,
                               'c12': 193.8,
                               'c44': 177.1,
                               'surf': 4.847717}),
                    ('WRe05', {'lat': 3.1684646,
                               'ugsf1': 2.0251,  # 211
                               'ugsf2': 1.944,
                               'c11': 537.460016,
                               'c12': 198.714877,
                               'c44': 181.790297,
                               'surf': 4.7490852}),
                    ('WRe10', {'lat': 3.16596693,
                               'ugsf1': 1.9301,
                               'ugsf2': 1.8667,
                               'c11': 536.324450,
                               'c12': 203.742842,
                               'c44': 184.820205,
                               'surf': 4.6426819}),
                    ('WRe15', {'lat': 3.16346922,
                               'ugsf1': 1.8322,
                               'ugsf2': 1.7677,
                               'c11': 531.897067,
                               'c12': 208.441515,
                               'c44': 187.130407,
                               'surf': 4.5335292}),
                    ('WRe20', {'lat': 3.160970015,
                               'ugsf1': 1.7375,
                               'ugsf2': 1.6873,
                               'c11': 527.776820,
                               'c12': 211.722869,
                               'c44': 189.347735,
                               'surf': 4.42611199}),
                    ('WRe25', {'lat': 3.158786,
                               'ugsf1': 1.6469,
                               'ugsf2': 1.6102,
                               'c11': 528.780194,
                               'c12': 217.275650,
                               'c44': 191.632656,
                               'surf': 4.323485376}),
                    ('WRe50', {'lat': 3.1485,
                               'ugsf1': 1.2428,
                               'ugsf2': 1.1315,
                               'c11': 511.892,
                               'c12': 242.872,
                               'c44': 197.935,
                               'surf': 3.904985515})
                    # ('Ta1', {'lat': 3.321,
                    #          'ugsf1': 1.000,
                    #          'ugsf2': 0.840,
                    #          'c11': 247,
                    #          'c12': 170,
                    #          'c44': 67,
                    #          'surf': 2.27}),
                    # ('Ta2', {'lat': 3.306,
                    #          'ugsf1': 0.947,
                    #          'ugsf2': 0.840,
                    #          'c11': 266,
                    #          'c12': 158,
                    #          'c44': 87,
                    #          'surf': 2.49}),
                    # ('Nb', {'lat': 3.3224040,
                    #         'ugsf1': 0.772,   # 211
                    #         'ugsf2': 0.677,   # 110
                    #         'c11': 246,
                    #         'c12': 137,
                    #         'c44': 20,
                    #         'surf': 2.255})
                    ])

# ('Nbcite', {'lat': 3.309,
#             'ugsf1': 0.823,
#             'ugsf2': 0.705,
#             'c11': 251,
#             'c12': 133,
#             'c44': 22,
#             'surf': 2.36})
# ('Mo', {'lat': 3.169,
#         'ugsf1': 0.835,  # 211
#         'ugsf2': 0.719,  # 110
#         'c11': 462,
#         'c12': 163,
#         'c44': 102,
#         'surf': 3.204})])


class cal_dis_emit(cal_md_dis_emit_curtin.cal_dis_emit_curtin,
                   cal_md_dis_emit_plt.cal_md_dis_emit_plt,
                   plt_drv.plt_drv):

    def __init__(self, pot=None):
        if pot is None:
            pot = md_pot_data.qe_pot.vca_W75Re25
        self.pot = pot
        self.mddis_drv = cal_md_dislocation.md_dislocation(self.pot)
        cal_md_dis_emit_curtin.cal_dis_emit_curtin.__init__(self)
        plt_drv.plt_drv.__init__(self)
        # self.theta = np.deg2rad(70.0)
        # x [1, 0, 0], y[0, 1, -1], z[0, 1, 1]
        # angle between [0, 1, -1] and [2, 1, -1]
        self.theta = 54.735610317245346
        # self.theta = 90 - self.theta
        return

    def loop_table(self):
        npts = len(vcaw.keys())
        data = np.ndarray([npts, 4])
        for key, i in zip(vcaw.keys(), range(npts)):
            print key
            data[i, 0] = i
            data[i, 1:] = self.get_bcc_w_result110(vcaw[key])
        np.savetxt('vcaw_110_Ke.txt', data)

        for key, i in zip(vcaw.keys(), range(npts)):
            print key
            data[i, 0] = i
            data[i, 1:] = self.get_bcc_w_result211(vcaw[key])
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

    def get_bcc_w_result211(self, param):
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
        burgers = param['lat'] / 2. * np.array([1., 1., -1.])

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
        surf = param['surf']
        print Gamma[1, 1]
        k1c = sqrt(2 * surf / abs(Gamma[1, 1] * 1e3))

        theta = np.deg2rad(54.735610317245346)
        omega = np.mat([[cos(theta), sin(theta), 0.0],
                        [-sin(theta), cos(theta), 0.0],
                        [0.0, 0.0, 1.0]])
        Gamma = omega * Gamma * omega.transpose()
        Gamma = np.abs(np.linalg.inv(Gamma))
        Gamma = Gamma * 1e9  # Pa

        # K1c
        # G11 = Gamma[1, 1]
        # k1c = sqrt(2 * surf / G11) * 1e6
        # print k1c, k1cn

        # Ke
        if axes is not None:
            T = axes_check.axes_check(axes)
            burgers = T.dot(burgers)
            c = c.transform(axes)

        (coeff, Kg) = self.cal_crack(param, c)
        usf = param['ugsf1']  # J/m^2

        # Gamma = (svect * Gamma * svect.transpose())[0, 0]  # in GPa
        G00 = Gamma[0, 0]
        k1e = np.sqrt(G00 * usf) * 1e-6
        k1e = k1e / coeff
        print k1e, k1c, Kg, k1e / k1c
        return (k1e, k1c, k1e / k1c)

    def get_bcc_w_result110(self, param):
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
        burgers = param['lat'] / 2. * np.array([1., 1., 1.])

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

        # k1c
        surf = param['surf']
        print Gamma[1, 1]
        k1c = sqrt(2 * surf / abs(Gamma[1, 1] * 1e3))

        theta = np.deg2rad(90.0)
        omega = np.mat([[cos(theta), sin(theta), 0.0],
                        [-sin(theta), cos(theta), 0.0],
                        [0.0, 0.0, 1.0]])

        Gamma = omega * Gamma * omega.transpose()
        Gamma = np.abs(np.linalg.inv(Gamma))
        Gamma = Gamma * 1e9  # Pa

        # # K1c
        # G11 = Gamma[1, 1]
        # surf = param['surf']
        # k1c = sqrt(2 * surf / G11) * 1e6

        # Ke
        phi = np.deg2rad(90 - 54.735610317245346)
        svect = np.mat(np.array([cos(phi), 0.0, sin(phi)]))
        Gamma = (svect * Gamma * svect.transpose())[0, 0]  # in GPa

        if axes is not None:
            T = axes_check.axes_check(axes)
            burgers = T.dot(burgers)
            c = c.transform(axes)

        (coeff, Kg) = self.cal_crack(param, c)
        usf = param['ugsf2']  # J/m^2

        G00 = Gamma
        k1e = np.sqrt(G00 * usf) * 1e-6 / svect[0, 0]
        k1e = k1e / coeff
        print k1e, k1c, Kg, k1e / k1c
        return (k1e, k1c, k1e / k1c)


if __name__ == '__main__':
    usage = "usage:%prog [options] arg1 [options] arg2"
    parser = OptionParser(usage=usage)
    parser.add_option('-t', "--mtype", action="store",
                      type="string", dest="mtype")
    parser.add_option('-p', "--param", action="store",
                      type='string', dest="fargs")
    (options, args) = parser.parse_args()
    drv = cal_dis_emit()
    dispatcher = {'run': drv.loop_table,
                  'plt': drv.plt_k1,
                  'loopcurtin': drv.loop_curtin}
    if options.fargs is not None:
        dispatcher[options.mtype.lower()](options.fargs)
    else:
        dispatcher[options.mtype.lower()]()

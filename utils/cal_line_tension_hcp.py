#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2018-01-31 16:09:47
# @Last Modified by:   chaomy
# @Last Modified time: 2018-03-02 10:08:24

from sys import stdout
from numpy import rad2deg
from numpy import sqrt, sin, cos
from numpy import pi
from numpy import linspace
GpaToEVA3 = 1. / 160.21766208
EVA3ToGAp = 160.21766208


class cal_line_tension(object):

    def cal_hcp(self):
        C11 = 63.5
        C12 = 25.9
        C33 = 66.4
        C13 = 21.7
        C44 = 18.42

        Ks = (1 / 2. * C44 * (C11 - C12))**(1 / 2.)

        C13p = sqrt(C11 * C33)

        A = C44 * (C13p - C13) / (C33 * (C13p + C13 + 2 * C44))
        Ke = (C13p + C13) * sqrt(A)

        stdout.write("Ks is {}\n".format(Ks))
        stdout.write("Ke is {}\n".format(Ke))

        # theta = 0.25 * pi
        theta = linspace(0, 0.5, 5) * pi
        lgR2r0 = 4.0
        b = 3.2
        Etheta = 3.2**2 / (4 * pi) * (Ke * sin(theta)**2 +
                                      Ks * cos(theta)**2 + 2 * (Ke - Ks) * cos(2 * theta)) * lgR2r0
        Etheta *= GpaToEVA3
        print(rad2deg(theta))
        print(Etheta)

        b = 3.2
        apb = 0.0157589
        # cr = 2 * apb / b * EVA3ToGAp * sqrt(0.017 * 750 / Etheta) * 750 / 5000
        f = 0.015
        M = apb / b * EVA3ToGAp
        N = sqrt(3 * pi * pi * apb * f * 380 / (32 * Etheta))
        print(M * N)

        # cr =  * sqrt(0.017 * 750 / Etheta) * 750 / 5000
        # print cr


if __name__ == '__main__':
    drv = cal_line_tension()
    drv.cal_hcp()

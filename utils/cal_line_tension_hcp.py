#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2018-01-31 16:09:47
# @Last Modified by:   chaomy
# @Last Modified time: 2018-08-28 13:06:21

from sys import stdout
from numpy import rad2deg
from numpy import sqrt, sin, cos
from numpy import pi
from numpy import linspace
from itertools import cycle
GpaToEVA3 = 1. / 160.21766208
EVA3ToGAp = 160.21766208
import plt_drv


# Mg3Nd AntiPhase boundary Energies
# APB No   = -345.71742306
# APB 1|b| = -345.09996314  DeltaE (0.61745992) -> 0.01588357695 eV/A^2
# APB 2|b| = -345.31738482  DeltaE (0.40003824) -> 0.01029060828 eV/A^2

# Area = 5.2429139999999999 * 7.4146000000000001 = 38.8741101444


class cal_line_tension(plt_drv.plt_drv):

    def __init__(self):
        plt_drv.plt_drv.__init__(self)

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
        # f = 0.015
        f = 0.013
        M = apb / b * EVA3ToGAp
        N = sqrt(3 * pi * pi * apb * f * 380 / (32 * Etheta))
        print(M * N)

        # cr =  * sqrt(0.017 * 750 / Etheta) * 750 / 5000
        # print cr

    def cal_hcp_apb4(self):
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

        b = 3.2
        apb1 = 0.0158
        apb2 = 0.0103
        apb12 = apb1 - apb2
        # cr = 2 * apb / b * EVA3ToGAp * sqrt(0.017 * 750 / Etheta) * 750 / 5000
        f = 0.013
        M = apb1 / b * EVA3ToGAp
        N = sqrt(3 * pi * pi * apb1 * f * 380 / (32 * Etheta)) - f

        M2 = apb12 / b * EVA3ToGAp
        N2 = sqrt(3 * pi * pi * apb12 * f * 380 / (32 * Etheta)) - f
        print(0.25 * (M * N + M2 * N2))

        # cr =  * sqrt(0.017 * 750 / Etheta) * 750 / 5000
        # print cr

    def plot_stress_to_radius(self):
        b = 3.209
        apb = 0.01747
        # cr = 2 * apb / b * EVA3ToGAp * sqrt(0.017 * 750 / Etheta) * 750 / 5000
        f = 0.015
        Etheta = 0.294
        r = linspace(20, 400, 10)  # A
        M = apb / b * EVA3ToGAp
        N = sqrt(3 * pi * pi * apb * f * r / (32 * Etheta))
        tau = (M * N)
        self.set_111plt()
        self.ax.plot(r / b, tau, label='tau-r', **next(self.keysiter))
        self.add_y_labels(
            cycle(["tau"]), *self.axls)
        self.add_x_labels(
            cycle(["<r>"]), *self.axls)
        self.add_legends(*self.axls)
        self.set_tick_size(*self.axls)
        self.fig.savefig('fig_theory.png', **self.figsave)


if __name__ == '__main__':
    drv = cal_line_tension()
    drv.cal_hcp()
    drv.cal_hcp_apb4()
    # drv.plot_stress_to_radius()

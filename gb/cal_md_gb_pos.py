#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-07-05 08:12:30
# @Last Modified by:   chaomy
# @Last Modified time: 2018-06-11 01:25:28

from numpy import loadtxt, linspace
import os
from math import cos, sin
from math import sqrt
from itertools import cycle
import numpy as np


class md_gb_pos(object):

    def loop_plt_angle(self):

        self.find_angles_1100()
        dd = np.zeros([len(self.ag) + 2, 2])
        dd[0, 0], dd[0, 1] = 0.0, 0.0
        dd[-1, 0], dd[-1, 1] = 90.0, 0.0

        for e, i in zip(self.ag, range(len(self.ag))):
            mdir = "1100_{:.2f}".format(e[0])
            dd[i + 1, 0] = e[0]
            dd[i + 1, 1] = np.loadtxt("{}/lmp.dat".format(mdir))
        print(dd)

        self.set_111plt((9, 6))
        self.set_keys()
        self.ax.plot(dd[:, 0], 1e3 * dd[:, 1], 'o--',
                     markersize=14, label='<100> GB')
        # ddp **next(self.keysiter)

        self.add_legends(self.ax)

        ylb = cycle(['GB[001] mJ/m^2'])
        xlb = cycle(['Angle (deg)'])

        self.add_y_labels(ylb, self.ax)
        self.add_x_labels(xlb, self.ax)
        self.set_tick_size(self.ax)
        self.fig.savefig('fig_gb.png', **self.figsave)

#!/usr/bin/env python
# encoding: utf-8
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-07-05 08:12:30
# @Last Modified by:   chaomy
# @Last Modified time: 2018-02-20 22:58:18

from numpy import loadtxt, linspace
import os
from math import cos, sin
from math import sqrt
from itertools import cycle


class md_gb_pos(object):

    def loop_plt(self):
        dd = loadtxt("dat.txt")
        self.set_111plt((9, 6))
        self.set_keys()
        self.ax.plot(dd[:, 0], dd[:, -1], 'o--', markersize=14,
                     label='<100> GB')
        self.add_legends(self.ax)

        ylb = cycle(['GB[001] mJ/m^2'])
        xlb = cycle(['Angle (deg)'])

        self.add_y_labels(ylb, self.ax)
        self.add_x_labels(xlb, self.ax)
        self.set_tick_size(self.ax)
        self.fig.savefig('fig_shift.png', **self.figsave)

    def loop_plt_angle(self):
        self.loop_angle(opt='read')
        print(self.gbdat)
        dd = self.gbdat

        self.set_111plt((9, 6))
        self.set_keys()
        self.ax.plot(dd[0, :], 1e3 * dd[1, :], 'o--', markersize=14,
                     label='<100> GB')
        # ddp **next(self.keysiter)

        self.add_legends(self.ax)

        ylb = cycle(['GB[001] mJ/m^2'])
        xlb = cycle(['Angle (deg)'])

        self.add_y_labels(ylb, self.ax)
        self.add_x_labels(xlb, self.ax)
        self.set_tick_size(self.ax)
        self.fig.savefig('fig_gb.png', **self.figsave)

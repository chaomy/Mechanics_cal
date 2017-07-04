#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-06-28 00:35:14
# @Last Modified by:   chaomy
# @Last Modified time: 2017-07-04 16:30:24

from itertools import cycle
import numpy as np
import plt_drv
import os


class cal_bcc_ideal_shear_plt(plt_drv.plt_drv):

    def __init__(self):
        plt_drv.plt_drv.__init__(self)
        return

    def plt_strain_vs_energy(self, infile='ishear.txt'):
        raw = np.loadtxt(infile)
        raw = raw[raw[:, 0].argsort()]
        print raw
        self.set_111plt()
        self.ax.plot(raw[:, 0], (raw[:, 1] - raw[0, 1]),
                     label='engy', **next(self.keysiter))
        self.fig.savefig("fig-engy.png", **self.figsave)
        return

    def plt_energy_stress(self, fname='stress.txt'):
        til = os.getcwd().split('/')[-1].split('_')[-2:]
        raw = np.loadtxt(fname)
        ylabeliter = cycle(['E [eV]', r'$\tau$ [Gpa]'])
        self.set_keys()
        self.set_211plt()
        axlist = [self.ax1, self.ax2]
        print raw[:, 0]
        print raw[:, -2]
        self.ax1.plot(raw[:, 0], (raw[:, 1] - raw[0, 1]),
                      label='engy', **next(self.keysiter))
        self.ax2.plot(raw[:, 0], raw[:, -2],
                      label='stress', **next(self.keysiter))
        self.add_legends(*axlist)
        self.set_tick_size(*axlist)
        self.add_y_labels(ylabeliter, *axlist)
        self.add_x_labels(cycle([r'$\epsilon$']), self.ax2)
        self.ax1.set_title('{} {}'.format(*til), fontsize=self.myfontsize)
        self.fig.savefig("fig-ishear.png", **self.figsave)
        return

    def plt_energy_stress_cell(self,
                               fname='ishear.txt'):
        self.set_311plt()
        raw = np.loadtxt(fname)
        raw = raw[raw[:, 0].argsort()]
        print raw
        ylabeliter = cycle(['dE', 'Sxx', 'Syy', 'Szz'])
        if fname == 'iten.txt':
            self.ax1.plot(raw[:, 0], raw[:, 1] - raw[0, 1],
                          label='engy', **next(self.keysiter))
            self.ax2.plot(raw[:, 0], -(raw[:, 4]),
                          label='sxx', **next(self.keysiter))
            syy = np.max([raw[:, 5], raw[:, 6]], axis=0)
            szz = np.min([raw[:, 5], raw[:, 6]], axis=0)
            #
            self.ax3.plot(raw[:, 0], -(syy),
                          label='syy', **next(self.keysiter))
            self.ax3.plot(raw[:, 0], -(szz),
                          label='szz', **next(self.keysiter))
        self.add_legends(*self.axlist)
        self.add_y_labels(ylabeliter, *self.axlist)
        self.set_tick_size(*self.axlist)
        self.fig.savefig('fig-engy-{}'.format(fname.split('.')[0]),
                         **self.figsave)
        return

    def cmp_plt(self):
        self.set_keys()
        self.set_211plt()
        self.plt_energy_stress(fname='stress_0.00.txt')
        self.plt_energy_stress(fname='stress_0.25.txt')
        self.fig.savefig("istress.png", **self.figsave)
        return

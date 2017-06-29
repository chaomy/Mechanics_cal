#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-06-28 00:35:14
# @Last Modified by:   chaomy
# @Last Modified time: 2017-06-29 00:10:38

import matplotlib.pylab as plt
from itertools import cycle
import numpy as np
import plt_drv


class cal_md_ideal_tensile_plt(plt_drv.plt_drv):

    def __init__(self):
        plt_drv.plt_drv.__init__(self)
        return

    def set_figs(self):

        return

    def plt_strain_vs_energy(self, infile='iten.txt'):
        raw = np.loadtxt(infile)
        raw = raw[raw[:, 0].argsort()]
        print raw
        self.set_keys()
        self.set_111plt()
        self.ax.plot(raw[:, 0], (raw[:, 1] - raw[0, 1]),
                     color=self.tableau[3],
                     label='engy',
                     **self.pltkwargs)
        self.fig.savefig("engy.png", **self.figsave)
        return

    def plt_energy_stress_tmp(self):
        raw = np.loadtxt("iten.txt")
        raw = raw[raw[:, 0].argsort()]
        self.set_keys()
        self.set_211plt()
        self.ax1.plot(raw[:, 0], (raw[:, 1] - raw[0, 1]),
                      label='engy',
                      color=self.tableau[0],
                      **self.pltkwargs)
        self.ax2.plot(raw[:, 0], (raw[:, 4] - raw[0, 4]),
                      label='stress',
                      color=self.tableau[3],
                      **self.pltkwargs)
        self.fig.savefig("iten.png", **self.figsave)
        return

    def plot_curv(self, opt='engy'):
        data = np.loadtxt("{}.txt".format(opt))
        if opt == 'engy':
            self.set_keys()
            self.set_311plt()
            self.ax1.plot(data[:, 0], data[:, 1],
                          marker='o', color=plt_drv.tableau20[1])
            self.ax2.plot(data[:, 0], data[:, 2],
                          marker='>', color=plt_drv.tableau20[3])
            self.ax3.plot(data[:, 0], data[:, 3],
                          marker='<', color=plt_drv.tableau20[5])
            self.ax3.plot(data[:, 0], data[:, 4],
                          marker='s', color=plt_drv.tableau20[7])
            self.fig.savefig("{}.png".format(opt), **self.figsave)

        elif opt == 'cell':
            self.set_keys()
            self.set_111plt()
            delta = np.linspace(0, 0.40, 41)
            self.ax.plot(delta, data[:, 0],
                         marker='o', color=plt_drv.tableau20[1])
            self.ax.plot(delta, data[:, 4],
                         marker='>', color=plt_drv.tableau20[3])
            self.ax.plot(delta, data[:, 8],
                         marker='<', color=plt_drv.tableau20[5])
            self.fig.savefig("{}.png".format(opt), **self.figsave)
        return

    def plt_energy_stress(self,
                          fname='ishear.txt'):
        self.set_keys()
        self.set_311plt()
        raw = np.loadtxt(fname)
        raw = raw[raw[:, 0].argsort()]
        print raw
        ylabeliter = cycle(['dE', 'Sxx', 'Syy', 'Szz'])
        if fname == 'iten.txt':
            self.ax1.plot(raw[:, 0], raw[:, 1] - raw[0, 1],
                          label='engy', **next(self.keysiter))
            self.ax2.plot(raw[:, 0], -(raw[:, 4] - raw[0, 4]) * 0.1,
                          label='sxx', **next(self.keysiter))
            self.ax3.plot(raw[:, 0], -(raw[:, 5] - raw[0, 5]) * 0.1,
                          label='syy', **next(self.keysiter))
            self.ax3.plot(raw[:, 0], -(raw[:, 6] - raw[0, 6]) * 0.1,
                          label='szz', **next(self.keysiter))
        self.add_legends(*self.axlist)
        self.add_y_labels(ylabeliter, *self.axlist)
        self.set_tick_size(*self.axlist)
        self.fig.savefig('fig-{}'.format(fname.split('.')[0]),
                         **self.figsave)
        return

    def plt_energy_stress_cmp(self):
        potlist = ['adp', 'pbe']
        self.set_211plt(mfigsize=(8.5, 4.3), lim=True)
        self.set_keys()
        self.set_pltkargs()
        plt.rc('xtick', labelsize=self.mlabelsize)
        plt.rc('ytick', labelsize=self.mlabelsize)
        for i in range(len(potlist)):
            pot = potlist[i]
            fname = 'stress.txt.{}'.format(pot)
            raw = np.loadtxt(fname)
            self.ax1.plot(raw[:, 0],
                          (raw[:, 1] - raw[0, 1]),
                          label=pot,
                          **self.keyslist[i])
            self.ax1.legend(**self.legendarg)
            self.ax2.plot(raw[:, 0],
                          (raw[:, -1] - raw[0, -1]),
                          label=pot,
                          **self.keyslist[i + 2])
            self.ax2.legend(**self.legendarg)
        plt.xlabel('strain', {'fontsize': self.myfontsize})
        self.fig.savefig("stress_cmp.png", **self.figsave)
        return

    def cmp_plt(self):
        self.set_keys()
        self.set_211plt()
        self.plt_energy_stress(fname='stress_0.00.txt')
        self.plt_energy_stress(fname='stress_0.25.txt')
        self.fig.savefig("istress.png", **self.figsave)
        return

#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-06-28 00:35:14
# @Last Modified by:   chaomy
# @Last Modified time: 2017-06-28 09:42:50

import matplotlib.pylab as plt
import numpy as np
import plt_drv


class cal_md_ideal_tensile_plt(plt_drv.plt_drv):

    def __init__(self):
        plt_drv.plt_drv.__init__(self)
        return

    def set_figs(self):

        return

    def plt_strain_vs_energy(self):
        raw = np.loadtxt("ishear.txt")
        raw = raw[raw[:, 0].argsort()]
        print raw
        self.set_keys()
        self.set_111plt()
        # energy
        self.ax.plot(raw[:, 0], (raw[:, 1] - raw[0, 1]),
                     color=self.tableau[3],
                     label='engy',
                     **self.pltkwargs)
        self.fig.savefig("engy.png", **self.figsave)
        return

    def plt_energy_stress(self, fname='ishear.txt', set=False):
        if set is True:
            self.set_keys()
            self.set_211plt()
        raw = np.loadtxt(fname)
        self.ax1.plot(raw[:, 0], (raw[:, 1] - raw[0, 1]),
                      label='engy',
                      **self.pltkwargs)
        self.ax2.plot(raw[:, 0], (raw[:, -1] - raw[0, -1]),
                      label='stress',
                      **self.pltkwargs)
        return

    def set_pltkargs(self):
        self.keyslist = [{'linestyle': '-.', 'color': 'b', 'linewidth': 2,
                          'marker': 'o'},
                         {'linestyle': '--', 'color': 'y', 'linewidth': 2,
                          'marker': '<'},
                         {'linestyle': '-.', 'color': 'g', 'linewidth': 2,
                          'marker': 'o'},
                         {'linestyle': '--', 'color': 'm', 'linewidth': 2,
                          'marker': '<'}]
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
            self.ax1.set_ylabel('energy [eV]', {'fontsize': self.myfontsize})

            self.ax2.plot(raw[:, 0],
                          (raw[:, -1] - raw[0, -1]),
                          label=pot,
                          **self.keyslist[i + 2])
            self.ax2.legend(**self.legendarg)
            self.ax2.set_ylabel('stress [Gpa]', {'fontsize': self.myfontsize})

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

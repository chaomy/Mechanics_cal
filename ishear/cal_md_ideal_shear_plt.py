#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-06-28 00:35:14
# @Last Modified by:   chaomy
# @Last Modified time: 2018-06-25 17:12:08


from itertools import cycle
from md_pot_data import unitconv
from matplotlib.ticker import FormatStrFormatter
from numpy import polyfit, polyval
import numpy as np
import plt_drv
import os


class cal_bcc_ideal_shear_plt(object):

    def __init__(self):
        plt_drv.plt_drv.__init__(self)

    def plt_vc(self):
        self.set_111plt()
        # tp_sxx = np.loadtxt('0.25Re-Tpath.txt')
        # sh_211 = np.loadtxt('stress.txt.0.25')
        # tp_sxx = np.loadtxt('0.50Re-Tpath.txt')
        # sh_211 = np.loadtxt('stress.txt.0.50')
        tp_sxx = np.loadtxt('0.50Re-Tpath.txt')
        sh_211 = np.loadtxt('stress.txt.0.50')
        self.ax.plot(tp_sxx[:, 0], tp_sxx[:, 1],
                     label='tp', **next(self.keysiter))
        self.ax.plot(sh_211[:, 0], sh_211[:, -3],
                     label='sh', **next(self.keysiter))
        self.add_legends(self.ax)
        self.fig.savefig('fig-stress.0.50.png')

    def plt_energy(self, infile='ishear.txt'):
        raw = np.loadtxt(infile)
        raw = raw[raw[:, 0].argsort()]
        print(raw)
        self.set_111plt()
        self.ax.plot(raw[:, 0], (raw[:, 1] - raw[0, 1]),
                     label='engy', **next(self.keysiter))
        self.fig.savefig("fig-engy.png", **self.figsave)

    def plt_energy_stress_lmp(self, fname='stress.txt'):
        # def plt_energy_stress_lmp(self, fname='ishear.txt'):
        til = os.getcwd().split('/')[-1].split('_')[-2:]
        raw = np.loadtxt(fname)
        ylabeliter = cycle(['E [eV]', r'$\tau$ [Gpa]'])
        self.set_keys()
        self.set_211plt()
        axlist = [self.ax1, self.ax2]
        self.ax1.plot(raw[:, 0], (raw[:, 1] - raw[0, 1]),
                      label='engy', **next(self.keysiter))
        self.ax2.plot(raw[:, 0], raw[:, -1],
                      label='stress', **next(self.keysiter))
        self.add_legends(*axlist)
        self.set_tick_size(*axlist)
        self.add_y_labels(ylabeliter, *axlist)
        self.add_x_labels(cycle([r'$\epsilon$']), self.ax2)
        self.fig.savefig("fig-ishear.png", **self.figsave)

    def plt_energy_stress(self, fname='stress.txt'):
        raw = np.loadtxt(fname)
        ylabeliter = cycle(['E [eV]', r'$\tau$ [Gpa]'])
        self.set_keys()
        self.set_211plt()
        axlist = [self.ax1, self.ax2]
        index = np.where(raw[:, 1] < -100.0)
        raw[:, 1][index] = raw[:, 1][index] / unitconv.uengy['rytoeV']
        self.ax1.plot(raw[:, 0], (raw[:, 1] - raw[0, 1]),
                      label='engy', **next(self.keysiter))
        yy = raw[:, -1]
        # call interp
        ply = polyfit(raw[:, 0], yy, 2)
        print(polyval(ply, [0.08]))
        self.ax2.plot(raw[:, 0], yy, label='stress', **next(self.keysiter))
        self.add_legends(*axlist)
        self.set_tick_size(*axlist)
        self.add_y_labels(ylabeliter, *axlist)
        self.add_x_labels(cycle([r'$\epsilon$']), self.ax2)
        self.fig.savefig("fig-ishear.png", **self.figsave)

    def plt_cmp_pth(self, tg='va'):
        self.set_211plt()
        self.set_keys()
        ylabeliter = cycle(['Energy per atom [eV]', r'Shear stress [Gpa]'])
        if tg in ['md']:
            self.ax1.set_ylim([-0.01, 0.170])
            self.ax2.set_ylim([-7.2, 8.0])
        elif tg in ['va']:
            self.ax1.set_ylim([-0.01, 0.170])
            self.ax2.set_ylim([-7.2, 8.0])
        self.ax1.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
        self.ax2.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

        cc = 1.0
        if tg in ['va']:
            cc = 0.1
        # vasp use Bar -> times 0.1 to be GPa
        lab = '{211}<111>'
        raw = np.loadtxt('stress.{}.211.txt'.format(tg))
        self.ax1.plot(raw[:, 0], raw[:, 1] - raw[0, 1],
                      label=lab, **next(self.keysiter))
        if tg in ['md']:
            self.ax2.plot(raw[:, 0], cc * raw[:, -3],
                          label=lab, **next(self.keysiter))
        if tg in ['va']:
            self.ax2.plot(raw[:, 0], cc * raw[:, -2],
                          label=lab, **next(self.keysiter))
        lab = '{110}<111>'
        raw = np.loadtxt('stress.{}.110.txt'.format(tg))
        self.ax1.plot(raw[:, 0], raw[:, 1] - raw[0, 1],
                      label=lab, **next(self.keysiter))
        self.ax2.plot(raw[:, 0], cc * raw[:, -1],
                      label=lab, **next(self.keysiter))
        self.add_legends(*self.axls)
        self.set_tick_size(*self.axls)
        self.add_y_labels(ylabeliter, *self.axls)
        self.add_x_labels(cycle([r'Shear strain']), self.ax2)
        self.fig.savefig("fig-cmp-{}.png".format(tg), **self.figsave)

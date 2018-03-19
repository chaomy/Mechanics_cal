#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-06-28 00:35:14
# @Last Modified by:   chaomy
# @Last Modified time: 2018-03-18 07:48:14


from itertools import cycle
from md_pot_data import unitconv
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
        print raw
        self.set_111plt()
        self.ax.plot(raw[:, 0], (raw[:, 1] - raw[0, 1]),
                     label='engy', **next(self.keysiter))
        self.fig.savefig("fig-engy.png", **self.figsave)

    def plt_energy_stress_lmp(self, fname='stress.txt'):
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

    def plt_energy_stress(self, ptype='211', fname='stress.txt'):
        til = os.getcwd().split('/')[-1].split('_')[-2:]
        raw = np.loadtxt(fname)
        ylabeliter = cycle(['E [eV]', r'$\tau$ [Gpa]'])
        self.set_keys()
        self.set_211plt()
        axlist = [self.ax1, self.ax2]
        print raw[:, 0]
        index = np.where(raw[:, 1] < -100.0)
        raw[:, 1][index] = raw[:, 1][index] / unitconv.uengy['rytoeV']
        self.ax1.plot(raw[:, 0], (raw[:, 1] - raw[0, 1]),
                      label='engy', **next(self.keysiter))

        if ptype in ['211']:
            # yy = raw[:, -3]
            yy = raw[:, -2]
        elif ptype in ['110']:
            yy = -raw[:, -1]

        # call interp
        ply = polyfit(raw[:, 0], yy, 2)
        print polyval(ply, [0.08])

        self.ax2.plot(raw[:, 0], yy,
                      label='stress', **next(self.keysiter))
        self.add_legends(*axlist)
        self.set_tick_size(*axlist)
        self.add_y_labels(ylabeliter, *axlist)
        self.add_x_labels(cycle([r'$\epsilon$']), self.ax2)
        self.ax1.set_title('{} {}'.format(*til), fontsize=self.myfontsize)
        self.fig.savefig("fig-ishear.png", **self.figsave)

    def plt_cmp(self):
        pln = '211'
        self.set_211plt()
        self.set_keys('upper left')
        axlist = [self.ax1, self.ax2]
        filelist = ['stress.adp.p{}'.format(pln), 'stress.pbe.p{}'.format(pln)]
        til = 'ideal shear along ({})'.format(pln)
        lablist = ['adp', 'pbe']
        ylabeliter = cycle(['E [eV]', r'$\tau$ [Gpa]'])
        for i in range(len(filelist)):
            raw = np.loadtxt(filelist[i])
            self.ax1.plot(raw[:, 0], (raw[:, 1] - raw[0, 1]),
                          label=lablist[i], **next(self.keysiter))
            self.ax2.plot(raw[:, 0], raw[:, -1],
                          label=lablist[i], **next(self.keysiter))
        self.add_legends(*axlist)
        self.set_tick_size(*axlist)
        self.add_y_labels(ylabeliter, *axlist)
        self.add_x_labels(cycle([r'$\epsilon$']), self.ax2)
        # self.ax1.set_title(til, fontsize=self.myfontsize)
        self.fig.savefig("fig-cmp-{}.png".format(pln), **self.figsave)

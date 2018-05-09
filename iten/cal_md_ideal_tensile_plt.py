#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-06-28 00:35:14
# @Last Modified by:   chaomy
# @Last Modified time: 2018-05-05 14:34:16


from optparse import OptionParser
from itertools import cycle
import matplotlib.pylab as plt
import numpy as np
import plt_drv


class cal_md_ideal_tensile_plt(plt_drv.plt_drv):

    def __init__(self):
        plt_drv.plt_drv.__init__(self)

    def plt_mesh(self):
        data = np.loadtxt("save.txt")
        xx = data[:, 0].reshape(20, 20)
        yy = data[:, 1].reshape(20, 20)
        zz = data[:, 2].reshape(20, 20)
        print(xx)
        self.set_111plt()
        cs = self.ax.contour(xx, yy, zz)
        self.ax.clabel(cs, inline=0.02, lw=5, fontsize=6)
        self.fig.savefig("fig_mesh.png")

    def plt_energy_stress(self, fname='iten.txt'):
        self.set_211plt()
        raw = np.loadtxt(fname)
        raw = raw[raw[:, 0].argsort()]
        ylabeliter = cycle(['dE', 'Sxx', 'Syy', 'Szz'])
        self.ax1.plot(raw[:, 0], (raw[:, 1] - raw[0, 1]),
                      label='engy', **next(self.keysiter))
        self.ax2.plot(raw[:, 0], raw[:, 4],
                      label='sxx', **next(self.keysiter))
        self.add_legends(*self.axls)
        self.add_y_labels(ylabeliter, *self.axls)
        self.add_x_labels(cycle([r'$\varepsilon_{xx}$']), self.ax2)
        self.set_tick_size(*self.axls)
        self.fig.savefig('fig-engy-{}'.format(fname.split('.')[0]),
                         **self.figsave)

    def cmp_pot_plt(self, tg='tp'):
        self.set_211plt()
        raw = np.loadtxt('iten.va.{}.txt'.format(tg))
        if tg in ['op']:
            raw[:, 1] /= 2.0
        ylabeliter = cycle(['dE', 'Sxx', 'Syy', 'Szz'])
        self.ax1.plot(raw[:, 0], raw[:, 1] - raw[0, 1],
                      label='PAW-PBE', **next(self.keysiter))
        self.ax2.plot(raw[:, 0], raw[:, 4],
                      label='PAW-PBE', **next(self.keysiter))
        raw = np.loadtxt('iten.md.{}.txt'.format(tg))
        if tg in ['op']:
            raw[:, 1] /= 2.0
        self.ax1.plot(raw[:, 0], raw[:, 1] - raw[0, 1],
                      label='MEAM', **next(self.keysiter))
        self.ax2.plot(raw[:, 0], raw[:, 4],
                      label='MEAM', **next(self.keysiter))
        self.add_legends(*self.axls)
        self.add_y_labels(ylabeliter, *self.axls)
        self.add_x_labels(cycle([r'$\varepsilon_{xx}$']), self.ax2)
        self.set_tick_size(*self.axls)
        self.fig.savefig('fig-iten-cmp-{}.png'.format(tg), **self.figsave)

    def cmp_path_plt(self, tg='md'):
        self.set_211plt()
        self.set_keys('upper left')
        raw = np.loadtxt('iten.{}.op.txt'.format(tg))
        raw = raw[:25, :]
        raw[:, 1] /= 2.0
        if tg in ['md']:
            self.ax1.set_ylim([-0.02, 0.63])
        elif tg in ['va']:
            self.ax1.set_ylim([-0.01, 0.61])
        self.ax2.set_ylim([-15, 19])
        ylabeliter = cycle(
            ['Energy per atom [eV]', 'Stress along tensile [Gpa]'])
        self.ax1.plot(raw[:, 0], (raw[:, 1] - raw[0, 1]),
                      label='OP', **next(self.keysiter))
        self.ax2.plot(raw[:, 0], raw[:, 4],
                      label='OP', **next(self.keysiter))
        raw = np.loadtxt('iten.{}.tp.txt'.format(tg))
        raw = raw[:25, :]
        self.ax1.plot(raw[:, 0], raw[:, 1] - raw[0, 1],
                      label='TP', **next(self.keysiter))
        self.ax2.plot(raw[:, 0], raw[:, 4],
                      label='TP', **next(self.keysiter))
        self.add_legends(*self.axls)
        self.add_y_labels(ylabeliter, *self.axls)
        self.add_x_labels(cycle([r'$\varepsilon_{xx}$']), self.ax2)
        self.set_tick_size(*self.axls)
        self.fig.savefig('fig-iten-path-{}.png'.format(tg), **self.figsave)

    def plt_cell(self, tg='op'):
        raw = np.loadtxt('iten.va.{}.txt'.format(tg))
        raw = raw[raw[:, 0].argsort()]
        self.set_111plt()
        ylabeliter = cycle(['lattice'])
        lyy = np.max([raw[:, 2], raw[:, 3]], axis=0)
        lzz = np.min([raw[:, 2], raw[:, 3]], axis=0)
        self.ax.plot(raw[:, 0], lyy,
                     'o--', label='vasp lyy', markersize=14, lw=4, c=self.colorlist[0])
        self.ax.plot(raw[:, 0], lzz,
                     'o--', label='vasp lzz', markersize=14, lw=4, c=self.colorlist[0])

        raw = np.loadtxt('iten.md.{}.txt'.format(tg))
        raw = raw[raw[:, 0].argsort()]
        lyy = np.max([raw[:, 2], raw[:, 3]], axis=0)
        lzz = np.min([raw[:, 2], raw[:, 3]], axis=0)
        self.ax.plot(raw[:, 0], lyy,
                     '<--', label='meam lyy', markersize=14, lw=4, c=self.colorlist[1])
        self.ax.plot(raw[:, 0], lzz,
                     '<--', label='meam lzz', markersize=14, lw=4, c=self.colorlist[1])
        self.add_legends(self.ax)
        self.add_y_labels(ylabeliter, self.ax)
        self.set_tick_size(self.ax)

        self.fig.savefig('fig-cell-cmp-{}.png'.format(tg), **self.figsave)


if __name__ == '__main__':
    usage = "usage:%prog [options] arg1 [options] arg2"
    parser = OptionParser(usage=usage)
    parser.add_option('-t', "--mtype", action="store",
                      type="string", dest="mtype")
    parser.add_option('-p', "--param", action="store",
                      type='string', dest="fargs")
    (options, args) = parser.parse_args()
    drv = cal_md_ideal_tensile_plt()
    dispatcher = {'plt': drv.plt_energy_stress,
                  'cell': drv.plt_cell,
                  'cmp': drv.cmp_pot_plt,
                  'pth': drv.cmp_path_plt,
                  'mesh': drv.plt_mesh}

    if options.fargs is not None:
        dispatcher[options.mtype.lower()](options.fargs)
    else:
        dispatcher[options.mtype.lower()]()

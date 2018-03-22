#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-06-28 00:35:14
# @Last Modified by:   chaomy
# @Last Modified time: 2018-03-20 21:00:32


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

    def plt_strain_vs_energy(self, infile='ishear.txt'):
        raw = np.loadtxt(infile)
        raw = raw[raw[:, 0].argsort()]
        print(raw)
        self.set_111plt()
        self.ax.plot(raw[:, 0], (raw[:, 1] - raw[0, 1]),
                     label='engy', **next(self.keysiter))
        self.fig.savefig("fig-engy.png", **self.figsave)

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

    def plt_energy_stress_two(self, fname='iten.txt'):
        self.set_111plt()
        raw = np.loadtxt("iten_tp.txt")
        raw = raw[raw[:, 0].argsort()]
        ylabeliter = cycle(['dE'])
        self.ax.plot(raw[:, 0], (raw[:, 1] - raw[0, 1]),
                     label='tp', **next(self.keysiter))
        raw = np.loadtxt("iten_op.txt")
        raw = raw[raw[:, 0].argsort()]
        self.ax.plot(raw[:, 0], (raw[:, 1] - raw[0, 1]) / 2.,
                     label='op', **next(self.keysiter))
        self.add_legends(*self.axls)
        self.add_y_labels(ylabeliter, *self.axls)
        self.add_x_labels(cycle([r'$\varepsilon_{xx}$']), self.ax)
        self.set_tick_size(*self.axls)
        self.fig.savefig('fig-engy-{}'.format(fname.split('.')[0]),
                         **self.figsave)

        ylabeliter = cycle(['Sxx'])
        self.set_111plt()
        raw = np.loadtxt("iten_tp.txt")
        raw = raw[raw[:, 0].argsort()]
        self.ax.plot(raw[:, 0], raw[:, 4],
                     label='tp', **next(self.keysiter))

        raw = np.loadtxt("iten_op.txt")
        raw = raw[raw[:, 0].argsort()]
        self.ax.plot(raw[:, 0], raw[:, 4],
                     label='op', **next(self.keysiter))
        self.add_legends(*self.axls)
        self.add_y_labels(ylabeliter, *self.axls)
        self.add_x_labels(cycle([r'$\varepsilon_{xx}$']), self.ax)
        self.set_tick_size(*self.axls)
        self.fig.savefig('fig-stss-{}'.format(fname.split('.')[0]),
                         **self.figsave)

    def plt_energy_stress_cell(self,
                               fname='ishear.txt'):
        self.set_311plt()
        raw = np.loadtxt(fname)
        raw = raw[raw[:, 0].argsort()]
        print(raw)
        ylabeliter = cycle(['dE', 'Sxx', 'Syy', 'Szz'])
        self.ax1.plot(raw[:, 0], raw[:, 1] - raw[0, 1],
                      label='engy', **next(self.keysiter))
        self.ax2.plot(raw[:, 0], -(raw[:, 4]),
                      label='sxx', **next(self.keysiter))
        syy = np.max([raw[:, 5], raw[:, 6]], axis=0)
        szz = np.min([raw[:, 5], raw[:, 6]], axis=0)

        self.ax3.plot(raw[:, 0], -(syy),
                      label='syy', **next(self.keysiter))
        self.ax3.plot(raw[:, 0], -(szz),
                      label='szz', **next(self.keysiter))
        self.add_legends(*self.axls)
        self.add_y_labels(ylabeliter, *self.axls)
        self.set_tick_size(*self.axls)
        self.fig.savefig('fig-engy-{}'.format(fname.split('.')[0]),
                         **self.figsave)

    def adjust_data_format(self, fname='iten.txt'):
        raw = np.loadtxt(fname)
        raw = raw[raw[:, 0].argsort()]
        raw[:, 1] = raw[:, 1] / 4.
        raw[:, 4] = 0.1 * raw[:, 4]
        np.savetxt('iten.save.txt', raw)

    def cmp_plt(self):
        self.set_211plt()
        self.set_keys('upper left')
        raw = np.loadtxt('iten.tp.txt')
        raw = raw[:25, :]
        ylabeliter = cycle(['Energy per atom[eV]', 'Sxx [Gpa]'])
        self.ax1.plot(raw[:, 0], raw[:, 1] - raw[0, 1],
                      label='tp', **next(self.keysiter))
        self.ax2.plot(raw[:, 0], -raw[:, 4],
                      label='tp', **next(self.keysiter))
        raw = np.loadtxt('iten.op.txt')
        raw = raw[:25, :]
        self.ax1.plot(raw[:, 0], raw[:, 1] - raw[0, 1],
                      label='op', **next(self.keysiter))
        self.ax2.plot(raw[:, 0], -raw[:, 4],
                      label='op', **next(self.keysiter))
        self.add_legends(*self.axls)
        self.add_y_labels(ylabeliter, *self.axls)
        self.add_x_labels(cycle([r'$\varepsilon_{xx}$']), self.ax2)
        self.set_tick_size(*self.axls)
        self.fig.savefig('fig-iten-cmp.png', **self.figsave)

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

    def cmp_pot_plt(self, tg='tp'):
        self.set_211plt()
        raw = np.loadtxt('iten.va.{}.txt'.format(tg))
        ylabeliter = cycle(['dE', 'Sxx', 'Syy', 'Szz'])
        self.ax1.plot(raw[:, 0], raw[:, 1] - raw[0, 1],
                      label='va', **next(self.keysiter))
        self.ax2.plot(raw[:, 0], raw[:, 4],
                      label='va', **next(self.keysiter))

        raw = np.loadtxt('iten.md.{}.txt'.format(tg))
        self.ax1.plot(raw[:, 0], raw[:, 1] - raw[0, 1],
                      label='md', **next(self.keysiter))
        self.ax2.plot(raw[:, 0], raw[:, 4],
                      label='md', **next(self.keysiter))
        self.add_legends(*self.axls)
        self.add_y_labels(ylabeliter, *self.axls)
        self.add_x_labels(cycle([r'$\varepsilon_{xx}$']), self.ax2)
        self.set_tick_size(*self.axls)
        self.fig.savefig('fig-iten-cmp-{}.png'.format(tg), **self.figsave)


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
                  'plt2': drv.plt_energy_stress_two,
                  'cell': drv.plt_cell,
                  'adj': drv.adjust_data_format,
                  'cmp': drv.cmp_pot_plt,
                  'mesh': drv.plt_mesh}

    if options.fargs is not None:
        dispatcher[options.mtype.lower()](options.fargs)
    else:
        dispatcher[options.mtype.lower()]()

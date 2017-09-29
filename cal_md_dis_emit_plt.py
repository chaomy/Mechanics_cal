#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-06-25 14:28:58
# @Last Modified by:   chaomy
# @Last Modified time: 2017-09-11 22:34:13

from numpy import loadtxt, min, array
from itertools import cycle


class cal_md_dis_emit_plt(object):

    def plt_k1(self):
        mlab = [r'WTa$_{0.50}$', r'W', r'WRe$_{0.05}$', r'WRe$_{0.10}$',
                r'WRe$_{0.15}$', r'WRe$_{0.20}$',
                r'WRe$_{0.25}$', r'WRe$_{0.50}$',
                r'Ta', 'Nb', 'Nb']

        self.set_111plt((9.2, 4.5))
        self.set_keys(loc='lower right')
        self.ax.set_ylim(1.4, 2.8)
        # self.ax.set_xlim(-0.12, 7.45)
        data = loadtxt('vcaw_112_Ke.txt')
        xlist = array([5.5, 6.0, 6.05, 6.10, 6.15, 6.20, 6.25, 6.5])
        self.ax.plot(xlist, data[:, 1],
                     label=r'K$_{1e}$ [11$\overline{1}$](2$\overline{1}$1)',
                     **next(self.keysiter))
        # for i, xx, yy, lab in zip(range(len(data[:, 0])),
        #                           data[:, 0],
        #                           min(data[:, 1:-1], axis=1), mlab):
        # for i, xx, yy, lab in zip(range(len(data[:, 0])),
        #                           xlist, min(data[:, 1:-1], axis=1), mlab):
        #     if i in [0]:
        #         self.ax.text(xx + 0.1, yy, lab,
        #                      fontsize=self.myfontsize - 4)
        #     elif i in [6, 7]:
        #         self.ax.text(xx - 0.33, yy - 0.13, lab,
        #                      fontsize=self.myfontsize - 4)
        #     else:
        #         self.ax.text(xx, yy - 0.0825, lab,
        #                      fontsize=self.myfontsize - 4)

        data = loadtxt('vcaw_110_Ke.txt')
        self.ax.plot(xlist, data[:, 1],
                     label=r'K$_{1e}$ [111](01$\overline{1}$)',
                     **next(self.keysiter))
        self.ax.plot(xlist, data[:, 2],
                     label=r'K$_{1c}$',
                     **next(self.keysiter))
        ylabiter = cycle([r"Stress intensity factor [MPa$\sqrt{m}$]"])
        xlabiter = cycle(["Number of valence electrons per atom"])

        self.add_legends(self.ax)
        self.set_tick_size(self.ax)
        # self.remove_xticks(self.ax)
        self.add_x_labels(xlabiter, self.ax)
        self.add_y_labels(ylabiter, self.ax)
        self.fig.savefig('k1e.png', **self.figsave)

        self.closefig()
        self.set_111plt((9.2, 4.5))
        self.ax.set_ylim(0.885, 1.285)
        # self.ax.set_xlim(-0.12, 7.45)
        self.set_keys(loc='lower left')
        data = loadtxt('vcaw_112_Ke.txt')
        next(self.keysiter)
        self.ax.plot(xlist, data[:, -1],
                     label=r'[11$\overline{1}$](2$\overline{1}$1)',
                     **next(self.keysiter))

        # for i, xx, yy, lab in zip(range(len(data[:, 0])),
        #                           data[:, 0],
        #                           min(data[:, 1:], axis=1), mlab):
        # for i, xx, yy, lab in zip(range(len(data[:, 0])), xlist,
        #                           min(data[:, 1:], axis=1), mlab):
        #     if i in [0]:
        #         self.ax.text(xx + 0.1, yy - 0.02, lab,
        #                      fontsize=self.myfontsize - 4)
        #     elif i in [6]:
        #         self.ax.text(xx - 0.42, yy - 0.04, lab,
        #                      fontsize=self.myfontsize - 4)
        #     elif i in [7]:
        #         self.ax.text(xx - 0.75, yy - 0.014, lab,
        #                      fontsize=self.myfontsize - 4)
        #     else:
        #         self.ax.text(xx - 0.1, yy - 0.0245, lab,
        #                      fontsize=self.myfontsize - 4)
        data = loadtxt('vcaw_110_Ke.txt')
        next(self.keysiter)
        next(self.keysiter)
        next(self.keysiter)
        self.ax.plot(xlist, data[:, -1],
                     label=r'[111](01$\overline{1}$)',
                     **next(self.keysiter))
        # for xx, yy, lab in zip(data[:, 0], data[:, -1], mlab):
        #     self.ax.text(xx, yy + 0.01, lab, fontsize=self.myfontsize - 3)

        ylabiter = cycle([r"K$_{Ie}$/K$_{Ic}$"])
        xlabiter = cycle(["Number of valence electrons per atom"])
        self.add_legends(self.ax)
        self.set_tick_size(self.ax)
        # self.remove_xticks(self.ax)
        self.add_x_labels(xlabiter, self.ax)
        self.add_y_labels(ylabiter, self.ax)
        self.fig.savefig('k1e2k1c.png', **self.figsave)
        return

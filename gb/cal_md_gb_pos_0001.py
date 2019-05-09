#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-07-05 08:12:30
# @Last Modified by:   chaomy
# @Last Modified time: 2018-11-19 22:16:21

from numpy import loadtxt, linspace
import os
from math import cos, sin
from math import sqrt
from itertools import cycle
import numpy as np


class md_gb_pos(object):
    
    def loop_aft_angle(self):
        self.find_angles_0001()
        total = np.ndarray([len(self.ag), 2])
        os.system("mkdir result_1")
        os.system("mkdir result_2")
        dd = np.zeros([len(self.ag) + 2, 2])
        dd[0, 0], dd[0, 1] = 0.0, 0.0
        dd[-1, 0], dd[-1, 1] = 60.0, 0.0
        grep="tail -n 1 log | awk \'{print $11}\' > final_gb_eng.txt"
        cc="cp data.20.txt ../result_1/"
        cc_2="cp data.20.txt ../result_2/"
        for e, i in zip(self.ag, range(len(self.ag))):
            e_rest=90.0-e[0]
            mdir = "0001_{:.2f}".format(e_rest)
            os.chdir(mdir)
            os.system(grep)
            os.system(cc+str(i)+".txt")
            os.system(cc_2+mdir+".lmp")
            dd[i + 1, 0] = e_rest*2
            dd[i + 1, 1] = np.loadtxt("final_gb_eng.txt")

            total[i, 0] = e_rest
            total[i, 1] = dd[i + 1, 1]
            os.chdir(os.pardir)

        np.savetxt('data_opt.txt', total, fmt='%1.8f')
        np.savetxt('data_opt_conv.txt', dd, fmt='%1.8f')
            
    def loop_plt_angle(self):
        self.find_angles_0001()

        dd = np.zeros([len(self.ag) + 2, 2])
        dd[0, 0], dd[0, 1] = 0.0, 0.0
        dd[-1, 0], dd[-1, 1] = 90.0, 0.0

        total = np.ndarray([len(self.ag), 2])

        for e, i in zip(self.ag, range(len(self.ag))):
            e_rest=90.0-e[0]
            mdir = "0001_{:.2f}".format(e_rest)
            dd[i + 1, 0] = e_rest
            dd[i + 1, 1] = np.loadtxt("{}/lmp.dat".format(mdir)) 

            total[i, 0] = e_rest
            total[i, 1] = dd[i + 1, 1]   

        np.savetxt('data_dirct.txt', total, fmt='%1.8f')
        print(dd)

        self.set_111plt((9, 6))
        self.set_keys()
        self.ax.plot(dd[1:-1, 0], 1e3 * dd[1:-1, 1], 'o--',
                    markersize=14, label='<0001> GB')
        # ddp **next(self.keysiter)

        self.add_legends(self.ax)

        ylb = cycle(['GB[0001] mJ/m^2'])
        xlb = cycle(['Angle (deg)'])

        self.add_y_labels(ylb, self.ax)
        self.add_x_labels(xlb, self.ax)
        self.set_tick_size(self.ax)
        self.fig.savefig('fig_gb.png', **self.figsave)

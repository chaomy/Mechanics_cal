#!/usr/bin/env python
# encoding: utf-8

###################################################################
#
# File Name : cal_vasp_itensile.py
#
###################################################################
#
# Purpose :
#
# Creation Date :
# Last Modified :
# Created By    : Chaoming Yang
#
###################################################################
import get_data
import md_pot_data
import os
import numpy as np
import plt_drv
import matplotlib.pylab as plt
from optparse import OptionParser


class cal_bcc_ideal_tensile(get_data.get_data,
                            plt_drv.plt_drv):

    def __init__(self):
        self._pot = md_pot_data.dft_data.Nb_pbe
        self.root_dir = os.getcwd()

        get_data.get_data.__init__(self)
        plt_drv.plt_drv.__init__(self)
        return

    def grab_engy(self):
        npts = 26
        data = np.ndarray([npts, 9])
        for i in range(26):
            dirname = "dir-{:03d}".format(i)
            print dirname
            os.chdir(dirname)
            delta = 0.01 * i
            engy, stress, vol = self.vasp_energy_stress_vol()
            (data[i, 0], data[i, 1], data[i, 2:8], data[i, -1]) = \
                delta, engy, stress.transpose(), vol
            os.chdir(self.root_dir)
        print data
        np.savetxt("iten.txt", data)
        return

    def plot_curv(self):
        data = np.loadtxt("iten.txt")
        self.set_keys()
        self.set_111plt()
        xlist = np.arange(0, 0.25, 0.01)
        self.ax.plot(xlist, data[:25, 0])
        self.fig.savefig("istress.png", **self.figsave)
        return


if __name__ == '__main__':
    usage = "usage:%prog [options] arg1 [options] arg2"
    parser = OptionParser(usage=usage)
    parser.add_option('-t', "--mtype",
                      action="store",
                      type="string",
                      dest="mtype", help="",
                      default="prp_r")
    (options, args) = parser.parse_args()

    drv = cal_bcc_ideal_tensile()

    if options.mtype.lower() == 'clc':
        drv.grab_engy()

    if options.mtype.lower() == 'plt':
        drv.plot_curv()

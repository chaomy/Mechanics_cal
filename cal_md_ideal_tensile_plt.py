#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-06-28 00:35:14
# @Last Modified by:   chaomy
# @Last Modified time: 2017-06-28 00:38:37

import sys
import matplotlib.pylab as plt
import pickle as pc
import glob
import numpy as np
import os
import plt_drv


class cal_md_ideal_tensile_plt(object):

    def __init__(self):

        return

    def set_figs(self):

        return

        def plt_energy_stress(self, fname='stress.txt', set=False):
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

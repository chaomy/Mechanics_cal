#!/usr/bin/env python
# encoding: utf-8
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-07-05 08:11:31
# @Last Modified by:   chaomy
# @Last Modified time: 2017-08-31 15:13:04


import md_pot_data
import cal_md_crack_uti
import cal_md_crack_pre
import cal_md_crack_run
import gn_config
from math import sqrt
from math import pi


class md_crack_ini(cal_md_crack_pre.md_crack_pre,
                   cal_md_crack_uti.md_crack_uti,
                   cal_md_crack_run.md_crack_run,
                   gn_config.bcc):

    def __init__(self):
        self.pot = md_pot_data.md_pot.Nb_eam
        cal_md_crack_pre.md_crack_pre.__init__(self)
        cal_md_crack_run.md_crack_run.__init__(self)
        cal_md_crack_uti.md_crack_uti.__init__(self)
        gn_config.bcc.__init__(self, self.pot)
        self.ckcoeff = cal_md_crack_uti.crack_coeff
        # self.set_params()
        return

    def set_params(self, elasticconstants):
        self.elastic = elasticconstants
        self.element = self.pot['element']
        self.surfe100, self.surfe110, self.surfe111 = \
            self.pot['surf100'], self.pot['surf110'], self.pot['surf111']
        self.c11, self.c12, self.c44 = \
            self.pot['c11'], self.pot['c12'], self.pot['c44']
        self.burger = sqrt(3) / 2. * self.pot['lattice']
        self.screw_coeff = self.burger / (2. * pi)
        self.edge_coeff = self.burger / (2. * pi)
        self.nu = 0.33
        return

#!/usr/bin/env python
# encoding: utf-8
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-07-05 08:11:31
# @Last Modified by:   chaomy
# @Last Modified time: 2018-03-27 16:46:07


from .cal_md_crack_uti import md_crack_uti
from .cal_md_crack_pre import md_crack_pre
from .cal_md_crack_run import md_crack_run
import gn_config
from math import sqrt
from math import pi


# class crack_coeff:
# Gg = None
# BB = None
# KK = None
# Kg = None
# p1 = None
# p2 = None
# q1 = None
# q2 = None
# u1 = None
# u2 = None

# class mtxs:
#     lij = None
#     sij = None

class crack_coeff:
    Gg = None
    BB = None
    K1 = None
    Kg = None
    p1 = None
    p2 = None
    q1 = None
    q2 = None
    u1 = None
    u2 = None

    lij = None
    sij = None
    bij_pstrain = {'b11': None, 'b12': None,
                   'b22': None, 'b16': None,
                   'b26': None, 'b66': None}


class md_crack_ini(md_crack_pre,
                   md_crack_uti,
                   md_crack_run,
                   gn_config.gnStructure):

    def __init__(self):
        md_crack_pre.__init__(self)
        md_crack_run.__init__(self)
        md_crack_uti.__init__(self)
        gn_config.gnStructure.__init__(self, self.pot)
        self.ckcoeff = crack_coeff
        # self.set_params()

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

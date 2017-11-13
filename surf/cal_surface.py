#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-06-28 00:35:14
# @Last Modified by:   chaomy
# @Last Modified time: 2017-11-09 16:14:52


from optparse import OptionParser
from md_pot_data import fluxdirs
import os
import gsf_data
import numpy as np
import md_pot_data
import ase.io
from gsf import cal_gsf
import cal_qe_surface
import cal_va_surface
import cal_md_surface


class cal_surface(cal_gsf.cal_gsf,
                  cal_qe_surface.cal_qe_surface,
                  cal_va_surface.cal_va_surface,
                  cal_md_surface.cal_md_surface):

    def __init__(self, pot=md_pot_data.va_pot.Nb_pbe,
                 msurf='x100z100'):
        self.pot = pot
        self.msurf = msurf
        cal_va_surface.cal_va_surface.__init__(self)
        cal_gsf.cal_gsf.__init__(self, self.pot, self.msurf)
        return

    def loop_prep_va(self):
        surlist = ['x100z100', 'x110z110', 'x112z111']
        for self.mgsf in surlist:
            print self.mgsf
            self.prep_va_surface()
        return

    def gn_surface_atoms(self):
        mgsf = self.mgsf
        atomss = self.set_bcc_convention(
            in_direction=gsf_data.gsfbase[mgsf],
            in_size=gsf_data.gsfsize[mgsf])
        for i in range(gsf_data.gsfpopn[mgsf]):
            atomss.pop()

        atomsb = self.set_bcc_convention(
            in_direction=gsf_data.gsfbase[mgsf],
            in_size=gsf_data.bulksize[mgsf])
        atomss.wrap()
        atomsb.wrap()
        return [atomss, atomsb]


if __name__ == '__main__':
    usage = "usage:%prog [options] arg1 [options] arg2"
    parser = OptionParser(usage=usage)
    parser.add_option("-t", "--mtype", action="store",
                      type="string", dest="mtype")
    parser.add_option('-p', "--param", action="store",
                      type='string', dest="fargs")
    (options, args) = parser.parse_args()
    drv = cal_surface()
    dispatcher = {'prepqe': drv.prep_qe_surface,
                  'prepva': drv.prep_va_surface,
                  'loopprepqe': drv.loop_pot_surf,
                  'loopsurfqe': drv.loop_cal_surf,
                  'loopprepva': drv.loop_prep_va,
                  'plt': drv.plt_surf,
                  'trans': drv.transdata}

    dispatcher[options.mtype.lower()]()
    if options.fargs is not None:
        dispatcher[options.mtype.lower()](options.fargs)
    else:
        dispatcher[options.mtype.lower()]()

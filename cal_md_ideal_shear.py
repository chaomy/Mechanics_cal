#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-06-25 14:28:58
# @Last Modified by:   chaomy
# @Last Modified time: 2017-07-09 14:05:58


from optparse import OptionParser
import md_pot_data as dat
import cal_md_ideal_shear_ini as init


if __name__ == '__main__':

    usage = "usage:%prog [options] arg1 [options] arg2"
    parser = OptionParser(usage=usage)

    parser.add_option('-t', "--mtype", action="store",
                      type="string", dest="mtype", help="",
                      default="prp_r")

    parser.add_option('-p', "--param", action="store",
                      type='string', dest="fargs",
                      default=None)

    (options, args) = parser.parse_args()
    drv = init.cal_bcc_ideal_shear(dat.qe_pot.vca_W75Re25,
                                   '211')

    dispatcher = {'qeone': drv.get_qe_stress,
                  'restart': drv.loop_prep_restart,
                  'twin': drv.gn_shear_twin_path,
                  'preprestart': drv.loop_prep_restart_from_log,
                  'loopre': drv.loop_prep_restart_from_log,
                  'clcqe': drv.qe_loop_stress,
                  'clcva': drv.va_loop_stress,
                  'clclmp': drv.convert_stress,
                  'prepva': drv.loop_prep_vasp,
                  'iva': drv.vasp_relax,
                  'ivasp': drv.vasp_relax,
                  'ilmp': drv.loop_shear_lmp,
                  'iqe': drv.qe_relax,
                  'gnprim': drv.gn_primitive_lmps,
                  'clctmp': drv.read_ofiles,
                  'pltengy': drv.plt_strain_vs_energy,
                  'pltstress': drv.plt_energy_stress,
                  'pltvc': drv.plt_vc}

    if options.fargs is not None:
        dispatcher[options.mtype.lower()](options.fargs)
    else:
        dispatcher[options.mtype.lower()]()

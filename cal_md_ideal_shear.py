#!/usr/bin/env python
# encoding: utf-8

###################################################################
#
# File Name : ./cal_bcc_ideal_shear.py
#
###################################################################
#
# Purpose :
#
# Creation Date :
# Last Modified : Sat Apr 22 21:22:17 2017
# Created By    : Chaoming Yang
#
###################################################################

from optparse import OptionParser
import cal_md_ideal_tensile_ini


if __name__ == '__main__':

    usage = "usage:%prog [options] arg1 [options] arg2"
    parser = OptionParser(usage=usage)

    parser.add_option('-t', "--mtype", action="store",
                      type="string", dest="mtype", help="",
                      default="prp_r")

    parser.add_option('-c', "--delta", action="store",
                      type='float', dest="delta",
                      default=0.02)

    (options, args) = parser.parse_args()
    drv = cal_md_ideal_tensile_ini.cal_bcc_ideal_shear()
    dispatcher = {'qe_stress': drv.get_qe_stress}
    dispatcher['qe_stress']()

    if options.mtype.lower() in ['clc_vasp', 'clc_lmp', 'clc_qe']:
        opt = options.mtype.lower().split('_')[-1]
        drv.md_ideal_shear('clc', opt)

    if options.mtype.lower() == 'vastress':
        drv.vasp_loop_stress()

    if options.mtype.lower() == 'qestress':
        drv.qe_loop_stress(opt='clc')

    if options.mtype.lower() == 'lmpstress':
        drv.convert_stress()

    if options.mtype.lower() == 'trans':
        drv.trans_stress()

    if options.mtype.lower() == 'twin':
        drv.shear_twin_path()

    if options.mtype.lower() in ['plt_engy', 'plt_cmp', 'plt_stress']:
        opt = options.mtype.lower().split('_')[-1]
        if opt in 'engy':
            drv.plt_strain_vs_energy()
        elif opt in 'stress':
            drv.plt_energy_stress()

    if options.mtype.lower() in ['vaspprep']:
        drv.loop_prep_vasp()

    if options.mtype.lower() in ['qeprep']:
        drv.loop_prep_qe()

    if options.mtype.lower() in ['ivasp', 'iva']:
        drv.vasp_relax()

    if options.mtype.lower() in ['ilmp']:
        drv.loop_shear_lmp()

    if options.mtype.lower() in ['iqe']:
        print drv.qe_relax()

    if options.mtype.lower() == 'cnt':
        drv.clc_data()
        # drv.prep_restart_from_log()

    if options.mtype.lower() == 'tmp':
        # drv.read_ofiles('clctmp')
        drv.read_ofiles('convert')

    if options.mtype.lower() == 'gnqe':
        drv.gn_primitive_lmps(tag='qe')

    if options.mtype.lower() in ['qe_restart', 'va_restart', 'cnt_restart']:
        opt = options.mtype.lower().split('_')[0]
        drv.loop_prep_restart(opt)

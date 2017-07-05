#!/usr/bin/env python
# encoding: utf-8

###################################################################
#
# File Name : ./cal_md_stat_crack.py
#
###################################################################
#
# Purpose :  introduce crack for md calculation (LAMMPS)
#
# Creation Date :
# Last Modified : Sun Apr  9 20:20:57 2016
# Created By    : Chaoming Yang
#
###################################################################

from optparse import OptionParser
import cal_md_crack_ini


if __name__ == "__main__":
    usage = "usage:%prog [options] arg1 [options] arg2"
    parser = OptionParser(usage=usage)

    parser.add_option('-t', "--mtype", action="store",
                      type="string", dest="mtype", help="",
                      default="prp_r")

    parser.add_option('-a', "--delta", action="store",
                      type='string', dest="fargs",
                      default=None)

    (options, args) = parser.parse_args()
    drv = cal_md_crack_ini.md_crack_ini()

    dispatcher = {'stat': drv.static_crack}
    if options.fargs is not None:
        dispatcher[options.mtype.lower()](options.fargs)
    else:
        dispatcher[options.mtype.lower()]()

#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-06-25 14:28:58
# @Last Modified by:   chaomy
# @Last Modified time: 2018-03-27 17:26:28

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
    dispatcher = {'stat': drv.static_crack,
                  'coeff': drv.cal_crack_anglecoeff,
                  'ani211': drv.aniso_211}

    if options.fargs is not None:
        dispatcher[options.mtype.lower()](options.fargs)
    else:
        dispatcher[options.mtype.lower()]()

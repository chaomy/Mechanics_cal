#!/usr/bin/env python
# encoding: utf-8
#
###################################################################
#
# File Name : ./cal_vasp_dislocation.py
#
###################################################################
#
# Purpose :
#
# Creation Date : Tue Apr 11 15:36:33 2017
# Last Modified : Sat Apr  1 23:15:41 2017
# Created By    : Chaoming Yang
#
###################################################################


import os
from optparse import OptionParser

usage = "usage:%prog [options] arg1 [options] arg2"
parser = OptionParser(usage=usage)
parser.add_option("-t", "--mtype", action="store",
                  type="string", dest="mtype", help="",
                  default="curv")
parser.add_option("-f", "--mfile", action="store",
                  type="string", dest="mfile",
                  default="./dummy.config.pair")

(options, args) = parser.parse_args()


def loop_sub_jobs():
    for i in range(20):
        dirname = 'dir-{:03d}'.format(i)
        print dirname
        os.chdir(dirname)
        os.system("qsub va.pbs")
        os.chdir(os.pardir)
    return

if __name__ == "__main__":
    loop_sub_jobs()

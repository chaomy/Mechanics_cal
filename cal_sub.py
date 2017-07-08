#!/usr/bin/env python
# encoding: utf-8
#
###################################################################
#
# File Name : ./cal_sub.py
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
import glob
from optparse import OptionParser

usage = "usage:%prog [options] arg1 [options] arg2"
parser = OptionParser(usage=usage)
parser.add_option("-t", "--mtype", action="store",
                  type="string", dest="mtype", help="",
                  default="curv")
(options, args) = parser.parse_args()


class subjobs(object):

    def __init__(self):
        self.diriter = None
        self.get_dirs()
        return

    def get_dirs(self):
        self.diriter = iter(glob.glob('dir-*'))
        return

    def loop_sub_jobs(self, opt):
        while True:
            try:
                mdir = next(self.diriter)
                print mdir
                self.gonadsub()
            except StopIteration:
                break
        return

    def gonadsub(self, mdir):
        os.chdir(mdir)
        os.system("qsub va.pbs")
        os.chdir(os.pardir)
        return

    def loop_shear_cnt(self):
        while True:
            try:
                mdir = next(self.diriter)
                os.chdir(mdir)
                if not os.path.isfile('ishear.txt'):
                    print mdir
                    os.system('qsub va.pbs')
                os.chdir(os.pardir)
            except StopIteration:
                break
        return

if __name__ == "__main__":
    drv = subjobs()
    dispatcher = {'sub': drv.loop_sub_jobs,
                  'shearcnt': drv.loop_shear_cnt}
    dispatcher[options.mtype.lower()]()

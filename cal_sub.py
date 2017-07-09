#!/usr/bin/env python
# encoding: utf-8
# -*- coding: utf-8 -*-
# @Author: yang37
# @Date:   2017-06-12 17:03:43
# @Last Modified by:   chaomy
# @Last Modified time: 2017-07-09 09:54:37


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

    def loop_sub_jobs(self):
        while True:
            try:
                mdir = next(self.diriter)
                print mdir
                self.gonadsub(mdir)
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

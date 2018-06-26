#!/usr/bin/env python
# encoding: utf-8
# -*- coding: utf-8 -*-
# @Author: yang37
# @Date:   2017-06-12 17:03:43
# @Last Modified by:   chaomy
# @Last Modified time: 2018-06-25 20:43:22


import os
import glob
from optparse import OptionParser


class subjobs(object):

    def __init__(self):
        self.diriter = None
        self.get_dirs()

    def get_dirs(self):
        self.diriter = iter(glob.glob('dir[-_]*'))

    def trans_to(self):
        pth = "/scratch/qiliang_flux/chaomy/MD/Nb/MEAMS/THERMO"
        fls = glob.glob("dir-*")
        for e in fls:
            os.system("scp {}/in.rst $FLUX:{}/{}".format(e, pth, e))

    def cnt_job(self):
        while True:
            try:
                mdir = next(self.diriter)
                print(mdir)
            except StopIteration:

    def loop_sub_jobs(self):
        while True:
            try:
                mdir = next(self.diriter)
                print(mdir)
                self.gonadsub(mdir)
            except StopIteration:
                break

    def gonadsub(self, mdir):
        os.chdir(mdir)
        os.system("qsub va.pbs")
        os.chdir(os.pardir)

    def mpbs(self):
        ml = glob.glob("va*.pbs")
        for ee in ml:
            os.system("qsub {}".format(ee))

    def loop_shear_cnt(self):
        while True:
            try:
                mdir = next(self.diriter)
                os.chdir(mdir)
                if not os.path.isfile('ishear.txt'):
                    print(mdir)
                    os.system('qsub va.pbs')
                os.chdir(os.pardir)
            except StopIteration:
                break


if __name__ == "__main__":
    usage = "usage:%prog [options] arg1 [options] arg2"
    parser = OptionParser(usage=usage)
    parser.add_option("-t", "--mtype", action="store",
                      type="string", dest="mtype", help="",
                      default="curv")
    (options, args) = parser.parse_args()
    drv = subjobs()
    dispatcher = {'sub': drv.loop_sub_jobs,
                  'shearcnt': drv.loop_shear_cnt,
                  'mpbs': drv.mpbs,
                  'to': drv.trans_to,
                  'cnt': drv.cnt_job}
    dispatcher[options.mtype.lower()]()

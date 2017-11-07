#!/usr/bin/env python
# encoding: utf-8

###################################################################
#
# File Name : cal_md_gb_MgNd.py
#
###################################################################
#
# Purpose :
#
# Creation Date : Thu Apr 20 01:51:02 2017
# Last Modified :
# Created By    : Chaoming Yang
#
###################################################################

import os
from optparse import OptionParser
import gn_lmp_infile
import gn_pbs
import numpy as np
import shutil
import Intro_vasp
import matplotlib.pylab as plt
from multiprocessing import Pool

try:
    import ase

except ImportError:
    print("error during import")


def unwrap_self_run_lammps(arg, **kwarg):
    return md_loop_tensile.lammps_job(*arg, **kwarg)


class cal_MgNd_driver(Intro_vasp.vasp_change_box,
                      gn_lmp_infile.gn_md_infile,
                      gn_pbs.gn_pbs):
    def __init__(self):
        self.hcp_a = 3.2019267694893
        self.hcp_c = 5.1969105399
        self.D03_a = 7.4662780330786
        gn_pbs.gn_pbs.__init__(self)
        gn_lmp_infile.gn_md_infile.__init__(self)
        self.root = os.getcwd()
        return

    def cal_mesh(self):
        inparam = gn_lmp_infile.gb_param
        inparam.ydisp = 0
        inparam.zdisp = 0

        # shift along x direction [1 1 -2 0] || [-1  1  -1]
        npts = 10
        delta_d = self.hcp_a / npts
        for i in range(npts):
            dirname = "dir-x-%03d" % (i)
            if not os.path.isdir(dirname):
                os.mkdir(dirname)
            if not os.path.isdir("%s/out" % (dirname)):
                os.mkdir("%s/out" % (dirname))
            inparam.xdisp = i * delta_d
            self.write_mgnd_infile("in.init", inparam)
            shutil.copy2("in.init", dirname)
            shutil.copy2("lib_MgNdPb.meam", dirname)
            shutil.copy2("MgNd_para.meam", dirname)
        return

    def run_lmp(self, jobid):
        dirname = "dir-x-%03d" % (jobid)
        os.chdir(dirname)
        os.system("mpirun lmp_linux -i in.init")
        os.chdir(self.root)
        return

    def loop_run_lmp(self, tag):
        npts = 10
        gbengylist = []
        deltalist = []
        for i in range(npts):
            dirname = "dir-x-%03d" % (i)
            deltalist.append(i * self.hcp_a / npts)
            os.chdir(dirname)
            if tag == 'run':
                os.system("mpirun lmp_linux -i in.init")
            if tag == 'grab':
                gbengylist.append(self.grab_data())
            os.chdir(self.root)
        print deltalist
        print gbengylist

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(deltalist, gbengylist)
        plt.show()
        return

    def grab_data(self):
        with open("log.lammps", 'r') as fid:
            raw = fid.readlines()
            gbengy = float(raw[-3].split()[-1])
        return gbengy


usage = "usage:%prog [options] arg1 [options] arg2"
parser = OptionParser(usage=usage)
parser.add_option("-t", "--mtype", action="store",
                  type="string", dest="mtype", help="",
                  default="prp_r")

(options, args) = parser.parse_args()

if __name__ == '__main__':
    drv = cal_MgNd_driver()
    if options.mtype.lower() == "prep":
        drv.cal_mesh()
    if options.mtype.lower() == "run":
        drv.run_lmp()
    if options.mtype.lower() == "loop":
        drv.loop_run_lmp()
    if options.mtype.lower() == "grab":
        drv.loop_run_lmp('grab')

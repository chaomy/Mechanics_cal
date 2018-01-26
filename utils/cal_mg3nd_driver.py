#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2018-01-15 23:30:10
# @Last Modified by:   chaomy
# @Last Modified time: 2018-01-16 18:12:56

import os
import glob
import gn_pbs
import get_data
from optparse import OptionParser


class cal_mg3nd_driver(gn_pbs.gn_pbs,
                       get_data.get_data):

    def __init__(self):
        gn_pbs.gn_pbs.__init__(self)
        get_data.get_data.__init__(self)

    def set_pbs(self, dirname, opt='va'):
        self.set_nnodes(1)
        self.set_ppn(12)
        self.set_job_title("{}".format(dirname))
        self.set_wall_time(2)
        self.set_main_job("""mpirun vasp""")
        self.write_pbs(od=True)
        return

    def cal_lat(self):
        for i in range(21):
            fname = "poscar.{:03d}".format(i)
            mdir = "dir_{:03d}".format(i)
            # os.system("cp INCAR {}".format(mdir))
            # os.system("cp KPOINTS {}".format(mdir))
            # os.system("cp POTCAR {}".format(mdir))
            # os.system("cp {}/{} {}/POSCAR".format(mdir, fname, mdir))
            # # self.mymkdir(mdir)
            os.chdir(mdir)
            os.system("qsub va.pbs")
            os.chdir(os.pardir)
            # self.set_pbs("{:03d}".format(i))
            # os.system("mv {} {}".format(fname, mdir))
            # os.system("mv va.pbs {}".format(mdir))
        return


if __name__ == '__main__':
    usage = "usage:%prog [options] arg1 [options] arg2"
    parser = OptionParser(usage=usage)
    parser.add_option('-t', "--mtype", action="store",
                      type="string", dest="mtype")
    parser.add_option('-p', "--param", action="store",
                      type='string', dest="fargs")
    (options, args) = parser.parse_args()
    drv = cal_mg3nd_driver()
    dispatcher = {'lat': drv.cal_lat}

    if options.fargs is not None:
        dispatcher[options.mtype.lower()](options.fargs)
    else:
        dispatcher[options.mtype.lower()]()

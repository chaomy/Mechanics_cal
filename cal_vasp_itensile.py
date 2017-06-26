#!/usr/bin/env python
# encoding: utf-8

###################################################################
#
# File Name : cal_vasp_itensile.py
#
###################################################################
#
# Purpose :
#
# Creation Date :
# Last Modified :
# Created By    : Chaoming Yang
#
###################################################################
import get_data
import md_pot_data
import os
import numpy as np
import gn_incar
import gn_pbs
import gn_kpoints
import plt_drv
import matplotlib.pylab as plt
from optparse import OptionParser


class cal_bcc_ideal_tensile(get_data.get_data,
                            plt_drv.plt_drv,
                            gn_incar.gn_incar,
                            gn_kpoints.gn_kpoints,
                            gn_pbs.gn_pbs):

    def __init__(self):
        self._pot = md_pot_data.dft_data.Nb_pbe
        self.root_dir = os.getcwd()
        gn_pbs.gn_pbs.__init__(self)
        gn_incar.gn_incar.__init__(self)
        gn_kpoints.gn_kpoints.__init__(self)
        get_data.get_data.__init__(self)
        plt_drv.plt_drv.__init__(self)
        return

    def set_pbs(self, dirname, delta, opt='vasp'):
        self.set_nnodes(2)
        self.set_ppn(12)
        self.set_job_title("%s" % (dirname))
        self.set_wall_time(8)
        self.set_main_job("""
            cal_md_ideal_shear.py  -t  i{}
                          """.format(opt))
        self.write_pbs(od=True)
        os.system("mv va.pbs %s" % (dirname))
        return

    def setup_recal(self, opt='prep'):
        self.set_intype('scf')
        for i in range(30):
            delta = 0.01 * i
            mdir = 'dir-{:4.3f}'.format(delta)
            if opt == 'prep':
                poscar = 'POSCAR{:4.3f}'.format(delta)
                self.mymkdir(mdir)
                self.copy_inputs(mdir, 'KPOINTS',
                                 'INCAR', 'POTCAR')
                os.system("cp {} {}/POSCAR".format(poscar, mdir))
            elif opt == 'run':
                os.chdir(mdir)
                os.system("mpirun vasp > vasp.log")
                os.chdir(os.pardir)
        return

    def grab_engy(self):
        npts = 26
        data = np.ndarray([npts, 9])
        for i in range(26):
            dirname = "dir-{:03d}".format(i)
            print dirname
            os.chdir(dirname)
            delta = 0.01 * i
            engy, stress, vol = self.vasp_energy_stress_vol()
            (data[i, 0], data[i, 1], data[i, 2:8], data[i, -1]) = \
                delta, engy, stress.transpose(), vol
            os.chdir(self.root_dir)
        print data
        np.savetxt("iten.txt", data)
        return

    def plot_curv(self):
        data = np.loadtxt("iten.txt")
        self.set_keys()
        self.set_111plt()
        self.ax.plot(data[:, 0], data[:, 1])
        self.fig.savefig("istress.png", **self.figsave)
        return


if __name__ == '__main__':
    usage = "usage:%prog [options] arg1 [options] arg2"
    parser = OptionParser(usage=usage)
    parser.add_option('-t', "--mtype",
                      action="store",
                      type="string",
                      dest="mtype", help="",
                      default="prp_r")
    (options, args) = parser.parse_args()

    drv = cal_bcc_ideal_tensile()

    if options.mtype.lower() == 'clc':
        drv.grab_engy()

    if options.mtype.lower() == 'plt':
        drv.plot_curv()

    if options.mtype.lower() == 'recal':
        drv.setup_recal()

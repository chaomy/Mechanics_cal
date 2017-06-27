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
import ase.io
import ase
from optparse import OptionParser


class cal_bcc_ideal_tensile(get_data.get_data,
                            plt_drv.plt_drv,
                            gn_incar.gn_incar,
                            gn_kpoints.gn_kpoints,
                            gn_pbs.gn_pbs):

    def __init__(self):
        self._pot = md_pot_data.dft_pot.Nb_pbe
        self.root_dir = os.getcwd()
        gn_pbs.gn_pbs.__init__(self)
        gn_incar.gn_incar.__init__(self)
        gn_kpoints.gn_kpoints.__init__(self)
        get_data.get_data.__init__(self)
        plt_drv.plt_drv.__init__(self)
        return

    def adjust(self):
        for i in range(0, 41):
            # dirname = 'dir-{:03}'.format(i)
            mdir = 'dir-{:4.3f}'.format(0.01 * i)
            os.system('cp {}/CONTCAR CONTCAR{:4.3f}'.format(mdir, 0.01 * i))
            # os.system('mv {} {}'.format(dirname, mdir))
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

    def grab_engy(self,  opt='engy'):
        npts = 41
        data = np.ndarray([npts, 9])
        for i in range(npts):
            # dirname = "dir-{:03d}".format(i)
            delta = 0.01 * i
            if opt == 'engy':
                dirname = 'dir-{:4.3f}'.format(delta)
                print dirname
                os.chdir(dirname)
                engy, stress, vol = self.vasp_energy_stress_vol()
                (data[i, 0], data[i, 1], data[i, 2:8], data[i, -1]) = \
                    delta, engy, stress.transpose(), vol
                os.chdir(self.root_dir)
            elif opt == 'cell':
                fname = 'POSCAR{:4.3f}'.format(delta)
                atoms = ase.io.read(fname,
                                    format='vasp')
                print atoms.get_cell()
        print data
        np.savetxt("iten.txt", data)
        return

    def plot_curv(self):
        data = np.loadtxt("iten.txt")
        self.set_keys()
        self.set_111plt()
        self.ax.plot(data[:, 0], data[:, 1], marker='o')
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

    if options.mtype.lower() in ['clc_engy', 'clc_cell']:
        tag = options.mtype.lower().split('_')[-1]
        drv.grab_engy(tag)

    if options.mtype.lower() in ['adj']:
        drv.adjust()

    if options.mtype.lower() == 'plt':
        drv.plot_curv()

    if options.mtype.lower() == 'recal_run':
        opt = options.mtype.lower().split('_')[-1]
        drv.setup_recal(opt)

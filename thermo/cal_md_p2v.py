#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: yang37
# @Date:   2017-06-21 18:42:47
# @Last Modified by:   chaomy
# @Last Modified time: 2018-04-29 03:44:48


from itertools import cycle
from optparse import OptionParser
from glob import glob
import matplotlib.pylab as plt
import ase
import ase.io
import os
import numpy as np
import ase.lattice
import shutil
import gn_config
import get_data
import plt_drv
import md_pot_data


class cal_md_thermo(gn_config.gnStructure,
                    get_data.get_data,
                    plt_drv.plt_drv):

    def __init__(self):
        self.pot = self.load_data("../BASICS/pot.dat")
        # self.pot = md_pot_data.va_pot.Nb_pbe
        # gn_lmp_infile.gn_md_infile.__init__(self)
        gn_config.gnStructure.__init__(self, self.pot)
        plt_drv.plt_drv.__init__(self)

    def pressure_vs_vol(self, opt='prep'):
        delta = -0.01
        npts = 30
        bas = np.mat([[-0.5, 0.5, 0.5],
                      [0.5, -0.5, 0.5],
                      [0.5, 0.5, -0.5]])
        lmp_bas = self.lmp_change_box(bas)

        if opt == 'clc':
            vol = np.zeros(npts)
            press = np.zeros(npts)

        for i in range(npts):
            rat = (1 + i * delta)**(1. / 3.)
            alat = rat * self.pot["lattice"]
            dirname = 'dir-%04d' % (i)
            if opt == 'prep':
                self.mymkdir(dirname)
                cell = alat * lmp_bas
                atoms = ase.Atoms(self.pot['element'],
                                  positions=[[0, 0, 0]],
                                  cell=cell,
                                  pbc=[1, 1, 1])
                lmp_bas = self.lmp_change_box(bas)
                self.write_lmp_config_data(atoms, 'init.txt')
                os.system("mv init.txt  {}".format(dirname))
                os.system("cp in.pv  {}".format(dirname))

            elif opt == 'run':
                os.chdir(dirname)
                os.system("lmp_mpi -i in.pv")
                os.chdir(os.pardir)

            elif opt == 'clc':
                os.chdir(dirname)
                data = np.loadtxt("out.txt")
                print(data)
                vol[i] = data[0]
                press[i] = data[1]
                os.chdir(os.pardir)
                shutil.rmtree(dirname)    # clean
        if opt == 'clc':
            np.savetxt("data.txt", (vol, press))

    def loop_pressure_vs_vol(self):
        dirlist = glob("dir-*")
        cnt = 0
        for mdir in dirlist[:]:
            print(mdir)
            if ((cnt % 1) == 0):
                if not os.path.isfile("fig-{}.png".format(mdir)):
                    self.pressure_vs_vol('prep')
                    self.pressure_vs_vol('run')
                    self.pressure_vs_vol('clc')
                    self.vasp_energy_stress_vol_plt(mdir, 30)
                    os.remove("dummy.lammps.ADP")
                    os.system("mv data.txt  {}".format(mdir))
                    os.system("cp p2v.png fig-{}.png".format(mdir))
            cnt += 1

    def vasp_energy_stress_vol_plt(self,
                                   inlabel='p-v',
                                   npt=30):
        self.set_keys()
        self.set_111plt()
        (vol, press) = np.loadtxt("data.txt")
        (dft_vol, dft_press) = np.loadtxt(
            "/Users/chaomingyang/src/Data_shares/DATA_DFT_PV.txt")
        vol = vol / vol[0]
        vol = vol**3
        dft_vol = dft_vol / dft_vol[0]

        self.ax.plot(vol[:npt], press[:npt], label='MEAM', **next(self.keysiter))
        self.ax.plot(dft_vol[:npt], dft_press[:npt], label='PAW-PBE',
                     **next(self.keysiter))

        self.add_legends(*self.axls)
        self.set_tick_size(*self.axls)
        self.add_y_labels(cycle(['Presssure (GPa)']), *self.axls)
        self.add_x_labels(cycle(['Relative volume (V / V$_0$)']), *self.axls)
        self.fig.savefig("fig-p2v.png", **self.figsave)

    def p2v_wrap(self):
        drv.pressure_vs_vol('prep')
        drv.pressure_vs_vol('run')
        drv.pressure_vs_vol('clc')
        drv.vasp_energy_stress_vol_plt(25)


if __name__ == '__main__':
    usage = "usage:%prog [options] arg1 [options] arg2"
    parser = OptionParser(usage=usage)
    parser.add_option("-t", "--mtype", action="store",
                      type="string", dest="mtype")
    parser.add_option('-p', "--param", action="store",
                      type='string', dest="fargs")

    (options, args) = parser.parse_args()
    drv = cal_md_thermo()
    dispatcher = {'auto': drv.p2v_wrap,
                  'loop': drv.loop_pressure_vs_vol,
                  'plt': drv.vasp_energy_stress_vol_plt}

    if options.fargs is not None:
        dispatcher[options.mtype.lower()](options.fargs)
    else:
        dispatcher[options.mtype.lower()]()

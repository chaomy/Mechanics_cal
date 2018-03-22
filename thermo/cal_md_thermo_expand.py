#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: yang37
# @Date:   2017-06-21 18:42:47
# @Last Modified by:   chaomy
# @Last Modified time: 2018-03-20 13:56:51


from optparse import OptionParser
import matplotlib.pylab as plt
import ase
import ase.io
import os
import numpy as np
import ase.lattice
import shutil
import gn_config
import get_data
import gn_lmp_infile
import gn_pbs
import plt_drv
import md_pot_data
from glob import glob


class cal_md_thermo(gn_config.gnStructure,
                    get_data.get_data,
                    gn_pbs.gn_pbs,
                    plt_drv.plt_drv,
                    gn_lmp_infile.gn_md_infile):

    def __init__(self):
        self.pot = self.load_data("../BASICS/pot.dat")
        # self.pot = md_pot_data.va_pot.Nb_pbe
        self.size = np.array([16, 16, 16])
        gn_lmp_infile.gn_md_infile.__init__(self, self.pot)
        gn_config.gnStructure.__init__(self, self.pot)
        plt_drv.plt_drv.__init__(self)

    def run_thermo(self, tag='run'):
        if tag == 'run':
            for i in range(1, 51):
                tend = 50 * i
                dirname = "dir-{:05.0f}".format(tend)
                os.chdir(dirname)
                os.system("mpirun -n 24 lmp_linux -i in.npt")
                os.chdir(os.pardir)
        elif tag == 'loopclc':
            temp_lx = []
            for i in range(1, 51):
                tend = 50 * i
                dirname = "dir-{:05.0f}".format(tend)
                os.chdir(dirname)
                temp_lx.append(self.get_temp_lat())
                os.chdir(os.pardir)
            print(temp_lx)
        temp_lx = np.array(temp_lx)
        np.savetxt("temp_lx.txt", temp_lx)

    def draw_temp_vol(self):
        self.set_figs()
        datan = 10
        lastline = 31
        for i in range(1, 51):
            tend = 50 * i
            dirname = "dir-{:05.0f}".format(tend)
            print(dirname)
            print(tend)
            if os.path.isfile("{}/log.lammps".format(dirname)):
                raw = self.mreadlines("{}/log.lammps".format(dirname))
                for j in range(lastline, lastline + datan):
                    print(raw[-j])

    def given_temp_prep(self):
        atoms = self.set_bcc_convention().repeat(([15, 15, 15]))
        self.write_lmp_config_data(atoms, "init.txt")
        oneatm = 1.01325
        tstart = 0.1
        initt = 50
        deltat = 50
        # initial
        dirname = "dir-{:05.0f}".format(0)
        self.mymkdir(dirname)
        self.write_md_thermo_expand('init')
        shutil.copy2("in.npt", dirname)
        # increase temp
        for i in range(50):
            # 2500  K
            tend = initt + deltat * i
            dirname = "dir-{:05.0f}".format(tend)
            self.mymkdir(dirname)
            if i == 0:
                self.write_md_thermo_expand(**{'tstart': tstart,
                                               'tend': tend,
                                               'pstart': 0.0,
                                               'pend': oneatm})
            else:
                self.write_md_thermo_expand(**{'tstart': tstart,
                                               'tend': tend,
                                               'pstart': oneatm,
                                               'pend': oneatm})
            tstart = tend
            shutil.copy2("in.npt", dirname)

    def get_temp_lat(self, filename="log.lammps"):
        raw = self.mreadlines(filename)
        temp, lx, cnt = 0, 0, 0
        for i in range(len(raw)):
            line = raw[i]
            seg = line.split()
            if len(seg) > 0:
                if seg[0] == 'Step':
                    if int(raw[i + 1].split()[0]) > 1:
                        for j in range(522, 602):          # 600 lines
                            segs = raw[i + j].split()
                            temp += float(segs[3])
                            lx += float(segs[-1])
                            cnt += 1
                        break
        temp /= cnt
        lx /= cnt
        print(temp, lx)
        return (temp, lx)

    def theormo_expand_plt(self):
        temp_lx = np.loadtxt("temp_lx.txt")
        print(temp_lx)

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
        self.set_keys("upper right")
        self.set_111plt((10, 6.5))
        (vol, press) = np.loadtxt("data.txt")
        (dft_vol, dft_press) = np.loadtxt(
            "/Users/chaomingyang/src/Data_shares/DATA_DFT_PV.txt")

        vol = vol / vol[0]
        vol = vol**3
        dft_vol = dft_vol / dft_vol[0]

        self.ax.plot(vol[:npt], press[:npt],
                     label='adp',
                     **next(self.keysiter))
        self.ax.plot(dft_vol[:npt], dft_press[:npt],
                     label='pbe',
                     **next(self.keysiter))
        plt.xlabel('relative volume (V / V$_0$)',
                   {'fontsize': self.myfontsize})
        plt.ylabel('presssure (GPa)', {'fontsize': self.myfontsize})
        self.add_legends(self.ax)
        self.fig.savefig("p2v.png", **self.figsave)

    def p2v_wrap(self):
        drv.pressure_vs_vol('prep')
        drv.pressure_vs_vol('run')
        drv.pressure_vs_vol('clc')
        drv.vasp_energy_stress_vol_plt(25)

    def add_vol_expan(self):
        atoms = ase.io.read("dump", format='lammps-dump')
        cell = atoms.get_cell()
        atoms.set_cell(1.0035 * cell, scale_atoms=True)
        self.write_lmp_config_data(atoms, "thermo.txt")
        return atoms


if __name__ == '__main__':
    usage = "usage:%prog [options] arg1 [options] arg2"
    parser = OptionParser(usage=usage)
    parser.add_option("-t", "--mtype", action="store",
                      type="string", dest="mtype")
    parser.add_option('-p', "--param", action="store",
                      type='string', dest="fargs")

    (options, args) = parser.parse_args()
    drv = cal_md_thermo()
    dispatcher = {'thermo': drv.given_temp_prep,
                  'runthermo': drv.run_thermo,
                  'clc': drv.run_thermo,
                  'p2vauto': drv.p2v_wrap,
                  'loopp2v': drv.loop_pressure_vs_vol,
                  'p2vplt': drv.vasp_energy_stress_vol_plt,
                  'plt': drv.draw_temp_vol,
                  'themoplt': drv.theormo_expand_plt,
                  'add': drv.add_vol_expan}

    if options.fargs is not None:
        dispatcher[options.mtype.lower()](options.fargs)
    else:
        dispatcher[options.mtype.lower()]()

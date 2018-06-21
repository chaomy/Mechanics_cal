#!/usr/bin/env python
#@Author : yang37
#@Date : 2017 - 06 - 21 18 : 42 : 47
#@Last Modified by : chaomy
#@Last Modified time : 2018 - 06 - 19 15 : 00 : 00

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
import gn_lmp_infile
import plt_drv
import md_pot_data


class cal_md_thermo(gn_config.gnStructure,
                    gn_lmp_infile.gn_md_infile,
                    get_data.get_data,
                    plt_drv.plt_drv):

    def __init__(self):
        #self.pot = self.load_data("../BASICS/pot.dat")
        self.pot = md_pot_data.md_pot.mg_Poco
#self.pot = md_pot_data.va_pot.Nb_pbe
        self.size = np.array([16, 16, 16])
        gn_lmp_infile.gn_md_infile.__init__(self)
        gn_config.gnStructure.__init__(self, self.pot)
        plt_drv.plt_drv.__init__(self)

    def run_thermo(self, tag='run'):
        if tag == 'run':
            for i in range(1, 51):
                tend = 50 * i
                mdir = "dir-{:05.0f}".format(tend)
                os.chdir(mdir)
                os.system("mpirun -n 24 lmp_linux -i in.npt")
                os.chdir(os.pardir)
        elif tag == 'loopclc':
            temp_lx = []
            for i in range(1, 51):
                tend = 50 * i
                mdir = "dir-{:05.0f}".format(tend)
                os.chdir(mdir)
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
            mdir = "dir-{:05.0f}".format(tend)
            print(mdir)
            print(tend)
            if os.path.isfile("{}/log.lammps".format(mdir)):
                raw = self.mreadlines("{}/log.lammps".format(mdir))
                for j in range(lastline, lastline + datan):
                    print(raw[-j])

    def given_temp_prep(self):
        # start with a guess thermoal coeffcient
        thermocoeff = 0.02 / 2000
        oneatm = 1.01325
        tstart = 0.1
        initt = 50
        deltat = 50
        lat = self.pot["lattice"]

# 50 to 2500 K
        for i in range(50):
            tend = initt + deltat * i
            self.pot["lattice"] = lat * (1 + tend * thermocoeff)
            atoms = self.set_bcc_convention().repeat(([15, 15, 15]))
            mdir = "dir-{:05.0f}".format(tend)
            self.mymkdir(mdir)
            self.write_lmp_config_data(atoms, mdir + "/lmp_init.txt")

            print(self.pot["lattice"], tend)
            self.write_md_thermo_expand(**{'tstart': tstart,
                                           'tend': tend,
                                           'pstart': 0.0,
                                           'pend': oneatm})
            tstart = tend
            shutil.copy2("in.npt", mdir)

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

    def get_lat_at_given_temp(self):
        data = np.loadtxt("log")
# Step TotEng Temp Lx Ly Lz v_S11 v_S22 v_S33 v_S12 v_S13 v_S23
        print(np.mean(data[-200:, 3]) / 7.0,
              np.mean(data[-200:, 4]) / 7.0,
              np.mean(data[-200:, 5]) / 7.0)

    def theormo_expand_plt(self):
        temp_lx = np.loadtxt("temp_lx.txt")
        print(temp_lx)

    def add_vol_expan(self):
        atoms = ase.io.read("dump", format='lammps-dump')
        cell = atoms.get_cell()
        atoms.set_cell(1.0040 * cell, scale_atoms=True)
        self.write_lmp_config_data(atoms)
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
    dispatcher = {'prep': drv.given_temp_prep,
                  'runthermo': drv.run_thermo,
                  'clc': drv.run_thermo,
                  'plt': drv.draw_temp_vol,
                  'themoplt': drv.theormo_expand_plt,
                  'add': drv.add_vol_expan,
                  'temp': drv.get_lat_at_given_temp}

    if options.fargs is not None:
        dispatcher[options.mtype.lower()](options.fargs)
    else:
        dispatcher[options.mtype.lower()]()

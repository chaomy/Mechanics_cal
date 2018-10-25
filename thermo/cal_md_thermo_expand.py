#!/usr/bin/env python
#@Author : yang37
#@Date : 2017 - 06 - 21 18 : 42 : 47
#@Last Modified by : chaomy
#@Last Modified time : 2018 - 06 - 19 15 : 00 : 00

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
import gn_lmp_infile
import plt_drv
import md_pot_data


class cal_md_thermo(gn_config.gnStructure, gn_lmp_infile.gn_md_infile,
                    get_data.get_data, plt_drv.plt_drv):

    def __init__(self):
        # self.pot = self.load_data("../BASICS/pot.dat")
        self.pot = md_pot_data.md_pot.mg_Poco
        # self.pot = md_pot_data.va_pot.Nb_pbe
        # self.size = np.array([16, 16, 16])
        self.size = np.array([8, 8, 8])
        gn_lmp_infile.gn_md_infile.__init__(self)
        gn_config.gnStructure.__init__(self, self.pot)
        plt_drv.plt_drv.__init__(self)

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
        # start with a guess thermal coeffcient
        thermocoeff = 0.02 / 2000
        oneatm = 1.01325
        tstart = 0.1
        initt = 100
        deltat = 100
        lat = self.pot["lattice"]
        for i in range(25):  # 50 to 2500 K
            tend = initt + deltat * i
            self.pot["lattice"] = lat * (1 + tend * thermocoeff)
            atoms = self.set_bcc_convention().repeat(([8, 8, 8]))
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
            shutil.copy2("in.rst", mdir)
            shutil.copy2("va.pbs", mdir)

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

    def get_lat_at_given_temp_HCP(self):
        data = np.loadtxt("log.out")
        print(len(data))
        for N in [1000]:
            I = 7
            lx, ly, lz = np.mean(
                data[-N:, 3]), np.mean(data[-N:, 4]), np.mean(data[-N:, 5])
            print(lx / I, ly / I, lz / I)

    def get_lat_at_given_temp(self):
        data = np.loadtxt("log.out")
        # Step TotEng Temp Lx Ly Lz v_S11 v_S22 v_S33 v_S12 v_S13 v_S23
        # Step TotEng PotEng Temp Press Pxx Pyy Pzz Lx Ly Lz. Nb
        Nodes = 200
        N = 12
        lat = np.mean([np.mean(data[-Nodes:, -3]) / N,
                       np.mean(data[-Nodes:, -2]) / N,
                       np.mean(data[-Nodes:, -1]) / N])
        lat0 = 3.32237981449018
        # print(lat, 100 * (lat / self.pot['lattice'] - 1.0))
        print(lat, 100 * (lat / lat0 - 1.0))
        return(lat, (lat / lat0 - 1.0))

    def add_vol_expan(self):
        atoms = ase.io.read("dump", format='lammps-dump')
        cell = atoms.get_cell()
        atoms.set_cell(1.0040 * cell, scale_atoms=True)
        self.write_lmp_config_data(atoms)
        return atoms

    def trans(self):
        pth = "$FLUX:/scratch/qiliang_flux/chaomy/MD/Nb/MEAMS/THERMO"
        initt = 50
        deltat = 50
        data = np.ndarray([51, 3])
        c = 1
        data[0] = 0.0, 3.32237981449018, 0.000000
        for i in range(50):  # 50 to 2500 K
            tend = initt + deltat * i
            mdir = "dir-{:05.0f}".format(tend)
            # os.system("scp {}/{}/log.lammps {}".format(pth, mdir, mdir))
            if os.path.isdir(mdir):
                print(mdir)
                os.chdir(mdir)
                # os.system("../a.out")
                data[c][0] = tend
                data[c][1:] = self.get_lat_at_given_temp()
                os.chdir(os.pardir)
                c += 1
        data = data[:c]
        np.savetxt("data.txt", data)

    def plot(self):
        self.set_keys()
        self.set_111plt()
        data = np.loadtxt("data.txt")
        self.ax.plot(data[:, 0], 1e2 * data[:, -1],
                     label='MEAM', **next(self.keysiter))
        next(self.keysiter)
        p = np.polyfit(data[:, 0] , data[:, -1], 2)
        print(p)
        self.ax.plot(data[:, 0], 1e2 * (p[0] * data[:, 0]**2 + p[1] * data[:, 0] + p[2]),
                     label='FITTING', linestyle='--', lw=2.5)
        # self.ax.plot(dft_vol[:npt], dft_press[:npt], label='PAW-PBE', **next(self.keysiter))
        self.add_legends(*self.axls)
        self.set_tick_size(*self.axls)
        self.add_x_labels(cycle(['Temperature [K]']), *self.axls)
        self.add_y_labels(cycle([r'$\Delta$ $a_0$ / $a_0$ [%]']), *self.axls)
        self.fig.savefig("FIG_THERMO_EXPAND.png", **self.figsave)

    def check_exp(self):
        self.set_keys()
        self.set_111plt()
        # this is Mo !!!
        data = np.array([[0.0, 291, 1126, 1162, 1327, 1510, 1639, 1819, 1968, 2073],
                         [3.142277, 3.1474, 3.1621, 3.1623, 3.1668, 3.1706, 3.1738, 3.1782, 3.1832, 3.1869]])
        # p = np.polyfit(data[0, :2], data[1, :2], 1)
        # print(p)
        coeff = (data[1, :] - data[1, 0]) / (data[1, 0])
        self.ax.plot(data[0, :], coeff, **next(self.keysiter))
        self.fig.savefig("FIG_EXP_EXPAND.png", **self.figsave)

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
                  'plt': drv.draw_temp_vol,
                  'add': drv.add_vol_expan,
                  'temp': drv.get_lat_at_given_temp,
                  'hcp': drv.get_lat_at_given_temp_HCP,
                  'trans': drv.trans,
                  'plt': drv.plot,
                  'exp': drv.check_exp}

    if options.fargs is not None:
        dispatcher[options.mtype.lower()](options.fargs)
    else:
        dispatcher[options.mtype.lower()]()

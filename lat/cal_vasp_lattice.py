#!/usr/bin/env python
# encoding: utf-8
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-07-05 08:12:30
# @Last Modified by:   chaomy
# @Last Modified time: 2018-03-12 22:04:04


import ase.lattice.hexagonal as Hexagonal
import ase.io
import numpy as np
import glob
from scipy.interpolate import interp1d
import matplotlib.pylab as plt
from optparse import OptionParser
import os
import ase.io
import ase
import gn_config
import get_data
import gn_kpoints
import gn_incar
import gn_pbs
import md_pot_data
import plt_drv


class cal_lattice(gn_config.bcc,
                  gn_config.fcc,
                  gn_config.hcp,
                  gn_pbs.gn_pbs,
                  get_data.get_data,
                  gn_kpoints.gn_kpoints,
                  gn_incar.gn_incar,
                  plt_drv.plt_drv):

    def __init__(self):
        self.figsize = (8, 6)
        self.npts = 20
        self.kpoints = [31, 31, 31]
        self.pot = md_pot_data.va_pot.Nb_pbe
        plt_drv.plt_drv.__init__(self)
        gn_kpoints.gn_kpoints.__init__(self)
        get_data.get_data.__init__(self)
        gn_incar.gn_incar.__init__(self)
        gn_pbs.gn_pbs.__init__(self)

    def interpolate_lattice(self, filename):
        data = np.loadtxt(filename)
        print data[:, 1]

        lattice = np.abs(data[:, 0])
        energy = data[:, 1]

        InterPoints = np.linspace(lattice[0], lattice[-1], 101)
        f = interp1d(lattice, energy)
        Ynew = f(InterPoints)
        i = np.argmin(Ynew)

        print "min energy ", np.min(energy)
        print "num", np.argmin(energy), "lattice", InterPoints[i]

        print(np.min(Ynew))

        self.set_111plt()
        self.ax.plot(lattice, energy, label="lat = {:5f}".format(
            2 * InterPoints[i]), **next(self.keysiter))
        self.add_legends(self.ax)

        plt.savefig("lattice.png")

    def set_pbs(self, mdir):
        self.set_pbs_type('va')
        self.set_wall_time(10)
        self.set_job_title(mdir)
        self.set_nnodes(1)
        self.set_ppn(12)
        self.set_main_job("mpirun vasp")
        self.write_pbs(od=True)

    def loop_kpts(self):
        files = glob.glob('./DATA*')
        for file in files:
            cal_lattice(file)

    def prepare_vasp_inputs(self, mdir):
        self.set_pbs(mdir)
        os.system("mv POSCAR {}".format(mdir))
        os.system("mv va.pbs {}".format(mdir))
        os.system("cp KPOINTS {}".format(mdir))
        os.system("cp INCAR {}".format(mdir))
        os.system("cp POTCAR {}".format(mdir))

    def gn_bcc(self):
        alat0 = self.pot['latbcc']
        delta = 0.001
        # rng = [-20, 20]
        # rng = [20, 100]
        rng = [-500, -400]
        # rng = [400, 500]
        gn_config.bcc.__init__(self, self.pot)
        for i in range(rng[0], rng[1]):
            self.pot["latbcc"] = alat0 + i * delta
            if i >= 0:
                mdir = "dir-p-{:03d}".format(i)
            else:
                mdir = "dir-n-{:03d}".format(abs(i))
            self.mymkdir(mdir)
            atoms = self.set_bcc_primitive((1, 1, 1))
            ase.io.write(filename="POSCAR", images=atoms, format='vasp')
            self.prepare_vasp_inputs(mdir)

    def gn_fcc(self):
        alat0 = self.pot["latfcc"]
        delta = 0.005
        rng = [-50, 50]
        gn_config.fcc.__init__(self, self.pot)
        for i in range(rng[0], rng[1]):
            self.pot["latfcc"] = alat0 + i * delta
            if i >= 0:
                mdir = "dir-p-{:03d}".format(i)
            else:
                mdir = "dir-n-{:03d}".format(abs(i))
            self.mymkdir(mdir)
            atoms = self.set_fcc_primitive((1, 1, 1))
            ase.io.write(filename="POSCAR", images=atoms, format='vasp')
            self.prepare_vasp_inputs(mdir)

    def gn_hcp_mesh(self):
        da = 0.01
        dc = 0.02
        shift = -da * (20 + 45)
        shift = -dc * 8
        for i in range(20):
            for j in range(16):
                mdir = "dir-{:03d}-{:03d}".format(i, j)

                aa = self.pot["ahcp"] + shift + da * i
                cc = self.pot["chcp"] + shift + dc * j
                atoms = Hexagonal.HexagonalClosedPacked(
                    latticeconstant={'a': aa, 'c': cc},
                    size=(1, 1, 1),
                    symbol=self.pot["element"],
                    pbc=(1, 1, 1))
                self.mymkdir(mdir)
                ase.io.write(filename="POSCAR", images=atoms, format='vasp')
                self.prepare_vasp_inputs(mdir)

    def gn_hcp(self):
        alat0 = self.pot['ahcp']
        for i in range(-15, 15):
            delta = 0.01
            alat = alat0 + i * delta
            if i >= 0:
                mdir = "dir-p-{:03d}".format(i)
            else:
                mdir = "dir-n-{:03d}".format(abs(i))
            atoms = Hexagonal.HexagonalClosedPacked(
                latticeconstant={'a': alat, 'c': self.pot["chcp"]},
                size=(1, 1, 1), symbol=self.pot["element"], pbc=(1, 1, 1))
            self.mymkdir(mdir)
            ase.io.write(filename="POSCAR", images=atoms, format='vasp')
            self.prepare_vasp_inputs(mdir)

    def collect_data(self, tag='hcp'):
        rng = [-50, 50]
        data = np.zeros([rng[1] - rng[0], 3])
        cnt = 0
        for i in range(rng[0], rng[1]):
            if i >= 0:
                mdir = "dir-p-{:03d}".format(i)
            else:
                mdir = "dir-n-{:03d}".format(abs(i))
            os.chdir(mdir)
            (energy, vol, atoms) = self.vasp_energy_stress_vol_quick()
            if (tag in 'fcc') or (tag == 'bcc'):
                data[cnt, 0] = atoms.get_cell()[0, 1]
                data[cnt, 1] = (energy)
                data[cnt, 2] = i
            if tag in 'hcp':
                data[cnt, 0] = atoms.get_cell()[0, 0]
                data[cnt, 1] = (energy)
                data[cnt, 2] = i
            cnt += 1
            os.chdir(os.pardir)
        np.savetxt('lat.dat', data)


if __name__ == '__main__':
    usage = "usage:%prog [options] arg1 [options] arg2"
    parser = OptionParser(usage=usage)
    parser.add_option("-t", "--mtype", action="store",
                      type="string", dest="mtype")
    parser.add_option('-p', "--param", action="store",
                      type='string', dest="fargs")
    (options, args) = parser.parse_args()
    drv = cal_lattice()
    dispatcher = {'plt': drv.interpolate_lattice,
                  'loop': drv.loop_kpts,
                  'bcc': drv.gn_bcc,
                  'fcc': drv.gn_fcc,
                  'hcp': drv.gn_hcp,
                  'hcp2': drv.gn_hcp_mesh,
                  'clc': drv.collect_data}

    if options.fargs is not None:
        dispatcher[options.mtype.lower()](options.fargs)
    else:
        dispatcher[options.mtype.lower()]()

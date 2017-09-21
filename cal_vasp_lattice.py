#!/usr/bin/env python
# encoding: utf-8
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-07-05 08:12:30
# @Last Modified by:   chaomy
# @Last Modified time: 2017-09-20 14:40:22


import numpy as np
import glob
from scipy.interpolate import interp1d
import matplotlib.pylab as plt
from optparse import OptionParser
import os
import ase.io
import ase


try:
    import gn_config
    import get_data
    import gn_kpoints
    import gn_incar
    import gn_pbs
    import md_pot_data

except ImportError:
    print "error in importing"


class cal_lattice(gn_config.bcc,
                  gn_config.fcc,
                  gn_config.hcp,
                  gn_pbs.gn_pbs,
                  get_data.get_data,
                  gn_kpoints.gn_kpoints,
                  gn_incar.gn_incar):

    def __init__(self):
        self.figsize = (8, 6)
        self.npts = 20
        self.kpoints = [31, 31, 31]
        self.pot = md_pot_data.va_pot.Mg_pbe
        gn_kpoints.gn_kpoints.__init__(self)
        get_data.get_data.__init__(self)
        gn_incar.gn_incar.__init__(self)
        gn_pbs.gn_pbs.__init__(self)
        return

    def interpolate_lattice(self, filename):
        data = np.loadtxt(filename)
        lattice = np.abs(data[:, 0])
        energy = data[:, 1]
        InterPoints = np.linspace(lattice[0], lattice[-1], 101)
        f = interp1d(lattice, energy)
        Ynew = f(InterPoints)

        i = np.argmin(Ynew)
        print "min energy ", np.min(energy)
        print "num", np.argmin(energy), "lattice", InterPoints[i]
        print(np.min(Ynew))

        fig = plt.figure(figsize=self.figsize)
        ax = fig.add_subplot(111)
        ax.plot(lattice, energy)
        plt.savefig("lattice.png")
        plt.show()
        return

    def set_pbs(self, mdir):
        self.set_pbs_type('va')
        self.set_wall_time(10)
        self.set_job_title(mdir)
        self.set_nnodes(1)
        self.set_ppn(12)
        self.set_main_job("mpirun vasp")
        self.write_pbs(od=None)
        return

    def loop_kpts(self):
        files = glob.glob('./DATA*')
        for file in files:
            cal_lattice(file)
        return

    def prepare_vasp_inputs(self, mdir):
        self.set_incar_type('dftunrelax')
        self.write_incar()
        self.set_diff_kpoints([31, 31, 31])
        self.set_intype('gamma')
        self.write_kpoints()
        self.set_pbs(mdir)
        os.system("cp ../../POTCAR .")
        return

    def gn_bcc(self):
        alat0 = self.pot['latbcc']
        delta = 0.001
        rng = [-15, 15]
        bcc_drv = gn_config.bcc(self.pot)

        for i in range(rng[0], rng[1]):
            alat = alat0 + i * delta

            if i >= 0:
                mdir = "dir-p-{:03d}".format(i)
            else:
                mdir = "dir-n-{:03d}".format(abs(i))

            self.mymkdir(mdir)
            os.chdir(mdir)

            bcc_drv.set_lattce_constant(alat)
            atoms = bcc_drv.set_bcc_primitive((1, 1, 1))
            ase.io.write(filename="POSCAR", images=atoms, format='vasp')

            self.prepare_vasp_inputs(mdir)
            os.chdir(os.pardir)
        return

    def gn_fcc(self):
        alat0 = 3.90
        delta = 0.005
        #  rng = [-20, 20]
        rng = [20, 50]
        fcc_drv = gn_config.fcc(self.pot)

        for i in range(rng[0], rng[1]):
            alat = alat0 + i * delta
            if i >= 0:
                mdir = "dir-p-{:03d}".format(i)
            else:
                mdir = "dir-n-{:03d}".format(abs(i))

            self.mymkdir(mdir)
            os.chdir(mdir)

            fcc_drv.set_lattce_constant(alat)
            atoms = fcc_drv.set_fcc_primitive((1, 1, 1))
            ase.io.write(filename="POSCAR", images=atoms, format='vasp')

            self.prepare_vasp_inputs(mdir)
            os.chdir(os.pardir)
        return

    def gn_hcp(self):
        alat0 = self.pot['ahcp']
        hcp_drv = gn_config.hcp(self.pot)
        for i in range(-15, 15):
            delta = 0.01
            alat = alat0 + i * delta
            if i >= 0:
                mdir = "dir-p-{:03d}".format(i)
            else:
                mdir = "dir-n-{:03d}".format(abs(i))
            self.mymkdir(mdir)
            os.chdir(mdir)
            hcp_drv.write_hcp_poscar(alat)

            self.set_incar_type("isif4")
            self.write_incar()

            self.set_diff_kpoints([36, 36, 19])
            self.set_intype('gamma')
            self.write_kpoints()
            self.set_pbs(mdir)
            os.system("cp ../POTCAR .")
            os.chdir(os.pardir)
        return

    def collect_data(self, tag='hcp'):
        rng = [-15, 15]
        rng = [-8, 8]
        rng = [-20, 50]

        data = np.zeros([rng[1] - rng[0], 2])
        cnt = 0
        for i in range(rng[0], rng[1]):
            if i >= 0:
                mdir = "dir-p-{:03d}".format(i)
            else:
                mdir = "dir-n-{:03d}".format(abs(i))
            os.chdir(mdir)

            (energy, vol, atoms) = self.vasp_energy_stress_vol_quick()
            if (tag == 'fcc') or (tag == 'bcc'):
                data[cnt, 0] = atoms.get_cell()[0, 1]
                data[cnt, 1] = (energy)

            if tag == 'hcp':
                data[cnt, 0] = atoms.get_cell()[0, 0]
                data[cnt, 1] = (energy)

            cnt += 1
            os.chdir(os.pardir)
        np.savetxt('lat.dat', data)
        return


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
                  'clc': drv.collect_data}

    if options.fargs is not None:
        dispatcher[options.mtype.lower()](options.fargs)
    else:
        dispatcher[options.mtype.lower()]()

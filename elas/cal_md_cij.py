#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: yang37
# @Date:   2017-06-21 18:42:47
# @Last Modified by:   chaomy
# @Last Modified time: 2018-05-08 22:38:44


import glob
import os
import ase.io
import numpy as np
import md_pot_data
from optparse import OptionParser
from scipy.optimize import leastsq
from utils import cal_add_strain
import ase.lattice.cubic as Cubic
import gn_config
import get_data
import gn_incar
import gn_pbs
import output_data

evA3toGpa = 160.21766208


class cal_cij(get_data.get_data,
              gn_incar.gn_incar,
              gn_pbs.gn_pbs,
              output_data.output_data,
              cal_add_strain.cal_add_strain,
              gn_config.gnStructure):

    def __init__(self):
        self.pot = md_pot_data.md_pot.Nb_eam
        # self.pot = md_pot_data.va_pot.Nb_pbe
        get_data.get_data.__init__(self)
        # self.pot = self.load_data('../BASICS/pot.dat')
        output_data.output_data.__init__(self)
        cal_add_strain.cal_add_strain.__init__(self)
        gn_config.gnStructure.__init__(self, self.pot)
        self.unit_delta = 1e-2 * 0.05
        self.rng = [-10, 11]

    def gn_primitive_lmps(self, strain=np.mat(np.identity(3)), fnm="lmp_init.txt"):
        cell = np.mat([[-0.5, 0.5, 0.5], [0.5, -0.5, 0.5], [0.5, 0.5, -0.5]])
        cell = strain * cell
        cell = self.pot['lattice'] * self.lmp_change_box(cell)
        atoms = ase.Atoms(self.pot['element'], positions=[
                          [0, 0, 0]], cell=cell, pbc=[1, 1, 1])
        with open(fnm, mode="w") as fout:
            fout.write("#lmp data config\n")
            fout.write("1 atoms\n")
            fout.write("1 atom types\n")
            fout.write("%f \t %.10f xlo xhi\n" % (0, cell[0, 0]))
            fout.write("%f \t %.10f ylo yhi\n" % (0, cell[1, 1]))
            fout.write("%f \t %.10f zlo zhi\n" % (0, cell[2, 2]))
            fout.write("{:.10f} {:.10f} {:.10f} xy xz yz\n".format(
                cell[1, 0], cell[2, 0], cell[2, 1]))
            fout.write("Atoms\n\n")
            fout.write("1 1 0.0 0.0 0.0\n")

    def obtain_cij(self):
        vol = 0.5 * self.pot["lattice"] ** 3
        v = np.zeros(3)

        data = np.loadtxt("vol.txt")
        res = np.polyfit(data[:, 0], data[:, 1] / (vol), deg=2)
        print(res)
        v[0] = (res[0])

        data = np.loadtxt("otho.txt")
        res = np.polyfit(data[:, 0], data[:, 1] / (vol), deg=2)
        print(res)
        v[1] = (res[0])

        data = np.loadtxt("mono.txt")
        res = np.polyfit(data[:, 0], data[:, 1] / (vol), deg=2)
        print(res)
        v[2] = (res[0])

        print(v)
        c = np.linalg.inv(
            np.mat([[1.5, 3., 0], [1, -1, 0], [0, 0, 2]])) * np.mat(v).transpose()
        print(c * evA3toGpa)

    # use one atom to calculate elastic constants
    def loop_prepare_cij(self):
        for kk in ['vol', 'otho', 'mono']:
            cnt = 0
            for j in range(self.rng[0], self.rng[1]):
                delta = self.unit_delta * j
                self.gn_primitive_lmps(self.getstrain[kk](delta),
                                       'lmp_{}_{:03}.txt'.format(kk, cnt))
                cnt += 1

    def run(self):
        for kk in ['vol', 'otho', 'mono']:
            data = np.zeros([self.rng[1] - self.rng[0], 3])
            cnt = 0
            for j in range(self.rng[0], self.rng[1]):
                data[cnt, 0] = self.unit_delta * j
                os.system("cp lmp_{}_{:03}.txt lmp_init.txt".format(kk, cnt))
                os.system('lmp_mpi -i in.init')
                data[cnt, 1] = np.loadtxt('out.txt')
                data[cnt, 2] = np.linalg.det(self.getstrain[kk](data[cnt, 0]))
                cnt += 1
            np.savetxt("{}.txt".format(kk), data)

    def auto(self):
        self.loop_prepare_cij()
        self.run()
        self.obtain_cij()

if __name__ == "__main__":
    usage = "usage:%prog [options] arg1 [options] arg2"
    parser = OptionParser(usage=usage)
    parser.add_option('-t', "--mtype", action="store",
                      type="string", dest="mtype")
    parser.add_option('-p', "--param", action="store",
                      type='string', dest="fargs")
    (options, args) = parser.parse_args()

    drv = cal_cij()

    dispatcher = {'prep': drv.loop_prepare_cij,
                  'run': drv.run,
                  'cal': drv.obtain_cij,
                  'auto': drv.auto}

    if options.fargs is not None:
        dispatcher[options.mtype.lower()](options.fargs)
    else:
        dispatcher[options.mtype.lower()]()

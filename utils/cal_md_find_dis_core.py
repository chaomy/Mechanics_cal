#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2018-05-01 21:43:58
# @Last Modified by:   chaomy
# @Last Modified time: 2018-05-02 22:29:05


import ase
import numpy as np
import ase.calculators.neighborlist
import ase.utils.geometry
import ase.io
from scipy.optimize import minimize

dismap = np.mat([[1, -2.36923792e-09,  1.33112098e-10],
                 [3., 0.0021926, 0.00196574]])


class md_find_core(object):

    def cost_method_find_core(self):
        data = []
        init = [0.1, -0.01]
        initlist = np.loadtxt("dis_traj.txt.start")

        for i in range(1, 21):
            init = initlist[i - 1]
            atoms_rlx = ase.io.read(
                "dump.all.{}".format(i), format="lammps-dump")

            # find 3 atom rows
            cut_tm = 2.7**2
            cidx = []
            rs = []
            centers = np.loadtxt("dis_center.txt")

            atoms_perf = ase.io.read("perf_poscar", format='vasp')
            atoms_perf.wrap(pbc=[1, 1, 1])
            prfpos = atoms_perf.get_positions()

            for i in range(len(prfpos)):
                dx = prfpos[i, 0] - centers[0]
                dy = prfpos[i, 2] - centers[1]
                if (dx**2 + dy**2) < cut_tm:
                    rs.append(np.array([dx, dy]))
                    cidx.append(i)

            for i in range(len(prfpos)):
                dx = prfpos[i, 0] - centers[2]
                dy = prfpos[i, 2] - centers[3]
                if (dx**2 + dy**2) < cut_tm:
                    rs.append(np.array([dx, dy]))
                    cidx.append(i)

            rs = rs[1:]
            cidx = cidx[1:]

            print(len(cidx))
            self.pos_rlx = atoms_rlx.get_positions()[cidx]
            self.pos_prf = prfpos[cidx]
            self.cidx = cidx

            print(rs)
            print(self.pos_prf)
            print(self.pos_rlx)

            res = minimize(self.cost_func, np.array(
                init), method='Nelder-Mead', tol=5e-4)
            print(res.x)
            data.append(res.x)
            # init = res.x
        np.savetxt("dis_traj.txt", np.array(data), fmt="%.7f")

    def mesh_cost(self):
        min_err = 1.0
        for dx in np.linspace(0.3, 0.7, 10):
            for dy in np.linspace(-1. / 3., 0.0, 40):
                err = self.cost_func([dx, dy])
                if err < min_err:
                    min_err = err
                    opt = [dx, dy]
        print(opt)

    def cost_func_ddmap(self, d):
        # pos_sln = self.dis_ref(d).get_positions()[self.cidx]
        drlx = self.pos_rlx[:, -1] - self.pos_prf[:, -1]
        # dz = pos_sln[:, -1] - self.pos_prf[:, -1]
        print(np.sum(drlx))

    def cost_func(self, d):
        pos_sln = self.dis_ref(d).get_positions()[self.cidx]
        e1 = np.linalg.norm(self.pos_rlx[:, 1] - pos_sln[:, 1])
        return e1

    def dis_ref(self, d):
        sx, sy, ix = 10.0, 5, 10.5
        unitx = np.sqrt(6) / 3. * self.pot['lattice']
        unity = np.sqrt(2) / 2. * self.pot['lattice']

        ci1 = [(sx + d[0]) * unitx, (sy + 1. / 3. + d[1]) * unity]
        atoms = self.set_dipole_box(1)
        atoms = self.bcc_screw_dipole_alongz_with_image(atoms, ci1)
        return self.convert_alongz_to_alongy(atoms)
        # ci2 = [(sx + ix) * unitx, (sy + 2. / 3.) * unity]

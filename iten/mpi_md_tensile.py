#!/usr/bin/env python
# encoding: utf-8
# -*- coding: utf-8 -*-
# @Author: yang37
# @Date:   2017-06-12 17:03:43
# @Last Modified by:   chaomy
# @Last Modified time: 2018-02-18 14:04:19


import copy
import os
import numpy as np
import ase
import ase.lattice
from multiprocessing import Pool
import shutil
import sys
import get_data
import output_data
import gn_lmp_infile
import md_pot_data
import gn_config

__author__ = 'Chaoming Yang'
__all__ = ['md_tensile']


def unwrap_self_run_lammps(arg, **kwarg):
    return md_loop_tensile.lammps_job(*arg, **kwarg)


class md_tensile(get_data.get_data,
                 output_data.output_data,
                 gn_lmp_infile.gn_md_infile,
                 gn_config.bcc):

    def __init__(self):
        self.pot = md_pot_data.md_pot.Nb_meamc
        get_data.get_data.__init__(self)
        output_data.output_data.__init__(self)
        gn_lmp_infile.gn_md_infile.__init__(self, self.pot)
        gn_config.bcc.__init__(self, self.pot)

        self.size = [1, 1, 1],
        self.orientation = np.mat([[1, 0, 0],
                                   [0, 1, 0],
                                   [0, 0, 1]])

        # self.M = self.get_M(self.orientation)
        self.M = self.orientation
        self.Sij = self.get_Sij()
        self.lx, self.ly, self.lz = 0, 0, 0
        self.addedstrain = np.mat([[0, 0, 0],
                                   [0, 0, 0],
                                   [0, 0, 0]], "float")
        self.stress_original = np.array([0, 0, 0, 0, 0, 0],
                                        dtype="float")

    def update_strain(self, Inputstrain,
                      stress, Correct_strain):
        addedstrain = copy.deepcopy(self.addedstrain)

        if abs(stress[1]) > self._stress_ThrValue:
            addedstrain[1, 1] = Inputstrain[1]
            Correct_strain[1, 1] += addedstrain[1, 1]
        if abs(stress[2]) > self._stress_ThrValue:
            addedstrain[2, 2] = Inputstrain[2]
            Correct_strain[2, 2] += addedstrain[2, 2]
        return Correct_strain

    def gn_bcc_tpath(self, delta, Correct_strain):
        # very important (vasp add strain is basis right time strain)
        bas = np.mat(np.identity(3))
        va_bas = bas * self.strainmtx
        cell = self.pot["lattice"] * va_bas
        atoms = ase.lattice.cubic.BodyCenteredCubic(directions=[[1., 0., 0.],
                                                                [0., 1., 0.],
                                                                [0., 0., 1.]],
                                                    latticeconstant=self.pot[
                                                        "lattice"],
                                                    size=(1, 1, 1),
                                                    symbol=self.pot['element'],
                                                    pbc=(1, 1, 1))
        pos = np.mat([[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]) * cell
        atoms.set_cell(cell)
        atoms.set_positions(pos)
        self.write_lmp_config_data(atoms, 'init.txt')
        return

    def gn_bcc_opath(self, delta, Correct_strain):
        alat = self.pot["lattice"]
        bas = np.mat(np.mat([[1., 0., 0.],
                             [0, np.sqrt(2), 0.],
                             [0., 0., np.sqrt(2)]]))
        va_bas = bas * self.strainmtx
        cell = alat * va_bas
        atoms = ase.lattice.cubic.FaceCenteredCubic(directions=[[1., 0., 0.],
                                                                [0., 1., 0.],
                                                                [0., 0., 1.]],
                                                    latticeconstant=alat,
                                                    size=(1, 1, 1),
                                                    symbol=self.pot['element'],
                                                    pbc=(1, 1, 1))
        pos = np.mat([[0.0, 0.0, 0.0], [0.5, 0.5, 0.0],
                      [0.5, 0.0, 0.5], [0.0, 0.5, 0.5]]) * cell
        atoms.set_cell(cell)
        atoms.set_positions(pos)
        self.write_lmp_config_data(atoms, 'init.txt')
        return


class md_loop_tensile(md_tensile):

    def __init__(self):
        self.pot = md_pot_data.md_pot.Nb_meamc
        self._flux_exe = 'lmp_mpi -i in.stat_tensile -screen no'
        self._looptime = 10
        self._increment = 0.02
        self._stress_ThrValue = 0.05
        self.root_dir = os.getcwd()

    def get_strain(self,
                   Sij,
                   stress,
                   stressPrime):
        tag = "constant"
        strain = np.zeros(6).transpose()
        if tag == "constant":
            coeff = 1.0
            strain[1] = Sij[0, 1] * stress[1] * coeff
            strain[2] = Sij[0, 2] * stress[2] * coeff
        else:
            stressPrime = float(stressPrime)
            coeff = 0.06 + 0.06 * stressPrime ** 0.20 + \
                0.03 * stressPrime ** 0.19 + \
                0.01 * stressPrime ** 0.1
            if stressPrime > 20:
                coeff = 0.08 + 0.01 * stressPrime ** 0.19 + \
                    0.09 * stressPrime ** 0.1
            strain = -Sij * stress * coeff
        return strain

    def lammps_job(self, option, delta):
        Correct_strain = np.mat([[0, 0, 0],
                                 [0, 0, 0],
                                 [0, 0, 0]], "float")
        md_tensile.__init__(self)
        count = 0
        coeff = 1.0
        self.strainmtx = np.matrix([[1 + delta, 0, 0],
                                    [0, 1 + 0.2 * delta, 0],
                                    [0, 0, 1 - 0.2 * delta]],
                                   "float") + Correct_strain
        strain = np.zeros(6).transpose()
        while True:
            if option in ['TP', 'tpath']:
                self.gn_bcc_tpath(delta, Correct_strain)
            elif option in ['OP', 'tpath']:
                self.gn_bcc_opath(delta, Correct_strain)

            os.system("%s" % (self._flux_exe))

            dat = np.loadtxt("len.txt")
            self.lx, self.ly, self.lz, self.eng = dat[
                0], dat[1], dat[2], dat[3]
            self.stress_original = np.loadtxt("stress.txt")

            stress = copy.deepcopy(self.stress_original)
            stress[0] = 0.0
            stress_abs = np.abs(stress)
            stressPrime = np.max(stress_abs)

            if stressPrime < self._stress_ThrValue:
                break
            elif count > 500:
                break
            else:
                if abs(stress[1]) > self._stress_ThrValue:
                    strain[1] = self.Sij[0, 1] * stress[1] * coeff
                    self.strainmtx[1, 1] += strain[1]

                if abs(stress[2]) > self._stress_ThrValue:
                    strain[2] = self.Sij[0, 2] * stress[2] * coeff
                    self.strainmtx[2, 2] += strain[2]

                with open('monitor.txt', 'a') as fid:
                    print >> fid, "delta ", delta
                    print >> fid, "Run times", count
                    print >> fid, "stress_original", self.stress_original
                    print >> fid, "stressPrime ", stressPrime
                    fid.close()
                count += 1

    def cal_md_tensile(self, options):
        # pool = Pool(processes=self._looptime)
        # List = np.arange(self._looptime)
        # results = pool.map(unwrap_self_run_lammps,
        #                    zip([self] * len(List),
        #                        [options] * len(List),
        #                        List))
        npts = 13
        dat = np.ndarray([npts, 10])
        for i in range(npts):
            delta = 0.02 * i
            self.lammps_job(options, delta)
            dat[i, 0], dat[i, 1] = delta, self.eng
            dat[i, 2], dat[i, 3] = self.ly, self.lz
            dat[i, 4:] = self.stress_original
        np.savetxt("iten.txt", dat, fmt='%6.5f')

if __name__ == "__main__":
    M = md_loop_tensile()
    M.cal_md_tensile('OP')
    # (strain2, stress2) = M.cal_md_tensile('OP')
    # M.output_log(strain, stress, stress2)

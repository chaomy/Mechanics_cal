#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-06-25 14:28:58
# @Last Modified by:   chaomy
# @Last Modified time: 2017-08-30 14:57:49

from numpy import cos, sin
from collections import OrderedDict
import ase
import ase.io
import numpy as np
import md_pot_data
import tool_elastic_constants
import stroh_solve
import ase.lattice
import cal_md_dislocation
import cal_md_crack_ini

matconsts = OrderedDict([('Al1', {'lat': 4.05,
                                  'ugsf': 0.167,
                                  'c11': 113.4,
                                  'c12': 61.5,
                                  'c44': 31.6}),
                         ('Al2', {'lat': 4.03,
                                  'ugsf': 0.119,
                                  'c11': 118.0,
                                  'c12': 62.2,
                                  'c44': 36.7}),
                         ('Gold', {'lat': 4.08,
                                   'ugsf': 0.097,
                                   'c11': 183.2,
                                   'c12': 158.7,
                                   'c44': 45.3}),
                         ('Silver', {'lat': 4.09,
                                     'ugsf': 0.119,
                                     'c11': 129.1,
                                     'c12': 91.7,
                                     'c44': 56.7})])


class cal_dis_dipole(object):

    def __init__(self, pot=None):
        if pot is None:
            pot = md_pot_data.qe_pot.vca_W75Re25
        self.pot = pot
        self.mddis_drv = cal_md_dislocation.md_dislocation(self.pot)
        return

    def set_dipole_box(self, sizen=1):
        n = 7 * sizen
        m = 11 * sizen
        t = 1 * sizen
        print self.pot
        alat = self.pot['lattice']
        atoms = ase.lattice.cubic.BodyCenteredCubic(
            directions=[[1., 1., -2.],
                        [-1., 1., 0],
                        [0.5, 0.5, 0.5]],
            latticeconstant=alat,
            size=(n, m, t),
            symbol=self.pot['element'])
        atoms = self.mddis_drv.cut_half_atoms_new(atoms, "cuty")
        supercell = atoms.get_cell()
        strain = np.mat([[1.0, 0.0, 0.0],
                         [0.5, 1.0, 0.5],
                         [0.0, 0.0, 1.0]])
        supercell = strain * supercell
        atoms.set_cell(supercell)
        atoms.wrap(pbc=[1, 1, 1])
        return atoms

    def set_dipole_triangular_box(self, sizen=1):
        u = 1. / 3. * np.array([1., -2., 1.])
        v = 1. / 3. * np.array([2., -1., -1.])
        z = 1. / 3. * np.array([1., 1., 1.])

        # n1u = 3 * n - 1
        # n1v = 0
        # n2u = 0
        # n2v = 3 * n - 1
        # c1z = 1. / 3.
        # c2z = -c1z

        # c1 = n1u * u + n1v * v + c1z * z
        # c2 = n2u * u + n2v * v + c2z * z

        # print c1; print c2
        atoms = ase.lattice.cubic.BodyCenteredCubic(
            directions=[[1., -2., 1.],
                        [2., -1., -1.],
                        [1., 1., 1.]],
            latticeconstant=self.pot['latbcc'],
            size=(4, 4, 1),
            symbol=self.pot['element'])
        supercell = atoms.get_cell()
        addstrain = False
        if addstrain is True:
            print len(atoms)
            strain = np.mat([[1.0, 0.0, 0.0],
                             [0.5, 1.0, 0.5],
                             [0.0, 0.0, 1.0]])
            supercell = strain * supercell
            atoms.set_cell(supercell)
            atoms.wrap(pbc=[1, 1, 1])
        ase.io.write('tri_perf_poscar.vasp', images=atoms, format='vasp')
        return atoms

    def loop_table(self):
        for key in matconsts.keys():
            print key
            self.get_cutin_result(matconsts[key])
        return

    def cal_crack(self):
        drv = cal_md_crack_ini.md_crack_ini()
        drv.cal_crack_anglecoeff()
        return

    def get_cutin_result(self, param):
        c = tool_elastic_constants.elastic_constants(
            C11=param['c11'],
            C12=param['c12'],
            C44=param['c44'])
        # A
        axes = np.array([[-1, -1, 2],
                         [1, 1, 1],
                         [-1, 1, 0]])
        burgers = param['lat'] * np.sqrt(2.) / 2. * np.array([-1., 1., 0])
        stroh = stroh_solve.Stroh(c, burgers, axes=axes)
        A = np.mat(np.zeros([3, 3]), dtype='complex')
        A[:, 0] = np.mat(stroh.A[0]).transpose()
        A[:, 1] = np.mat(stroh.A[2]).transpose()
        A[:, 2] = np.mat(stroh.A[4]).transpose()

        B = np.mat(np.zeros([3, 3]), dtype='complex')
        B[:, 0] = np.mat(stroh.L[0]).transpose()
        B[:, 1] = np.mat(stroh.L[2]).transpose()
        B[:, 2] = np.mat(stroh.L[4]).transpose()

        Gamma = 0.5 * np.real(np.complex(0, 1) * A * np.linalg.inv(B))
        theta = np.deg2rad(70.0)
        omega = np.mat([[cos(theta), sin(theta), 0.0],
                        [-sin(theta), cos(theta), 0.0],
                        [0.0, 0.0, 1.0]])
        phi = 0
        svect = np.mat(np.array([cos(phi), 0.0, sin(phi)]))
        usf = param['ugsf']  # J/m^2
        Gamma = np.abs(np.linalg.inv(Gamma))
        # Gamma = omega * Gamma * omega

        Gamma = (svect * Gamma * svect.transpose())[0, 0]  # in GPa
        print Gamma

        Gamma = Gamma * 1e9  # Pa
        ke1 = np.sqrt(Gamma * usf)
        print ke1 * 1e-6
        return

        # A = np.mat(np.zeros([3, 3]), dtype='complex')
        # A[:, 0] = np.mat(stroh.A[1]).transpose()
        # A[:, 1] = np.mat(stroh.A[3]).transpose()
        # A[:, 2] = np.mat(stroh.A[5]).transpose()

        # B = np.mat(np.zeros([3, 3]), dtype='complex')
        # B[:, 0] = np.mat(stroh.L[1]).transpose()
        # B[:, 1] = np.mat(stroh.L[3]).transpose()
        # B[:, 2] = np.mat(stroh.L[5]).transpose()
        # Gamma = np.real(np.complex(0, 1) * A * np.linalg.inv(B))

    def print_dis_constants(self):
        c = tool_elastic_constants.elastic_constants(
            C11=self.pot['c11'],
            C12=self.pot['c12'],
            C44=self.pot['c44'])
        # A
        axes = np.array([[1, -1, 1],
                         [2, 1, -1],
                         [0, 1, 1]])

        burgers = self.pot['lattice'] / 2 * np.array([1., 1., 1.])
        stroh = stroh_solve.Stroh(c, burgers, axes=axes)
        #
        print stroh.A[0]
        # print stroh.L
        # print "K tensor", stroh.K_tensor
        # print "K (biKijbj)", stroh.K_coeff, "eV/A"
        # print "pre-ln alpha = biKijbj/4pi", stroh.preln, "ev/A"
        return

    def bcc_screw_dipole_triangular_atoms(self, atoms=None, fname='qe.in'):
        c = tool_elastic_constants.elastic_constants(
            C11=self.pot['c11'],
            C12=self.pot['c12'],
            C44=self.pot['c44'])
        axes = np.array([[1, 1, -2],
                         [-1, 1, 0],
                         [1, 1, 1]])
        burgers = self.pot['lattice'] / 2 * np.array([1., 1., 1.])
        stroh = stroh_solve.Stroh(c, burgers, axes=axes)

        atoms = self.set_dipole_triangular_box()
        pos = atoms.get_positions()
        cell = atoms.get_cell()

        c1 = [0.51 * cell[0, 0], 1. / 3. * cell[1, 1]]
        c2 = [2 * cell[1, 0], 2. / 3. * cell[1, 1]]
        print cell, c1, c2

        shiftc1 = np.ones(np.shape(pos)) * np.array([c1[0], c1[1], 0.0])
        shiftc2 = np.ones(np.shape(pos)) * np.array([c2[0], c2[1], 0.0])

        disp1 = stroh.displacement(pos - shiftc1)
        disp2 = stroh.displacement(pos - shiftc2)

        atoms.set_positions(pos + np.real(disp1) - np.real(disp2))
        ase.io.write('tri_dis_poscar.vasp', atoms, format='vasp')
        return

    def bcc_screw_dipole_configs_alongz(self, sizen=1):
        c = tool_elastic_constants.elastic_constants(
            C11=self.pot['c11'],
            C12=self.pot['c12'],
            C44=self.pot['c44'])

        axes = np.array([[1, 1, -2],
                         [-1, 1, 0],
                         [1, 1, 1]])

        burgers = self.pot['lattice'] / 2 * np.array([1., 1., 1.])
        stroh = stroh_solve.Stroh(c, burgers, axes=axes)

        atoms = self.set_dipole_box()
        atoms_perf = atoms.copy()
        pos = atoms.get_positions()

        unitx = np.sqrt(6) / 3. * self.pot['lattice']
        unity = np.sqrt(2) / 2. * self.pot['lattice']
        sx = 10.0 * sizen
        sy = 5 * sizen
        ix = 10.5 * sizen

        # c1 = 1. / 3. * np.sum(self.pot['core1'], axis=0)
        # c2 = 1. / 3. * np.sum(self.pot['core2'], axis=0)
        # shiftc1 = \
        # np.ones(np.shape(pos)) * np.array([c1[0, 0], c1[0, 1], 0.0])
        # shiftc2 = \
        # np.ones(np.shape(pos)) * np.array([c2[0, 0], c2[0, 1], 0.0])

        opt = 'split'
        if opt == 'split':
            c1 = self.pot['posleft'] + \
                np.array([0.0, 0.21 * self.pot['yunit']])
            c2 = self.pot['posrigh'] + \
                np.array([0.0, -0.21 * self.pot['yunit']])
        else:
            c1 = [(sx) * unitx, (sy + 1. / 3.) * unity]
            c2 = [(sx + ix) * unitx, (sy + 2. / 3.) * unity]
        shiftc1 = np.ones(np.shape(pos)) * np.array([c1[0], c1[1], 0.0])
        shiftc2 = np.ones(np.shape(pos)) * np.array([c2[0], c2[1], 0.0])

        disp1 = stroh.displacement(pos - shiftc1)
        disp2 = stroh.displacement(pos - shiftc2)

        atoms.set_positions(pos + np.real(disp1) - np.real(disp2))
        # ase.io.write('POSCAR', atoms, format='vasp')
        return (atoms, atoms_perf)


if __name__ == '__main__':
    drv = cal_dis_dipole()
    # drv.bcc_screw_dipole_configs_alongz()
    # drv.bcc_screw_dipole_triangular_atoms()
    # drv.print_dis_constants()
    # drv.get_cutin_result()
    # drv.loop_table()
    drv.cal_crack()

#!/usr/bin/env python
# encoding: utf-8
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-07-05 08:12:30
# @Last Modified by:   chaomy
# @Last Modified time: 2018-05-01 22:52:15

import numpy as np
import atomman as am
import atomman.unitconvert as uc
import matplotlib.pyplot as plt
import stroh_solve
import tool_elastic_constants


class cubic_cij:
    c11 = None
    c12 = None
    c44 = None


class cal_intro_ani_dis(object):

    def cal_nye(self):
        torient = 'z'
        if torient == 'y':
            e1 = 1. / 3. * np.array([1., 1., -2.])
            e2 = np.array([0.5, 0.5, 0.5])
            e3 = 1. / 2. * np.array([1, -1, 0])

        e1 = 1. / 3. * np.array([1., 1., -2.])
        e2 = 1. / 2. * np.array([-1., 1., 0])
        e3 = np.array([0.5, 0.5, 0.5])

        # unit cell
        unit_cell = np.array([e1, e2, e3])

        sizen = 1
        n = 7 * sizen
        m = 11 * sizen
        t = 3

        atoms = self.set_bcc_convention(directions=[[e1[0], e1[1], e1[2]],
                                                    [e2[0], e2[1], e2[2]],
                                                    [e3[0], e3[1], e3[2]]],
                                        size=(n, m, t))

        # perfect positions
        p = self.pot['lattice'] * self.bccneigh

        # atoms = self.intro_dipole_screw_atoms(atoms, self.pot['lattice'])
        atoms = self.intro_single_screw_atoms(atoms)

        ase.io.write("lmp_dis.cfg", atoms, "cfg")

        # system
        system, elements = am.convert.ase_Atoms.load(atoms)

        r = (3**0.5 / 2. + 1) / 2.
        system.nlist(self.pot['lattice'] * r)
        nye_rlt = am.defect.nye_tensor(system, p, axes=unit_cell)

        print(np.max(nye_rlt['strain'][:]))
        # for key, value in nye_rlt.iteritems():
        #     print key, value
        #     system.atoms_prop(key=key, value=value)

        # int_sum = np.empty(3)
        # for i in range(3):
        #     int_sum[i] = am.plot.interpolate_contour(system, 'Nye_tensor',
        #                                              index=[1, i],
        #                                              cmap='bwr')[0]
        # print('burgers vector estimate = ', int_sum)

        ############################################################
        # since lammmps can only use  xy xz yz
        # new x is old y ; new y is old z ; new z is old x
        # x  1  1 -2
        # y  1  1  1
        # z  1 -1  0
        ############################################################

    def print_dis_constants(self):
        struct = "hex"
        if struct in ["cubic"]:
            # Cubic
            c = tool_elastic_constants.elastic_constants(
                C11=self.pot['c11'],
                C12=self.pot['c12'],
                C44=self.pot['c44'])
            axes = np.array([[1, -1, 1],
                             [2, 1, -1],
                             [0, 1, 1]])
            burgers = self.pot['lattice'] / 2 * np.array([1., 1., 1.])
            stroh = stroh_solve.Stroh(c, burgers, axes=axes)
            print(stroh.A[0])

        # hexagonal
        if struct in ["hex"]:
            print(self.pot["lattice"])
            axes = np.array([[1, 0, 0],
                             [0, 1, 0],
                             [0, 0, 1]])

        burgers = self.pot['lattice'] / 2 * np.array([1., 1, 0])
        c = am.ElasticConstants()
        c.hexagonal(C11=326.08, C33=357.50, C12=129.56, C13=119.48, C44=92.54)
        stroh = stroh_solve.Stroh(c, burgers, axes=axes)
        print(stroh.A)
        print(stroh.L)

        # print(c)
        # print stroh.L
        # print "K tensor", stroh.K_tensor
        # print "K (biKijbj)", stroh.K_coeff, "eV/A"
        # print "pre-ln alpha = biKijbj/4pi", stroh.preln, "ev/A"

    def calQ(self):
        c11, c12, c44 = self.cij.c11, self.cij.c12, self.cij.c44
        e = c44 / (c12 + c44)
        f = (c11 - c44) / (c12 + c44)
        n1sqr = n1 * n1
        n2sqr = n2 * n2
        n3sqr = n3 * n3
        fsqr = f * f

        D = e * e * (e + f) + e * (fsqr - 1) * (n1sqr + n2sqr + n3sqr) \
            + (f - 1) * (f - 1) * (f + 2) * n1sqr * n2sqr * n3sqr
        Dp = (c12 + c44) * D

        nn11 = e * (e + f) - e * f * n1 * n1 + \
            (f * f - 1) * (n2 * n3) * (n2 * n3)
        nn11 /= Dp

        nn12 = -n1 * n2 * ((f - 1) * n3sqr + e)
        nn12 /= Dp
        Q11 = self.integrate(nn11)

    def fcc_edge(self):
        axes = np.array([[1, 0, -1],
                         [1, 1, 1],
                         [1, -2, 1]])

        alat = uc.set_in_units(4.0248, 'angstrom')

        C11 = uc.set_in_units(113.76, 'GPa')
        C12 = uc.set_in_units(61.71, 'GPa')
        C44 = uc.set_in_units(31.25, 'GPa')

        c = am.ElasticConstants(C11=C11, C12=C12, C44=C44)
        burgers = alat / 2 * np.array([1., 0., -1.])

        # initializing a new Stroh object using the data
        stroh = am.defect.Stroh(c, burgers, axes=axes)

        pos_test = uc.set_in_units(np.array([12.4, 13.5, -10.6]), 'angstrom')

        disp = stroh.displacement(pos_test)

        print("displacement =", uc.get_in_units(disp, 'angstrom'), 'angstrom')

        # monopole system
        box = am.Box(a=alat, b=alat, c=alat)
        atoms = am.Atoms(natoms=4, prop={'atype': 1, 'pos': [[0.0, 0.0, 0.0],
                                                             [0.5, 0.5, 0.0],
                                                             [0.0, 0.5, 0.5],
                                                             [0.5, 0.0, 0.5]]})
        ucell = am.System(atoms=atoms, box=box, scale=True)
        system = am.rotate_cubic(ucell, axes)

        shift = np.array(
            [0.12500000000000, 0.50000000000000, 0.00000000000000])
        new_pos = system.atoms_prop(key='pos', scale=True) + shift
        system.atoms_prop(key='pos', value=new_pos, scale=True)

        system.supersize((-7, 7), (-6, 6), (0, 1))
        disp = stroh.displacement(system.atoms_prop(key='pos'))

        system.atoms_prop(key='pos',
                          value=system.atoms_prop(key='pos') + disp)

        system.pbc = (False, False, True)
        system.wrap()

        pos = system.atoms_prop(key='pos')
        x = uc.get_in_units(pos[:, 0], 'angstrom')
        y = uc.get_in_units(pos[:, 1], 'angstrom')

        plt.figure(figsize=(8, 8))
        plt.scatter(x, y, s=30)
        plt.xlim(min(x), max(x))
        plt.ylim(min(y), max(y))
        plt.xlabel('x-position (Angstrom)', fontsize='large')
        plt.ylabel('y-position (Angstrom)', fontsize='large')
        plt.show()

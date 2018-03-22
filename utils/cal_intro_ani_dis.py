#!/usr/bin/env python
# encoding: utf-8
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-07-05 08:12:30
# @Last Modified by:   chaomy
# @Last Modified time: 2017-09-20 23:25:21

try:
    import numpy as np
    import atomman as am
    import atomman.unitconvert as uc
    import matplotlib.pyplot as plt

except ImportError:
    print("error during import")


class cubic_cij:
    c11 = None
    c12 = None
    c44 = None


class cal_intro_ani_dis(object):
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
        return

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
        return

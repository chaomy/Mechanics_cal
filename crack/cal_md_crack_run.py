#!/usr/bin/env python
# encoding: utf-8
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-07-05 08:10:22
# @Last Modified by:   chaomy
# @Last Modified time: 2018-03-27 17:16:08


import os


class md_crack_run(object):

    def cal_crack_anglecoeff(self):
        self.set_plane_strain_bij()
        self.get_scalarB()
        self.get_coeffs()
        # s1 = self.ckcoeff.u1
        # s2 = self.ckcoeff.u2

    def static_crack(self):
        self.set_plane_strain_bij()
        self.get_scalarB(self.pot["surf110"])
        self.get_coeffs()
        atoms = self.intro_crack_k1(atoms=self.gn_perf_plate())
        self.write_lmp_config_data(atoms, 'crack.txt')

    def static_crack_continue(self):
        execuable = "mpirun lmp_linux -in"
        p1, p2, q1, q2, u1, u2 = self.get_coeffs()
        Krestart = 2.4589812
        for i in range(1, 120):
            K = Krestart + 0.007 * i
            self.intro_crack_k1(K, p1, p2, q1, q2, u1, u2)
            os.system("rm crackxyz/*")
            os.system("%s in.read" % (execuable))
            self.rename_cfg(K)

    def Test_init(self):
        Kg = self.get_scalarB()
        self.get_coeffs()
        for i in range(50):
            K = Kg + 0.01 * i
            self.intro_crack_k1(self.crackcoeff)
            os.system("cp ./Crack.txt Init/Crack_%g" % (0.01 * K))

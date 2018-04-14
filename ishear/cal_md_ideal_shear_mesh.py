# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2018-03-04 21:16:49
# @Last Modified by:   chaomy
# @Last Modified time: 2018-03-07 12:52:13

import os
import numpy as np


class cal_bcc_ideal_shear_mesh(object):

    def mesh(self):
        basis = self.basis
        cnt = 0
        for i in list(range(-4, -1)) + list(range(1, 4)):
            for j in list(range(-4, -1)) + list(range(1, 4)):
                for k in list(range(-4, -1)) + list(range(1, 4)):
                    for d in range(1, 6):
                        strain = np.mat([[1 + 0.02 * i, 0.0, 0.0],
                                         [-0.02 * d, 1 + 0.02 * j, 0.0],
                                         [0.0, 0.0, 1 + 0.02 * k]])

                        cnt += 1
                        mdir = "dir_{:04d}".format(cnt)
                        new_strain = basis.transpose() * strain * basis
                        self.gn_primitive_lmps(new_strain, 'vasp')

                        self.mymkdir(mdir)
                        os.system("cp va.pbs {}".format(mdir))
                        os.system("cp POTCAR {}".format(mdir))
                        os.system("cp POSCAR {}".format(mdir))
                        os.system("cp KPOINTS {}".format(mdir))
                        os.system("cp INCAR {}".format(mdir))

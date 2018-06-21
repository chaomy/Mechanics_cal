# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2018-06-11 02:16:45
# @Last Modified by:   chaomy
# @Last Modified time: 2018-06-11 12:45:20

# add or delete some atoms

import ase
import numpy as np


class md_gb_del(object):

    def analysize_atomic_strain(self):
        atoms = ase.io.read("dump.final", format='lammps-dump')
        print(len(atoms))
        stss = np.loadtxt("dump.stss")

        mxval = np.max([np.max(stss[:, -6]), -np.min(stss[:, -6]),
                        np.max(stss[:, -5]), -np.min(stss[:, -5]),
                        np.max(stss[:, -4]), -np.min(stss[:, -4])])
        mxid = np.argmax([np.max(stss[:, -6]), -np.min(stss[:, -6]),
                          np.max(stss[:, -5]), -np.min(stss[:, -5]),
                          np.max(stss[:, -4]), -np.min(stss[:, -4])])

        print(mxval, mxid)  
        a = np.where(stss[:, -4] < -mxval + 400)
        print(a[0])

        del atoms[a[0]]
        self.write_lmp_config_data(atoms, "lmp_init.txt")

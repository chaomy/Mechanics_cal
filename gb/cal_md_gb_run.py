#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-07-05 08:12:30
# @Last Modified by:   chaomy
# @Last Modified time: 2019-05-13 11:56:48

import os
from glob import glob
import ase.io
import numpy as np
from ase import Atom


class md_gb_run(object):
    
    def loop_gb_run(self):  
        dlist = glob("1210_*")
        for mdir in dlist:
            os.system("cp in.min {}".format(mdir))
            os.chdir(mdir)
            self.mymkdir("dump")
            os.system("mpirun -n 1 lmp_mpi -i in.min")
            os.chdir(os.pardir)

    def loop_set_usp_run(self):
        fls = glob("coo.*")
        for i in range(len(fls)):
            mdir = "dir_{:03d}".format(i + 1)
            self.mymkdir(mdir)
            os.system("cp coo.{} {}/coo".format(i + 1, mdir))
            os.system("cp in.npt {}".format(mdir))
            os.system("mkdir {}/dump".format(mdir))

    def loop_grand(self):
        for o in ['n', 'p']:
            for i in range(0, 6):
                mdir = 'dir_{}_{:02d}'.format(o, i)
                self.mymkdir(mdir)
                self.make_gb_grand_canonical(o, i)
                os.system("cp in.npt va.pbs {}".format(mdir))
                os.system(
                    "cp init.{}.{:02d}.txt {}/init.txt".format(o, i, mdir))
                self.mymkdir("{}/dump".format(mdir))

    def make_gb_grand_canonical(self, opt='p', n=4):
        atoms = ase.io.read("dump/dump.00000", format='lammps-dump')
        cell = atoms.get_cell()
        usp = 8.0
        idx = []
        miny = cell[1, 1]
        maxy = 0.0

        for atom in atoms:
            if (atom.position[1] >= 0.5 * cell[1, 1] - 0.5 * usp) and (atom.position[1] <= 0.5 * cell[1, 1] + 0.5 * usp):
                atom.symbol = 'Re'
                idx.append(atom.index)

        if (opt == 'n'):
            del atoms[[idx[t] for t in np.random.randint(0, len(idx), n)]]
        elif (opt == 'p'):
            [atoms.append(Atom('Re', [np.random.rand() * cell[0, 0], np.random.rand() * 8.0 +
                                      0.5 * cell[1, 1] - 4.0, np.random.rand() * cell[2, 2]])) for i in range(n)]
        self.write_lmp_config_data(atoms, "init.{}.{:02d}.txt".format(opt, n))

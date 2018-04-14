# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-12-03 11:07:29
# @Last Modified by:   chaomy
# @Last Modified time: 2018-02-21 00:29:35

import os
from numpy import sqrt, arccos, arcsin, rad2deg
from numpy import ndarray, loadtxt
from numpy import mat, unique


class md_gb_loop(object):

    def find_angles(self):
        ags = []
        res = []

        # npts = 7

        ux = self.pot["ahcp"]
        uy = sqrt(3) / 2. * self.pot["ahcp"]

        cnt = 0

        jl = [1, 2]
        il = [[1, 2, 3, 4, 5, 6, 7], [1, 2, 3, 4, 5, 6, 7, 8]]

        for j in jl:
            for i in il[j - 1]:
                if j % 2 == 1:
                    idx = i + 0.5
                else:
                    idx = i
                idy = j
                lx = idx * ux
                ly = idy * uy

                ll = sqrt(lx**2 + ly**2)
                ang = rad2deg(arcsin(ly / ll))

                if (ang <= 30):
                    cnt += 1
                    ags.append(ang)
                    res.append([ang, ll, i, j])

        print(len(ags))
        (a, b) = unique(ags, return_index=True)

        print(len(a), len(b))
        self.ag = []

        for ee in b:
            self.ag.append(res[ee])

        print(self.ag)

    def loop_angle(self, opt='prep'):

        self.find_angles()

        npts = len(self.ag)

        j = 0
        self.gbdat = ndarray([2, npts + 1])
        self.gbdat[0, 0], self.gbdat[0, 1] = 0, 0

        for ee in self.ag:
            mdir = "dir_{:02.03f}_{:02d}_{:02d}".format(ee[0], ee[2], ee[3])

            if opt in ['prep']:
                self.mymkdir(mdir)
                self.build_hcp_gb_lmp(ee[0], ee[1], (30, 40), (0.0, 0.0, 0.0))
                os.system('mv in.gb {}'.format(mdir))
                self.mymkdir('{}/out'.format(mdir))

            elif opt in ['run']:
                os.chdir(mdir)
                os.system('mpirun -n 4 lmp_mpi -i in.gb')
                os.chdir(os.pardir)

            elif opt in ['read']:
                os.chdir(mdir)
                j += 1
                self.gbdat[0, j] = ee[0]
                self.gbdat[1, j] = loadtxt('lmp.dat')
                os.chdir(os.pardir)

        if (1 == 0):  # this is not correct
            for i, j in zip(list(range(1, npts)), list(range(0, npts - 1))):
                idx = i + 0.5
                idy = 1

                lx = idx * ux
                ly = idy * uy

                ll = sqrt(lx**2 + ly**2)
                # ll = sqrt(idx**2 + 3./4.) * self.pot["ahcp"]

                ang = rad2deg(arcsin(ly / ll))

                mdir = "dir_{:03}_{:02.3f}".format(i, ang)

                if opt in ['prep']:
                    self.mymkdir(mdir)
                    self.build_hcp_gb_lmp(ang, ll, (30, 40), (0.0, 0.0, 0.0))
                    os.system('mv in.gb {}'.format(mdir))
                    self.mymkdir('{}/out'.format(mdir))

                elif opt in ['run']:
                    os.chdir(mdir)
                    os.system('mpirun -n 4 lmp_mpi -i in.gb')
                    os.chdir(os.pardir)

                elif opt in ['read']:
                    os.chdir(mdir)
                    self.gbdat[0, j] = ang
                    self.gbdat[1, j] = loadtxt('lmp.dat')
                    os.chdir(os.pardir)

    def loop_thickness(self, opt='prep'):
        for i in range(30, 100, 10):
            thk = (i, i + 10)
            mdir = 'thk_{:4.3f}'.format(thk[0])
            if opt in ['prep']:
                self.mymkdir(mdir)
                self.build_hcp_gb(15, thk, (0.0, 0.0, 0.0))
                os.system('mv in.gb {}'.format(mdir))
                self.mymkdir('{}/out'.format(mdir))

            elif opt in ['run']:
                os.chdir(mdir)
                os.system('mpirun -n 4 lmp_mpi -i in.gb')
                os.chdir(os.pardir)

    def loop_dispx_hcp(self, opt='prep'):
        for i in range(20):
            dispx = 0.05 * i
            mdir = 'dispx-{:4.3f}'.format(dispx)
            if opt in ['prep']:
                self.mymkdir(mdir)
                self.build_hcp_gb(0., (30, 40), (dispx, 0.0, 0.0))
                os.system('mv in.gb {}'.format(mdir))
                self.mymkdir('{}/out'.format(mdir))

            elif opt in ['run']:
                os.chdir(mdir)
                os.system('mpirun -n 4 lmp_mpi -i in.gb')
                os.chdir(os.pardir)

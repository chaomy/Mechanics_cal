# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-12-03 11:07:29
# @Last Modified by:   chaomy
# @Last Modified time: 2018-02-21 03:02:26

import os
from numpy import sqrt, arccos, arcsin, rad2deg
from numpy import ndarray, loadtxt
from numpy import mat, unique


class md_gb_hcp_1100(object):

    def find_angles_1100(self):
        ags = []
        res = []

        # npts = 7
        ux = 0.5 * self.pot["chcp"]
        uy = self.pot["ahcp"]

        cnt = 0

        il = [[1, 2, 3, 4, 5, 6, 7, 8],
              [1, 2, 3, 4, 5, 6, 7, 8],
              [1, 2, 3, 4, 5, 6, 7, 8],
              [1, 2, 3, 4, 5, 6, 7, 8]]
        jl = [1, 2, 3, 4]

        for j in jl:
            for i in il[j - 1]:
                idx = i
                idy = j
                lx = idx * ux
                ly = idy * uy

                ll = sqrt(lx**2 + ly**2)
                ang = rad2deg(arcsin(ly / ll))

                if (ang <= 90):
                    cnt += 1
                    ags.append(ang)
                    while (ll < 0):
                        ll *= 2
                    res.append([ang, ll, i, j])

        print(len(ags))
        (a, b) = unique(ags, return_index=True)

        print(len(a), len(b))
        self.ag = []

        for ee in b:
            self.ag.append(res[ee])
        print(self.ag)

    def loop_angle_1100(self, opt='prep'):

        self.find_angles()

        npts = len(self.ag)

        j = 0
        self.gbdat = ndarray([2, npts + 1])
        self.gbdat[0, 0], self.gbdat[0, 1] = 0, 0

        for ee in self.ag:
            mdir = "dir_{:02.03f}_{:02d}_{:02d}".format(ee[0], ee[2], ee[3])

            if opt in ['prep']:
                self.mymkdir(mdir)
                self.build_hcp_gb_lmp(ee[0], ee[1], (40, 50), (0.0, 0.0, 0.0))
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

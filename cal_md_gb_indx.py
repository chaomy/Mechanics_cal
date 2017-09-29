#!/usr/bin/env python
# encoding: utf-8
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-07-05 08:12:30
# @Last Modified by:   chaomy
# @Last Modified time: 2017-09-29 16:03:54

import numpy as np
import os
from cal_md_gb_hcp import inhcp
from math import cos, sin
from math import sqrt


class md_gb_indx(object):

    def hcp_tilt_index(self, thetarange=[0., 60.]):
        lo = np.deg2rad(thetarange[0])
        hi = np.deg2rad(thetarange[1])
        v0 = np.array([1, 0, 0])
        vn = np.array([0, 0, 1])
        SQRT3 = sqrt(3)

        cnt = 0
        ratiolist = []
        rot = {'x': None,
               'y': None,
               'z': None,
               'theta': None} 
        for xidx in range(1, 100):
            for yidx in range(1, 80):
                vec = np.array([xidx, yidx * SQRT3, 0])
                theta = np.arccos(np.dot(vec, v0) / np.linalg.norm(vec))
                if theta < hi:
                    for mulx in range(1, 10):
                        v2 = np.cross(vn, vec)
                        y2x = (v2[1] / SQRT3) / (mulx * v2[0])
                        y2xrnd = np.round(y2x)
                        if np.isclose(y2x, y2xrnd, 1e-3):
                            ratio = y2xrnd / mulx
                            if ratio not in ratiolist:
                                ratiolist.append(ratio)
                                rot['x'] = np.array([xidx, yidx, 0.])
                                rot['y'] = np.array([mulx, y2xrnd, 0])
                                rot['z'] = vn
                                rot['theta'] = np.rad2deg(theta)
                                cnt += 1
                                print rot
            # v1indx = np.array([2 * xidx, yidx, 0])
            # v2indx = np.round(v2)
        print cnt
        return

    def hcp_til_mtx(self, angdeg=30.):
        crat = self.pot['chcp'] / self.pot['ahcp']
        theta = np.deg2rad(angdeg)
        base = np.mat([[1.0, 0.0, 0.0],
                       [0.0, sqrt(3.), 0.0],
                       [0.0, 0.0, crat]])
        matx = np.mat([[cos(theta), sin(theta), 0.0],
                       [-sin(theta), cos(theta), 0.0],
                       [0., 0., 1.]])
        nmat = (matx * base.transpose()).transpose()
        return nmat

    def loop_dispx_hcp(self, opt='prep'):
    	for i in range(20):
    		dispx = 0.05 * i 
    		mdir = 'dispx-{:4.3f}'.format(dispx) 
    		if opt in ['prep']:
	    		self.mymkdir(mdir)
	    		self.build_hcp_gb(20, (dispx, 0.0, 0.0))
	    		os.system('mv in.gb {}'.format(mdir))
	    		self.mymkdir('{}/out'.format(mdir))
	    	elif opt in ['run']:
	    		os.chdir(mdir)
	    		os.system('mpirun -n 8 lmp_mpi -i in.gb')
	    		os.chdir(os.pardir)
    	return

    def build_hcp_gb(self, ang=20., disp=(0, 0, 0)):
        fid = open("in.gb", 'w')
        pmat = self.hcp_til_mtx(ang)
        nmat = self.hcp_til_mtx(-ang)
        fid.write(inhcp % (self.pot['ahcp'], self.pot['chcp'],
                           pmat[0, 0], pmat[0, 1], pmat[0, 2],
                           pmat[1, 0], pmat[1, 1], pmat[1, 2],
                           pmat[2, 0], pmat[2, 1], pmat[2, 2],
                           nmat[0, 0], nmat[0, 1], nmat[0, 2],
                           nmat[1, 0], nmat[1, 1], nmat[1, 2],
                           nmat[2, 0], nmat[2, 1], nmat[2, 2],
                           self.pot['pair_style'], self.pot['file'],
                           disp[0], self.pot['ehcp']))
        return
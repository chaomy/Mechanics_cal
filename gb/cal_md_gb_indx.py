#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-07-05 08:12:30
# @Last Modified by:   chaomy
# @Last Modified time: 2018-05-10 17:01:25

import numpy as np
import os
from math import cos, sin
from math import sqrt


class md_gb_indx(object):

    def hcp_tilt_index(self, thetarange=[0., 60.]):
        # lo = np.deg2rad(thetarange[0])
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
                                print(rot)
            # v1indx = np.array([2 * xidx, yidx, 0])
            # v2indx = np.round(v2)
        print(cnt)

    def hcp_til_mtx_z1120(self, angdeg=15.):
        crat = self.pot['chcp'] / self.pot['ahcp']
        theta = np.deg2rad(angdeg)
        base = np.mat([[sqrt(3.), 0.0, 0.0],
                       [0.0, crat, 0.0],
                       [0.0, 0.0, 1.0]])
        matx = np.mat([[cos(theta), sin(theta), 0.0],
                       [-sin(theta), cos(theta), 0.0],
                       [0.0, 0.0, 1.0]])
        nmat = (matx * base.transpose()).transpose()
        # basis z || 1100
        self.bs = """basis   0.0   0.0   0.0    &
basis   0.5    0.0   0.5    &
basis   ${b1}  0.5   0.0    &
basis   ${b2}  0.5   0.5 

"""
        return nmat

    def hcp_til_mtx_z1100(self, angdeg=15.):
        crat = self.pot['chcp'] / self.pot['ahcp']
        theta = np.deg2rad(angdeg)
        base = np.mat([[1.0, 0.0, 0.0],
                       [0.0, crat, 0.0],
                       [0.0, 0.0, sqrt(3.)]])
        matx = np.mat([[cos(theta), sin(theta), 0.0],
                       [-sin(theta), cos(theta), 0.0],
                       [0.0, 0.0, 1.0]])
        nmat = (matx * base.transpose()).transpose()
        # x || [12-10]; y || 0001;   basis z || 1100
        self.bs = """basis   0.0   0.0   0.0    &
basis   0.5   0.0   0.5    &
basis   0.0   0.5   ${b1}  &
basis   0.5   0.5   ${b2}

"""
        return nmat

    def hcp_til_mtx_z0001(self, angdeg=15.):
        crat = self.pot['chcp'] / self.pot['ahcp']
        theta = np.deg2rad(angdeg)
        base = np.mat([[1.0, 0.0, 0.0],
                       [0.0, sqrt(3.), 0.0],
                       [0.0, 0.0, crat]])
        matx = np.mat([[cos(theta), sin(theta), 0.0],
                       [-sin(theta), cos(theta), 0.0],
                       [0., 0., 1.]])
        nmat = (matx * base.transpose()).transpose()
        # basis z || 0001
        self.bs = """basis   0.0   0.0   0.0    &
basis   0.5   0.5   0.0    &
basis   0.0   ${b1} 0.5    &
basis   0.5   ${b2} 0.5

"""
        return nmat

#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2018-06-03 16:28:53
# @Last Modified by:   chaomy
# @Last Modified time: 2018-06-03 22:42:23

import numpy as np
from numpy import cos, sin, tan, sqrt, arccos, dot, sum


class latConverter():

    # 		output = zeros(6,1);
    # 		output(1,1) = sqrt(sum(M(1,:).^2));
    # 		output(2,1) = sqrt(sum(M(2,:).^2));
    # 		output(3,1) = sqrt(sum(M(3,:).^2));
    #       output(6,1) = acos(sum(M(1,:).*M(2,:))/(output(1,1)*output(2,1)));
    # 		output(5,1) = acos(sum(M(1,:).*M(3,:))/(output(1,1)*output(3,1)));
    #       output(4,1) = acos(sum(M(2,:).*M(3,:))/(output(2,1)*output(3,1)));

    def matToVec(self, M):
        v = np.zeros([6, 1])
        v[0] = sqrt(sum(dot(M[0, :], M[0, :].transpose())))
        v[1] = sqrt(sum(dot(M[1, :], M[1, :].transpose())))
        v[2] = sqrt(sum(dot(M[2, :], M[2, :].transpose())))
        v[5] = arccos(sum(dot(M[0, :], M[1, :].transpose())) / (v[0] * v[1]))
        v[4] = arccos(sum(dot(M[0, :], M[2, :].transpose())) / (v[0] * v[2]))
        v[3] = arccos(sum(dot(M[1, :], M[2, :].transpose())) / (v[1] * v[2]))
        return v

        # output = zeros(3);
        # output(1,1) = v(1);
        # output(2,1) = v(2)*cos(v(6));
        # output(2,2) = v(2)*sin(v(6));
        # output(3,1) = v(3)*cos(v(5));
        # output(3,2) = v(3)*cos(v(4))*sin(v(6))-((v(3)*cos(v(5))-v(3)*cos(v(4))*cos(v(6)))/tan(v(6)));
        # output(3,3) = sqrt(v(3)^2 -output(3,1)^2 - output(3,2)^2);

    def vecToMat(self, v):
        M = np.zeros([3, 3])
        M[0, 0] = v[0]

        M[1, 0] = v[1] * cos(v[5])
        M[1, 1] = v[1] * sin(v[5])
        M[2, 0] = v[2] * cos(v[4])
        M[2, 1] = v[2] * cos(v[3]) * sin(v[5]) - ((v[2] *
                                                   cos(v[4]) - v[2] * cos(v[3]) * cos(v[5])) / tan(v[5]))
        M[2, 2] = sqrt(v[2]**2 - M[2, 0]**2 - M[2, 1]**2)
        return M

if __name__ == '__main__':
    a = latConverter()
    v = a.matToVec(np.mat([[300.,   99.0000,        32],
                           [-1.2644,    4.7391,       0.5],
                           [-5,       8,    2.6935]]))
    m = a.vecToMat(v)
    print(m)

#!/usr/bin/env python
# encoding: utf-8
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-07-05 08:12:30
# @Last Modified by:   chaomy
# @Last Modified time: 2018-02-20 22:50:43

import shutil
import re
from math import sqrt


class md_gb_pre(object):

    def loop_gb_list(self):
        bicrystals = self.read_gb_list()
        for bicrystal in bicrystals:
            self.build_tiltgb_bcc(bicrystal)
        return

    def read_gb_list(self):
        with open('BuildGBlist.lst_Marinica110', 'r') as f:
            bicrystals = []
            raw = f.readlines()
            for lines in raw:
                angle = lines.split()[0]
                angle_phi = lines.split()[1]
                xupper = lines.split()[3] + ' ' + \
                    lines.split()[4] + ' ' + lines.split()[5]
                yupper = lines.split()[6] + ' ' + \
                    lines.split()[7] + ' ' + lines.split()[8]
                zupper = lines.split()[9] + ' ' + \
                    lines.split()[10] + ' ' + lines.split()[11]

                xlower = lines.split()[13] + ' ' + \
                    lines.split()[14] + ' ' + lines.split()[15]
                ylower = lines.split()[16] + ' ' + \
                    lines.split()[17] + ' ' + lines.split()[18]
                zlower = lines.split()[19] + ' ' + \
                    lines.split()[20] + ' ' + lines.split()[21]
                u_mult = lines.split()[22]
                l_mult = lines.split()[23]
                xsh = lines.split()[24]
                zsh = lines.split()[25]
                # print xlower+' '+ylower+' '+zlower
                bicrystals.append([angle, angle_phi,
                                   xupper, yupper, zupper,
                                   xlower, ylower, zlower,
                                   u_mult, l_mult, xsh, zsh])
        f.close()
        return bicrystals

    def build_tiltgb_bcc(self, bicrystal):
        # initialize bicrystal: angle and misorientations
        angle = bicrystal[0]
        angle_phi = bicrystal[1]
        xupper = bicrystal[2]
        yupper = bicrystal[3]
        zupper = bicrystal[4]
        xlower = bicrystal[5]
        ylower = bicrystal[6]
        zlower = bicrystal[7]
        u_mult = int(bicrystal[8])
        l_mult = int(bicrystal[9])
        xsh = bicrystal[10]
        zsh = bicrystal[11]
        alat = 3.1433632
        # symmetric for now but change for asymmetric
        Lx = alat * sqrt(float(xupper.split()[0])**2 + float(
            xupper.split()[1])**2 + float(xupper.split()[2])**2) * u_mult
        Lz = alat * sqrt(float(zupper.split()[0])**2 + float(
            zupper.split()[1])**2 + float(zupper.split()[2])**2)

        # Lz used to be 2
        print(Lx)
        print(Lz)
        stepx = 0.05
        max_x = int(Lx / stepx)
        stepz = 0.05
        max_z = int(Lz / stepz)
        for i in range(1):  # (0,max_x):
            for j in range(1):  # (0,max_z):
                x = float(xsh) + stepx * (i - 0)
                z = float(zsh) + stepz * (j - 0)
                filename = './gb-' + \
                    str(angle) + '-' + str(angle_phi) + \
                    '/mesh-' + str(x) + '-' + str(z)
                try:
                    # shutil.rmtree('filename')
                    shutil.copytree('Base_Marinica', filename)
                    o = open(filename + '/gb4.in', "r")
                    w = open(filename + '/gb4.in_final', "w")
                    data = open(filename + '/gb4.in').read()

                    data = re.sub("latticeparam", str(alat), data)
                    data = re.sub("xsize", str(Lx), data)
                    data = re.sub("zsize", str(Lz), data)

                    data = re.sub("xmove", str(x), data)
                    data = re.sub("zmove", str(z), data)

                    data = re.sub("xorientup", xupper, data)
                    data = re.sub("yorientup", yupper, data)
                    data = re.sub("zorientup", zupper, data)

                    data = re.sub("xorientlow", xlower, data)
                    data = re.sub("yorientlow", ylower, data)
                    data = re.sub("zorientlow", zlower, data)

                    w.write(data)
                    w.close()
                    o.close()
                except OSError:
                    print("Hmmmmm \n")
                else:
                    print(filename)

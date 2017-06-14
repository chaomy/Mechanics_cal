#!/usr/bin/env python

import re
import shutil
import math

def build_gb_dir(bicrystal):
    # initialize bicrystal: angle and misorientations
    angle = bicrystal[0]      # misorientations
    angle_phi = bicrystal[1]  # inclination
    xupper = bicrystal[2]
    yupper = bicrystal[3]
    zupper = bicrystal[4]

    xlower = bicrystal[5]
    ylower = bicrystal[6]
    zlower = bicrystal[7]
    u_mult = int(bicrystal[8])  #
    l_mult = int(bicrystal[9])  #

    alat = 3.615  # BCC  W!!!

    # symmetric for now but change for asymmetric
    Lx = alat * math.sqrt(float(xupper.split()[0])**2 +
                          float(xupper.split()[1])**2 +
                          float(xupper.split()[2])**2) * u_mult

    Lz = alat * math.sqrt(float(zupper.split()[0])**2 +
                          float(zupper.split()[1])**2 +
                          float(zupper.split()[2])**2)

    # Lz used to be 2
    print Lx
    print Lz

    stepx = 0.1
    max_x = int(Lx / stepx)
    stepz = 3.614 / 2
    max_z = int(Lz / stepz)
    for i in range(0, max_x):
        for j in range(1, max_z):
            x = stepx * i
            z = stepz * j
            filename = './gb-' + str(angle) + '-' + str(angle_phi) + \
                '/mesh-' + str(x) + '-' + str(z)
            try:
                # shutil.rmtree('filename')
                shutil.copytree('Base_Cu', filename)

                o = open(filename + '/gb4.in', "r")
                w = open(filename + '/gb4.in_final', "w")
                data = open(filename + '/gb4.in').read()

 #              datax=re.sub("xmove",str(x),data)
 #              dataz=re.sub("zmove",str(z),datax)

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
                print "Hmmmmm \n"
#                msublist.append(dir)
            else:
                print "There \n"
                print filename

# build_gb_dir(2)

def read_gb_list():

    with open('GBlist.lst', 'r') as f:
        answer = []
        bicrystals = []
        read_data = f.readlines()
        for lines in read_data:
            # print lines.split()[5]
            angle = lines.split()[0]
            angle_phi = lines.split()[1]
            xupper = lines.split()[3] + ' ' + lines.split()[4] + ' ' + lines.split()[5]
            yupper = lines.split()[6] + ' ' + lines.split()[7] + ' ' + lines.split()[8]
            zupper = lines.split()[9] + ' ' + lines.split()[10] + ' ' + lines.split()[11]

            xlower = lines.split()[13] + ' ' + lines.split()[14] + ' ' + lines.split()[15]
            ylower = lines.split()[16] + ' ' + lines.split()[17] + ' ' + lines.split()[18]
            zlower = lines.split()[19] + ' ' + lines.split()[20] + ' ' + lines.split()[21]
            u_mult = lines.split()[22]
            l_mult = lines.split()[23]

            # print xlower+' '+ylower+' '+zlower
            bicrystals.append([angle, angle_phi, xupper, yupper, zupper,
                               xlower, ylower, zlower, u_mult, l_mult])
    f.closed
    return bicrystals


if __name__ == '__main__':
    # list of bicrystals: anfle, orient xyz upper, orient xyz lower
    bicrystals = read_gb_list()
    for bicrystal in bicrystals:
        build_gb_dir(bicrystal)

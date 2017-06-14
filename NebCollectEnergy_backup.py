#!/usr/local/bin/python
# encoding: utf-8
#
import os, glob
import re
import matplotlib.pylab as plt
import numpy as np
import pickle as pc

def Nebcollect():
    LogFiles = glob.glob("./screen.*")
    NebEnergy = []
    mEnergy1 = "\s*Energy initial, next-to-last, final =\s*"
    mEnergy2 = "(-?\d*\.\d*)\s*(-?\d*\.\d*)\s*(-?\d*\.\d*)"
    mEnergy = mEnergy1 + mEnergy2
    FindEnergy = re.compile(mEnergy, re.DOTALL)

    #  lz="\striclinic box = \(-?\d*\s*-?\d*\s*-?\d*\)\s*to\s*\((\d*\.\d*)\s*(\d*\.\d*)\s*(\d*\.\d*)\)"
    #  findz = re.compile(lz, re.DOTALL)

    #  with open(LogFiles[0], 'r') as fid:
    #  raw = fid.read()
    #  #  lz = float(findz.findall(raw)[-1][-1])  #  read lz , unit A
    count = 0;

    for file in LogFiles:
        print file
        with open(file, 'r') as fid:
            Raw = fid.read()
            fid.close()

        NebEnergy.append(float(FindEnergy.findall(Raw)[-1][-1]))
        print "%d %f\n"%(count, NebEnergy[count]);
        count += 1

    NebEnergy = np.array(NebEnergy);
    ########################################################
    # for change the unit to be meV/ burger
    ########################################################
    #  NebEnergy = 1000 * NebEnergy * 0.5  #  meV / burger
    #  NebEnergy = np.delete(NebEnergy, np.argmin(NebEnergy));

    print "after delete", len(NebEnergy);
    NebEnergy = NebEnergy - np.min(NebEnergy)

    # plot the results
    fig=plt.figure(figsize=(8,6))
    ax = fig.add_subplot(111)
    ax.get_xaxis().get_major_formatter().set_useOffset(False)

    x = np.linspace(0,1,len(NebEnergy))

    ax.plot(x, NebEnergy, linestyle = '--',
            marker = 'o', markersize = 12, label = 'energy')

    plt.legend(bbox_to_anchor=(0.1,0.1),
                mode='expand',
                borderaxespad=0.01,
                fontsize = 19)

    plt.xlabel("Normalized reaction coordinate",
            {'fontsize': 19})   # (110): -110  (11-2) -110
    plt.ylabel("Energy per length [meV/A]",
            {'fontsize': 19})   #

    plt.yticks(size = 19)
    plt.xticks(size = 19)
    plt.savefig("neb.png",
                bbox_inches='tight', pad_inches=0.03)

#   plt.show()
    fid = open("pickle_data",'w')
    A = pc.Pickler(fid)
    A.dump([x, NebEnergy])
    return

if __name__ == '__main__':
    Nebcollect()

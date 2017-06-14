#!/usr/bin/env python
# encoding: utf-8

###################################################################
#
# File Name : cal_vasp_thermo.py
#
###################################################################
#
# Purpose :
#
# Creation Date :
# Last Modified :
# Created By    : Chaoming Yang
#
###################################################################

import os
import numpy as np
import matplotlib.pylab as plt
from optparse import OptionParser
import ase.io
import ase

try:
    import gn_config
    import get_data
    import gn_kpoints
    import gn_incar
    import gn_pbs
    import Intro_vasp

except ImportError:
    print "error during import"


class cal_vasp_thermo(gn_config.bcc,
                      gn_config.fcc,
                      gn_config.hcp,
                      get_data.get_data,
                      gn_kpoints.gn_kpoints,
                      gn_incar.gn_incar,
                      gn_pbs.gn_pbs,
                      Intro_vasp.vasp_change_box):

    def __init__(self,
                 in_structure='bcc',
                 lattice_constant=3.3224040,
                 in_element='Nb',
                 in_kpoints=[31, 31, 31]):

        self._strct = in_structure
        self._lat0 = lattice_constant
        self._elem = in_element
        self._root = os.getcwd()
        self._kpts = in_kpoints
        return

    def pressure_vs_vol(self, opt='prep'):
        lat0 = self._lat0

        npts = 15
        delta = -0.02

        # pbs
        self.set_nnodes(1)
        self.set_ppn(12)
        self.set_wall_time(30)
        self.set_pbs_type("va")

        stresslist = []
        vollist = []

        for i in range(npts):
            dirname = 'dir-%04d' % (i)
            self.mymkdir(dirname)

            rat = (1 + i * delta) ** (1. / 3.)
            lat = rat * lat0

            if opt == 'prep':
                # build cell and write POSCAR
                bcc_drv = gn_config.bcc(self._elem, lat)
                atoms = bcc_drv.set_bcc_primitive((1, 1, 1))
                ase.io.write(filename="POSCAR",
                             images=atoms,
                             format='vasp')

                # prepare pbs
                self.set_job_title("%s" % (dirname))
                self.set_main_job("mpirun  vasp")
                self.write_pbs(od=False)

                os.system("mv  va.pbs   {}".format(dirname))
                os.system("mv  POSCAR   {}".format(dirname))
                os.system("cp  KPOINTS  {}".format(dirname))
                os.system("cp  POTCAR   {}".format(dirname))
                os.system("cp  INCAR    {}".format(dirname))

            elif opt == 'clc':
                os.chdir(dirname)
                (engy, stress, vol) = self.vasp_energy_stress_vol()
                stresslist.append(0.1 * stress[0])  # Gpa
                print stress
                vollist.append(vol)
                os.chdir(self._root)

        if opt == 'clc':
            plt.plot(vollist, stresslist)
            plt.savefig("v2p.png")
            plt.show()
            np.savetxt("dft_p_v.txt", [vollist, stresslist])
        return


usage = "usage:%prog [options] arg1 [options] arg2"
parser = OptionParser(usage=usage)
parser.add_option("-t",
                  "--mtype",
                  action="store",
                  type="string",
                  dest="mtype",
                  help="",
                  default="prp_r")

(options, args) = parser.parse_args()

if __name__ == '__main__':
    drv = cal_vasp_thermo()

    if options.mtype.lower() == 'p2vprep':
        drv.pressure_vs_vol('prep')

    if options.mtype.lower() == 'p2vclc':
        drv.pressure_vs_vol('clc')

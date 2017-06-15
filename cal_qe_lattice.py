#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: yangchaoming
# @Date:   2017-06-13 15:37:47
# @Last Modified by:   chaomy
# @Last Modified time: 2017-06-15 15:50:06

import os
import numpy as np
import gn_qe_inputs
import md_pot_data
import gn_config
import get_data
import output_data
from optparse import OptionParser
from scipy import interpolate


class cal_lattice(gn_config.bcc,
                  gn_config.fcc,
                  gn_config.hcp,
                  get_data.get_data,
                  output_data.output_data,
                  gn_qe_inputs.gn_qe_infile):

    def __init__(self, inpot=None):
        self.pot = md_pot_data.qe_pot.vca_W75Re25
        output_data.output_data.__init__(self)
        get_data.get_data.__init__(self)
        gn_qe_inputs.gn_qe_infile.__init__(self, self.pot)
        self.alat0 = self.pot['lattice']
        self.element = self.pot['element']
        self.mass = self.pot['mass']
        self.kpnts = [43, 43, 43]
        return

    def cal_bcc_lattice(self):
        bcc_drv = gn_config.bcc(self.pot)
        for i in range(-10, 10):
            if i >= 0:
                dirname = "dir-p-{:03d}".format(i)
            else:
                dirname = "dir-n-{:03d}".format(abs(i))
            self.mymkdir(dirname)
            alat = self.alat0 + i * 0.006
            bcc_drv.set_lattce_constant(alat)
            atoms = bcc_drv.set_bcc_primitive((1, 1, 1))
            self.gn_qe_bcc_lattice_infile(atoms)
            os.system("mv qe.in {}".format(dirname))
        return

    def loop_pots(self):
        count = 0
        potlist = None
        for pot in potlist:
            count += 1
            latticeList = []
            energy = []
            stress = []
            for i in range(-10, 10):
                alat = self.alat0 + i * 0.005
                latticeList.append(alat)
                self.gnInfile(alat,
                              self._kpoint,
                              pot)
                data = self.qe_get_energy_stress()
                energy.append(data[0])
                stress.append(data[1])
                self.output_delta_energy(delta=latticeList,
                                         energy=energy,
                                         stress=stress,
                                         file_name='lat.txt')
        return

    # to be changed
    def loop_kpoints(self):
        for kpoint in range(33, 36):
            latticeList = []
            energy = []
            stress = []
            for i in range(-15, 15):
                alat = self.alat0 + i * 0.004
                latticeList.append(alat)
                self.gnInfile(alat, kpoint)
                data = self.get_data()
                energy.append(data[0])
                stress.append(data[1])
            self.output(latticeList,
                        kpoint,
                        energy,
                        stress)
        return

    def loop_ecut(self):
        bcc_drv = gn_config.bcc(self.pot)
        bcc_drv.set_lattce_constant(self.alat0)
        for ecut in range(30, 50):
            self.set_ecut('{}'.format(ecut))
            dirname = 'dir-ecut-{}'.format(ecut)
            self.mymkdir(dirname)
            atoms = bcc_drv.set_bcc_primitive((1, 1, 1))
            self.gn_qe_bcc_lattice_infile(atoms)
            os.system('mv qe.in {}'.format(dirname))
            os.system('cp $POTDIR/{} {}'.format(self.pot['file'],
                                                dirname))
        return

    def gn_qe_bcc_lattice_infile(self, atoms):
        self.set_thr('1.0D-6')
        with open('qe.in', 'w') as fid:
            fid = self.qe_write_control(fid, atoms)
            fid = self.qe_write_system(fid, atoms)
            fid = self.qe_write_electrons(fid)
            fid = self.qe_write_cell(fid, atoms.get_cell())
            fid = self.qe_write_species(fid, atoms, self.pot)
            fid = self.qe_write_pos(fid, atoms)
            fid = self.qe_write_kpts(fid, self.kpnts)
            fid.close()
        return


def find_lattice(self,
                 Lattice_column,
                 energy):
    Leftpoint = Lattice_column[0]
    Rightpoint = Lattice_column[-1]
    InterPoints = np.linspace(Leftpoint,
                              Rightpoint,
                              201)
    f = interpolate.UnivariateSpline(Lattice_column[:],
                                     energy[:], s=0)(InterPoints)
    i = np.argmin(f)
    return InterPoints[i - 1]


if __name__ == '__main__':
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
    drv = cal_lattice()

    if options.mtype.lower() == 'kpoints':
        drv.loop_kpoints()

    elif options.mtype.lower() == 'ecut':
        drv.loop_ecut()

    elif options.mtype.lower() == 'pot':
        drv.loop_pots()

    elif options.mtype.lower() == 'bcc':
        drv.cal_bcc_lattice()

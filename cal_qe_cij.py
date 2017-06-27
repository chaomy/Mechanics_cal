#!/usr/bin/env python

###################################################################
#
# File Name : ./cal_qe_cij.py
#
###################################################################
#
# Purpose : cal cij by vasp
#
# Creation Date :
# Last Modified : Sat Apr  1 23:15:50 2017
# Created By    : Chaoming Yang
#
###################################################################

import glob
import os
import numpy as np
import md_pot_data
from optparse import OptionParser
from scipy.optimize import leastsq

try:
    import gn_config
    import gn_qe_inputs
    import get_data
    import gn_kpoints
    import gn_incar
    import gn_pbs
    import output_data

except ImportError:
    print "error during import"


class cal_cij(gn_config.bcc,
              gn_config.fcc,
              gn_config.hcp,
              get_data.get_data,
              gn_kpoints.gn_kpoints,
              gn_incar.gn_incar,
              gn_pbs.gn_pbs,
              output_data.output_data,
              gn_qe_inputs.gn_qe_infile):

    def __init__(self):
        self.unit_delta = 0.005
        self.looptimes = 10
        self.volume = None
        self.energy0 = None

        self.pot = md_pot_data.qe_pot.vca_W75Re25
        self.alat = self.pot['lattice']
        self.elem = self.pot['element']
        self.cij_type_list = ['c11', 'c12', 'c44']
        self.cij_type = 'c11'
        self.struct = self.pot['structure']
        self.root = os.getcwd()

        gn_kpoints.gn_kpoints.__init__(self)
        get_data.get_data.__init__(self)
        gn_incar.gn_incar.__init__(self)
        gn_pbs.gn_pbs.__init__(self)
        output_data.output_data.__init__(self)
        gn_qe_inputs.gn_qe_infile.__init__(self, self.pot)

        if self.struct == 'bcc':
            gn_config.bcc.__init__(self, self.pot)
        elif self.struct == 'fcc':
            gn_config.fcc.__init__(self, self.pot)
        elif self.struct == 'hcp':
            gn_config.hcp.__init__(self, self.pot)
        return

    def fit_para(self, delta_list, energy_list):
        xdata = delta_list
        ydata = energy_list - self.energy0

        def residuals(p):
            a, c = p
            return ydata - (a * (0.5 * xdata**2) + c)
        r = leastsq(residuals, [1, 0])
        a, c = r[0]
        a = a / self.volume * self.ev_angstrom3_to_GPa
        return a

    def obtain_cij(self, opt='np'):
        self.set_volume_energy0()
        if opt is 'np':
            raw1 = np.loadtxt("data_c11.txt")
            raw2 = np.loadtxt("data_c12.txt")
            raw3 = np.loadtxt("data_c44.txt")
        else:
            delta_list, energy_list = self.get_delta_energy("data_c11.txt")

        delta_list, energy_list = raw1[:, 0], raw1[:, 1]
        c11_plus_c12 = self.fit_para(delta_list, energy_list)

        delta_list, energy_list = raw2[:, 0], raw2[:, 1]
        c11_minus_c12 = self.fit_para(delta_list, energy_list)

        delta_list, energy_list = raw3[:, 0], raw3[:, 1]
        c44 = self.fit_para(delta_list, energy_list)

        c11 = (0.111111111111111 * c11_plus_c12 +
               0.333333333333333 * c11_minus_c12)
        c12 = (0.111111111111111 * c11_plus_c12 -
               0.166666666666667 * c11_minus_c12)
        c44 = 0.25 * c44

        with open("cij.dat", 'w') as fout:
            fout.write("C11\t%f\t\nC12\t%f\t\nC44\t%f\t\n" % (c11, c12, c44))
            fout.close()
        return

    def set_volume_energy0(self):
        raw = np.loadtxt("equilibrium.txt")
        self.volume = np.average(raw[:, 0])
        self.energy0 = np.average(raw[:, 1])
        return

    def set_cij_type(self, cij_type):
        self.cij_type = cij_type
        return

    def gn_qe_cij_infile(self, atoms):
        self.set_thr('1.0D-6')
        self.set_ecut('48')
        self.set_kpnts((44, 44, 44))
        self.set_degauss('0.03D0')
        with open('qe.in', 'w') as fid:
            fid = self.qe_write_control(fid, atoms)
            fid = self.qe_write_system(fid, atoms)
            fid = self.qe_write_electrons(fid)
            fid = self.qe_write_cell(fid, atoms.get_cell())
            fid = self.qe_write_species(fid, atoms, self.pot)
            fid = self.qe_write_pos(fid, atoms)
            fid = self.qe_write_kpts(fid)
            fid.close()
        return

    def loop_prepare_cij(self):
        for i in range(len(self.cij_type_list)):
            mtype = self.cij_type_list[i]
            self.set_cij_type(mtype)

            for j in range(-self.looptimes, self.looptimes):
                delta = self.unit_delta * j
                if j >= 0:
                    dirname = "dir-%s-p%03d" % (mtype, j)
                else:
                    dirname = "dir-%s-n%03d" % (mtype, np.abs(j))
                self.mymkdir(dirname)
                os.chdir(dirname)
                self.set_pbs(dirname)
                atoms = self.write_bcc_primitive_with_strain(delta=delta,
                                                             in_tag=mtype,
                                                             write=False)
                self.gn_qe_cij_infile(atoms)
                os.system("cp $POTDIR/{} .".format(self.pot['file']))
                os.chdir(self.root)
        return

    def collect_data_cij(self):
        for i in range(len(self.cij_type_list)):
            mtype = self.cij_type_list[i]
            self.set_cij_type(mtype)
            out_file_name = "data_%s.txt" % (mtype)

            for j in range(-self.looptimes, self.looptimes):
                delta = self.unit_delta * j
                if j >= 0:
                    dirname = "dir-%s-p%03d" % (mtype, j)
                else:
                    dirname = "dir-%s-n%03d" % (mtype, np.abs(j))
                os.chdir(dirname)
                print "i am in ", dirname
                (energy, vol, stress) = self.qe_get_energy_stress()
                os.chdir(self.root)
                self.output_delta_energy(delta, energy,
                                         file_name=out_file_name)
                if j == 0:
                    self.output_equilibrium(energy=energy, volume=volume)
        return

    def set_pbs(self, dirname):
        self.set_wall_time(2)
        self.set_job_title(dirname)
        self.set_nnodes(1)
        self.set_ppn(12)
        self.set_main_job("mpirun pw.x < qe.in > qe.out")
        self.write_pbs(od=True)
        return

    def loop_sub_jobs(self):
        dir_list = glob.glob("dir-*")
        for i in range(len(dir_list)):
            os.chdir(dir_list[i])
            os.system("qsub va.pbs")
            os.chdir(self.root)
        return

if __name__ == "__main__":
    usage = "usage:%prog [options] arg1 [options] arg2"
    parser = OptionParser(usage=usage)
    parser.add_option("-t",
                      "--mtype",
                      action="store",
                      type="string",
                      dest="mtype",
                      help="",
                      default="prep")
    (options, args) = parser.parse_args()
    Job = cal_cij()
    if options.mtype.lower() in ['prep']:
        Job.loop_prepare_cij()

    if options.mtype.lower() == 'sub':
        Job.loop_sub_jobs()

    if options.mtype.lower() == 'cij':
        Job.obtain_cij()

    if options.mtype.lower() == 'cnt':
        Job.cal_cij_continue()

    if options.mtype.lower() == 'clc':
        Job.collect_data_cij()

    #  A.set_volume_energy0()
    #  A.obtain_cij_old()
    #  c12 = 0.25 * C11plus2C12 + 0.75 * C11minusC12
    #  c11 = 0.25 * C11plus2C12 - 0.25 * C11minusC12
    #  c44 = float(get_cxx("c44",Lattice_constant,delta_increment,nt,E0)) * (1./4.) * inv_V0  * unit_conversion1

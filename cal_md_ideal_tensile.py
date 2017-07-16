#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: yang37
# @Date:   2017-06-12 17:03:43
# @Last Modified by:   chaomy
# @Last Modified time: 2017-07-16 12:16:07


import os
import ase
import ase.io
import ase.lattice
import glob
import gn_pbs
import numpy as np
import gn_config
import get_data
import plt_drv
import md_pot_data
from scipy.optimize import minimize
from md_pot_data import unitconv
from optparse import OptionParser


class cal_bcc_ideal_tensile(get_data.get_data,
                            gn_pbs.gn_pbs,
                            plt_drv.plt_drv):

    def __init__(self):
        # self.pot = self.load_data('../pot.dat')
        self.pot = md_pot_data.md_pot.Nb_adp
        gn_pbs.gn_pbs.__init__(self)
        plt_drv.plt_drv.__init__(self)
        self.alat = self.pot['lattice']
        self.npts = 25
        self.delta = 0.02

        e1 = np.array([1., 0., 0.])
        e2 = np.array([0., 1., 0.])
        e3 = np.array([0., 0., 1.])

        e1 = e1 / np.linalg.norm(e1)
        e2 = e2 / np.linalg.norm(e2)
        e3 = e3 / np.linalg.norm(e3)

        self.basis = np.mat([e1, e2, e3])
        get_data.get_data.__init__(self)

        self.va_prim = np.mat([[-0.5, 0.5, 0.5],
                               [0.5, -0.5, 0.5],
                               [0.5, 0.5, -0.5]])

        self.configdrv = gn_config.bcc(self.pot)
        self.lm_prim = self.configdrv.lmp_change_box(self.va_prim)
        self.root = os.getcwd()
        return

    def read_stress(self):
        (engy, stress, vol) = self.vasp_energy_stress_vol()
        stress = stress * 0.1
        print "engy {}".format(engy)
        print "{}, {}, {}, {}, {}, {}".format(
            stress[0][0], stress[1][0], stress[2][0],
            stress[3][0], stress[4][0], stress[5][0])
        return

    def lmp_relax(self):
        delta = 0.04
        x0 = np.array([1., 1.])
        res = minimize(self.runlmp, x0, delta,
                       method='Nelder-Mead',
                       options={'xtol': 1e-3, 'disp': True})
        print res.fun
        print res.x
        return

    def loop_collect_vasp(self):
        dirlist = glob.glob("dir-*")
        npts = len(dirlist)
        data = np.ndarray([npts, 10])
        for i in range(npts):
            dirname = dirlist[i]
            print dirname
            os.chdir(dirname)
            atoms = ase.io.read('CONTCAR', format='vasp')
            (engy, stress, vol) = self.vasp_energy_stress_vol()
            cell = atoms.get_cell()
            os.chdir(self.root)
            data[i, 0:4] = cell[0, 0] - 1.0, engy, cell[1, 1], cell[2, 2]
            data[i, 4:] = stress.transpose()
        np.savetxt("iten.txt", data)
        return

    def loop_tensile_lmp(self):
        x0 = np.array([0.91, 1.11])
        npts = self.npts
        data = np.ndarray([npts, 4])
        for i in range(npts):
            delta = self.delta * i
            res = minimize(self.runlmp, x0, delta,
                           method='Nelder-Mead',
                           options={'disp': True})
            x0 = res.x
            print res
            data[i][0] = (delta)
            data[i][1] = res.fun
            data[i][2:] = res.x
        np.savetxt("ishear.txt", data)
        return

    def runvasp(self, x, delta):
        basis = self.basis
        strain = np.mat([[1.0 + delta, 0.0,  0.0],
                         [0.0,  x[0],  0.0],
                         [0.0,  0.0,  x[1]]])
        new_strain = basis.transpose() * strain * basis
        self.gn_primitive_lmps(new_strain, 'vasp')
        os.system("mpirun vasp > vasp.log")
        (engy, stress, vol) = self.vasp_energy_stress_vol()
        print engy
        return engy

    def runlmp(self, x, delta):
        basis = self.basis
        strain = np.mat([[1.0 + delta, 0.0,  0.0],
                         [0.0,  x[0],  0.0],
                         [0.0,  0.0,  x[1]]])
        new_strain = basis.transpose() * strain * basis
        self.gn_primitive_lmps(new_strain, 'lmp')
        os.system("lmp_mpi -i in.init -screen  no")
        raw = np.loadtxt("out.txt")
        engy = raw[0]
        print engy
        return engy

    def gn_primitive_lmps(self,
                          strain=np.mat([[1., 0., 0.],
                                         [0., 1., 0.],
                                         [0., 0., 1.]]),
                          tag='lmp'):
        alat = self.alat
        bas = np.mat([[-0.5, 0.5, 0.5],
                      [0.5, -0.5, 0.5],
                      [0.5, 0.5, -0.5]])
        ##########################################################
        # very important (vasp add strain is basis right time strain)
        ##########################################################
        if tag == 'vasp':
            va_bas = bas * strain
            # poscar input type
            cell = alat * va_bas
            atoms = ase.Atoms('Nb',
                              positions=[[0, 0, 0]],
                              cell=cell,
                              pbc=[1, 1, 1])
            ase.io.write("POSCAR", images=atoms, format='vasp')
        if tag == 'lmp':
            # convert to lammps data style
            lmp_bas = bas * strain
            lmp_bas = self.configdrv.lmp_change_box(lmp_bas)
            cell = alat * lmp_bas
            pos = np.array([[0, 0, 0]])
            file_name = 'init.txt'
            atom_num = 1
            with open(file_name, mode="w") as fout:
                fout.write("#lmp data config")
                fout.write("\n")
                fout.write("%d atoms\n" % (1))
                fout.write("1 atom types\n")
                fout.write("%f\t%f xlo xhi\n" % (0, cell[0, 0]))
                fout.write("%f\t%f ylo yhi\n" % (0, cell[1, 1]))
                fout.write("%f\t%f zlo zhi\n" % (0, cell[2, 2]))
                fout.write("%f  %f  %f xy xz yz\n"
                           % (cell[1, 0],
                              cell[2, 0],
                              cell[2, 1]))
                fout.write("Atoms\n")
                fout.write("\n")
                for i in range(atom_num):
                    fout.write("%d  1  %12.7f %12.7f %12.7f\n"
                               % (i + 1, pos[i, 0], pos[i, 1], pos[i, 2]))
            fout.close()
        return

    def vasp_relax(self, given=True):
        data = np.loadtxt("strain.txt")
        if given is True:
            delta = data[0]
            x0 = data[-2:]
            print delta
            print x0
        else:
            delta = data
            x0 = np.array([1., 1.])
        data = np.zeros(4)
        res = minimize(self.runvasp, x0, delta,
                       method='Nelder-Mead',
                       options={'fatol': 5e-4, 'disp': True})
        print res
        data[0] = delta
        data[1] = res.fun
        data[2:] = res.x
        np.savetxt("ishear.txt", data)
        return

    def set_pbs(self, dirname, delta):
        self.set_nnodes(1)
        self.set_ppn(12)
        self.set_job_title("%s" % (dirname))
        self.set_wall_time(50)
        self.set_main_job("""
        python ~/My_cal/Mechnical_cal/cal_md_ideal_tensile.py  -t  ivasp
                          """)
        self.write_pbs(od=False)
        os.system("mv va.pbs %s" % (dirname))
        return

    def loop_prep_vasp(self):
        npts = self.npts
        for i in range(npts):
            delta = self.delta * i
            dirname = "dir-{:03d}".format(i)
            self.mymkdir(dirname)

            os.system("echo {} > strain.txt".format(delta))
            os.system("mv strain.txt {}".format(dirname))
            self.set_pbs(dirname, delta)
        return

    def loop_prep_vasp_given(self):
        data = np.loadtxt("ishear.txt")
        for i in range(len(data)):
            dirname = "dir-{:03d}".format(i)
            self.mymkdir(dirname)
            os.system("cp KPOINTS {}".format(dirname))
            os.system("cp INCAR   {}".format(dirname))
            os.system("cp POTCAR  {}".format(dirname))
            np.savetxt("strain.txt", data[i])
            os.system("mv strain.txt {}".format(dirname))
            self.set_pbs(dirname, data[i][0])
        return

    def convert_stress(self):
        from scipy.interpolate import InterpolatedUnivariateSpline
        raw = np.loadtxt("ishear.txt")
        data = np.zeros((len(raw),
                         len(raw[0]) + 1))
        data[:, :-1] = raw
        convunit = unitconv.ustress['evA3toGpa']
        vol = np.zeros(len(raw))
        vperf = 0.5 * self.alat**3
        strmat = np.zeros([3, 3])
        for i in range(len(raw)):
            strmat[0, 0], strmat[1, 1], strmat[2, 2] = \
                1 + raw[i, 0], raw[i, 2], raw[i, 3]
            strmat = np.mat(strmat)
            vol[i] = vperf * np.linalg.det(strmat)
        tag = 'interp'
        if tag == 'interp':
            # interpolate
            spl = InterpolatedUnivariateSpline(raw[:, 0], raw[:, 1], k=1)
            splder1 = spl.derivative()
            for i in range(len(raw)):
                # append the stress to the last column
                print "vol", vol[i]
                data[i, -1] = splder1(raw[i, 0]) * convunit / vol[i]
        print data
        np.savetxt("istress.txt", data)
        return


if __name__ == '__main__':
    usage = "usage:%prog [options] arg1 [options] arg2"
    parser = OptionParser(usage=usage)
    parser.add_option('-t', "--mtype",
                      action="store",
                      type="string",
                      dest="mtype", help="",
                      default="prp_r")

    parser.add_option('-c', "--delta",
                      action="store",
                      type='float', dest="delta",
                      default=0.02)

    (options, args) = parser.parse_args()
    drv = cal_bcc_ideal_tensile()
    if options.mtype.lower() == 'ilmp':
        drv.loop_tensile_lmp()

    if options.mtype.lower() == 'lmprelax':
        drv.lmp_relax()

    if options.mtype.lower() == 'vaspprep':
        if os.path.isfile("ishear.txt"):
            drv.loop_prep_vasp_given()
        else:
            drv.loop_prep_vasp()

    if options.mtype.lower() == 'ivasp':
        drv.vasp_relax()

    if options.mtype.lower() == 'clcvasp':
        drv.loop_collect_vasp()

    if options.mtype.lower() == 'plt':
        drv.convert_stress()
        drv.plt_energy_stress()

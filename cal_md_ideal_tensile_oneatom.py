#!/usr/bin/env python
# encoding: utf-8

###################################################################
#
# File Name : cal_md_ideal_tensile_oneatom.py
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
from optparse import OptionParser


class cal_bcc_ideal_tensile(get_data.get_data,
                            gn_pbs.gn_pbs,
                            plt_drv.plt_drv,
                            gn_config.bcc):

    def __init__(self):
        # self.pot = self.load_data('../pot.dat')
        self.pot = md_pot_data.dft_pot.Nb_pbe
        gn_pbs.gn_pbs.__init__(self)
        plt_drv.plt_drv.__init__(self)
        get_data.get_data.__init__(self)
        gn_config.bcc.__init__(self, self.pot)

        self.alat = self.pot['lattice']
        self.npts = 5
        self.delta = 0.02

        e1 = np.array([1., 0., 0.])
        e2 = np.array([0., 1., 0.])
        e3 = np.array([0., 0., 1.])

        e1 = e1 / np.linalg.norm(e1)
        e2 = e2 / np.linalg.norm(e2)
        e3 = e3 / np.linalg.norm(e3)

        self.basis = np.mat([e1, e2, e3])
        self.va_prim = np.mat([[-0.5, 0.5, 0.5],
                               [0.5, -0.5, 0.5],
                               [0.5, 0.5, -0.5]])
        self.root = os.getcwd()
        return

    def loop_collect_vasp(self):
        dirlist = glob.glob("dir-*")
        npts = len(dirlist)
        data = np.ndarray([npts, 10])
        for i in range(npts):
            #  dirname = "dir-{:03d}".format(i)
            dirname = dirlist[i]
            print dirname
            os.chdir(dirname)
            raw = np.loadtxt("iten.txt")
            (engy, stress, vol) = self.vasp_energy_stress_vol()
            os.chdir(self.root)
            data[i, 0:4] = raw
            data[i, 4:] = stress.transpose()
        np.savetxt("istress.txt", data)
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
            data[i][0] = delta
            data[i][1] = res.fun
            data[i][2:] = res.x
        np.savetxt("iten.txt", data)
        return

    def runvasp(self, x, delta):
        basis = self.basis
        strain = np.mat([[1.0 + delta, 0.0, 0.0],
                         [0.0, x[0], 0.0],
                         [0.0, 0.0, x[1]]])
        new_strain = basis.transpose() * strain * basis
        self.gn_primitive_lmps(new_strain, 'vasp')
        os.system("mpirun vasp > vasp.log")
        (engy, stress, vol) = self.vasp_energy_stress_vol()
        self.recordstrain(delta, x, engy)
        print engy
        return engy

    def runlmp(self, x, delta):
        basis = self.basis
        strain = np.mat([[1.0 + delta, 0.0, 0.0],
                         [0.0, x[0], 0.0],
                         [0.0, 0.0, x[1]]])
        new_strain = basis.transpose() * strain * basis
        self.gn_primitive_lmps(new_strain, 'lmp')
        os.system("lmp_mpi -i in.init -screen  no")
        raw = np.loadtxt("out.txt")
        print raw
        self.recordstrain(delta, x, raw)
        return raw[0]

    def gn_primitive_lmps(self,
                          strain=np.mat(np.identity(3)),
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
            lmp_bas = self.lmp_change_box(lmp_bas)
            cell = alat * lmp_bas
            atoms = ase.Atoms('Nb',
                              positions=[[0., 0., 0.]],
                              cell=cell,
                              pbc=[1, 1, 1])
            self.write_lmp_config_data(atoms, 'init.txt')
        return

    def load_input_params(self):
        if os.path.isfile('restart.txt'):
            data = np.loadtxt("restart.txt")
            delta = data[0]
            x0 = data[2:4]
            print delta
        else:
            data = np.loadtxt("strain.txt")
            delta = data
            x0 = np.array([1., 1.])
        return (delta, x0)

    def vasp_relax(self):
        (delta, x0) = self.load_input_params()
        data = np.zeros(4)
        res = minimize(self.runvasp, x0, delta,
                       method='Nelder-Mead',
                       options={'fatol': 2e-3, 'disp': True})
        print res
        data[0] = delta
        data[1] = res.fun
        data[2:] = res.x
        np.savetxt("iten.txt", data)
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

    def recordstrain(self, delta, x, fval):
        fid = open("s{:4.3f}.txt".format(delta), "a")
        formatstr = '{:6.5f} ' * (len(x) + len(fval))
        formatstr += '\n'
        fid.write(formatstr.format(x[0], x[1], *fval))
        fid.close()
        return

    def loop_prep_restart(self, opt='va'):
        raw = np.mat(np.loadtxt("iten.txt"))
        for i in range(len(raw)):
            dirname = "dir-{:03d}".format(i)
            self.mymkdir(dirname)
            np.savetxt("restart.txt", raw[i])
            if opt in ['va', 'vasp']:
                self.copy_inputs(dirname, 'KPOINTS',
                                 'INCAR', 'POTCAR', 'restart.txt')
            elif opt in ['qe']:
                os.system("mv restart.txt {}".format(dirname))
                os.system('cp $POTDIR/{}  {}'.format(self.pot['file'],
                                                     dirname))
            self.set_pbs(dirname, raw[i][0])
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

    if options.mtype.lower() == 'ivasp':
        drv.vasp_relax()

    if options.mtype.lower() == 'clcvasp':
        drv.loop_collect_vasp()

    if options.mtype.lower() == 'plt':
        drv.convert_stress()

    if options.mtype.lower() in ['qe_restart', 'va_restart']:
        opt = options.mtype.lower().split('_')[0]
        drv.loop_prep_restart(opt)

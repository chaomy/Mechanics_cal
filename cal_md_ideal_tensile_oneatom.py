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
        self.npts = 20
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
        self.stress = np.zeros(6)
        return

    def runvasp_tp(self, x, delta):
        basis = self.basis
        strain = np.mat([[1.0 + delta, 0.0, 0.0],
                         [0.0, x, 0.0],
                         [0.0, 0.0, x]])
        new_strain = basis.transpose() * strain * basis
        self.gn_primitive_lmps(new_strain, 'vasp')
        os.system("mpirun vasp > vasp.log")
        (engy, stress, vol) = self.vasp_energy_stress_vol()
        self.stress = stress.flatten()
        print engy, self.stress, vol
        self.recordstrain(delta, [x[0], x[0]], engy)
        return engy

    def runlmp_tp(self, x, delta):
        basis = self.basis
        strain = np.mat([[1.0 + delta, 0.0, 0.0],
                         [0.0, x, 0.0],
                         [0.0, 0.0, x]])
        new_strain = basis.transpose() * strain * basis
        self.gn_primitive_lmps(new_strain, 'lmp')
        os.system("lmp_mpi -i in.init -screen  no")
        raw = np.loadtxt("out.txt")
        engy = raw[0]
        self.stress = raw[1:]
        self.recordstrain(delta, [x[0], x[0]], engy)
        return engy

    def runvasp_op(self, x, delta):
        basis = self.basis
        strain = np.mat([[1.0 + delta, 0.0, 0.0],
                         [0.0, x[0], 0.0],
                         [0.0, 0.0, x[1]]])
        new_strain = basis.transpose() * strain * basis
        self.gn_primitive_lmps(new_strain, 'vasp')
        os.system("mpirun vasp > vasp.log")
        (engy, stress, vol) = self.vasp_energy_stress_vol()
        print engy, stress, vol
        self.stress = stress.flatten()
        self.recordstrain(delta, x, engy)
        return engy

    def runlmp_op(self, x, delta):
        basis = self.basis
        strain = np.mat([[1.0 + delta, 0.0, 0.0],
                         [0.0, x[0], 0.0],
                         [0.0, 0.0, x[1]]])
        new_strain = basis.transpose() * strain * basis
        self.gn_primitive_lmps(new_strain, 'lmp')
        os.system("lmp_mpi -i in.init -screen  no")
        raw = np.loadtxt("out.txt")
        engy = raw[0]
        self.stress = raw[1:]
        self.recordstrain(delta, x, engy)
        return engy

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

    def loop_tensile_lmp(self, opt='op'):
        npts = self.npts
        data = np.ndarray([npts, 10])
        if opt == 'op':
            x0 = np.array([0.91, 1.11])
        elif opt == 'tp':
            x0 = 0.91
        for i in range(npts):
            delta = self.delta * i
            if opt == 'op':
                res = minimize(self.runlmp_op, x0, delta,
                               method='Nelder-Mead',
                               options={'disp': True})
                data[i][2], data[i][3] = res.x[0], res.x[1]
            elif opt == 'tp':
                res = minimize(self.runlmp_tp, x0, delta,
                               method='Nelder-Mead',
                               options={'disp': True})
                data[i][2], data[i][3] = res.x, res.x
            x0 = res.x
            print res
            data[i][0] = delta
            data[i][1] = res.fun
            data[i][4:] = self.stress
        np.savetxt("iten.txt", data)
        return

    def vasp_relax(self, opt='op'):
        (delta, x0) = self.load_input_params()
        data = np.zeros(8 + len(x0))
        if opt == 'op':
            res = minimize(self.runvasp_op, x0, delta,
                           method='Nelder-Mead',
                           options={'fatol': 5e-4, 'disp': True})
            data[2], data[3] = res.x[0], res.x[1]

        elif opt == 'tp':
            x0 = x0[0]
            res = minimize(self.runvasp_tp, x0, delta,
                           method='Nelder-Mead',
                           options={'fatol': 5e-4, 'disp': True})
            data[2], data[3] = res.x, res.x
        data[0] = delta
        data[1] = res.fun
        data[-6:] = self.stress
        print res
        np.savetxt("iten.txt", data)
        return

    def set_pbs(self, dirname, delta, opt='tp'):
        self.set_nnodes(1)
        self.set_ppn(12)
        self.set_job_title("iva_{}_{}".format(opt, dirname))
        self.set_wall_time(50)
        self.set_main_job("""
        ../cal_md_ideal_tensile_oneatom.py  -t ivasp_{}
                          """.format(opt))
        self.write_pbs(od=True)
        os.system("mv va.pbs %s" % (dirname))
        return

    def recordstrain(self, delta, x, fval):
        fid = open("s{:4.3f}.txt".format(delta), "a")
        formatstr = '{:6.5f} ' * (1 + 2 + 6)
        formatstr += '\n'
        fid.write(formatstr.format(fval, x[0], x[1], *self.stress))
        fid.close()
        return

    def loop_prep_restart(self, opt1='va', opt2='tp'):
        raw = np.mat(np.loadtxt("iten.txt"))
        for i in range(len(raw)):
            dirname = "dir-{:03d}".format(i)
            self.mymkdir(dirname)
            np.savetxt("restart.txt", raw[i])
            if opt1 in ['va', 'vasp']:
                self.copy_inputs(dirname, 'KPOINTS',
                                 'INCAR', 'POTCAR', 'restart.txt')
            elif opt1 in ['qe']:
                os.system("mv restart.txt {}".format(dirname))
                os.system('cp $POTDIR/{}  {}'.format(self.pot['file'],
                                                     dirname))
            self.set_pbs(dirname, raw[i][0], opt2)
        return

    def loop_collect(self, opt='va'):
        dirlist = glob.glob("dir-*")
        npts = len(dirlist)
        data = np.ndarray([npts, 10])
        # delta, engy, x, stress
        for i in range(npts):
            dirname = dirlist[i]
            print dirname
            if os.path.isfile('{}/iten.txt'.format(dirname)):
                data[i, :] = np.loadtxt('{}/iten.txt'.format(dirname))
        np.savetxt("iten.txt", data)
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
    if options.mtype.lower() in ['ilmp_tp', 'ilmp_op']:
        opt = options.mtype.lower().split('_')[-1]
        drv.loop_tensile_lmp(opt)

    if options.mtype.lower() in ['ivasp_op', 'ivasp_tp']:
        opt = options.mtype.lower().split('_')[-1]
        drv.vasp_relax(opt)

    if options.mtype.lower() in ['clc_va', 'clc_qe']:
        opt = options.mtype.lower().split('_')[-1]
        drv.loop_collect(opt)

    if options.mtype.lower() in ['restart_qe_op', 'restart_qe_tp',
                                 'restart_va_op', 'restart_va_tp']:
        opt1 = options.mtype.lower().split('_')[-2]
        opt2 = options.mtype.lower().split('_')[-1]
        drv.loop_prep_restart(opt1, opt2)

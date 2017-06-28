#!/usr/bin/env python
# encoding: utf-8

###################################################################
#
# File Name : cal_md_ideal_tensile_tp.py
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
import gn_pbs
import numpy as np
import gn_config
import get_data
import plt_drv
import md_pot_data
import gn_config
from scipy.optimize import minimize
from optparse import OptionParser


class cal_bcc_ideal_tensile_tp(get_data.get_data,
                               gn_pbs.gn_pbs,
                               plt_drv.plt_drv,
                               gn_config.bcc):

    def __init__(self):
        self.pot = self.load_data('../pot.dat')
        # self.pot = md_pot_data.md_pot.Nb_adp
        gn_pbs.gn_pbs.__init__(self)
        plt_drv.plt_drv.__init__(self)
        get_data.get_data.__init__(self)
        gn_config.bcc.__init__(self, self.pot)

        self.alat = self.pot['lattice']
        self.range = (16, 21)
        self.npts = self.range[1] - self.range[0]
        self.delta = 0.02
        e1 = [1., 0., 0.]
        e2 = [0., 1., 0.]
        e3 = [0., 0., 1.]
        self.basis = np.mat([e1, e2, e3])
        self.stress = None
        self.root = os.getcwd()
        return

    def loop_tensile_lmp(self):
        x0 = np.array([0.93, 1.06])
        npts = self.range[1] - self.range[0]
        data = np.ndarray([npts, (8 + len(x0))])
        for i in range(npts):
            delta = self.delta * (self.range[0] + i)
            res = minimize(self.runlmp, x0, delta,
                           method='Nelder-Mead',
                           options={'fatol': 2e-3, 'disp': True})
            x0 = res.x
            print res
            data[i][0] = delta
            data[i][1] = res.fun
            data[i][2:(2 + len(res.x))] = res.x
            data[i][-6:] = self.stress
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
        self.stress = stress.flatten()
        print engy, self.stress
        return engy

    def runlmp(self, x, delta):
        basis = self.basis
        strain = np.mat([[1.0 + delta, 0.0, 0.0],
                         [0.0, x[0], 0.0],
                         [0.0, 0.0, x[1]]])
        new_strain = basis.transpose() * strain * basis
        self.gn_convention_lmps(new_strain, 'lmp')
        os.system("lmp_mpi -i in.init -screen  no")
        raw = np.loadtxt("out.txt")
        self.recordstrain(delta, x, raw)
        self.stress = raw[1:]
        print raw
        return raw[0]

    def gn_convention_lmps(self,
                           strain=np.mat(np.identity(3)),
                           tag='lmp'):
        alat = self.alat
        bas = np.mat(np.identity(3))
        ##########################################################
        # very important (vasp add strain is basis right time strain)
        ##########################################################
        va_bas = bas * strain
        cell = alat * va_bas
        atoms = ase.lattice.cubic.BodyCenteredCubic(directions=[[1., 0., 0.],
                                                                [0., 1., 0.],
                                                                [0., 0., 1.]],
                                                    latticeconstant=self.alat,
                                                    size=(1, 1, 1),
                                                    symbol=self.pot['element'],
                                                    pbc=(1, 1, 1))
        pos = np.mat([[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]) * cell
        atoms.set_cell(cell)
        atoms.set_positions(pos)
        if tag == 'vasp':
            ase.io.write("POSCAR", images=atoms, format='vasp')
        if tag == 'lmp':
            self.write_lmp_config_data(atoms, 'init.txt')
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

    def plt_energy_stress(self):
        raw = np.loadtxt("iten.txt")
        raw = raw[raw[:, 0].argsort()]
        print raw
        self.set_keys()
        self.set_211plt()
        self.ax1.plot(raw[:, 0], (raw[:, 1] - raw[0, 1]),
                      label='engy',
                      color=self.tableau[0],
                      **self.pltkwargs)
        self.ax2.plot(raw[:, 0], (raw[:, 4] - raw[0, 4]),
                      label='stress',
                      color=self.tableau[3],
                      **self.pltkwargs)
        self.fig.savefig("iten.png", **self.figsave)
        return

    def loop_prep_restart(self, opt='va'):
        raw = np.mat(np.loadtxt("iten.txt"))
        for i in range(len(raw)):
            dirname = "dir-{:03d}".format(i)
            self.mymkdir(dirname)
            np.savetxt("restart.txt", raw[i])
            if opt in ['va', 'vasp']:
                self.copy_inputs(dirname, 'KPOINTS',
                                 'INCAR', 'POTCAR',
                                 'restart.txt')
            elif opt in ['qe']:
                os.system("mv restart.txt {}".format(dirname))
                os.system('cp $POTDIR/{}  {}'.format(self.pot['file'],
                                                     dirname))
            self.set_pbs(dirname, raw[i][0])
        return

    def recordstrain(self, delta, x, fval):
        fid = open("s{:4.3f}.txt".format(delta), "a")
        formatstr = '{:6.5f} ' * (len(x) + len(fval))
        formatstr += '\n'
        fid.write(formatstr.format(x[0], x[1], *fval))
        fid.close()
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
    drv = cal_bcc_ideal_tensile_tp()
    if options.mtype.lower() == 'ilmp':
        drv.loop_tensile_lmp()

    if options.mtype.lower() == 'ivasp':
        drv.vasp_relax()

    if options.mtype.lower() == 'clcvasp':
        drv.loop_collect_vasp()

    if options.mtype.lower() == 'plt':
        drv.plt_energy_stress()

    if options.mtype.lower() in ['qe_restart', 'va_restart']:
        opt = options.mtype.lower().split('_')[0]
        drv.loop_prep_restart(opt)

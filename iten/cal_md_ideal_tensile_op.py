#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: yang37
# @Date:   2017-06-12 17:03:43
# @Last Modified by:   chaomy
# @Last Modified time: 2018-06-25 17:53:51


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
from md_pot_data import fluxdirs
from scipy.optimize import minimize
from optparse import OptionParser

sqrt2 = 2**(0.5)


class cal_bcc_ideal_tensile_op(get_data.get_data,
                               gn_pbs.gn_pbs,
                               gn_config.gnStructure,
                               plt_drv.plt_drv):

    def __init__(self):
        self.pot = self.load_data('../BASICS/pot.dat')
        gn_pbs.gn_pbs.__init__(self)
        plt_drv.plt_drv.__init__(self)
        get_data.get_data.__init__(self)
        gn_config.gnStructure.__init__(self, self.pot)
        self.alat = self.pot['lattice']
        self.npts = 13
        self.range = (0, 13)
        self.delta = 0.02
        self.basis = np.mat(np.identity(3))
        self.stress = None

    def loop_collect_vasp(self):
        npts = 30
        data = np.ndarray([npts, 10])
        lat = 3.322404
        for i in range(npts):
            mdir = "dir-{:4.3f}".format(0.01 * i)
            print(mdir)
            os.chdir(mdir)
            (engy, stress, vol) = self.vasp_energy_stress_vol()
            atoms = ase.io.read('CONTCAR', format='vasp')
            os.chdir(os.pardir)
            cell = atoms.get_cell() / lat
            data[i, 0], data[i, 1], data[i, 2], data[i, 3] = \
                cell[0, 0] - 1.0, engy, cell[1, 1] / sqrt2, cell[2, 2] / sqrt2
            data[i, 4:] = stress.transpose()
        np.savetxt("iten.txt", data)

    def load_input_params(self):
        if os.path.isfile('restart.txt'):
            data = np.loadtxt("restart.txt")
            delta = data[0]
            x0 = data[2:4]
            print(delta)
        else:
            data = np.loadtxt("strain.txt")
            delta = data
            x0 = np.array([1., 1.])
        return (delta, x0)

    def loop_tensile_lmp(self):
        npts = self.npts
        data = np.ndarray([npts, 10])
        x0 = np.array([0.80, 1.20])
        npts = self.range[1] - self.range[0]
        data = np.ndarray([npts, 10])
        for i in range(self.range[0], self.range[1]):
            delta = self.delta * (self.range[0] + i)
            res = minimize(self.runlmp, x0, delta,
                           tol=1e-4, method='Nelder-Mead')
            data[i][2], data[i][3] = res.x[0], res.x[1]
            print(res)
            data[i][0] = delta
            data[i][1] = res.fun
            data[i][4:] = self.stress
            print(self.stress)
        np.savetxt("iten.txt", data, fmt='%6.5f')
        np.savetxt("iten.md.op.txt", data, fmt='%6.5f')

    def runvasp(self, x, delta):
        basis = self.basis
        strain = np.mat([[1.0 + delta, 0.0, 0.0],
                         [0.0, x[0], 0.0],
                         [0.0, 0.0, x[1]]])
        new_strain = basis.transpose() * strain * basis
        self.gn_op_convention_lmp(new_strain, 'vasp')
        os.system("mpirun vasp > vasp.log")
        (engy, stress, vol) = self.vasp_energy_stress_vol()
        print(engy)
        return engy

    def runlmp(self, x, delta):
        basis = self.basis
        strain = np.mat([[1.0 + delta, 0.0, 0.0],
                         [0.0, x[0], 0.0],
                         [0.0, 0.0, x[1]]])
        new_strain = basis.transpose() * strain * basis
        self.gn_op_convention_lmp(new_strain, 'lmp')
        os.system("lmp_mpi -i in.init -screen  no")
        raw = np.loadtxt("out.txt")
        engy = raw[0]
        self.stress = raw[1:]
        return engy

    def gn_op_convention_lmp(self, strain=np.mat(np.identity(3)), tag='lmp'):
        alat = self.alat
        bas = np.mat(np.mat([[1., 0., 0.],
                             [0, np.sqrt(2), 0.],
                             [0., 0., np.sqrt(2)]]))
        va_bas = bas * strain
        cell = alat * va_bas

        atoms = ase.lattice.cubic.FaceCenteredCubic(directions=[[1., 0., 0.],
                                                                [0., 1., 0.],
                                                                [0., 0., 1.]],
                                                    latticeconstant=self.alat,
                                                    size=(1, 1, 1),
                                                    symbol=self.pot['element'],
                                                    pbc=(1, 1, 1))
        pos = np.mat([[0.0, 0.0, 0.0], [0.5, 0.5, 0.0],
                      [0.5, 0.0, 0.5], [0.0, 0.5, 0.5]]) * cell
        atoms.set_cell(cell)
        atoms.set_positions(pos)
        if tag == 'vasp':
            ase.io.write("POSCAR", images=atoms, format='vasp')
        if tag == 'lmp':
            self.write_lmp_config_data(atoms, 'init.txt')

    def gn_convention_lmps(self,
                           strain=np.mat(np.identity(3)),
                           tag='lmp'):
        alat = self.alat
        bas = np.mat(np.identity(3))
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
        (delta, x0) = self.load_input_params()
        data = np.zeros(10)
        res = minimize(self.runvasp, x0, delta, method='Nelder-Mead')
        data[2], data[3] = res.x[0], res.x[1]
        data[0] = delta
        data[1] = res.fun
        data[-6:] = self.stress
        print(res)
        np.savetxt("iten.txt", data)

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
            self.set_pbs(dirname, raw[i][0])

    def set_pbs(self, dirname, delta):
        self.set_nnodes(1)
        self.set_ppn(12)
        self.set_job_title("iva_{}".format(dirname))
        self.set_wall_time(50)
        self.set_main_job("""
        ../cal_md_ideal_tensile_op.py  -t ivasp
                          """.format(opt))
        self.write_pbs(od=True)
        os.system("mv va.pbs %s" % (dirname))

    def trans(self):
        for i in range(30):
            mdir = "dir-{:4.3f}".format(0.01 * i)
            self.mymkdir(mdir)
            fdir = fluxdirs['VA'] + \
                'Nb/Tensile/OneTensileOpath100/{}'.format(mdir)
            os.system('scp {}/OUTCAR {}'.format(fdir, mdir))
            os.system('scp {}/CONTCAR {}'.format(fdir, mdir))


if __name__ == '__main__':
    usage = "usage:%prog [options] arg1 [options] arg2"
    parser = OptionParser(usage=usage)
    parser.add_option('-t', "--mtype", action="store",
                      type="string", dest="mtype")
    parser.add_option('-p', "--param", action="store",
                      type='float', dest="fargs")

    (options, args) = parser.parse_args()
    drv = cal_bcc_ideal_tensile_op()
    dispatcher = {'ilmp': drv.loop_tensile_lmp,
                  'ivasp': drv.vasp_relax,
                  'clcva': drv.loop_collect_vasp,
                  'restart': drv.loop_prep_restart,
                  'trans': drv.trans}

    if options.fargs is not None:
        dispatcher[options.mtype.lower()](options.fargs)
    else:
        dispatcher[options.mtype.lower()]()

#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: yang37
# @Date:   2017-06-12 17:03:43
# @Last Modified by:   chaomy
# @Last Modified time: 2018-05-02 13:59:18


import os
import ase
import ase.io
import ase.lattice
import gn_pbs
import numpy as np
import gn_config
import get_data
from scipy.optimize import minimize
from optparse import OptionParser
from md_pot_data import fluxdirs
import cal_md_ideal_tensile_pre
import md_pot_data
from glob import glob


class cal_bcc_ideal_tensile_tp(get_data.get_data,
                               cal_md_ideal_tensile_pre.iten_pre,
                               gn_pbs.gn_pbs,
                               gn_config.gnStructure):

    def __init__(self, pot=md_pot_data.va_pot.Nb_pbe):
        self.pot = self.load_data('../BASICS/pot.dat')
        gn_pbs.gn_pbs.__init__(self)
        get_data.get_data.__init__(self)
        gn_config.gnStructure.__init__(self, self.pot)
        cal_md_ideal_tensile_pre.iten_pre.__init__(self)
        self.alat = self.pot['lattice']
        self.range = (0, 13)
        self.npts = self.range[1] - self.range[0]
        self.delta = 0.02
        self.basis = np.mat(np.identity(3))

    def loop_tensile_lmp(self):
        x0 = 0.91
        npts = self.range[1] - self.range[0]
        data = np.ndarray([npts, 10])  # delta + engy + yy, zz, + stress
        # inpath = "/Users/yangchaoming/src/LMPinfiles/in.init"
        # os.system("cp {} .".format(inpath))
        for i in range(npts):
            delta = self.delta * (self.range[0] + i)
            res = minimize(self.runlmp, x0, delta,
                           tol=5e-5, method='Nelder-Mead')
            x0 = res.x
            print(res)
            data[i][0] = delta
            data[i][1] = res.fun
            data[i][2:2 + 2] = res.x
            data[i][-6:] = self.stress
        np.savetxt("iten.txt", data, fmt='%6.5f')
        np.savetxt("iten.md.tp.txt", data, fmt='%6.5f')

    def mesh(self):
        cn = 0
        basis = self.basis
        op = 'lmp'
        nps = 20
        raw = np.ndarray([nps * nps, 9])
        for i in range(nps):
            for j in range(nps):
                strain = np.mat([[1 + 0.01 * j, 0.0, 0.0],
                                 [0.0, 1 - i * 0.001, 0.0],
                                 [0.0, 0.0, 1 - i * 0.001 * 2]])
                mdir = "dir_{:04d}".format(cn)
                new_strain = basis.transpose() * strain * basis
                if op in ['vasp']:
                    self.gn_primitive_lmps(new_strain, 'vasp')
                    self.mymkdir(mdir)
                    os.system("cp va.pbs {}".format(mdir))
                    os.system("cp POTCAR {}".format(mdir))
                    os.system("cp POSCAR {}".format(mdir))
                    os.system("cp KPOINTS {}".format(mdir))
                    os.system("cp INCAR {}".format(mdir))
                elif op in ['lmp']:
                    self.gn_convention_lmps(new_strain, 'lmp')
                    os.system("lmp_mpi -i in.init -screen no")
                    raw[cn][:2] = 0.01 * j, 0.001 * i
                    raw[cn][2:] = np.loadtxt("out.txt")
                cn += 1
        np.savetxt("save.txt", raw, fmt='%10.6f')

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
        print(engy, self.stress)
        return engy

    def runlmp(self, x, delta):
        basis = self.basis
        strain = np.mat([[1.0 + delta, 0.0, 0.0],
                         [0.0, x, 0.0],
                         [0.0, 0.0, x]])
        new_strain = basis.transpose() * strain * basis
        self.gn_convention_lmps(new_strain, 'lmp')
        os.system("lmp_mpi -i in.init -screen  no")
        raw = np.loadtxt("out.txt")
        engy = raw[0]
        self.stress = raw[1:]
        return engy

    def gn_convention_lmps(self, strain=np.mat(np.identity(3)), tag='lmp'):
        alat = self.alat
        bas = np.mat(np.identity(3))
        # very important (vasp add strain is basis right time strain)
        va_bas = bas * strain
        cell = alat * va_bas
        atoms = self.set_bcc_convention()
        pos = np.mat([[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]) * cell
        atoms.set_cell(cell)
        atoms.set_positions(pos)
        if tag == 'vasp':
            ase.io.write("POSCAR", images=atoms, format='vasp')
        if tag == 'lmp':
            self.write_lmp_config_data(atoms, 'init.txt')

    def vasp_relax(self, given=True):
        data = np.loadtxt("strain.txt")
        if given is True:
            delta = data[0]
            x0 = data[-2:]
            print(delta)
            print(x0)
        else:
            delta = data
            x0 = np.array([1., 1.])
        data = np.zeros(4)
        res = minimize(self.runvasp, x0, delta, method='Nelder-Mead',
                       options={'fatol': 5e-4, 'disp': True})
        print(res)
        data[0] = delta
        data[1] = res.fun
        data[2:] = res.x
        np.savetxt("iten.txt", data)

    def set_pbs(self, mdir, delta):
        self.set_nnodes(1)
        self.set_ppn(12)
        self.set_job_title("%s" % (mdir))
        self.set_wall_time(50)
        self.set_main_job("""
        python ~/My_cal/Mechnical_cal/cal_md_ideal_tensile.py  -t  ivasp
                          """)
        self.write_pbs(od=False)
        os.system("mv va.pbs %s" % (mdir))

    def loop_prep_restart(self, opt='va'):
        raw = np.mat(np.loadtxt("iten.txt"))
        for i in range(len(raw)):
            mdir = "dir-{:03d}".format(i)
            self.mymkdir(mdir)
            np.savetxt("restart.txt", raw[i])
            if opt in ['va', 'vasp']:
                self.copy_inputs(mdir, 'KPOINTS',
                                 'INCAR', 'POTCAR',
                                 'restart.txt')
            elif opt in ['qe']:
                os.system("mv restart.txt {}".format(mdir))
                os.system('cp $POTDIR/{}  {}'.format(self.pot['file'], mdir))
            self.set_pbs(mdir, raw[i][0])

    def recordstrain(self, delta, x, fval):
        fid = open("s{:4.3f}.txt".format(delta), "a")
        formatstr = '{:6.5f} ' * (1 + 2 + 6)
        formatstr += '\n'
        fid.write(formatstr.format(fval, x[0], x[1], *self.stress))
        fid.close()

    def loop_collect_vasp(self):
        npts = 40
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
                cell[0, 0] - 1.0, engy, cell[1, 1], cell[2, 2]
            data[i, 4:] = stress.transpose()
        np.savetxt("iten.txt", data)

    def loop_clc_qe(self):
        flist = glob("WRe.out.*")
        data = np.ndarray([len(flist), 7])
        for i in range(len(flist)):
            file = flist[i]
            strain = float(file[-5:])
            (engy, vol, stress) = self.qe_get_energy_stress(file)
            print(stress)
            data[i, 0] = strain
            data[i, 1] = engy / 4.
            data[i, 4], data[i, 5], data[i, 6] = \
                stress[0, 0], stress[1, 1], stress[2, 2]
        np.savetxt('iten.txt', data)


if __name__ == '__main__':
    usage = "usage:%prog [options] arg1 [options] arg2"
    parser = OptionParser(usage=usage)
    parser.add_option('-t', "--mtype", action="store",
                      type="string", dest="mtype")
    parser.add_option('-p', "--param", action="store",
                      type='float', dest="fargs")

    (options, args) = parser.parse_args()
    drv = cal_bcc_ideal_tensile_tp()
    dispatcher = {'ilmp': drv.loop_tensile_lmp,
                  'ivasp': drv.vasp_relax,
                  'clcqe': drv.loop_clc_qe,
                  'clcva': drv.loop_collect_vasp,
                  'restart': drv.loop_prep_restart,
                  'mesh': drv.mesh}

    if options.fargs is not None:
        dispatcher[options.mtype.lower()](options.fargs)
    else:
        dispatcher[options.mtype.lower()]()

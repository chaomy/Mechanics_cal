#!/usr/bin/env python
# encoding: utf-8

###################################################################
#
# File Name : ./cal_bcc_ideal_shear.py
#
###################################################################
#
# Purpose :
#
# Creation Date :
# Last Modified : Sat Apr 22 21:22:17 2017
# Created By    : Chaoming Yang
#
###################################################################

from md_pot_data import unitconv
from optparse import OptionParser
from scipy.optimize import minimize
from scipy.interpolate import InterpolatedUnivariateSpline
import cal_md_ideal_tensile_plt
import os
import ase
import ase.io
import ase.lattice
import ase.lattice.cubic as cubic
import numpy as np
import gn_config
import gn_pbs
import get_data
import plt_drv
import md_pot_data
import gn_qe_inputs


class cal_bcc_ideal_shear(get_data.get_data,
                          gn_config.bcc,
                          gn_pbs.gn_pbs,
                          plt_drv.plt_drv):

    def __init__(self,
                 shtype='110'):
        # self.pot = self.load_data('../pot.dat')
        self.pot = md_pot_data.qe_pot.vca_W50Re50
        gn_pbs.gn_pbs.__init__(self)
        plt_drv.plt_drv.__init__(self)
        get_data.get_data.__init__(self)
        gn_config.bcc.__init__(self, self.pot)
        self.alat = self.pot['lattice']
        self.npts = 20
        self.delta = 0.02
        shd111p211 = {'e1': np.array([1., 1., 1.]),
                      'e2': np.array([1., 1., -2.]),
                      'e3': np.array([-1., 1., 0.])}

        shd111p110 = {'e1': np.array([1., 1., 1.]),
                      'e2': np.array([1., -1, 0]),
                      'e3': np.array([1, 1., -2])}

        if shtype == '211':
            e1 = shd111p211['e1']
            e2 = shd111p211['e2']
            e3 = shd111p211['e3']
        elif shtype == '110':
            e1 = shd111p110['e1']
            e2 = shd111p110['e2']
            e3 = shd111p110['e3']

        e1 = e1 / np.linalg.norm(e1)
        e2 = e2 / np.linalg.norm(e2)
        e3 = e3 / np.linalg.norm(e3)

        self.basis = np.mat([e1, e2, e3])
        self.va_prim = np.mat([[-0.5, 0.5, 0.5],
                               [0.5, -0.5, 0.5],
                               [0.5, 0.5, -0.5]])
        self.lm_prim = self.lmp_change_box(self.va_prim)
        self.qedrv = gn_qe_inputs.gn_qe_infile(self.pot)
        # set qe simulation setup
        self.qedrv.set_degauss('0.03D0')
        self.qedrv.set_ecut('45')
        self.qedrv.set_kpnts((33, 33, 33))
        self.root = os.getcwd()
        return

    def gn_convention(self,
                      strain=np.mat([[1., 0., 0.],
                                     [0., 1., 0.],
                                     [0., 0., 1.]])):
        alat = self.alat
        atoms = cubic.BodyCenteredCubic(directions=[[1, 0, 0],
                                                    [0, 1, 0],
                                                    [0, 0, 1]],
                                        latticeconstant=alat,
                                        size=(1, 1, 1),
                                        symbol=self.pot['element'],
                                        pbc=(1, 1, 1))

        ase.io.write("POSCAR_perf", images=atoms, format='vasp')
        # add strain
        atoms.set_cell((strain * atoms.get_cell()))
        pos = np.mat(atoms.get_positions())
        pos = pos * strain
        atoms.set_positions(pos)
        ase.io.write("POSCAR_c", images=atoms, format='vasp')
        self.write_lmp_config_data(atoms, "con.txt")
        return

    def gn_primitive_lmps(self,
                          strain=np.mat(np.identity(3)),
                          tag='lmp'):

        alat = self.alat
        bas = np.mat([[-0.5, 0.5, 0.5],
                      [0.5, -0.5, 0.5],
                      [0.5, 0.5, -0.5]])

        if tag == 'vasp':
            va_bas = bas * strain
            cell = alat * va_bas
            atoms = ase.Atoms(self.pot['element'],
                              positions=[[0, 0, 0]],
                              cell=cell,
                              pbc=[1, 1, 1])
            ase.io.write("POSCAR", images=atoms, format='vasp')

        if tag == 'qe':
            qe_bas = bas * strain
            cell = alat * qe_bas
            atoms = ase.Atoms(self.pot['element'],
                              positions=[[0, 0, 0]],
                              cell=cell,
                              pbc=[1, 1, 1])
            self.qedrv.gn_infile_dipole_ideal_shear(atoms)

        if tag == 'lmp':
            lmp_bas = bas * strain
            lmp_bas = self.lmp_change_box(lmp_bas)
            cell = alat * lmp_bas
            atoms = ase.Atoms(self.pot['element'],
                              positions=[[0., 0., 0.]],
                              cell=cell,
                              pbc=[1, 1, 1])
            self.write_lmp_config_data(atoms, 'init.txt')
        return

    def shear_twin_path(self):
        e1 = 0.5 * np.array([-1, 1, 1])
        e2 = 0.5 * np.array([1, -1, 1])
        e3 = 0.5 * np.array([1, 1, -1])

        et = np.array([-1, -1, 1])
        npts = 10
        delta = 1. / npts
        oneo6 = 1. / 6.
        for i in range(npts + 1):
            xtos = i * delta
            a1 = e1 + oneo6 * xtos * et
            a2 = e2 + oneo6 * xtos * et
            a3 = e3
            cell = np.mat([a1, a2, a3])
            print cell
            atoms = ase.Atoms('Nb',
                              positions=[[0, 0, 0]],
                              cell=cell,
                              pbc=[1, 1, 1])
            ase.io.write("POSCAR_{:03d}".format(i),
                         images=atoms,
                         format='vasp')
        return

    def set_pbs(self, dirname, delta, opt='vasp'):
        self.set_nnodes(1)
        self.set_ppn(12)
        self.set_job_title("%s" % (dirname))
        self.set_wall_time(70)
        self.set_main_job("""../cal_md_ideal_shear.py  -t  i{}
                          """.format(opt))
        self.write_pbs(od=True)
        os.system("mv va.pbs %s" % (dirname))
        return

    def loop_prep_vasp(self):
        npts = self.npts
        for i in range(npts):
            delta = self.delta * i
            dirname = "dir-{:03d}".format(i)
            self.mymkdir(dirname)
            os.system("echo {} > strain.txt".format(delta))
            self.copy_inputs(dirname, 'KPOINTS',
                             'INCAR', 'POTCAR', 'strain.txt')
            self.set_pbs(dirname, delta)
        return

    def loop_prep_restart(self, opt='va'):
        raw = np.mat(np.loadtxt("ishear.txt"))
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
            self.set_pbs(dirname, raw[i][0], opt)
        return

    def loop_sub(self):
        npts = self.npts
        for i in range(npts):
            dirname = "dir-{:03d}".format(i)
            os.chdir(dirname)
            os.system("qsub va.pbs")
            os.chdir(self.root)
        return

    def load_input_params(self):
        if os.path.isfile('restart.txt'):
            data = np.loadtxt("restart.txt")
            delta = data[0]
            x0 = data[-5:]
            print delta
            print x0
        else:
            data = np.loadtxt("strain.txt")
            delta = data
            x0 = np.array([1., 1., 1., 0.0, 0.0])
        return (delta, x0)

    def qe_relax(self):
        (delta, x0) = self.load_input_params()
        data = np.zeros(7)
        res = minimize(self.runqe, x0, delta,
                       method='Nelder-Mead',
                       options={'fatol': 2e-3, 'disp': True})

        print res
        data[0] = delta
        data[1] = res.fun
        data[2:] = res.x
        np.savetxt("ishear.txt", data)
        return

    def vasp_relax(self):
        (delta, x0) = self.load_input_params()
        data = np.zeros(7)
        res = minimize(self.runvasp, x0, delta,
                       method='Nelder-Mead',
                       options={'fatol': 2e-3, 'disp': True})

        print res
        data[0] = delta
        data[1] = res.fun
        data[2:] = res.x
        np.savetxt("ishear.txt", data)
        return

    def loop_shear_lmp(self):
        x0 = np.array([1., 1., 1., 0., 0.])
        npts = self.npts
        data = np.ndarray([npts, 7])
        for i in range(npts):
            delta = self.delta * i
            res = minimize(self.runlmp, x0, delta,
                           method='Nelder-Mead',
                           options={'xtol': 1e-3, 'disp': True})
            x0 = res.x
            print res
            data[i][0] = (delta)
            data[i][1] = res.fun
            data[i][2:] = res.x
        np.savetxt("ishear.txt", data)
        return

    def recordstrain(self, delta, x, fval):
        fid = open("s{:4.3f}.txt".format(delta), "a")
        fid.write('{} {} {} {} {} {}\n'.format(x[0], x[1], x[2],
                                               x[3], x[4], fval))
        fid.close()
        return

    def runlmp(self, x, delta):
        basis = self.basis
        # y shear toward x direction
        strain = np.mat([[x[0], 0.0, 0.0],
                         [-delta, x[1], 0.0],
                         [x[3], x[4], x[2]]])
        new_strain = basis.transpose() * strain * basis
        self.gn_primitive_lmps(new_strain, 'lmp')
        os.system("lmp_mpi -i in.init -screen  no")
        raw = np.loadtxt("out.txt")
        engy = raw[0]
        self.recordstrain(delta, x, engy)
        return engy

    def runqe(self, x, delta):
        basis = self.basis
        strain = np.mat([[x[0], 0.0, 0.0],
                         [-delta, x[1], 0.0],
                         [x[3], x[4], x[2]]])
        new_strain = basis.transpose() * strain * basis
        self.gn_primitive_lmps(new_strain, 'qe')
        os.system("mpirun pw.x < qe.in > qe.out")
        (engy, vol, stress) = self.qe_get_energy_stress('qe.out')
        self.recordstrain(delta, x, engy)
        return engy

    def runvasp(self, x, delta):
        basis = self.basis
        strain = np.mat([[x[0], 0.0, 0.0],
                         [-delta, x[1], 0.0],
                         [x[3], x[4], x[2]]])
        new_strain = basis.transpose() * strain * basis
        self.gn_primitive_lmps(new_strain, 'vasp')
        os.system("mpirun vasp > vasp.log")
        (engy, stress, vol) = self.vasp_energy_stress_vol()
        self.recordstrain(delta, x, engy)
        return engy

    def qe_loop_stress(self, opt='clc'):
        npts = self.npts
        convunit = unitconv.ustress['evA3toGpa']
        # conveng = unitconv.uengy['rytoeV']
        data = np.ndarray([npts, 3])
        if opt == 'clc':
            for i in range(npts):
                dirname = "dir-{:03d}".format(i)
                print dirname
                if os.path.isdir(dirname):
                    os.chdir(dirname)
                    (engy, vol, stress) = self.qe_get_energy_stress('qe.out')
                    raw = np.loadtxt("ishear.txt")
                    os.chdir(self.root)
                    vol = vol * (unitconv.ulength['BohrtoA']**3)
                    data[i, 0] = raw[0]
                    data[i, 1] = raw[1]
                    data[i, 2] = vol
            np.savetxt('stress.txt', data)
        elif opt is 'stress':
            data = np.loadtxt('stress.txt')
            spl = InterpolatedUnivariateSpline(data[:, 0], data[:, 1], k=3)
            splder1 = spl.derivative()
            for i in range(len(data)):
                data[i, -1] = splder1(data[i, 0]) * convunit / data[i, 2]
            print data
        return

    def vasp_loop_stress(self):
        npts = self.npts
        data = np.ndarray([npts, 4])
        convunit = unitconv.ustress['evA3toGpa']
        for i in range(npts):
            dirname = "dir-{:03d}".format(i)
            os.chdir(dirname)
            atoms = ase.io.read("CONTCAR", format='vasp')
            raw = np.loadtxt("ishear.txt")
            os.chdir(self.root)
            data[i, 0] = raw[0]
            data[i, 1] = raw[1]
            data[i, 2] = np.linalg.det(atoms.get_cell())

        spl = InterpolatedUnivariateSpline(data[:, 0], data[:, 1], k=3)
        splder1 = spl.derivative()
        for i in range(len(data)):
            # append the stress to the last column
            data[i, -1] = splder1(data[i, 0]) * convunit / data[i, 2]
        print data
        np.savetxt("stress.txt", data)
        return

    ##########################################################
    # calculate derivative of energy to get stress
    ##########################################################
    def convert_stress(self):
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
                raw[i, 2], raw[i, 3], raw[i, 4]
            strmat[1, 0], strmat[2, 0], strmat[2, 1] = \
                raw[i, 0], raw[i, 5], raw[i, 6]
            strmat = np.mat(strmat)
            vol[i] = vperf * np.linalg.det(strmat)
        tag = 'interp'
        if tag == 'interp':
            # interpolate
            spl = InterpolatedUnivariateSpline(raw[:, 0], raw[:, 1], k=1)
            splder1 = spl.derivative()

            for i in range(len(raw)):
                # append the stress to the last column
                print "coeff", convunit / vol[i]
                data[i, -1] = splder1(raw[i, 0]) * convunit / vol[i]
        print data
        np.savetxt("stress.txt", data)
        return

    def trans_stress_to_cartesian(self, stssvect, opt='vasp'):
        basis = self.basis
        stssmtx = np.mat(np.zeros([3, 3]))

        stssmtx[0, 0] = stssvect[0]
        stssmtx[1, 1] = stssvect[1]
        stssmtx[2, 2] = stssvect[2]

        stssmtx[1, 0] = stssvect[3]
        stssmtx[2, 0] = stssvect[4]
        stssmtx[2, 1] = stssvect[5]

        stssmtx[0, 1] = stssvect[3]
        stssmtx[0, 2] = stssvect[4]
        stssmtx[1, 2] = stssvect[5]

        stssmtx = basis * stssmtx * basis.transpose()

        stssvect[0] = stssmtx[0, 0]
        stssvect[1] = stssmtx[1, 1]
        stssvect[2] = stssmtx[2, 2]

        stssvect[3] = stssmtx[1, 0]
        stssvect[4] = stssmtx[2, 0]
        stssvect[5] = stssmtx[2, 1]
        return stssvect

    def clc_data(self):
        npts = self.npts
        data = np.ndarray([npts, 7])
        for i in range(npts):
            dirname = "dir-{:03d}".format(i)
            if os.path.isdir(dirname):
                raw = np.loadtxt("{}/ishear.txt".format(dirname))
            data[i, :] = raw
        np.savetxt('ishear.txt', data)
        return

    def read_ofiles(self, opt='clccell'):
        import glob
        flist = glob.glob('dir-*')
        if opt == 'clcengy':
            data = np.ndarray([2, len(flist)])
            for i in range(len(flist)):
                file = flist[i]
                cnt = int(file[4:7])
                fid = open(file, 'r')
                raw = fid.readlines()
                fid.close()
                for line in raw[:-1]:
                    if len(line.split()) == 1:
                        dat = line.split()[0]
                if dat[0] == '[':
                    dat = dat.split('\'')[1]
                data[0, i] = 0.02 * cnt
                data[1, i] = dat
            print data
            np.savetxt('ishear.txt', data)
        elif opt == 'clccell':
            for i in range(len(flist)):
                mdir = flist[i]
                print self.qe_get_cell('{}/qe.in'.format(mdir))
        return


if __name__ == '__main__':
    usage = "usage:%prog [options] arg1 [options] arg2"
    parser = OptionParser(usage=usage)

    parser.add_option('-t', "--mtype", action="store",
                      type="string", dest="mtype", help="",
                      default="prp_r")

    parser.add_option('-c', "--delta", action="store",
                      type='float', dest="delta",
                      default=0.02)

    (options, args) = parser.parse_args()

    drv = cal_bcc_ideal_shear()
    pltdrv = cal_md_ideal_tensile_plt.cal_md_ideal_tensile_plt()

    if options.mtype.lower() in ['run_vasp', 'run_lmp']:
        opt = options.mtype.lower().split('_')[-1]
        drv.md_ideal_shear('run', opt)

    if options.mtype.lower() in ['clc_vasp', 'clc_lmp']:
        opt = options.mtype.lower().split('_')[-1]
        drv.md_ideal_shear('clc', opt)

    if options.mtype.lower() == 'clcqe' or \
            options.mtype.lower() == 'qeclc':
        drv.clc_data()

    if options.mtype.lower() == 'vastress':
        drv.vasp_loop_stress()

    if options.mtype.lower() == 'qestress':
        drv.qe_loop_stress(opt='stress')

    if options.mtype.lower() == 'lmpstress':
        drv.convert_stress()

    if options.mtype.lower() == 'trans':
        drv.trans_stress()

    if options.mtype.lower() == 'twin':
        drv.shear_twin_path()

    if options.mtype.lower() == 'plt_engy':
        pltdrv.plt_strain_vs_energy()

    if options.mtype.lower() == 'cmpplt':
        drv.cmp_plt()

    if options.mtype.lower() == 'cmp':
        drv.plt_energy_stress_cmp()

    if options.mtype.lower() == 'vaspprep':
        drv.loop_prep_vasp()

    if options.mtype.lower() == 'qeprep':
        drv.loop_prep_qe()

    if options.mtype.lower() in ['ivasp', 'iva']:
        drv.vasp_relax()

    if options.mtype.lower() == 'ilmp':
        drv.loop_shear_lmp()

    if options.mtype.lower() == 'iqe':
        print drv.qe_relax()

    if options.mtype.lower() == 'sub':
        drv.loop_sub()

    if options.mtype.lower() == 'tmp':
        drv.read_ofiles()

    if options.mtype.lower() == 'gnqe':
        drv.gn_primitive_lmps(tag='qe')

    if options.mtype.lower() in ['qe_restart', 'va_restart']:
        opt = options.mtype.lower().split('_')[0]
        drv.loop_prep_restart(opt)

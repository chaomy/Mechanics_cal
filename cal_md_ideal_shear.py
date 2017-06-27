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
import matplotlib.pylab as plt
import md_pot_data
import gn_qe_inputs
from optparse import OptionParser
from scipy.optimize import minimize
from scipy.interpolate import InterpolatedUnivariateSpline
from md_pot_data import unitconv


class cal_bcc_ideal_shear(get_data.get_data,
                          gn_pbs.gn_pbs,
                          plt_drv.plt_drv):

    def __init__(self):
        self.pot = md_pot_data.qe_pot.vca_W75Re25
        gn_pbs.gn_pbs.__init__(self)
        plt_drv.plt_drv.__init__(self)
        self.alat = self.pot['lattice']
        self.npts = 20
        self.delta = 0.02
        shd111p211 = {'e1': np.array([1., 1., 1.]),
                      'e2': np.array([1., 1., -2.]),
                      'e3': np.array([-1., 1., 0.])}

        shd111p110 = {'e1': np.array([1., 1., 1.]),
                      'e2': np.array([1., -1, 0]),
                      'e3': np.array([1, 1., -2])}
        shtype = '110'
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
        get_data.get_data.__init__(self)
        self.va_prim = np.mat([[-0.5, 0.5, 0.5],
                               [0.5, -0.5, 0.5],
                               [0.5, 0.5, -0.5]])

        self.configdrv = gn_config.bcc(self.pot)
        self.lm_prim = self.configdrv.lmp_change_box(self.va_prim)
        self.qedrv = gn_qe_inputs.gn_qe_infile(self.pot)
        self.qedrv.set_degauss('0.03D0')
        self.qedrv.set_ecut('45')
        self.qedrv.set_kpnts('35')
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
        self.configdrv.write_lmp_config_data(atoms, "con.txt")
        return

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

    def read_stress(self):
        (engy, stress, vol) = self.vasp_energy_stress_vol()
        stress = stress * 0.1
        print "engy {}".format(engy)
        print "{}, {}, {}, {}, {}, {}".format(
            stress[0][0], stress[1][0], stress[2][0],
            stress[3][0], stress[4][0], stress[5][0])
        return

    def lmp_relax(self):
        delta = 0.4
        x0 = np.array([1., 1., 1., 0., 0.0])
        res = minimize(self.runlmp, x0, delta,
                       method='Nelder-Mead',
                       options={'xtol': 1e-3, 'disp': True})
        print res.fun
        print res.x
        return

    def set_pbs(self, dirname, delta, opt='vasp'):
        self.set_nnodes(1)
        self.set_ppn(12)
        self.set_job_title("%s" % (dirname))
        self.set_wall_time(70)
        self.set_main_job("""cal_md_ideal_shear.py  -t  i{}
                          """.format(opt))
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
            self.copy_inputs(dirname, 'KPOINTS',
                             'INCAR', 'POTCAR', 'strain.txt')
            self.set_pbs(dirname, delta)
        return

    def loop_prep_vasp_restart(self):
        data = np.loadtxt("ishear.txt")
        for i in range(len(data)):
            dirname = "dir-{:03d}".format(i)
            self.mymkdir(dirname)
            np.savetxt("restart.txt", data[i])
            self.copy_inputs(dirname, 'KPOINTS',
                             'INCAR', 'POTCAR', 'restart.txt')
            self.set_pbs(dirname, data[i][0])
        return

    def loop_prep_qe_restart(self):
        raw = np.mat(np.loadtxt("ishear.txt"))
        for i in range(len(raw)):
            dirname = "dir-{:03d}".format(i)
            self.mymkdir(dirname)
            np.savetxt("restart.txt", raw[i])
            os.system("mv restart.txt {}".format(dirname))
            os.system('cp $POTDIR/{}  {}'.format(self.pot['file'],
                                                 dirname))
            self.set_pbs(dirname, raw[i][0], opt='qe')
        return

    def loop_prep_qe(self):
        npts = self.npts
        for i in range(npts):
            delta = self.delta * i
            dirname = "dir-{:03d}".format(i)
            self.mymkdir(dirname)
            os.system("echo {} > strain.txt".format(delta))
            os.system("mv strain.txt {}".format(dirname))
            os.system('cp $POTDIR/{} {}'.format(self.pot['file'],
                                                dirname))
            self.set_pbs(dirname, delta, opt='qe')
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
                       options={'disp': True})
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
                       options={'xtol': 1e-3, 'disp': True})
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
                           options={'fatol': 2e-4, 'disp': True})
            x0 = res.x
            print res
            data[i][0] = (delta)
            data[i][1] = res.fun
            data[i][2:] = res.x
        np.savetxt("ishear.txt", data)
        return

    def runlmp(self, x, delta):
        basis = self.basis
        # y shear toward x direction
        # 211 -delta
        # 110  delta
        strain = np.mat([[x[0], 0.0, 0.0],
                         [-delta, x[1], 0.0],
                         [x[3], x[4], x[2]]])

        new_strain = basis.transpose() * strain * basis
        self.gn_primitive_lmps(new_strain, 'lmp')
        os.system("lmp_mpi -i in.init -screen  no")
        raw = np.loadtxt("out.txt")
        engy = raw[0]
        print engy
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
        print engy
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
        print engy
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

    def gn_infile_dipole_ideal_shear(self,
                                     atoms=None):
        with open('qe.in', 'w') as fid:
            fid = self.qe_write_control(fid, atoms)
            fid = self.qe_write_system(fid, atoms)
            fid = self.qe_write_electrons(fid)
            fid = self.qe_write_cell(fid, atoms.get_cell())
            fid = self.qe_write_species(fid, atoms, self.pot)
            fid = self.qe_write_pos(fid, atoms)
            fid = self.qe_write_kpts(fid, (33, 33, 33))
            fid.close()
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

    def cmp_plt(self):
        self.set_keys()
        self.set_211plt()
        self.plt_energy_stress(fname='stress_0.00.txt')
        self.plt_energy_stress(fname='stress_0.25.txt')
        self.fig.savefig("istress.png", **self.figsave)
        return

    def plt_energy_stress(self, fname='stress.txt', set=False):
        if set is True:
            self.set_keys()
            self.set_211plt()
        raw = np.loadtxt(fname)
        self.ax1.plot(raw[:, 0], (raw[:, 1] - raw[0, 1]),
                      label='engy',
                      **self.pltkwargs)
        self.ax2.plot(raw[:, 0], (raw[:, -1] - raw[0, -1]),
                      label='stress',
                      **self.pltkwargs)
        return

    def set_pltkargs(self):
        self.keyslist = [{'linestyle': '-.', 'color': 'b', 'linewidth': 2,
                          'marker': 'o'},
                         {'linestyle': '--', 'color': 'y', 'linewidth': 2,
                          'marker': '<'},
                         {'linestyle': '-.', 'color': 'g', 'linewidth': 2,
                          'marker': 'o'},
                         {'linestyle': '--', 'color': 'm', 'linewidth': 2,
                          'marker': '<'}]
        return

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

    def plt_energy_stress_cmp(self):
        potlist = ['adp', 'pbe']
        self.set_211plt(mfigsize=(8.5, 4.3), lim=True)
        self.set_keys()
        self.set_pltkargs()
        plt.rc('xtick', labelsize=self.mlabelsize)
        plt.rc('ytick', labelsize=self.mlabelsize)
        for i in range(len(potlist)):
            pot = potlist[i]
            fname = 'stress.txt.{}'.format(pot)
            raw = np.loadtxt(fname)
            self.ax1.plot(raw[:, 0],
                          (raw[:, 1] - raw[0, 1]),
                          label=pot,
                          **self.keyslist[i])
            self.ax1.legend(**self.legendarg)
            self.ax1.set_ylabel('energy [eV]', {'fontsize': self.myfontsize})

            self.ax2.plot(raw[:, 0],
                          (raw[:, -1] - raw[0, -1]),
                          label=pot,
                          **self.keyslist[i + 2])
            self.ax2.legend(**self.legendarg)
            self.ax2.set_ylabel('stress [Gpa]', {'fontsize': self.myfontsize})

        plt.xlabel('strain', {'fontsize': self.myfontsize})
        self.fig.savefig("stress_cmp.png", **self.figsave)
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

    if options.mtype.lower() == 'lmprelax':
        drv.lmp_relax()

    if options.mtype.lower() == 'plt':
        if not os.path.isfile("stress.txt"):
            drv.convert_stress()
        drv.plt_energy_stress()

    if options.mtype.lower() == 'cmpplt':
        drv.cmp_plt()

    if options.mtype.lower() == 'cmp':
        drv.plt_energy_stress_cmp()

    if options.mtype.lower() == 'restart':
        drv.loop_prep_vasp_restart()

    if options.mtype.lower() == 'qerestart':
        drv.loop_prep_qe_restart()

    if options.mtype.lower() == 'vaspprep':
        drv.loop_prep_vasp()

    if options.mtype.lower() == 'qeprep':
        drv.loop_prep_qe()

    if options.mtype.lower() == 'ivasp':
        drv.vasp_relax()

    if options.mtype.lower() == 'ilmp':
        drv.loop_shear_lmp()

    if options.mtype.lower() == 'iqe':
        print drv.qe_relax()

    if options.mtype.lower() == 'sub':
        drv.loop_sub()

    if options.mtype.lower() == 'gnqe':
        drv.gn_primitive_lmps(tag='qe')

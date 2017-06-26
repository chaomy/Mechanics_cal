#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: yangchaoming
# @Date:   2017-06-13 15:37:47
# @Last Modified by:   chaomy
# @Last Modified time: 2017-06-25 19:30:33

import os
import numpy as np
import gn_qe_inputs
import md_pot_data
import gn_config
import get_data
import gn_pbs
import output_data
import glob
import plt_drv
from optparse import OptionParser
from scipy.interpolate import InterpolatedUnivariateSpline


class cal_lattice(gn_config.bcc,
                  gn_config.fcc,
                  gn_config.hcp,
                  get_data.get_data,
                  output_data.output_data,
                  gn_qe_inputs.gn_qe_infile,
                  plt_drv.plt_drv,
                  gn_pbs.gn_pbs):

    def __init__(self, inpot=None):
        self.pot = md_pot_data.qe_pot.pbe_w
        output_data.output_data.__init__(self)
        get_data.get_data.__init__(self)
        gn_pbs.gn_pbs.__init__(self)
        gn_qe_inputs.gn_qe_infile.__init__(self, self.pot)
        plt_drv.plt_drv.__init__(self)
        self.alat0 = self.pot['lattice']
        self.element = self.pot['element']
        self.mass = self.pot['mass']
        self.kpnts = [42, 42, 42]
        self.root = os.getcwd()
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
            os.system('cp $POTDIR/{} {}'.format(self.pot['file'],
                                                dirname))
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

    def goandsub(self, dirname, rootdir):
        os.chdir(dirname)
        os.system('qsub va.pbs')
        os.chdir(rootdir)
        return

    def loop_degauss(self, opt='prep'):
        degauss0 = 0.02
        for i in range(7):
            degauss = degauss0 + 0.005 * i
            mdir = 'degauss{:4.3f}'.format(degauss)
            if opt == 'prep':
                self.set_degauss('{}D0'.format(degauss))
                self.mymkdir(mdir)
                self.loop_kpoints()
                self.set_pbs(mdir)
                os.system('mv dir-* {}'.format(mdir))
            elif opt == 'sub':
                self.goandsub(mdir, self.root)
            elif opt == 'clc':
                os.system("cp {}/kpts.txt kpts_{:4.3f}.txt".format(mdir,
                                                                   degauss))
        return

    def loop_kpoints(self):
        bcc_drv = gn_config.bcc(self.pot)
        bcc_drv.set_lattce_constant(self.alat0)
        self.set_ecut('{}'.format(48))
        self.set_disk_io('none')
        for kpts in range(32, 50):
            self.set_kpnts((kpts, kpts, kpts))
            dirname = 'dir-kpt-{}'.format(kpts)
            self.mymkdir(dirname)
            atoms = bcc_drv.set_bcc_primitive((1, 1, 1))
            self.gn_qe_bcc_lattice_infile(atoms)
            os.system('mv qe.in {}'.format(dirname))
            os.system('cp $POTDIR/{} {}'.format(self.pot['file'],
                                                dirname))
        return

    def loop_kpoints_ecut(self):
        for i in range(43, 50):
            kpts = i
            mdir = 'kpts_{:02}'.format(kpts)
            self.mymkdir(mdir)
            os.chdir(mdir)
            self.set_kpnts((i, i, i))
            self.loop_ecut()
            self.set_pbs(mdir)
            os.chdir(os.pardir())
        return

    def loop_ecut(self):
        bcc_drv = gn_config.bcc(self.pot)
        bcc_drv.set_lattce_constant(self.alat0)
        for ecut in range(28, 52):
            self.set_ecut('{}'.format(ecut))
            dirname = 'dir-ecut-{}'.format(ecut)
            self.mymkdir(dirname)
            atoms = bcc_drv.set_bcc_primitive((1, 1, 1))
            self.gn_qe_bcc_lattice_infile(atoms)
            os.system('mv qe.in {}'.format(dirname))
            os.system('cp $POTDIR/{} {}'.format(self.pot['file'],
                                                dirname))
        return

    def clc_data(self, tag='ecut'):
        dirlist = glob.glob("dir-*")
        npts = len(dirlist)
        ecutlist = np.zeros(npts)
        engylist = np.zeros(npts)
        for i in range(npts):
            print dirlist[i]
            os.chdir(dirlist[i])
            ecutlist[i] = int(dirlist[i][-2:])
            (engylist[i], vol, stress) = self.qe_get_energy_stress()
            os.chdir(self.root)
        np.savetxt('{}.txt'.format(tag), [ecutlist, engylist])
        return

    def clc_lattice(self, tag='bcc'):
        rng = [-10, 10]
        cnt = 0
        data = np.zeros([rng[1] - rng[0], 2])
        for i in range(rng[0], rng[1]):
            if i >= 0:
                dirname = "dir-p-{:03d}".format(i)
            else:
                dirname = "dir-n-{:03d}".format(abs(i))
            os.chdir(dirname)
            (energy, vol, stress) = self.qe_get_energy_stress('qe.out')
            cellmtx = self.qe_get_cell('qe.in')
            if (tag == 'fcc') or (tag == 'bcc'):
                data[cnt, 0] = 2 * cellmtx[0, 1]
                data[cnt, 1] = (energy)
            cnt += 1
            os.chdir(self.root)
        print data
        np.savetxt('lat.dat', data)
        return

    def plt_data(self, tag='ecut', data=None):
        if data is None:
            [val, engy] = np.loadtxt('{}.txt'.format(tag))
        else:
            [val, engy] = data[:, 0], data[:, 1]
        engy = engy[np.argsort(val)]
        self.set_111plt()
        self.set_keys()
        self.ax.plot(np.sort(val), engy)
        self.fig.savefig('{}.png'.format(tag))
        return

    def loop_clc_data(self, opt):
        degauss0 = 0.02
        for i in range(7):
            degauss = degauss0 + 0.005 * i
            mdir = 'degauss{:4.3f}'.format(degauss)
            print mdir
            os.chdir(mdir)
            self.clc_data(tag='kpts')
            os.system("cp kpts.txt ../kpts_{:4.3f}.txt".format(degauss))
            os.chdir(self.root_dir)
        return

    def loop_plt_data(self):
        filelist = glob.glob('kpts_*')
        self.set_111plt()
        self.set_keys()
        for mfile in filelist:
            data = np.loadtxt(mfile)
            tag = mfile[5:-4]
            print tag
            data[1] = data[1][np.argsort(data[0])]
            self.ax.plot(np.sort(data[0]), data[1], label=tag)
            self.ax.legend()
        self.fig.savefig('fig-default.png')
        return

    def loop_run(self):
        dirlist = glob.glob("dir-*")
        for dirname in dirlist:
            print dirname
            os.chdir(dirname)
            os.system("mpirun -n 24 pw.x < qe.in > qe.out")
            os.chdir(self.root)
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
            fid = self.qe_write_kpts(fid)
            fid.close()
        return

    def find_lattice(self, data=None):
        if data is None:
            data = np.loadtxt('lat.txt')
            data[:, 0] = np.abs(data[:, 0])  # in case
        interps = np.linspace(data[0, 0], data[-1, 0], 201)
        print interps
        spl = InterpolatedUnivariateSpline(data[:, 0], data[:, 1])
        print "min lat", interps[np.argmin(spl(interps))]
        return data

    def set_pbs(self, dirname, od=True):
        self.set_nnodes(2)
        self.set_ppn(12)
        self.set_job_title("%s" % (dirname))
        self.set_wall_time(7)
        self.set_main_job("""
    cal_qe_lattice.py -t run
                        """)
        self.write_pbs(od=od)
        os.system("mv va.pbs %s" % (dirname))
        return

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

    if options.mtype.lower() in ['kpts', 'prepkpts', 'loopkpts']:
        drv.loop_kpoints()

    if options.mtype.lower() in ['degauss_prep', 'degauss_sub', 'degauss_clc']:
        drv.loop_degauss(opt=options.mtype.lower().split('_')[-1])

    elif options.mtype.lower() in ['ecut', 'loopecut']:
        drv.loop_ecut()

    elif options.mtype.lower() in ['kptsecut']:
        drv.loop_kpoints_ecut()

    elif options.mtype.lower() == 'pot':
        drv.loop_pots()

    elif options.mtype.lower() == 'bcc':
        drv.cal_bcc_lattice()

    elif options.mtype.lower() == 'run':
        drv.loop_run()

    elif options.mtype.lower() in ['clcecut', 'clckpts']:
        tag = options.mtype.lower()[3:]
        drv.clc_data(tag=tag)

    elif options.mtype.lower() == 'clclat':
        drv.clc_lattice()

    elif options.mtype.lower() == 'pltecut':
        drv.plt_data(tag='ecut')

    elif options.mtype.lower() == 'pltlat':
        data = drv.find_lattice()
        drv.plt_data(tag='lat', data=data)

    elif options.mtype.lower() in ['pltkpts', 'pltecut']:
        tag = options.mtype.lower()[3:]
        drv.plt_data(tag=tag)

    elif options.mtype.lower() in ['loopclckpts', 'loopclcecut']:
        tag = options.mtype.lower()[7:]
        drv.loop_clc_data(opt=tag)

    elif options.mtype.lower() in ['loopplt']:
        drv.loop_plt_data()

#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: yangchaoming
# @Date:   2017-06-13 15:37:47
# @Last Modified by:   chaomy
# @Last Modified time: 2017-11-02 00:43:52

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
        if inpot is None: 
            self.pot = md_pot_data.qe_pot.vca_W75Ta25
        output_data.output_data.__init__(self)
        get_data.get_data.__init__(self)
        gn_pbs.gn_pbs.__init__(self)
        gn_qe_inputs.gn_qe_infile.__init__(self, self.pot)
        plt_drv.plt_drv.__init__(self)
        gn_config.bcc.__init__(self, self.pot)
        self.alat0 = self.pot['lattice']
        self.element = self.pot['element']
        self.mass = self.pot['mass']
        self.kpnts = [44, 44, 44]
        self.root = os.getcwd()
        return

    def cal_bcc_lattice(self):
        for i in range(-15, 15):
            if i >= 0:
                mdir = "dir-p-{:03d}".format(i)
            else:
                mdir = "dir-n-{:03d}".format(abs(i))
            self.mymkdir(mdir)
            self.pot['latbcc'] = self.alat0 + i * 0.006
            atoms = self.set_bcc_primitive((1, 1, 1))
            self.gn_qe_bcc_lattice_infile(atoms)
            self.set_my_pbs(mdir)  
            os.system("mv qe.in {}".format(mdir))
            os.system('cp $POTDIR/{} {}'.format(self.pot['file'],
                                                mdir))
        return

    def loop_pots(self):
        potlist = {'WTa0.25': md_pot_data.qe_pot.vca_W75Ta25,
                   'WTa0.20': md_pot_data.qe_pot.vca_W80Ta20,
                   'WTa0.15': md_pot_data.qe_pot.vca_W85Ta15,
                   'WTa0.10': md_pot_data.qe_pot.vca_W90Ta10,
                   'WTa0.05': md_pot_data.qe_pot.vca_W95Ta05}
        for key in potlist.keys():
            self.mymkdir(key)
            self.__init__(potlist[key])           
            os.chdir(key)
            self.cal_bcc_lattice()
            os.chdir(os.pardir)
        return 

    def goandsub(self, mdir, rootdir):
        os.chdir(mdir)
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
                os.system('mv va.pbs {}'.format(mdir))
                os.system('mv dir-* {}'.format(mdir))
            elif opt == 'sub':
                self.goandsub(mdir, self.root)
            elif opt == 'clc':
                os.system("cp {}/kpts.txt kpts_{:4.3f}.txt".format(mdir,
                                                                   degauss))
        return

    def loop_kpoints(self):
        self.set_ecut('{}'.format(48))
        self.set_disk_io('none')
        for kpts in range(32, 50):
            self.set_kpnts((kpts, kpts, kpts))
            mdir = 'dir-kpt-{}'.format(kpts)
            self.mymkdir(mdir)
            atoms = self.set_bcc_primitive((1, 1, 1))
            self.gn_qe_bcc_lattice_infile(atoms)
            os.system('mv qe.in {}'.format(mdir))
            os.system('cp $POTDIR/{} {}'.format(self.pot['file'],
                                                mdir))
        return

    def loop_kpoints_ecut(self, opt='clc'):
        for i in range(43, 50):
            kpts = i
            mdir = 'kpts_{:02}'.format(kpts)
            self.mymkdir(mdir)
            os.chdir(mdir)
            if opt == 'prep':
                self.set_kpnts((i, i, i))
                self.loop_ecut()
                self.set_pbs(mdir)
            elif opt == 'clc':
                self.clc_data(tag='ecut')
            elif opt == 'dat':
                os.system("cp ecut.txt ../ecut_{:02}.txt".format(kpts))
            os.chdir(self.root)
        return

    def loop_ecut(self):
        for ecut in range(28, 52):
            self.set_ecut('{}'.format(ecut))
            mdir = 'dir-ecut-{}'.format(ecut)
            self.mymkdir(mdir)
            atoms = self.set_bcc_primitive((1, 1, 1))
            self.gn_qe_bcc_lattice_infile(atoms)
            os.system('mv qe.in {}'.format(mdir))
            os.system('cp $POTDIR/{} {}'.format(self.pot['file'],
                                                mdir))
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
            os.chdir(os.pardir)
        np.savetxt('{}.txt'.format(tag), [ecutlist, engylist])
        return

    def clc_lattice(self, tag='bcc'):
        rng = [-15, 15]
        cnt = 0
        data = np.zeros([rng[1] - rng[0], 2])
        for i in range(rng[0], rng[1]):
            if i >= 0: 
                mdir = "dir-p-{:03d}".format(i)
            else: 
                mdir = "dir-n-{:03d}".format(abs(i))
            os.chdir(mdir)
            (energy, vol, stress) = self.qe_get_energy_stress('qe.out')
            cellmtx = self.qe_get_cell('qe.in')
            if tag in ['fcc', 'bcc']:
                data[cnt, 0] = 2 * cellmtx[0, 1]
                data[cnt, 1] = (energy)
            cnt += 1
            os.chdir(self.root)
        print data
        np.savetxt('lat.txt', data)
        return

    def plt_data(self, tag='ecut', data=None):
        if tag in ['lat']:
            data = self.find_lattice()
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

    def loop_plt_data(self, tag='kpts'):
        filelist = glob.glob('{}_*'.format(tag))
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
        for mdir in dirlist:
            print mdir
            os.chdir(mdir)
            os.system("mpirun pw.x < qe.in > qe.out")
            os.chdir(self.root)
        return

    def gn_qe_bcc_lattice_infile(self, atoms):
        self.set_thr('1.0D-6')
        self.set_ecut('40')
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

    def find_lattice(self, data=None):
        if data is None:
            data = np.loadtxt('lat.txt')
            data[:, 0] = np.abs(data[:, 0])  # in case
        interps = np.linspace(data[0, 0], data[-1, 0], 201)
        print interps
        spl = InterpolatedUnivariateSpline(data[:, 0], data[:, 1])
        print "min lat", interps[np.argmin(spl(interps))]
        return data

    def set_my_pbs(self, mdir, od=True):
        self.set_nnodes(1)
        self.set_ppn(4)
        self.set_mem(1)
        self.set_job_title("%s" % (mdir))
        self.set_wall_time(8)
        self.set_main_job("""
mpirun pw.x < qe.in > qe.out
                         """)
        self.write_pbs(od=od)
        os.system("mv va.pbs %s" % (mdir))
        return

if __name__ == '__main__':
    usage = "usage:%prog [options] arg1 [options] arg2"
    parser = OptionParser(usage=usage)
    parser.add_option("-t", "--mtype", action="store",
                      type="string", dest="mtype")
    parser.add_option('-p', "--param", action="store",
                      type='string', dest="fargs")
    (options, args) = parser.parse_args()
    drv = cal_lattice()

    dispatcher = {'kpts': drv.loop_kpoints,
                  'degauss': drv.loop_degauss,
                  'ecut': drv.loop_ecut,
                  'kpt_ecut': drv.loop_kpoints_ecut,
                  'pots': drv.loop_pots,
                  'bcc': drv.cal_bcc_lattice,
                  'run': drv.loop_run,
                  'clcdata': drv.clc_data,
                  'clclat': drv.clc_lattice,
                  'plt': drv.plt_data,
                  'loop': drv.loop_clc_data}

    if options.fargs is not None:
        dispatcher[options.mtype.lower()](options.fargs)
    else:
        dispatcher[options.mtype.lower()]()

#!/usr/bin/env Gython
# encoding: utf-8
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-07-05 08:12:30
# @Last Modified by:   chaomy
# @Last Modified time: 2018-09-11 12:28:06


from optparse import OptionParser
from itertools import cycle
import matplotlib.ticker as ticker
import os
import re
import glob
import matplotlib.pylab as plt
import numpy as np
import pickle as pc
import ase
import ase.io
import get_data
import gn_config
import plt_drv


class lmps_neb_tools(get_data.get_data, gn_config.bcc, plt_drv.plt_drv):

    def __init__(self):
        get_data.get_data.__init__(self)
        plt_drv.plt_drv.__init__(self)

    def loop_mk_files(self):
        for i in range(8):
            Strain = 0.5 * i

            self.change_structure(0.01 * Strain)
            mdir = 'mdir-%2.2f' % (Strain)
            os.mkdir(mdir)
            os.system("cp Strain_Neb/*     %s" % (mdir))
            os.system("cp gnStructure.py   %s" % (mdir))

            os.chdir(mdir)
            os.mkdir('Init_cfg')
            os.mkdir('Final_cfg')
            os.mkdir('Init_Restart')
            os.mkdir('Final_custom')

            os.system("python gnStructure.py")
            os.system("python Intro_init.py")
            os.system("python Intro_final.py")

            gn_config.bcc.__init__(self)
            os.chdir(os.pardir)

    def loop_run_init_final(self):
        dir_list = glob.glob("z*")
        for i in range(len(dir_list)):
            mdir = dir_list[i]
            print(mdir)
            os.chdir(mdir)
            # os.system("sh clean.sh")
            os.system("lmp_linux -in in.init  > log.init  &")
            # os.system("lmp_linux -in in.final > log.final &")
            os.chdir(os.pardir)

    def loop_prepare_Neb(self):
        root_dir = os.getcwd()
        dir_list = glob.glob("z*")
        for i in range(len(dir_list)):
            mdir = dir_list[i]
            print(mdir)
            os.chdir(mdir)
            self.change_in_Neb()
            self.create_final_screw()
            os.chdir(root_dir)

    def change_structure(self, delta):
        with open("./SamplegnStructure.py", 'r') as fid:
            raw = fid.readlines()
            fid.close()

        print(raw[55])
        raw[55] = '                    [%5.4f, 1,0],\n' % (delta)
        with open("./gnStructure.py", 'w') as fid:
            for i in range(len(raw)):
                fid.write(raw[i])
            fid.close()

    def change_in_Neb(self):
        fileList = glob.glob("./Init_Restart/*")
        print(fileList[-1])
        mfile = fileList[-1]

        with open("./in.neb_dislocation", 'r') as fid:
            raw = fid.readlines()
        print(raw[14])
        raw[14] = "read_restart  %s\n" % (mfile)
        with open("in.new", 'w') as fid:
            fid.writelines(raw)
            fid.close()
        os.system("cp in.new in.neb_dislocation")

    def create_final_screw(self):
        # os.system("lmp_mpi -i in.init_dipole")
        # os.system("lmp_mpi -i in.final_dipole")
        fileList = glob.glob("./Final_custom/dump.custom.*")
        mydir = os.getcwd().split('/')[-1]
        os.system("cp  %s  ." % (fileList[-1]))
        with open(fileList[-1], 'r') as fid:
            raw = fid.readlines()
        print(raw[3])
        print(raw[9])
        with open("final.coord", 'w') as fid:
            fid.write(raw[3])
            fid.writelines(raw[9:])

        # sshdir = "$FLUX:/home/chaomy/{}".format(mydir)
        # os.system("scp init_restart  final.coord  ../dummy.lammps.ADP {}".format(sshdir))
        # os.system("rm ./Final_custom/dump.custom.*")
        # os.system("rm dummp.custom.*")

    def plot_neb_energy(self, neb_energy, figname='neb.png'):
        self.set_111plt()
        next(self.keysiter)
        next(self.keysiter)
        x = np.linspace(0, 1, len(neb_energy))
        self.ax.plot(x, 1e3 * neb_energy, label='MEAM', **next(self.keysiter))
        self.add_legends(*self.axls)
        # (110): -110  (11-2) -110
        xlabeliter = cycle(["Normalized reaction coordinate"])
        # ylabeliter = cycle(['Energy per length [meV/|b|]'])
        ylabeliter = cycle(['Energy increment [meV]'])
        # ylabeliter = cycle(['Energy [meV]'])
        self.add_x_labels(xlabeliter, *self.axls)
        self.add_y_labels(ylabeliter, *self.axls)
        self.set_tick_size(*self.axls)
        self.fig.savefig(figname, **self.figsave)

    def plt_kink_stress(self):
        self.set_111plt()
        data = np.loadtxt("d00.txt")
        # self.ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%1.f'))
        self.ax.plot(data[:, 0], data[:, 1],
                     label="NEB", **next(self.keysiter))
        self.add_x_labels(
            cycle(["Normalized reaction coordinate"]), *self.axls)
        self.add_y_labels(cycle(['Energy [eV]']), *self.axls)
        self.add_legends(*self.axls)
        self.set_tick_size(*self.axls)
        self.fig.savefig("FIG_NEB.png", **self.figsave)

    def read_lmp_log_file_dis_kink(self, figname='neb.png'):
        file_list = glob.glob("log.lammps.*")
        nlogs = len(file_list)
        neb_energy = []
        for i in range(nlogs):
            mfile = "log.lammps.{:02d}".format(i)
            neb_energy.append(self.md_get_final_energy(mfile))
        print(np.argmax(neb_energy))
        neb_energy = np.array(neb_energy)
        neb_energy -= np.min(neb_energy)
        self.plot_neb_energy(neb_energy, figname)
        data = np.ndarray([len(neb_energy), 2])
        data[:, 0] = np.linspace(0, 1, len(neb_energy))
        data[:, 1] = neb_energy
        print("max_engy", np.max(data[:, 1]))
        np.savetxt('data.txt', data)

    def plot_two(self):
        Nb_neb = np.loadtxt("data.nb.txt")
        Mo_neb = np.loadtxt("data.mo.txt")

        self.set_111plt()
        next(self.keysiter)
        next(self.keysiter)
        self.ax.plot(Nb_neb[:, 0], 1e3 * Nb_neb[:, 1],
                     label='Nb', **next(self.keysiter))
        self.ax.plot(Mo_neb[:, 0], 1e3 * Mo_neb[:, 1],
                     label='Mo', **next(self.keysiter))

        self.add_legends(*self.axls)
        # (110): -110  (11-2) -110
        xlabeliter = cycle(["Normalized reaction coordinate"])
        ylabeliter = cycle(['Energy per length [meV/|b|]'])
        # ylabeliter = cycle(['Energy [meV]'])
        self.add_x_labels(xlabeliter, *self.axls)
        self.add_y_labels(ylabeliter, *self.axls)
        self.set_tick_size(*self.axls)
        self.fig.savefig("FIG_NEB.png", **self.figsave)

    def read_lmp_log_file(self, figname='neb.png'):
        # mydir = os.getcwd().split('/')[-1]
        # sshdir = "$FLUX:/home/chaomy/{}".format(mydir)
        # os.system("scp {}/log.lammps.* .".format(sshdir))
        file_list = glob.glob("log.lammps.*")
        nlogs = len(file_list)
        neb_energy = []
        for i in range(1, nlogs):
            mfile = "log.lammps.%d" % (i)
            print(mfile)
            neb_energy.append(self.md_get_final_energy(mfile))
        print(np.argmax(neb_energy))
        neb_energy = np.array(neb_energy)
        neb_energy -= np.min(neb_energy)

        # optional !!!
        # area = 68.2737796696
        # neb_energy /= area

        self.plot_neb_energy(neb_energy, figname)
        data = np.ndarray([len(neb_energy), 2])
        data[:, 0] = np.linspace(0, 1, len(neb_energy))
        data[:, 1] = neb_energy
        print("max_engy", np.max(data[:, 1]))
        np.savetxt('data.txt', data)

    def read_screen(self):
        log_files = glob.glob("./screen.*")
        neb_energy = []

        find_energy1 = "\s*Energy initial, next-to-last, final =\s*"
        find_energy2 = "(-?\d*\.\d*)\s*(-?\d*\.\d*)\s*(-?\d*\.\d*)"
        find_energy = find_energy1 + find_energy2

        re_energy = re.compile(find_energy, re.DOTALL)
        count = 0

        for file in log_files:
            print(file)
            with open(file, 'r') as fid:
                raw = fid.read()
                fid.close()

            neb_energy.append(float(re_energy.findall(raw)[-1][-1]))
            print("%d %f\n" % (count, neb_energy[count]))
            count += 1

        neb_energy = np.array(neb_energy)
        ########################################################
        # for change the unit to be meV/ burger
        ########################################################
        #  neb_energy = 1000 * neb_energy * 0.5  #  meV / burger
        #  neb_energy = np.delete(neb_energy, np.argmin(neb_energy));

        print("after delete", len(neb_energy))
        neb_energy = neb_energy - np.min(neb_energy)
        return neb_energy

    def restart_neb(self, cut=False):
        cut = True
        numlines = 739200 + 9
        if cut is True:
            dumpfilelist = glob.glob("dump.all.*")
            for i in range(len(dumpfilelist)):
                file = "dump.all.%d" % (i + 1)
                os.system("tail -n %d %s  > dump.%d" % (numlines, file, i + 1))
        else:
            atoms = ase.io.read("./dump.1", format="lammps-dump")
            print(atoms)
            self.write_lmp_config_data(atoms, "dump.init.data")

    def interp(self):
        atomsi = ase.io.read("contcar.init", format='vasp')
        print(atomsi.get_cell())

        atomsf = ase.io.read("contcar.final", format='vasp')
        print(atomsf.get_cell())

        npts = 6
        delta = 1 / (npts - 1)
        for i in range(npts):
            r = (i) * delta
            print(r)
            pos = (1 - r) * atomsi.get_positions() + r * atomsf.get_positions()
            atoms = atomsi.copy()
            atoms.set_positions(pos)
            mdir = "{:02d}".format(i)
            self.mymkdir(mdir)
            # os.system("cp INCAR KPOINTS POTCAR {}".format())
            ase.io.write("{}/POSCAR".format(mdir), images=atoms, format='vasp')

if __name__ == '__main__':
    usage = "usage:%prog [options] arg1 [options] arg2"
    parser = OptionParser(usage=usage)
    parser.add_option("-t", "--mtype", action="store",
                      type="string", dest="mtype")
    parser.add_option('-p', "--param", action="store",
                      type='string', dest="fargs")

    (options, args) = parser.parse_args()

    drv = lmps_neb_tools()
    dispatcher = {'plt': drv.read_lmp_log_file,
                  'plt2': drv.plot_two,
                  'screen': drv.read_screen,
                  'rst': drv.restart_neb,
                  'adj': drv.create_final_screw,
                  'inter': drv.interp,
                  'kink': drv.plt_kink_stress}

    if options.fargs is not None:
        dispatcher[options.mtype.lower()](options.fargs)
    else:
        dispatcher[options.mtype.lower()]()

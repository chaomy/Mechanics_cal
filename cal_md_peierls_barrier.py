#!/usr/bin/env python
# encoding: utf-8

###################################################################
#
# File Name : ./cal_md_peierls_barrier.py
#
###################################################################
#
# Purpose :     # use the NEB method in LAMMPS
#
# Creation Date :
# Last Modified : Wed Mar 29 15:06:51 2017
# Created By    : Chaoming Yang
#
###################################################################

import os
import ase
import glob
import numpy as np
from multiprocessing import Pool
import matplotlib.pylab as plt
from scipy.interpolate import InterpolatedUnivariateSpline
from optparse import OptionParser
import ase.lattice
import gn_config
import get_data
import Intro_vasp
import gn_lmp_infile
import cal_md_dis_schmid
import cal_md_dislocation
import gn_pbs
import plt_drv

#  import atomman as am
#  import atomman.lammps as lmp
#  import atomman.unitconvert as uc


def unwrap_self_run_lammps(arg, **kwarg):
    return md_reaction_coordinate.run_lmp_minimize(*arg, **kwarg)


class md_reaction_coordinate(gn_config.bcc,
                             gn_config.fcc,
                             gn_config.hcp,
                             get_data.get_data,
                             gn_pbs.gn_pbs,
                             Intro_vasp.vasp_change_box,
                             gn_lmp_infile.gn_md_infile,
                             plt_drv.plt_drv):

    def __init__(self):
        self.pot = self.load_data('../pot.dat')
        #  self.pot = md_pot_data.md_pot.Nb_adp_tmp
        gn_pbs.gn_pbs.__init__(self)

        self._alat = self.pot['lattice']
        self._burger = self._alat * np.sqrt(3.) / 2.
        self.potential_file = self.pot['file']
        self.structure = self.pot['structure']

        Intro_vasp.vasp_change_box.__init__(self, self.pot)
        gn_lmp_infile.gn_md_infile.__init__(self, self.pot)

        self.shmid_drv = cal_md_dis_schmid.cal_bcc_schmid(self.pot)
        self.mddis_drv = cal_md_dislocation.md_dislocation(self.pot)

        if self.structure == 'bcc':
            gn_config.bcc.__init__(self, self.pot)

        elif self.structure == 'fcc':
            gn_config.fcc.__init__(self, self.pot)

        elif self.structure == 'hcp':
            gn_config.hcp.__init__(self, self.pot)

        self.set_config_file_format('lmp')
        self.num_configures = 15
        plt_drv.plt_drv.__init__(self)
        self.set_keys()
        self.root_dir = os.getcwd()
        self.mfigsize = (8.5, 4.3)
        return

    def get_init_final_state(self):
        # state initial #
        poscar_atoms = ase.io.read("./state_init/lmp_init.cfg", format='cfg')
        atoms_tobe_changed = ase.io.read("./state_init/W.0.cfg", format='cfg')

        map_list = self.map_atoms_list(poscar_atoms, atoms_tobe_changed)
        new_atoms = self.sort_atoms(atoms_tobe_changed, map_list)

        ase.io.write("peierls_init.cfg",
                     new_atoms,
                     format='cfg')

        # state final #
        poscar_atoms = ase.io.read("./state_final/lmp_init.cfg",
                                   format='cfg')
        atoms_tobe_changed2 = ase.io.read("./state_final/W.0.cfg",
                                          format='cfg')
        map_list = self.map_atoms_list(poscar_atoms,
                                       atoms_tobe_changed2)

        new_atoms2 = self.sort_atoms(atoms_tobe_changed2,
                                     map_list)

        ase.io.write("peierls_final.cfg",
                     new_atoms2,
                     format='cfg')
        return

        ############################################################
        # use dipole method to calcualte the peierls barrier
        ############################################################
    def dipole_peierls_barrier(self):
        e1 = 1. / 3. * np.array([1., 1., -2.])
        e2 = np.array([0.5, 0.5, 0.5])
        e3 = 1. / 2. * np.array([1, -1, 0])

        sizen = 1
        n = 7 * sizen
        m = 11 * sizen
        t = 1

        atoms = ase.lattice.cubic.BodyCenteredCubic(directions=[e1, e2, e3],
                                                    latticeconstant=self._alat,
                                                    size=(n,  t,  m),
                                                    symbol=self.pot['element'],
                                                    pbc=(1, 1, 1))
        # add shiftment to the supercell #
        atoms = self.mddis_drv.cut_half_atoms_new(atoms, "cutz")
        supercell = atoms.get_cell()
        strain = np.mat([[1.0, 0.0, 0.0],
                         [0.0, 1.0, 0.0],
                         [0.5, 0.5, 1.0]])

        supercell = strain * supercell
        atoms.set_cell(supercell)

        atoms2 = atoms.copy()
        unitx = np.sqrt(6) / 3. * self._alat
        unity = np.sqrt(2) / 2. * self._alat

        sx = 10.0 * sizen
        sy = 5 * sizen
        ix = 10.5 * sizen

        c1 = [(sx) * unitx, (sy + 1. / 3.) * unity]
        c2 = [(sx + ix) * unitx, (sy + 2. / 3.) * unity]
        center = [c1, c2]
        atoms = self.mddis_drv.intro_dipole_screw_atoms_LMP(atoms,
                                                            center=center,
                                                            lattice=self._alat)
        self.write_lmp_config_data(atoms, "init.txt")
        movex = 1.0
        c1 = [(sx + movex) * unitx, (sy + 1. / 3.) * unity]
        c2 = [(sx + ix + movex) * unitx, (sy + 2. / 3.) * unity]
        center = [c1, c2]
        atoms = self.mddis_drv.intro_dipole_screw_atoms_LMP(atoms2,
                                                            center=center,
                                                            lattice=self._alat)
        self.write_lmp_config_data(atoms, "final.txt")
        # self.write_lmp_coords(atoms, "final.coord")
        return

        ############################################################
        # cluster method  (large supercell)
        ############################################################
    def cluster_peierls_barrier(self):
        alat = self._alat
        move1 = [-alat * np.sqrt(6.) / 3., 0, 0]
        self.shmid_drv.make_screw_plate(size=[200, 260, 2],
                                        rad=[340, 380],
                                        move=move1,
                                        tag='[211]',
                                        filename="final.txt",
                                        opt='neb')

        os.system("cp init.txt  init0.txt")
        os.system("cp final.txt init1.txt")

        # Mid is 150 160
        # Mid2   200 210
        # Lar    400 410
        return

    def interp_peierls(self):
        init_file = glob.glob("./state_init/*")
        final_file = glob.glob("./state_final/*")

        atoms_init = ase.io.read(init_file[-1], format='cfg')
        atoms_final = ase.io.read(final_file[-1], format='cfg')

        self.write_lmp_config_data(atoms_init, "init.data")
        self.write_lmp_config_data(atoms_final, "final.data")

        atoms_interp = atoms_init.copy()

        position_init = atoms_init.get_positions()
        position_final = atoms_final.get_positions()

        delta = position_final - position_init

        num_list = [47, 78, 79, 49, 50, 181, 180, 178, 179, 149]

        delta_copy = np.zeros(np.shape(delta))
        for i in range(len(num_list)):
            delta_copy[num_list[i]] = delta[num_list[i]]

        print delta_copy
        delta = delta_copy

        inter_num = self.num_configures
        delta = delta / inter_num

        for i in range(inter_num):
            dir_name = 'dir-%03d' % (i)

            if not os.path.isdir(dir_name):
                os.mkdir(dir_name)

            os.system("cp  ./w_eam4.fs  %s" % (dir_name))
            os.chdir(dir_name)

            filename_data = "lmp_init_%d.txt" % (i)

            atoms_interp = position_init + i * delta

            atoms_init.set_positions(atoms_interp)

            self.gn_md_minimize_cfg(filename_data,
                                    "w_eam4.fs",
                                    "W",
                                    fix=True)

            self.write_lmp_config_data(atoms_init, filename_data)
            if not os.path.isdir("cfg"):
                os.mkdir("cfg")
            os.chdir(self.root_dir)
        return

    #  reaction coordinate method #
    def peierls(self):
        lattice = 3.143390
        e1 = np.array([1.,   1.,  -2.])
        e2 = np.array([-1.,  1.,   0])
        e3 = np.array([0.5,  0.5,  0.5])

        # 3, 5;   5, 9;  7, 11
        r, s = 7,  11
        v1 = r * e1
        v2 = 0.5 * (r * e1 + s * e2) + 0.5 * e3
        v3 = e3

        self.set_lattce_constant(self._alat)
        self.set_element('W')

        inter_num = self.num_configures
        delta = 1. / inter_num

        for i in range(inter_num):
            dir_name = 'dir-%03d' % (i)

            if not os.path.isdir(dir_name):
                os.mkdir(dir_name)

            os.system("cp  ./w_eam4.fs  %s" % (dir_name))
            os.chdir(dir_name)

            # generate config #
            atoms = self.set_bcc_convention([v1, v2, v3],
                                            (r, 1, 1))  # z periodic 12
            atoms2 = atoms.copy()
            atoms2 = self.cut_half_atoms(atoms2)

            #  movex = np.sqrt(6.)/3. * lattice * delta * i;

            ase.io.write("lmp_perf.cfg",
                         atoms2,
                         "cfg")
            s = i * delta
            atoms = self.intro_dipole_screw_atoms(
                atoms, lattice, move_x=None, input_s=s)
            atoms = self.cut_half_atoms(atoms)

            ase.io.write("lmp_init.cfg",
                         atoms,
                         "cfg")

            filename_data = "lmp_init_%d.txt" % (i)

            self.gn_md_minimize_cfg(filename_data,
                                    "w_eam4.fs",
                                    "W",
                                    fix=True)

            self.write_lmp_config_data(atoms, filename_data)

            if not os.path.isdir("cfg"):
                os.mkdir("cfg")

            os.chdir(self.root_dir)
        return

    def run_lmp_minimize(self, loc_dir):
        os.chdir(loc_dir)
        os.system("lmp_mpi -in in.minimize")
        os.chdir(self.root_dir)
        return

    def multi_thread_minimize(self):
        dir_list = glob.glob("dir-*")
        num_threads = len(dir_list)
        pool = Pool(processes=num_threads)
        pool.map(unwrap_self_run_lammps,
                 zip([self] * num_threads, dir_list))
        return

    def collect_peierls_energy(self):
        num = self.num_configures
        energy_list = []
        count_list = []
        for i in range(0, num):
            dir_name = 'dir-%03d' % (i)
            config_name = 'config-%03d.cfg' % (i)
            config_png = 'init-%03d.png' % (i)
            os.chdir(dir_name)

            energy_list.append(self.md_get_final_energy())
            count_list.append(i)

            filelist = glob.glob("./cfg/*")
            file = filelist[-1]

            os.system("python ../../DD_map_dis.py")
            os.system("cp  %s  ../%s" % (file, config_name))
            os.system("cp  init.png  ../%s" % (config_png))

            os.chdir(self.root_dir)

        energy_list = np.array(energy_list)
        energy_list = energy_list - np.min(energy_list)

        with open("peierls.txt", 'w') as fid:
            for i in range(num):
                fid.write("%d %7.6f\n" % (count_list[i], energy_list[i]))
        print energy_list
        self.plot()
        return

    def plot_engy(self, **kwargs):
        if 'engy' not in kwargs.keys():
            (distances, engy) = np.loadtxt("./peierls.txt")
        else:
            distances, engy = kwargs['distances'], kwargs['engy']

        if 'ax' not in kwargs.keys():
            self.set_111plt(self.mfigsize)
            ax = self.ax
        else:
            ax = kwargs['ax']

        if 'pltkeys' in kwargs.keys():
            pltkeys = kwargs['pltkeys']
        else:
            pltkeys = self.pltkwargs

        ax.plot(distances, engy, label=kwargs['label'], **pltkeys)
        ax.legend(**self.legendarg)
        return ax

    def read_peierls_barrier_neb(self):
        file_list = glob.glob("log.lammps.*")
        nlogs = len(file_list)
        neb_energy = []
        for i in range(nlogs):
            mfile = "log.lammps.%d" % (i)
            neb_energy.append(self.md_get_final_energy_e(mfile))
        neb_energy = np.array(neb_energy)
        neb_energy -= np.min(neb_energy)
        # neb_energy *= (1. / 4.)     # dislocation dipole (count how many
        # burgers vector)
        disp = np.linspace(0.0, self._burger, len(neb_energy))
        data = np.array([disp, neb_energy])
        np.savetxt('pengy.txt', data)
        return data

    def set_pltkargs(self):
        keyslist = [{'linestyle': '-', 'color': 'b', 'linewidth': 2,
                     'marker': 'o'},
                    {'linestyle': '--', 'color': 'r', 'linewidth': 2,
                        'marker': '<'},
                    {'linestyle': ':', 'color': 'g', 'linewidth': 2,
                        'marker': '>'}]
        return keyslist

    def plot_multi_peierls_barrier(self):
        configlist = ['pengy.txt.eam', 'pengy.txt.adp']
        labellist = ['eam', 'adp']
        self.set_111plt(self.mfigsize)
        keyslist = self.set_pltkargs()
        cnt = 0

        for config in configlist:
            (distances, engy) = np.loadtxt(config)
            inargs = {'distances': distances,
                      'engy': engy,
                      'ax': self.ax,
                      'pltkeys': keyslist[cnt],
                      'label': labellist[cnt]}
            cnt += 1
            self.ax = self.plot_engy(**inargs)
        self.ax.set_xlabel('displacement [A]', {'fontsize': self.mlabelsize})
        self.ax.set_ylabel('energy [eV]', {'fontsize': self.mlabelsize})
        plt.yticks(size=self.mlabelsize)
        plt.xticks(size=self.mlabelsize)
        plt.savefig("mpeierls.png", **self.figsave)
        return

    def cal_peierls_stress(self):
        if os.path.isfile('pengy.txt'):
            data = np.loadtxt('pengy.txt')
        else:
            data = self.read_peierls_barrier_neb()
        spl = InterpolatedUnivariateSpline(data[0], data[1], k=3)
        splder1 = spl.derivative()
        interpnts = np.linspace(0, self._burger, 51)
        stress = splder1(interpnts)
        plt.savefig("tmp.png")
        print np.max(stress)
        return

    def loop_rcut(self):
        import cal_md_neb
        drvneb = cal_md_neb.lmps_neb_tools()
        npt = 7
        for i in range(npt):
            rcut = 5.17 + 0.015 * i
            dirname = 'dir-%5.4f' % (rcut)
            os.system("cp ../looprcut/{}/dummy.lamm*  .".format(dirname))
            potname = 'pot_%5.4f_lat' % (rcut)
            self.pot = self.load_data('../{}'.format(potname))
            drv.dipole_peierls_barrier()
            drvneb.create_final_screw()
            os.system("mpirun -n 16 lmp_mpi -i in.neb_dislocation_dipole -p 16x1")
            drvneb.read_lmp_log_file(figname='fig.%5.4f.png' % (rcut))
        return

    def record(self):
        init_nofix = [-2047.73032095]
        init_fix = [-2047.73032095]
        init_fix = [-2031.95617241]
        final_fix = [-2034.69835008]
        init_nofix = [-2031.95617241]
        final_nofix = [-2038.13741398]
        return

    def convert_dump_to_poscar(self):
        initfile = glob.glob("bcc.init.*")[-1]
        finalfile = glob.glob("bcc.final.*")[-1]
        print initfile, finalfile
        initatoms = ase.io.read(initfile, format='lammps-dump')
        finalatoms = ase.io.read(finalfile, format='lammps-dump')

        ase.io.write("POSCAR00", images=initatoms, format='vasp')
        ase.io.write("POSCAR01", images=finalatoms, format='vasp')
        return


if __name__ == "__main__":
    usage = "usage:%prog [options] arg1 [options] arg2"
    parser = OptionParser(usage=usage)
    parser.add_option("-t", "--mtype",
                      action="store",
                      type="string",
                      dest="mtype",
                      help="",
                      default="prp_r")

    (options, args) = parser.parse_args()
    drv = md_reaction_coordinate()

    if options.mtype.lower() == 'dipole':
        drv.dipole_peierls_barrier()

    elif options.mtype.lower() == 'cluster':
        drv.cluster_peierls_barrier()

    elif options.mtype.lower() == 'vasp':
        drv.convert_dump_to_poscar()

    elif options.mtype.lower() == 'stress':
        drv.cal_peierls_stress()

    elif options.mtype.lower() == 'cmp':
        drv.plot_multi_peierls_barrier()

    elif options.mtype.lower() == 'looprcut':
        drv.loop_rcut()

    #  job.collect_peierls_energy()
    #  drv.peierls()
    #  job.multi_thread_minimize()

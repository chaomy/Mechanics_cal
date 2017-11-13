#!/usr/bin/env python
# encoding: utf-8
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-07-05 08:12:30
# @Last Modified by:   chaomy
# @Last Modified time: 2017-11-09 16:05:31


from optparse import OptionParser
import os
import ase
import shutil
import numpy as np
import ase.lattice.cubic as Cubic
import md_pot_data
import gn_config
import get_data
from utils import Intro_vasp
import gn_pbs
import gn_lmp_infile
import atomman as am
import cal_md_dislocation_hcp
import cal_md_dislocation_bcc
import cal_md_dislocation_fcc


class md_dislocation(gn_config.bcc,
                     gn_config.fcc, gn_config.hcp,
                     get_data.get_data, gn_pbs.gn_pbs,
                     Intro_vasp.vasp_change_box,
                     gn_lmp_infile.gn_md_infile,
                     cal_md_dislocation_hcp.md_dislocation_hcp,
                     cal_md_dislocation_bcc.md_dislocation_bcc,
                     cal_md_dislocation_fcc.md_dislocation_fcc):

    def __init__(self, pot=None):
        if pot is None:
            self.pot = md_pot_data.md_pot.feyo
            # self.pot = md_pot_data.md_pot.Nb_eam
        else:
            self.pot = pot
        gn_pbs.gn_pbs.__init__(self)
        Intro_vasp.vasp_change_box.__init__(self, self.pot)
        gn_lmp_infile.gn_md_infile.__init__(self, self.pot)
        if self.pot['structure'] in 'bcc':
            gn_config.bcc.__init__(self, self.pot)
            cal_md_dislocation_bcc.md_dislocation_bcc.__init__(self)

        elif self.pot['structure'] in 'fcc':
            gn_config.fcc.__init__(self, self.pot)
            cal_md_dislocation_fcc.md_dislocation_fcc.__init__(self)

        elif self.pot['structure'] in 'hcp':
            gn_config.hcp.__init__(self, self.pot)
            cal_md_dislocation_hcp.md_dislocation_hcp.__init__(self)

        self.set_config_file_format('lmp')
        return

    def intro_kink_pair(self):
        e1 = 1. / 3. * np.array([1., 1., -2.])
        e2 = 1. / 2. * np.array([-1., 1., 0])
        e3 = np.array([0.5, 0.5, 0.5])

        self.set_lattce_constant(self.pot['lattice'])
        self.set_element('W')

        atoms = self.set_bcc_convention([e1, e2, e3], (30, 30, 60))

        xc1 = (0.0 + (-2.56656)) / 2. + 45 * \
            np.sqrt(6.) / 3. * self.pot['lattice']
        yc1 = (0.0 + (2.22271)) / 2. + 15 * np.sqrt(2.) * self.pot['lattice']
        H = np.sqrt(2. / 3.0) * self.pot['lattice']

        h = 0.0 * H
        atoms = self.intro_kink_screw_dislocations(
            atoms, (xc1, yc1), (xc1 + H, yc1), h, 1. / 4.)

        ase.io.write("lmp_init.cfg",
                     atoms,
                     "cfg")

        fname = "init.data"
        self.write_lmp_config_data(atoms, fname)
        self.gn_md_minimize_cfg("init.data",
                                "./w_eam4.fs",
                                "W")
        return

    def cal_nye(self):
        torient = 'z'
        if torient == 'y':
            e1 = 1. / 3. * np.array([1., 1., -2.])
            e2 = np.array([0.5, 0.5, 0.5])
            e3 = 1. / 2. * np.array([1, -1, 0])

        e1 = 1. / 3. * np.array([1., 1., -2.])
        e2 = 1. / 2. * np.array([-1., 1., 0])
        e3 = np.array([0.5, 0.5, 0.5])

        # unit cell
        unit_cell = np.array([e1, e2, e3])

        sizen = 1
        n = 7 * sizen
        m = 11 * sizen
        t = 3

        atoms = Cubic.BodyCenteredCubic(directions=[e1, e2, e3],
                                        latticeconstant=self.pot['lattice'],
                                        size=(n, m, t),
                                        symbol='W',
                                        pbc=(1, 1, 1))

        # perfect positions
        p = self.pot['lattice'] * self.bccneigh

        atoms = self.intro_dipole_screw_atoms(atoms,
                                              self.pot['lattice'])

        ase.io.write("lmp_dis.cfg",
                     atoms,
                     "cfg")
        # system
        system, elements = am.convert.ase_Atoms.load(atoms)

        r = (3**0.5 / 2. + 1) / 2.
        system.nlist(self.pot['lattice'] * r)

        nye_rlt = am.defect.nye_tensor(system, p, axes=unit_cell)

        for key, value in nye_rlt.iteritems():
            print key, value
            system.atoms_prop(key=key, value=value)

        print atoms.get_positions()

        int_sum = np.empty(3)
        for i in range(3):
            int_sum[i] = am.plot.interpolate_contour(system,
                                                     'Nye_tensor',
                                                     index=[1, i],
                                                     cmap='bwr')[0]
        print('burgers vector estimate = ', int_sum)
        return

        ############################################################
        # since lammmps can only use  xy xz yz
        # new x is old y ; new y is old z ; new z is old x
        # x  1  1 -2
        # y  1  1  1
        # z  1 -1  0
        ############################################################
    def modified_cal_disp_dipo(self, movex=0.0, torient='y'):
        if torient == 'y':
            e1 = 1. / 3. * np.array([1., 1., -2.])
            e2 = np.array([0.5, 0.5, 0.5])
            e3 = 1. / 2. * np.array([1, -1, 0])

        self.set_lattce_constant(self.pot['lattice'])
        self.set_element('W')
        sizen = 1
        n = 7 * sizen
        m = 11 * sizen
        t = 1
        atoms = Cubic.BodyCenteredCubic(directions=[e1, e2, e3],
                                        latticeconstant=self.pot['lattice'],
                                        size=(n, t, m),
                                        symbol=self._element,
                                        pbc=(1, 1, 1))
        atoms = self.cut_half_atoms_new(atoms, "cutz")
        supercell = atoms.get_cell()
        strain = np.mat([[1.0, 0.0, 0.0],
                         [0.0, 1.0, 0.0],
                         [0.5, 0.5, 1.0]])
        supercell = strain * supercell
        atoms.set_cell(supercell)
        atoms2 = atoms.copy()
        unitx = np.sqrt(6) / 3. * self.pot['lattice']
        unity = np.sqrt(2) / 2. * self.pot['lattice']
        c1 = [(10 + movex) * unitx, (5 + 1. / 3.) * unity]
        c2 = [(10 + 10.5 + movex) * unitx, (5 + 2. / 3.) * unity]
        center = [c1, c2]
        atoms = self.intro_dipole_screw_atoms_LMP(atoms, center=center,
                                                  lattice=self.pot['lattice'])
        self.write_lmp_config_data(atoms)
        ase.io.write("lmp_perf_modified.cfg", atoms2, "cfg")
        ase.io.write("lmp_dis.cfg", atoms, "cfg")
        os.system("mv lmp_perf_modified.cfg output")
        os.system("mv lmp_dis.cfg output")
        return atoms

    def cal_disp_dipo_lisa(self):
        e1 = 1. / 3. * np.array([-1., -1., 2.])
        e2 = 1. / 2. * np.array([1., -1., 0])
        e3 = 1. / 2. * np.array([1., 1., 1.])
        n, m = 21, 13

        v1 = n * e1 - 1. / (3. * m) * e3
        v2 = 0.5 * n * e1 + m * e2 + (0.5 - 1. / (6 * m)) * e3
        v3 = e3

        print v1, v2, v3

        v1 = np.round(v1, decimals=0)
        v2 = np.round(v2, decimals=0)
        v3 = np.round(v3, decimals=1)

        self.set_lattce_constant(self.pot['lattice'])
        self.set_element('Nb')
        print v1, v2, v3

        atoms = Cubic.BodyCenteredCubic(directions=[v1, v2, v3],
                                        latticeconstant=3.0,
                                        size=(n, 1, 1),
                                        symbol='Nb',
                                        pbc=(1, 1, 1))
        ase.io.write("lmp_init.xyz", atoms, "xyz")
        return

    def add_strain(self, atoms, delta):
        cell = atoms.get_cell()
        strain = np.mat([[1, 0, 0],
                         [delta, 1, 0],
                         [0, 0, 1]], 'float')

        positions = np.mat(atoms.get_positions())
        positions = positions * strain
        cell = strain * cell
        atoms.set_positions(positions)
        atoms.set_cell(cell)
        return atoms

    def static_add_stress(self):
        atoms = self.prep_relaxed_dislocation()
        self.gn_md_tensile("Nb.eam.alloy.webarchive", "Nb")
        strain_dir = "add_strain"
        if not os.path.isdir(strain_dir):
            os.mkdir(strain_dir)
        os.chdir(strain_dir)

        for i in range(0, 20):
            delta = 0.05 * i
            filename = "delta%4.3f.txt" % (delta)
            atoms_new = self.add_strain(atoms.copy(), delta)
            self.write_lmp_config_data(atoms_new, filename)
        return

    def prep_relaxed_dislocation(self):
        # atoms = ase.io.read('bcc.init.213', format='lammps-dump')
        atoms = ase.io.read('bcc.init.383', format='lammps-dump')
        self.write_lmp_config_data(atoms, "relaxed.txt")
        return atoms

    def cal_non_periodic_screw(self):
        e1 = np.array([1., 1., -2.])
        e2 = np.array([-1., 1., 0.])
        e3 = np.array([1., 1., 1.])
        atoms = self.set_bcc_convention(
            [e1, e2, e3], (80, 50, 40))  # z periodic 12
        atoms = self.intro_single_screw_atoms(atoms)
        self.write_lmp_config_data(atoms)

        #  if not os.path.isdir("restart"):
        #  os.mkdir("restart")
        #  os.mkdir("cfg")
        return

    # calculate the dislocation velocity under shear by LMP #
    def loop_write_pbs(self):
        #  templist = [300, 600, 900, 1200, 1500, 1800];
        #  templist = [400, 500, 700, 800];

        templist = [1000, 1100]
        stress = '0.05'

        for i in range(len(templist)):
            temp = templist[i]
            dirname = 'dir-T%d' % (temp)

            if not os.path.isdir(dirname):
                os.mkdir(dirname)

            shutil.copy("relaxed.txt", dirname)
            shutil.copy("W.set.txt", dirname)

            os.chdir(dirname)

            self.mymkdir('restart')
            self.mymkdir('force_cfg')

            self.set_nnodes(1)
            self.set_wall_time(90)
            self.set_job_title("W-%s-T%d" % (stress, temp))
            self.set_main_job(
                "mpirun lmp_linux -in in.md_addforce  > screen.log ")
            self.write_pbs()
            self.gn_md_add_force(temp, stress)
            os.chdir(os.pardir)
        return


if __name__ == "__main__":
    usage = "usage:%prog [options] arg1 [options] arg2"
    parser = OptionParser(usage=usage)
    parser.add_option("-t", "--mtype", action="store",
                      type="string", dest="mtype")
    parser.add_option('-p', "--param", action="store",
                      type='string', dest="fargs")
    (options, args) = parser.parse_args()
    drv = md_dislocation()
    dispatcher = {'kink': drv.intro_kink_pair,
                  # ouptut the rotated dipole configuration #
                  'modify': drv.modified_cal_disp_dipo,
                  'hcpedge': drv.hcp_edge_dislocation,
                  'bccedge': drv.cal_single_edge_dislocations,
                  'bccscrew': drv.cal_single_screw_dislocations,
                  'nye': drv.cal_nye,
                  'cuau': drv.cal_cu3Au_dis,
                  'ani': drv.intro_ani_edge_fcc,
                  'looppbs': drv.loop_write_pbs,
                  'static': drv.static_add_stress,
                  'prep': drv.prep_relaxed_dislocation,
                  'gnedge': drv.cal_single_edge_dislocations_read,
                  'gnscrew': drv.cal_single_screw_dislocatoins_read}

    if options.fargs is not None:
        dispatcher[options.mtype.lower()](options.fargs)
    else:
        dispatcher[options.mtype.lower()]()

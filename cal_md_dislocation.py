#!/usr/bin/env python
# encoding: utf-8

###################################################################
#
# File Name : ./cal_md_dislocation.py
#
###################################################################
#
# Purpose :  introduce dislocation for md calculation (LAMMPS)
#
# Creation Date :
# Last Modified : Sun Apr  9 20:20:57 2017
# Created By    : Chaoming Yang
#
###################################################################

import os
import ase
import glob
import shutil
import numpy as np
import ase.lattice.cubic as Cubic
from optparse import OptionParser
import md_pot_data
import gn_config
import get_data
import Intro_vasp
import gn_pbs
import gn_lmp_infile
import atomman as am
import atomman.lammps as lmp


class md_dislocation(gn_config.bcc,
                     gn_config.fcc,
                     gn_config.hcp,
                     get_data.get_data,
                     gn_pbs.gn_pbs,
                     Intro_vasp.vasp_change_box,
                     gn_lmp_infile.gn_md_infile):

    def __init__(self, pot=None):
        if pot is None:
            self.pot = md_pot_data.md_pot.Nb_eam
        else:
            self.pot = pot

        gn_pbs.gn_pbs.__init__(self)
        #  self._alat = 3.1648492  # W
        #  self._alat = 3.1711
        self._alat = self.pot.get('lattice')   # w_eam4
        self.structure = self.pot['structure']

        #  self._alat = 3.304   # Ta
        #  self._alat = 3.16741543   # Mo

        # 0 K   300 K
        # Nb    3.307   3.3102
        # W     3.1652  3.165258
        # W.set 3.1648492

        print self.pot
        Intro_vasp.vasp_change_box.__init__(self,
                                            self.pot)
        gn_lmp_infile.gn_md_infile.__init__(self,
                                            self.pot)

        if self.structure == 'bcc':
            gn_config.bcc.__init__(self, self.pot)

        elif self.structure == 'fcc':
            gn_config.fcc.__init__(self, self.pot)

        elif self.structure == 'hcp':
            gn_config.hcp.__init__(self, self.pot)

        self.set_config_file_format('lmp')
        self.bccneigh = np.array([[0.5, 0.5, 0.5],
                                  [-0.5, 0.5, 0.5],
                                  [0.5, -0.5, 0.5],
                                  [0.5, 0.5, -0.5],
                                  [-0.5, -0.5, 0.5],
                                  [0.5, -0.5, -0.5],
                                  [-0.5, 0.5, -0.5],
                                  [-0.5, -0.5, -0.5]])

        self.root_dir = os.getcwd()
        return

    def fcc_screw_dipo(self):
        e1 = np.array([-2, -1., 1.])
        e2 = np.array([1., -1., 1.0])
        e3 = np.array([0.0, 1.0, 1.0])

        self.set_lattce_constant(4.05000)
        self.set_element('Al')

        atoms = Cubic.FaceCenteredCubic(directions=[[-2, -1, 1],
                                                    [1, -1, 1],
                                                    [0,  1, 1]],
                                        size=(25, 10, 1),
                                        latticeconstant=4.050,
                                        symbol='Al',
                                        pbc=(1, 1, 1))

        movex = 0.0  # -np.sqrt(6.)/3. * self._alat
        ase.io.write("lmp_perf.cfg",
                     atoms,
                     "cfg")

        atoms = self.intro_dipole_screw_atoms(atoms)
        ase.io.write("lmp_init.cfg",
                     atoms,
                     "cfg")

        self.write_lmp_config_data(atoms)
        return

    def hcp_edge_dislocation(self):
        lata = 3.2019267694893
        latc = 5.1969105399
        atoms = ase.io.read("../MgNd_2NNmeam/Mgcfg/Mg.0.cfg")

        # cut a layers to generate free surface
        atoms = self.cut_y_normal_atoms(atoms, gp_n=2)
        atoms = self.intro_single_edge_atoms(atoms)

        # cut a layer normal the burger direction
        atoms = self.cut_x_normal_atoms(atoms,
                                        lata,
                                        1,
                                        np.sqrt(3) / 4.0)
        ase.io.write("edge.cfg",
                     atoms,
                     "cfg")

        self.write_lmp_config_data(atoms)
        os.system("cp ./lmp_init.txt ../MgNd_2NNmeam")
        return

        ############################################################
        # introduce dislocation kink pair
        ############################################################
    def intro_kink_pair(self):
        e1 = 1. / 3. * np.array([1.,  1., -2.])
        e2 = 1. / 2. * np.array([-1., 1.,  0])
        e3 = np.array([0.5,  0.5,  0.5])

        self.set_lattce_constant(self._alat)
        self.set_element('W')

        atoms = self.set_bcc_convention([e1, e2, e3], (30, 30, 60))

        xc1 = (0.0 + (-2.56656)) / 2. + 45 * np.sqrt(6.) / 3. * self._alat
        yc1 = (0.0 + (2.22271)) / 2. + 15 * np.sqrt(2.) * self._alat
        H = np.sqrt(2. / 3.0) * self._alat

        h = 0.0 * H
        atoms = self.intro_kink_screw_dislocations(
            atoms, (xc1, yc1), (xc1 + H, yc1), h, 1. / 4.)

        ase.io.write("lmp_init.cfg",
                     atoms,
                     "cfg")

        fname = "init.data"
        self.write_lmp_config_data(atoms, file_name=fname)
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
                                        latticeconstant=self._alat,
                                        size=(n, m, t),
                                        symbol='W',
                                        pbc=(1, 1, 1))

        # perfect positions
        p = self._alat * self.bccneigh

        atoms = self.intro_dipole_screw_atoms(atoms,
                                              self._alat)

        ase.io.write("lmp_dis.cfg",
                     atoms,
                     "cfg")
        # system
        system, elements = am.convert.ase_Atoms.load(atoms)

        r = (3**0.5 / 2. + 1) / 2.
        system.nlist(self._alat * r)

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

        self.set_lattce_constant(self._alat)
        self.set_element('W')
        sizen = 1
        n = 7 * sizen
        m = 11 * sizen
        t = 1
        atoms = Cubic.BodyCenteredCubic(directions=[e1, e2, e3],
                                        latticeconstant=self._alat,
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
        unitx = np.sqrt(6) / 3. * self._alat
        unity = np.sqrt(2) / 2. * self._alat
        c1 = [(10 + movex) * unitx, (5 + 1. / 3.) * unity]
        c2 = [(10 + 10.5 + movex) * unitx, (5 + 2. / 3.) * unity]
        center = [c1, c2]
        atoms = self.intro_dipole_screw_atoms_LMP(atoms, center=center,
                                                  lattice=self._alat)
        self.write_lmp_config_data(atoms)
        ase.io.write("lmp_perf_modified.cfg",
                     atoms2,
                     "cfg")
        ase.io.write("lmp_dis.cfg",
                     atoms,
                     "cfg")
        os.system("mv lmp_perf_modified.cfg output")
        os.system("mv lmp_dis.cfg output")
        return atoms

    def cal_disp_dipo_lisa(self):
        e1 = 1. / 3. * np.array([-1., -1.,  2.])
        e2 = 1. / 2. * np.array([1.,  -1.,   0])
        e3 = 1. / 2. * np.array([1.,   1.,  1.])
        n, m = 21, 13

        v1 = n * e1 - 1. / (3. * m) * e3
        v2 = 0.5 * n * e1 + m * e2 + (0.5 - 1. / (6 * m)) * e3
        v3 = e3

        print v1, v2, v3

        v1 = np.round(v1, decimals=0)
        v2 = np.round(v2, decimals=0)
        v3 = np.round(v3, decimals=1)

        self.set_lattce_constant(self._alat)
        self.set_element('Nb')
        print v1, v2, v3

        atoms = Cubic.BodyCenteredCubic(directions=[v1, v2, v3],
                                        latticeconstant=3.0,
                                        size=(n, 1, 1),
                                        symbol='Nb',
                                        pbc=(1, 1, 1))
        ase.io.write("lmp_init.xyz",
                     atoms,
                     "xyz")

        return

    def cal_single_edge_dislocations(self):
        e1 = np.array([1., 1., 1.])
        e2 = np.array([-1., 1., 0.])
        e3 = np.array([-1., -1., 2.])

        self.set_lattce_constant(self._alat)
        self.set_element('W')

        atoms = self.set_bcc_convention([e1, e2, e3], (80, 60, 5))  # 5

        atoms = self.cut_y_normal_atoms(atoms)
        atoms = self.intro_dipole_edge_with_image_atoms(atoms)

        # if we don't cut a layer of atoms, will generate two dislocations
        # cut a layer normal the burger direction
        atoms = self.cut_x_normal_atoms(atoms,
                                        self._alat)
        self.write_lmp_config_data(atoms)

        self.gn_md_minimize_cfg("lmp_init.txt",
                                "W.set.txt",    # "Nb.eam.alloy.webarchive",
                                "W")

        os.system("rm cfg/*; mpirun -n 12 lmp_mpi -in in.minimize")
        return

    def cal_single_screw_dislocatoins(self):
        e1 = np.array([-1., 1., 0])
        e2 = np.array([1., 1., -2])
        e3 = np.array([0.5, 0.5, 0.5])

        self.set_lattce_constant(self._alat)
        self.set_element('Nb')

        atoms = self.set_bcc_convention(
            [e1, e2, e3], (20, 20, 10))  # z periodic 12

        atoms = self.cut_y_normal_atoms(atoms)
        ase.io.write("lmp_perf.xyz",
                     atoms,
                     "xyz")

        atoms = self.intro_single_screw_atoms(atoms)
        #  if we don't cut a layer of atoms,
        #  it will generate two screw dislocations
        #  atoms = self.cut_z_normal_atoms(atoms, self._alat);

        self.write_lmp_config_data(atoms)
        ase.io.write("lmp_init.xyz",
                     atoms,
                     "xyz")

        #  os.system("rm cfg/*; mpirun -n 12 lmp_mpi -in in.minimize");
        return

    def add_strain(self, atoms, delta):
        cell = atoms.get_cell()
        strain = np.mat([[1, 0,  0],
                         [delta, 1,  0],
                         [0, 0,  1]], 'float')

        positions = np.mat(atoms.get_positions())
        positions = positions * strain
        cell = strain * cell
        atoms.set_positions(positions)
        atoms.set_cell(cell)
        return atoms

    def static_add_stress(self):
        atoms = self.prepare_md_dislocation()
        self.gn_md_tensile("Nb.eam.alloy.webarchive",
                           "Nb")
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

    def prepare_md_dislocation(self):
        file_list = glob.glob("./cfg/W.398.cfg")
        print file_list[-1]
        atoms = ase.io.read(file_list[-1],
                            format='cfg')
        self.write_lmp_config_data(atoms, file_name="relaxed.txt")
        return atoms

    def cal_non_periodic_screw(self):
        e1 = np.array([1.,   1., -2.])
        e2 = np.array([-1.,  1.,  0.])
        e3 = np.array([1.,   1.,  1.])

        self.set_lattce_constant(self._alat)
        self.set_element('Nb')
        atoms = self.set_bcc_convention(
            [e1, e2, e3], (80, 50, 12))  # z periodic 12

        atoms = self.intro_single_screw_atoms(atoms)
        self.write_lmp_config_data(atoms)

        #  if not os.path.isdir("restart"):
        #  os.mkdir("restart")
        #  os.mkdir("cfg")
        return

    # calculate the nano particle  #
    def cal_non_periodic_screw_xdislocation(self):
        e1 = np.array([-1.,   1.,  0])
        e2 = np.array([-1.,  -1.,  2])
        e3 = np.array([1.,   1.,  1.])

        self.set_lattce_constant(self._alat)
        self.set_element('W')

        # atoms = self.set_bcc_convention([e1, e2, e3], (60, 40, 3))  # z peri
        # 18
        atoms = self.set_bcc_convention(
            [e1, e2, e3], (60, 40, 80))  # z peri 18

        atoms = self.intro_single_screw_atoms(atoms)

        #  atoms = self.intro_dipole_screw_with_image_atoms(atoms);

        atoms = self.cut_y_normal_atoms(atoms)

        self.write_lmp_config_data(atoms)

        if not os.path.isdir("restart"):
            os.mkdir("restart")
            os.mkdir("cfg")

        #  self.gn_md_nano_tensile(
            #  potential_file = "Nb.eam.alloy.webarchive",
            #  element = "Nb",
            #  temperature = "20",
            #  deform_direction = "xz",
            #  deform_rate = "5e-7")

        self.gn_md_minimize_cfg("lmp_init.txt",
                                "W.set.txt",    # "Nb.eam.alloy.webarchive",
                                "W")
        os.system("rm cfg/* ; mpirun -n 4 lmp_mpi -in in.minimize")
        # write pbs
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

            if not os.path.isdir("restart"):
                os.mkdir("restart")
                os.mkdir("force_cfg")

            self.set_nnodes(1)
            self.set_wall_time(90)
            self.set_job_title("W-%s-T%d" % (stress, temp))
            self.set_main_job(
                "mpirun lmp_linux -in in.md_addforce  > screen.log ")
            self.write_pbs()
            self.gn_md_add_force(temp, stress)
            os.chdir(self.root_dir)
        return

    def cal_cu3Au_dis(self):
        atoms = ase.io.read("./custom/perf.00000.dump",
                            format='lammps-dump')

        # delete some Cu atoms to create missing rows #
        unity = 3.63902984582228 * np.sqrt(2.) / 2.
        index_list = []
        for atom in atoms:
            if atom.position[1] < (8 * unity) and atom.symbol is 'H':
                index_list.append(atom.index)
        del atoms[index_list]

        # introduce edge dislocation #
        cell = atoms.get_cell()
        zhigh = cell[2, 2]
        center = [0.5 * cell[0, 0], 10. * unity, 0.0]
        atoms = self.intro_single_edge_atoms(atoms,
                                             center, 1,
                                             [0, 1, 2])

        #  atoms = self.intro_edge_nuclei(atoms, center, -1,
        #  [0, 1, 2],
        #  [0.0*zhigh,
        #  zhigh])

        # cut extra half plane #
        #  atoms = self.cut_x_normal_atoms(atoms, self._alat,
        #  0, np.sqrt(2)/3.)

        ase.io.write("ase.cfg", atoms, format='cfg')
        system, elements = am.convert.ase_Atoms.load(atoms)
        lmp.atom_data.dump(system, "lmp_init.txt")
        return

    def cal_edge_nucleate_neb(self):
        return

    def cal_ani_dis(self):
        self.intro_ani_edge_fcc()
        return

    def loop(self):
        for i in range(1, 24):
            index = 24 - i
            os.system("cp coord_%d coord_%d" % (index, index + 1))
        return


if __name__ == "__main__":
    usage = "usage:%prog [options] arg1 [options] arg2"
    parser = OptionParser(usage=usage)
    parser.add_option("-t", "--mtype", action="store",
                      type="string", dest="mtype", help="",
                      default="prp_r")

    (options, args) = parser.parse_args()

    N = md_dislocation()

    # static calculate the kink pair migration barrier by neb method #
    if options.mtype == "nebkink":
        # nano particle only Peierdic along z #
        N.cal_kink_pair_neb_pre()

    if options.mtype == "kink":
        N.intro_kink_pair()

    if options.mtype == "modify":
        # ouptut the rotated dipole configuration #
        N.modified_cal_disp_dipo()

    if options.mtype.lower() == "hcpedge":
        N.hcp_edge_dislocation()

    if options.mtype.lower() == 'bccedge':
        N.cal_single_edge_dislocations()

    if options.mtype.lower() == 'bccscerw':
        N.cal_single_screw_dislocatoins()

    if options.mtype.lower() == 'nye':
        N.cal_nye()

    if options.mtype.lower() == 'cuau':
        N.cal_cu3Au_dis()

    if options.mtype.lower() == 'sub':
        N.loop_sub_jobs()

    if options.mtype.lower() == 'ani':
        N.cal_ani_dis()

    #  calculate the mobility of screw dislocation by MD
    #  N.cal_non_periodic_screw_xdislocation()
    #  N.static_add_stress();
    #  prepare md input files
    #  N.prepare_md_dislocation()
    #  prepare md input files
    #  N.loop_write_pbs()

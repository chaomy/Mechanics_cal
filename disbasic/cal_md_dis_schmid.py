#!/usr/bin/env python
# encoding: utf-8
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2018-02-06 14:17:35
# @Last Modified by:   chaomy
# @Last Modified time: 2018-03-18 22:11:32


from optparse import OptionParser
import get_data
import ase
import ase.io
import copy
import os
import numpy as np
import glob
import md_pot_data
import atomman as am
import atomman.lammps as lmp
import atomman.unitconvert as uc


class cal_bcc_schmid(get_data.get_data):

    def __init__(self, pot=md_pot_data.md_pot.Nb_eam):
        get_data.get_data.__init__(self)
        self.pot = self.load_data("../BASICS/pot.dat")
        self._range = (0, 25)
        #  self._delta = 0.02
        self._delta = 0.002   # totla is  5  percent
        #  self._delta = 0.0025   # totla is  5  percent

    # calculate elastic constant of screw dislocation
    def cal_screw_const(self, tag='intro'):
        axes = np.array([[1, 1, -2],
                         [-1, 1, 0],
                         [1, 1, 1]])

        alat = uc.set_in_units(self.pot['lattice'], 'angstrom')
        C11 = uc.set_in_units(self.pot['c11'], 'GPa')
        C12 = uc.set_in_units(self.pot['c12'], 'GPa')
        C44 = uc.set_in_units(self.pot['c44'], 'GPa')

        c = am.ElasticConstants(C11=C11, C12=C12, C44=C44)
        burgers = alat / 2 * np.array([1., 1., 1.])

        stroh = am.defect.Stroh(c, burgers, axes=axes)
        print "K tensor", stroh.K_tensor
        print "K (biKijbj)", stroh.K_coeff, "eV/A"
        print "pre-ln alpha = biKijbj/4pi", stroh.preln, "ev/A"

    #     drv.make_screw_plate(size=[80, 120, 2], rad=[200, 230],
    #                          move=[0., 0., 0.], tag='[211]',
    #                          filename="lmp_init.txt")

    def make_screw_plate(self, size=[40, 60, 2], rad=[100, 115],
                         move=[0., 0., 0.], tag='[211]', filename="lmp_init.txt", opt=None):

        alat = uc.set_in_units(self.pot['lattice'], 'angstrom')
        C11 = uc.set_in_units(self.pot['c11'], 'GPa')
        C12 = uc.set_in_units(self.pot['c12'], 'GPa')
        C44 = uc.set_in_units(self.pot['c44'], 'GPa')

        if tag == '[110]':
            axes = np.array([[-1, 0, 1],
                             [1, -2, 1],
                             [1, 1, 1]])
            unitx = alat * np.sqrt(2)
            unity = alat * np.sqrt(6)
            sizex = size[1]
            sizey = size[0]
            shftx = 0.5 * alat * np.sqrt(2) / 2.
            shfty = 0.5 * alat * np.sqrt(6.) / 2.

        elif tag == '[211]':
            axes = np.array([[1, -2, 1],
                             [1, 0, -1],
                             [1, 1, 1]])
            unitx = alat * np.sqrt(6)
            unity = alat * np.sqrt(2)
            sizex = size[0]
            sizey = size[1]
            shftx = -0.5 * alat * np.sqrt(6.) / 3.
            shfty = -1. / 3. * alat * np.sqrt(2) / 2.

        sizez = size[2]
        c = am.ElasticConstants(C11=C11, C12=C12, C44=C44)
        burgers = 0.5 * alat * np.array([1., 1., 1.])

        # initializing a new Stroh object using the data
        stroh = am.defect.Stroh(c, burgers, axes=axes)

        # monopole system
        box = am.Box(a=alat, b=alat, c=alat)
        atoms = am.Atoms(natoms=2, prop={'atype': 2, 'pos': [[0.0, 0.0, 0.0],
                                                             [0.5, 0.5, 0.5]]})

        ucell = am.System(atoms=atoms, box=box, scale=True)
        system = am.rotate_cubic(ucell, axes)

        center = [0.5 * unitx * sizex, 0.5 * unity * sizey]
        shift = np.array([shftx, shfty, 0.00000000000000])
        new_pos = system.atoms_prop(key='pos') + shift
        system.atoms_prop(key='pos', value=new_pos)

        system.supersize((0, sizex), (0, sizey), (0, sizez))

        # to make a plate #
        radius = rad[0]
        radiusout = rad[1]
        radius2 = radius * radius
        radiusout2 = radiusout * radiusout

        elements = []
        for atom in system.atoms:
            elements.append('Nb')

        ase_atoms = am.convert.ase_Atoms.dump(system, elements)
        pos = ase_atoms.get_positions()
        delindex = []

        for i in range(len(pos)):
            atom = ase_atoms[i]
            dx = pos[i, 0] - center[0]
            dy = pos[i, 1] - center[1]
            r = dx * dx + dy * dy

            if r > radiusout2:
                delindex.append(atom.index)
            if r < radius2:
                atom.symbol = 'W'

        del ase_atoms[delindex]
        (system, elements) = am.convert.ase_Atoms.load(ase_atoms)

        ############################################################
        # use neb, it's to generate init configuration
        ############################################################
        if opt == 'neb':
            system_init = copy.deepcopy(system)

            shift = np.array([-0.50000000000,
                              -0.500000000000,
                              0.00000000000000])
            new_pos = system_init.atoms_prop(key='pos',
                                             scale=True) + shift
            system_init.atoms_prop(key='pos',
                                   value=new_pos,
                                   scale=True)

            disp = stroh.displacement(system_init.atoms_prop(key='pos'))
            system_init.atoms_prop(key='pos',
                                   value=system_init.atoms_prop(key='pos') + disp)

            shift = np.array([0.50000, 0.50000, 0.00000])
            new_pos = system_init.atoms_prop(key='pos', scale=True) + shift
            system_init.atoms_prop(key='pos', value=new_pos, scale=True)

            # for lammps read structure
            lmp.atom_data.dump(system_init, "init.txt")

        # for dd map plot
        ase.io.write("lmp_perf.cfg", images=ase_atoms, format='cfg')
        lmp.atom_data.dump(system, "lmp_perf.txt")

        shift = np.array([-0.50000, -0.500, 0.0000])
        new_pos = system.atoms_prop(key='pos', scale=True) + shift
        system.atoms_prop(key='pos',
                          value=new_pos,
                          scale=True)

        new_pos = system.atoms_prop(key='pos') + move
        system.atoms_prop(key='pos', value=new_pos)

        disp = stroh.displacement(system.atoms_prop(key='pos'))
        system.atoms_prop(key='pos', value=system.atoms_prop(key='pos') + disp)

        new_pos = system.atoms_prop(key='pos') - move
        system.atoms_prop(key='pos', value=new_pos)

        shift = np.array([0.500000, 0.500000, 0.000000])
        new_pos = system.atoms_prop(key='pos', scale=True) + shift
        system.atoms_prop(key='pos', value=new_pos, scale=True)

        # for lammps read structure
        lmp.atom_data.dump(system, filename)

    #############################################################
    # change orient such that dislocation line is along x
    #############################################################
    def change_orient(self, ase_atoms=None, outfile="new_dis.cfg"):
        if ase_atoms is None:
            files = glob.glob("out/*")
            ase_atoms = ase.io.read(files[-1], format='lammps-dump')

        # we add the shear on the relaxed structure #
        cell = ase_atoms.get_cell()
        pos = ase_atoms.get_positions()
        print pos[:, 0]

        #  ase.io.write("old_dis.cfg", images=ase_atoms,
        #  format='cfg')

        #  here we rotate 90 degrees along the y axis
        #  and switch the x and z coordinates
        newcell = copy.deepcopy(cell)
        newcell[0, 0] = cell[2, 2]
        newcell[2, 2] = cell[0, 0]

        newpos = np.zeros([len(pos), 3])
        newpos[:, 0] = pos[:, 2]
        newpos[:, 1] = pos[:, 1]
        newpos[:, 2] = pos[:, 0]

        ase_atoms.set_cell(newcell)
        ase_atoms.set_positions(newpos)

        system, elements = am.convert.ase_Atoms.load(ase_atoms)
        lmp.atom_data.dump(system, "lmp_init_along_x.txt")
        return ase_atoms

    def add_shear(self, atoms, delta):
        #  system, elements = am.convert.ase_Atoms.load(atoms)
        #  lmp.atom_data.dump(system, "lmp_tmp.txt")
        #  os.system("mv lmp_tmp.txt {}/lmp_init_0.txt".format(dirname))

        cell = np.mat(atoms.get_cell())
        pos = np.mat(atoms.get_positions())
        # add along  xz  which is  kai = 90
        #  strain = np.mat([[1.0, 0.0, 0.0],
        #  [0.0, 1.0, 0.0],
        #  [delta, 0.0, 1.0]])

        strain = np.mat([[1.0, 0.0, 0.0],
                         [delta, 1.0, 0.0],
                         [0.0, 0.0, 1.0]])

        new_cell = strain * cell
        pos = pos * strain
        atoms.set_cell(new_cell)
        atoms.set_positions(pos)
        return atoms

    def prep_mrss(self):
        files = glob.glob("out/*")
        self._range = (25, 40)
        self._delta = 0.002   # totla is  5  percent

        for i in range(self._range[0], self._range[1]):
            delta = self._delta * i

            # add espsilon xz direction  (z direction shear along x)
            dirname = 'dir-%.4f' % (delta)
            if not os.path.isdir(dirname):
                os.mkdir(dirname)

            dis_atoms = ase.io.read(files[-1], format='lammps-dump')
            dis_atoms.wrap(pbc=[1, 1, 1])

            perf_atoms = ase.io.read('lmp_perf.cfg', format='cfg')
            perf_atoms.wrap(pbc=[1, 1, 1])

            # change orient #
            dis_atoms = self.change_orient(dis_atoms)
            perf_atoms = self.change_orient(perf_atoms)

            # add shear #
            dis_atoms = self.add_shear(dis_atoms, delta)
            perf_atoms = self.add_shear(perf_atoms, delta)

            # output #
            system, elements = am.convert.ase_Atoms.load(dis_atoms)
            lmp.atom_data.dump(system, "lmp_tmp.txt")
            ase.io.write("lmp_tmp.cfg", images=perf_atoms, format='cfg')

            os.system("mv lmp_tmp.cfg {}/lmp_perf.cfg".format(dirname))
            os.system("mv lmp_tmp.txt {}/lmp_init.txt".format(dirname))
            os.system("cp in.minimize {}".format(dirname))

    def mvdir(self):
        self._range = (0, 30)
        self._delta = 0.0025   # totla is  5  percent
        for i in range(self._range[0], self._range[1]):
            delta = self._delta * i
            dirname = 'dir-%.4f' % (delta)
            newname = 'dir-%.4f' % (delta)
            os.system("mv %s  %s " % (dirname, newname))

    def mvdump(self):
        for i in range(0, 20):
            delta = 2 * self._delta * i
            dirname = 'dir-%.4f' % (delta)
            filelist = glob.glob("%s/bcc.*.dump" % (dirname))
            print filelist[-1]
            os.system("cp %s dump_%03d" % (filelist[-1], i))
            #  ase_atoms = ase.io.read(filelist[-1], format='lammps-dump')
            #  cfgname = "dis_z_%.4f" % (delta)
            #  ase.io.write(cfgname, images=ase_atoms,
            #  format='cfg')

    def draw_ddmap(self):
        #  for i in range(self._range[0], self._range[1]):
        self._delta = 0.004
        for i in range(0, 10):
            #  delta = self._delta * i
            delta = self._delta * i
            dirname = 'dir-%.4f' % (delta)
            filelist = glob.glob("%s/bcc.*.dump" % (dirname))
            #  dumpname = "dump%.4f" % (delta)
            dumpname = filelist[-1]
            ase_atoms = ase.io.read(dumpname, format='lammps-dump')

            # change the orientation back
            strain = np.mat([[1.0, 0.0, 0.0],
                             [delta, 1.0, 0.0],
                             [0.0, 0.0, 1.0]])

            invstrain = np.linalg.inv(strain)

            cell = ase_atoms.get_cell()
            print "before add strain", cell

            cell = invstrain * cell
            print "after add strain ", cell

            pos = np.mat(ase_atoms.get_positions())
            pos = pos * invstrain

            newcell = copy.deepcopy(cell)
            newcell[0, 0] = cell[2, 2]
            newcell[2, 2] = cell[0, 0]

            newpos = copy.deepcopy(pos)

            newpos[:, 0] = pos[:, 2]
            newpos[:, 1] = pos[:, 1]
            newpos[:, 2] = pos[:, 0]

            ase_atoms.set_cell(newcell)
            ase_atoms.set_positions(newpos)

            cfgname = "dis_z_%.4f" % (delta)
            print cfgname
            ase.io.write(cfgname, images=ase_atoms,
                         format='cfg')

    def run_lmp(self, tag='lmp'):
        #  for i in range(self._range[0], self._range[1]):
        for i in range(10, 20):
            delta = 2 * self._delta * i
            dirname = 'dir-%.4f' % (delta)
            os.chdir(dirname)
            if tag == 'lmp':
                os.system("cp ../in.minimize .")
                os.system("lmp_mpi -i in.minimize")

            elif tag == 'show':
                filelist = glob.glob("bcc.*")
                print filelist[-1]
                os.system("cp  %s ../bcc.%03d" % (filelist[-1], i))
            os.chdir(os.pardir)

    def static_shear(self):
        # add shear based on relaxed configurations
        # the orientation must be change to x along [111] before
        delta = 0.004
        lastdir = 'out'
        for i in range(0, 10):
            # add espsilon xz direction  (z direction shear along x)
            dirname = 'dir-%.4f' % (i * delta)
            if not os.path.isdir(dirname):
                os.mkdir(dirname)

            files = glob.glob("%s/bcc.*" % (lastdir))
            print files
            dis_atoms = ase.io.read(files[-1], format='lammps-dump')
            dis_atoms.wrap(pbc=[1, 1, 1])

            # add shear #
            dis_atoms = self.add_shear(dis_atoms, delta)

            # output #
            system, elements = am.convert.ase_Atoms.load(dis_atoms)
            lmp.atom_data.dump(system, "lmp_tmp.txt")

            os.system("mv lmp_tmp.txt {}/lmp_init.txt".format(dirname))
            os.system("cp in.minimize {}".format(dirname))

            # run #
            os.chdir(dirname)
            os.system("lmp_mpi -i in.minimize")
            os.system("rm in.minimize")
            os.chdir(os.pardir)
            lastdir = dirname

    def wrap_auto(self):
        self.make_screw_plate()
        self.prep_mrss()
        self.run_lmp('lmp')


if __name__ == '__main__':
    usage = "usage:%prog [options] arg1 [options] arg2"
    parser = OptionParser(usage=usage)
    parser.add_option("-t", "--mtype", action="store",
                      type="string", dest="mtype")
    parser.add_option('-p', "--param", action="store",
                      type='string', dest="fargs")

    (options, args) = parser.parse_args()
    drv = cal_bcc_schmid()

    dispatcher = {'plate': drv.make_screw_plate,
                  'change': drv.change_orient,
                  'mrss': drv.prep_mrss,
                  'run': drv.run_lmp,
                  'ddmap': drv.draw_ddmap,
                  'static': drv.static_shear,
                  'const': drv.cal_screw_const}

    if options.fargs is not None:
        dispatcher[options.mtype.lower()](options.fargs)
    else:
        dispatcher[options.mtype.lower()]()

    # if options.mtype == 'plate':
    #     drv.make_screw_plate(size=[80, 120, 2], rad=[200, 230],
    #                          move=[0., 0., 0.], tag='[211]',
    #                          filename="lmp_init.txt")

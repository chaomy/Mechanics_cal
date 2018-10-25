#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2018-02-06 14:17:35
# @Last Modified by:   chaomy
# @Last Modified time: 2018-08-30 13:13:39


import ase
import ase.io
import copy
import os
import numpy as np
import glob
import atomman as am
import atomman.lammps as lmp
import atomman.unitconvert as uc
from utils import stroh_solve


# make plate -> change orient -> add shear


class cal_bcc_schmid(object):

    def __init__(self):
        self._range = (0, 25)
        self._delta = 0.002   # total is  5  percent

    # calculate elastic constant of screw dislocation
    def cal_screw_const(self, tag='intro'):
        axes = np.array([[1, 1, -2], [-1, 1, 0], [1, 1, 1]])

        alat = uc.set_in_units(self.pot['lattice'], 'angstrom')
        C11 = uc.set_in_units(self.pot['c11'], 'GPa')
        C12 = uc.set_in_units(self.pot['c12'], 'GPa')
        C44 = uc.set_in_units(self.pot['c44'], 'GPa')

        c = am.ElasticConstants(C11=C11, C12=C12, C44=C44)
        burgers = alat / 2 * np.array([1., 1., 1.])

        stroh = am.defect.Stroh(c, burgers, axes=axes)
        print("K tensor", stroh.K_tensor)
        print("K (biKijbj)", stroh.K_coeff, "eV/A")
        print("pre-ln alpha = biKijbj/4pi", stroh.preln, "ev/A")

    # def make_screw_plate(self, size=[40, 60, 3], rad=[90, 100], move=[0.,
    # 0., 0.], filename="lmp_init.txt", opt=None):
    def make_screw_plate(self, size=[70, 90, 3], rad=[150, 160], move=[0., 0., 0.], filename="lmp_init.txt", opt=None):
        e1 = [1, -2, 1]
        e2 = [1, 0, -1]
        e3 = [1, 1, 1]
        axes = np.array([e1, e2, e3])
        alat = self.pot['lattice']
        c = am.ElasticConstants(C11=self.pot['c11'],
                                C12=self.pot['c12'],
                                C44=self.pot['c44'])
        ux = np.sqrt(6) / 3. * self.pot['lattice']
        uy = np.sqrt(2) / 2. * self.pot['lattice']

        burgers = 0.5 * alat * np.array([1., 1., 1.])
        stroh = am.defect.Stroh(c, burgers, axes=axes)

        atoms = self.set_bcc_convention(
            [e1, e2, e3], (size[0], size[1], size[2]))
        pos = atoms.get_positions()

        center = np.array([3 * 0.5 * size[0] * ux, size[1] * uy])

        delindex = []
        radius2 = rad[0] * rad[0]
        radiusout2 = rad[1] * rad[1]
        for i in range(len(pos)):
            atom = atoms[i]
            dx = pos[i, 0] - center[0]
            dy = pos[i, 1] - center[1]
            r = dx * dx + dy * dy
            if r > radiusout2:
                delindex.append(atom.index)
            if r < radius2:
                # atom.symbol = 'W'
                continue
        del atoms[delindex]

        pos = atoms.get_positions()
        discenter = np.array(
            [center[0] + 0.5 * ux, center[1] + 1. / 3. * uy, 0.0])
        shf = np.ones(pos.shape) * discenter
        print(pos - shf)
        d1 = stroh.displacement(pos - shf)

        # before displace generate perfect atoms
        self.write_lmp_config_data(atoms, 'lmp_perf.txt')
        atoms.set_positions(pos + np.real(d1))
        self.write_lmp_config_data(atoms, 'lmp_init.txt')
        np.savetxt("dis_center.txt", discenter)

    def make_screw_plate_old(self, size=[40, 60, 3], rad=[100, 115], move=[0., 0., 0.], filename="lmp_init.txt", opt=None):
        alat = uc.set_in_units(self.pot['lattice'], 'angstrom')
        C11 = uc.set_in_units(self.pot['c11'], 'GPa')
        C12 = uc.set_in_units(self.pot['c12'], 'GPa')
        C44 = uc.set_in_units(self.pot['c44'], 'GPa')

        axes = np.array([[1, -2, 1], [1, 0, -1], [1, 1, 1]])

        unitx = alat * np.sqrt(6)
        unity = alat * np.sqrt(2)
        sizex = size[0]
        sizey = size[1]

        sizez = size[2]
        c = am.ElasticConstants(C11=C11, C12=C12, C44=C44)
        burgers = 0.5 * alat * -np.array([1., 1., 1.])

        # initializing a new Stroh object using the data
        stroh = am.defect.Stroh(c, burgers, axes=axes)

        # monopole system
        box = am.Box(a=alat, b=alat, c=alat)
        atoms = am.Atoms(natoms=2, prop={'atype': 2, 'pos': [
                         [0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]})
        ucell = am.System(atoms=atoms, box=box, scale=True)
        system = am.rotate_cubic(ucell, axes)

        # shftx = 0.5 * alat * np.sqrt(6.) / 3.
        shftx = 0.0
        shfty = -1. / 3. * alat * np.sqrt(2) / 2.
        # shfty = 2. / 3. * alat * np.sqrt(2) / 2.

        center = [0.5 * unitx * sizex, 0.5 * unity * sizey]
        new_pos = system.atoms_prop(key='pos') + np.array([shftx, shfty, 0.0])
        system.atoms_prop(key='pos', value=new_pos)
        system.supersize((0, sizex), (0, sizey), (0, sizez))

        # to make a plate #
        radius2 = rad[0] * rad[0]
        radiusout2 = rad[1] * rad[1]

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

        # use neb, it's to generate init configuration
        if opt in ['neb']:
            system_init = copy.deepcopy(system)

            shift = np.array([-0.5, -0.5, 0.0])
            new_pos = system_init.atoms_prop(key='pos', scale=True) + shift
            system_init.atoms_prop(key='pos', value=new_pos, scale=True)

            disp = stroh.displacement(system_init.atoms_prop(key='pos'))
            system_init.atoms_prop(
                key='pos', value=system_init.atoms_prop(key='pos') + disp)

            shift = np.array([0.5, 0.50, 0.0])
            new_pos = system_init.atoms_prop(key='pos', scale=True) + shift
            system_init.atoms_prop(key='pos', value=new_pos, scale=True)

            # for lammps read structure
            lmp.atom_data.dump(system_init, "init.txt")

        # for dd map plot
        ase.io.write("lmp_perf.cfg", images=ase_atoms, format='cfg')
        lmp.atom_data.dump(system, "lmp_perf.txt")

        shift = np.array([-0.5, -0.5, 0.0])
        new_pos = system.atoms_prop(key='pos', scale=True) + shift
        system.atoms_prop(key='pos', value=new_pos, scale=True)

        new_pos = system.atoms_prop(key='pos') + move
        system.atoms_prop(key='pos', value=new_pos)

        disp = stroh.displacement(system.atoms_prop(key='pos'))

        # pull
        pull = False
        if pull is True:
            core_rows = [disp[:, 2].argsort()[-3:][::-1]]
            print(disp[core_rows])
            exclus = np.arange(len(disp), dtype=int)
            unitburger = np.mean(disp[core_rows][:, 2])
            print(unitburger)
            exclus = np.delete(exclus, core_rows)
            disp[core_rows] -= 1. / 3. * unitburger
            # disp[exclus] -= 1. / 3. * unitburger

        system.atoms_prop(key='pos', value=system.atoms_prop(key='pos') + disp)

        new_pos = system.atoms_prop(key='pos') - move
        system.atoms_prop(key='pos', value=new_pos)

        shift = np.array([0.500000, 0.500000, 0.000000])
        new_pos = system.atoms_prop(key='pos', scale=True) + shift
        system.atoms_prop(key='pos', value=new_pos, scale=True)

        new_pos = system.atoms_prop(key='pos', scale=False)

        # for lammps read structure
        lmp.atom_data.dump(system, filename)

        dis_atoms = am.convert.ase_Atoms.dump(system, elements)
        return (ase_atoms, dis_atoms)

    # change orient such that dislocation line is along x
    def change_orient(self, ase_atoms=None, outfile="new_dis.cfg"):
        if ase_atoms is None:
            files = glob.glob("out/*")
            ase_atoms = ase.io.read(files[-1], format='lammps-dump')

        # we add the shear on the relaxed structure #
        cell = ase_atoms.get_cell()
        pos = ase_atoms.get_positions()
        print(pos[:, 0])

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

        #  add along  xz  which is  kai = 90
        #  strain = np.mat([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [delta, 0.0, 1.0]])

        strain = np.mat([[1.0, 0.0, 0.0], [delta, 1.0, 0.0], [0.0, 0.0, 1.0]])
        new_cell = strain * np.mat(atoms.get_cell())
        pos = np.mat(atoms.get_positions()) * strain
        atoms.set_cell(new_cell)
        atoms.set_positions(pos)
        return atoms

    def add_shear_yz(self, atoms, delta):
        # add shear
        strain = np.mat([[1.0, 0.0, 0.0],
                         [0.0, 1.0, 0.0],
                         [0.0, delta, 1.0]])
        new_cell = strain * np.mat(atoms.get_cell())
        pos = np.mat(atoms.get_positions()) * strain
        atoms.set_cell(new_cell)
        atoms.set_positions(pos)
        return atoms

    def prep_mrss(self):
        self._range = (1, 5)
        self._delta = 0.0001
        perf_atoms, dis_atoms = self.make_screw_plate()
        for i in range(self._range[0], self._range[1]):
            delta = self._delta * i
            mdir = "dir_{:04f}".format(delta)
            self.mymkdir(mdir)
            dis_atoms = self.add_shear_yz(dis_atoms.copy(), delta)
            perf_atoms = self.add_shear_yz(perf_atoms.copy(), delta)
            self.write_lmp_config_data(
                dis_atoms, "{}/lmp_init.txt".format(mdir))
            self.write_lmp_config_data(
                perf_atoms, "{}/lmp_perf.txt".format(mdir))
            os.system("cp in.minimize {}".format(mdir))

    def prep_mrss_change_orient(self):
        files = glob.glob("out/*")
        self._range = (10, 12)
        self._delta = 0.002   # total is  5  percent

        for i in range(self._range[0], self._range[1]):
            delta = self._delta * i

            # add espsilon xz direction  (z direction shear along x)
            dirname = 'dir-%.4f' % (delta)
            self.mymkdir(dirname)

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
            print(filelist[-1])
            os.system("cp %s dump_%03d" % (filelist[-1], i))
            #  ase_atoms = ase.io.read(filelist[-1], format='lammps-dump')
            #  cfgname = "dis_z_%.4f" % (delta)
            #  ase.io.write(cfgname, images=ase_atoms,
            #  format='cfg')

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
                print(filelist[-1])
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
            self.mymkdir(dirname)

            files = glob.glob("%s/bcc.*" % (lastdir))
            print(files)
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


# if __name__ == '__main__':
#     usage = "usage:%prog [options] arg1 [options] arg2"
#     parser = OptionParser(usage=usage)
#     parser.add_option("-t", "--mtype", action="store",
#                       type="string", dest="mtype")
#     parser.add_option('-p', "--param", action="store",
#                       type='string', dest="fargs")

#     (options, args) = parser.parse_args()
#     drv = cal_bcc_schmid()

#     dispatcher = {'change': drv.change_orient,
#                   'run': drv.run_lmp,
#                   'static': drv.static_shear,
#                   'const': drv.cal_screw_const}

#     if options.fargs is not None:
#         dispatcher[options.mtype.lower()](options.fargs)
#     else:
#         dispatcher[options.mtype.lower()]()

    # if options.mtype == 'plate':
    #     drv.make_screw_plate(size=[80, 120, 2], rad=[200, 230],
    #                          move=[0., 0., 0.], tag='[211]',
    #                          filename="lmp_init.txt")

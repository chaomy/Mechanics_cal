#!/usr/bin/env python
# encoding: utf-8
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-07-05 08:12:30
# @Last Modified by:   chaomy
# @Last Modified time: 2018-05-01 22:35:40


import os
import ase
import numpy as np
import md_pot_data
import gn_config
import get_data
import gn_pbs
import gn_lmp_infile
import atomman as am
from optparse import OptionParser
from utils import Intro_vasp
from utils import cal_md_find_dis_core
from prec import cal_md_prec
from gb import cal_md_gb_hcp_dis
import cal_md_dis_dipole
import cal_md_dislocation_hcp
import cal_md_dislocation_fcc
import cal_md_dislocation_bcc
import cal_md_peierls_barrier
import cal_md_dis_schmid
import plt_drv


class md_dislocation(gn_config.gnStructure,
                     get_data.get_data, gn_pbs.gn_pbs,
                     plt_drv.plt_drv,
                     Intro_vasp.vasp_change_box,
                     gn_lmp_infile.gn_md_infile,
                     cal_md_peierls_barrier.cal_barrier,
                     cal_md_dislocation_hcp.md_dislocation_hcp,
                     cal_md_dislocation_bcc.md_dislocation_bcc,
                     cal_md_dislocation_fcc.md_dislocation_fcc,
                     cal_md_prec.md_prec,
                     cal_md_dis_dipole.cal_dis_dipole,
                     cal_md_gb_hcp_dis.gb_hcp_dis,
                     cal_md_dis_schmid.cal_bcc_schmid,
                     cal_md_find_dis_core.md_find_core):

    def __init__(self, pot=md_pot_data.md_pot.mg_kim):
        # self.pot = pot
        # self.pot = md_pot_data.md_pot.mg_Poco
        self.pot = self.load_data('../BASICS/pot.dat')
        # self.pot = self.load_data('../BASICS_MO/pot.dat')
        plt_drv.plt_drv.__init__(self)
        gn_config.gnStructure.__init__(self, self.pot)
        gn_pbs.gn_pbs.__init__(self)
        Intro_vasp.vasp_change_box.__init__(self)
        cal_md_prec.md_prec.__init__(self)
        cal_md_gb_hcp_dis.gb_hcp_dis.__init__(self)
        cal_md_peierls_barrier.cal_barrier.__init__(self)
        gn_lmp_infile.gn_md_infile.__init__(self)
        cal_md_dis_schmid.cal_bcc_schmid.__init__(self)
        cal_md_dis_dipole.cal_dis_dipole.__init__(self)
        cal_md_dislocation_bcc.md_dislocation_bcc.__init__(self)
        cal_md_dislocation_fcc.md_dislocation_fcc.__init__(self)
        cal_md_dislocation_hcp.md_dislocation_hcp.__init__(self)
        cal_md_find_dis_core.md_find_core.__init__(self)

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

        ase.io.write("lmp_init.cfg", atoms, "cfg")
        fname = "init.data"
        self.write_lmp_config_data(atoms, fname)
        self.gn_md_minimize_cfg("init.data", "./w_eam4.fs", "W")

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

        atoms = self.set_bcc_convention(directions=[[e1[0], e1[1], e1[2]],
                                                    [e2[0], e2[1], e2[2]],
                                                    [e3[0], e3[1], e3[2]]],
                                        size=(n, m, t))

        # perfect positions
        p = self.pot['lattice'] * self.bccneigh

        # atoms = self.intro_dipole_screw_atoms(atoms, self.pot['lattice'])
        atoms = self.intro_single_screw_atoms(atoms)

        ase.io.write("lmp_dis.cfg", atoms, "cfg")

        # system
        system, elements = am.convert.ase_Atoms.load(atoms)

        r = (3**0.5 / 2. + 1) / 2.
        system.nlist(self.pot['lattice'] * r)
        nye_rlt = am.defect.nye_tensor(system, p, axes=unit_cell)

        print(np.max(nye_rlt['strain'][:]))
        # for key, value in nye_rlt.iteritems():
        #     print key, value
        #     system.atoms_prop(key=key, value=value)

        # int_sum = np.empty(3)
        # for i in range(3):
        #     int_sum[i] = am.plot.interpolate_contour(system, 'Nye_tensor',
        #                                              index=[1, i],
        #                                              cmap='bwr')[0]
        # print('burgers vector estimate = ', int_sum)

        ############################################################
        # since lammmps can only use  xy xz yz
        # new x is old y ; new y is old z ; new z is old x
        # x  1  1 -2
        # y  1  1  1
        # z  1 -1  0
        ############################################################

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

    def cnv_poscar_lmp(self):
        atoms = ase.io.read("POSCAR", format='vasp')
        self.write_lmp_config_data(atoms)

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
                  'bccedge': drv.cal_single_edge_dislocations,
                  'bccscrew': drv.cal_single_screw_dislocations,
                  'nye': drv.cal_nye,
                  'cuau': drv.cal_cu3Au_dis,
                  'ani': drv.intro_ani_edge_fcc,
                  'static': drv.static_add_stress,
                  'hcpedge': drv.hcp_edge_dislocation,  # isotropic
                  'bedge': drv.build_edge_basal_hcp,  # hcp basal
                  'bscrew': drv.build_screw_basal_hcp,  # hcp basal
                  'thermo': drv.cal_thermo,
                  'prec': drv.make_prec,
                  'sprec': drv.make_screw_prec,
                  'd03': drv.buildd03small,
                  'hcp': drv.buildHCP,
                  'gb': drv.make_gb,
                  'peierls': drv.dipole_peierls_barrier,
                  'plate': drv.make_screw_plate,
                  'cnv': drv.cnv_poscar_lmp,
                  'dipole': drv.bcc_screw_dipole_configs_alongz,
                  'mrss': drv.prep_mrss,
                  'find': drv.cost_method_find_core}

    if options.fargs is not None:
        dispatcher[options.mtype.lower()](options.fargs)
    else:
        dispatcher[options.mtype.lower()]()

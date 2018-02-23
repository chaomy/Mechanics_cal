#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-08-07 20:35:25
# @Last Modified by:   chaomy
# @Last Modified time: 2018-02-20 19:28:26

from optparse import OptionParser
from numpy import array, max, min
from numpy import real
from numpy import ones
from numpy import shape
from numpy import sqrt
from numpy import dot
from numpy.linalg import norm
from numpy import sin, cos, mat
import numpy as np
import md_pot_data
import stroh_solve
import ase.io
import tool_elastic_constants
import ase.lattice
import cal_md_dislocation
import gn_config
import Intro_vasp
import cal_intro_iso_dis


class cal_dis_edge(gn_config.bcc,
                   cal_intro_iso_dis.cal_intro_iso_dis,
                   Intro_vasp.vasp_change_box):

    def __init__(self, pot=None):
        self.pot = md_pot_data.md_pot.W_eam
        cal_intro_iso_dis.cal_intro_iso_dis.__init__(self)
        gn_config.bcc.__init__(self, self.pot)
        Intro_vasp.vasp_change_box.__init__(self, self.pot)
        self.b2e = 74.2068309515  # x[110]
        self.e2v = 35.2643896828  # x[111]
        self.rot = 109.471220
        # self.rot = self.b2e

    def rotate_cell(self, atoms, phi):
        pos = atoms.get_positions()
        # rotate x, y
        # rotmat = np.mat([[cos(phi), sin(phi), 0.],
        #                  [-sin(phi), cos(phi), 0.],
        #                  [0., 0., 1.]])
        # rotate x, z
        rotmat = np.mat([[cos(phi), 0.0, sin(phi)],
                         [0.0, 1.0, 0.0],
                         [-sin(phi), 0., cos(phi)]])
        pos = (rotmat * pos.transpose()).transpose()
        atoms.set_positions(pos)
        return atoms

    def duplic(self, fname='1096_unit'):
        atoms = ase.io.read(fname, format='lammps-dump')
        mx = max(atoms.get_positions())
        mi = min(atoms.get_positions())
        for atom in atoms:
            if atom.position[2] > 24.5:
                if atom.symbol is not 'H':
                    atom.symbol = 'Mo'
                else:
                    atom.symbol = 'W'
        ase.io.write(filename='1096_symbol.aims', images=atoms, format='aims')

    def rotateto111(self, fname='1096_30x02y.aims'):
        # atoms = ase.io.read(fname, format='lammps-dump')
        atoms = ase.io.read(fname, format='aims')
        # self.rotate_cell(atoms, phi=self.rot)
        atoms.rotate(a=self.rot, v='y')
        # ase.io.write(filename='x111_30x02y.aims', images=atoms, format='aims')
        return atoms

    def gn_edge(self):
        e1 = [1., 1., 1.]
        e2 = [-1., -1., 2.]
        e3 = [1., -1., 0]

        # first rotate
        atoms = ase.io.read('1096_30x02y.aims', format='aims')
        atoms.rotate(a=self.rot, v='y')

        # atoms = ase.io.read('x111_30x02y.aims', format='aims')
        pos1 = array([272.249, 0.000, -338.769])
        pos2 = array([274.791, 0.000, -338.769 + np.sqrt(2) / 2.])
        pos = 0.5 * (pos1 + pos2)
        # center1 = [pos[0], 0.0, 10]
        # center2 = [pos[0], 0.0, -10]

        atoms = self._intro_single_edge_atoms(atoms, center=pos)
        # atoms = self._intro_dipole_edge_atoms(atoms, c1=center1, c2=center2)
        # atoms = self.cut_x_normal_atoms(atoms,self.pot['lattice'], ratio=1./sqrt(5.))

        # rotate back
        atoms.rotate(a=-self.rot, v='y')

        ase.io.write('single_edge_30x02y.aims', images=atoms, format='aims')
        self.write_lmp_config_data(atoms)

    def cal_trans_ang(self):
        self.calang(array([1, 1, 0]), array([1, 1, 1]))
        self.calang(array([0, 0, 1]), array([-1, -1, 2]))

    def calang(self, a, b):
        x = dot(a, b) / (norm(a) * norm(b))
        rotangle = np.rad2deg(np.arccos(x))
        print rotangle
        return rotangle


if __name__ == '__main__':
    usage = "usage:%prog [options] arg1 [options] arg2"
    parser = OptionParser(usage=usage)
    parser.add_option('-t', "--mtype", action="store",
                      type="string", dest="mtype", help="",
                      default="prp_r")
    parser.add_option('-p', "--param", action="store",
                      type='string', dest="fargs",
                      default=None)
    (options, args) = parser.parse_args()

    drv = cal_dis_edge()
    dispatcher = {'dup': drv.duplic,
                  'rotate1': drv.rotateto111,
                  'edge': drv.gn_edge}

    if options.fargs is not None:
        dispatcher[options.mtype.lower()](options.fargs)
    else:
        dispatcher[options.mtype.lower()]()
    # drv.read_tim()
    # drv.gn_edge()
    # drv.duplic()

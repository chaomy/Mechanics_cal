#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2018-04-25 23:53:39
# @Last Modified by:   chaomy
# @Last Modified time: 2018-04-25 23:53:57


import os
import ase
from optparse import OptionParser
import ase.lattice.hexagonal as Hexagonal
import ase.lattice.cubic as Cubic
import gn_config
import gn_pbs
import get_data
import Intro_vasp
import gn_lmp_infile


class md_precipitate(gn_config.bcc,
                     gn_config.fcc,
                     gn_config.hcp,
                     get_data.get_data,
                     gn_pbs.gn_pbs,
                     Intro_vasp.vasp_change_box,
                     gn_lmp_infile.gn_md_infile):

    def __init__(self,
                 element,
                 structure='hcp'):
        gn_pbs.gn_pbs.__init__(self)

        self._element = element
        self.lattice_a = 3.18421469227
        self.lattice_c = 5.18442354563

        gn_config.hcp.__init__(self, self.lattice_a, self.lattice_c)
        self._default_direction = [[1, 0, 0],
                                   [0, 1, 0],
                                   [0, 0, 1]]
        self.root_dir = os.getcwd()

    def bulk(self):
        size = (40, 40, 30)
        atoms = \
            Hexagonal.HexagonalClosedPacked(latticeconstant={'a': self.lattice_a,
                                                             'c': self.lattice_c},
                                            size=size,
                                            symbol='Mg',
                                            pbc=(1, 1, 1))

        self.write_lmp_config_data(atoms)

    def interface(self):
        l_size = (1, 1, 1)
        l_direction = self._default_direction
        atoms = Cubic.FaceCenteredCubic(directions=l_direction,
                                        latticeconstant=self.lattice_a,
                                        size=l_size,
                                        symbol='Mg',
                                        pbc=(1, 1, 1))

        print(atoms.get_cell())
        ase.io.write(filename="init.cfg",
                     images=atoms,
                     format='cfg')

usage = "usage:%prog [options] arg1 [options] arg2"
parser = OptionParser(usage=usage)
parser.add_option("-t", "--mtype", action="store",
                  type="string", dest="mtype", help="",
                  default="prp_r")
(options, args) = parser.parse_args()

if __name__ == "__main__":
    #  MgNd  D03  lattice   7.4662780330786   (Nd at fcc center)
    #  Mg    a  3.209065    c   5.196911
    #  3.2019267694893   5.1969105399
    Job = md_precipitate(element='Mg', structure='hcp')
    if options.mtype == "bulk":
        Job.bulk()

    if options.mtype.lower() == "inter":
        Job.interface()

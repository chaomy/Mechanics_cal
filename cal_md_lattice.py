#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: yang37
# @Date:   2017-06-21 18:42:47
# @Last Modified by:   chaomy
# @Last Modified time: 2017-09-24 10:44:41


import os
import ase.io
from optparse import OptionParser

try:
    import get_data
    import gn_config
    import gn_lmp_infile

except ImportError:
    print "error during import"


class cal_md_lattice(gn_config.bcc,
                     gn_config.fcc,
                     gn_config.hcp,
                     get_data.get_data,
                     gn_lmp_infile.gn_md_infile):

    def __init__(self,
                 in_structure='bcc',
                 lattice_constant=3.31109116862,
                 in_element='W',
                 in_potential='dummy.lammps.EAM'):

        get_data.get_data.__init__(self)
        gn_lmp_infile.gn_md_infile.__init__(self)

        self._lattice_element = in_element
        self._in_lattice_constant = lattice_constant
        self._structure = in_structure
        self._lattice_potential = in_potential

        if self._structure == 'bcc':
            gn_config.bcc.__init__(self)
        elif self._structure == 'fcc':
            gn_config.fcc.__init__(self)
        elif self._structure == 'hcp':
            gn_config.hcp.__init__(self)

        # set for gn_config
        self.set_element(self._lattice_element)
        self.set_lattce_constant(self._in_lattice_constant)
        self.root_dir = os.getcwd()
        return

    def gn_lattice_atoms(self):
        self.set_lattce_constant(self._in_lattice_constant)

        lattice_direction = [[1, 0, 0],
                             [0, 1, 0],
                             [0, 0, 1]]
        atoms = self.set_bcc_convention(in_direction=lattice_direction,
                                        in_size=(1, 1, 1))
        return atoms

    def gn_temp_atoms(self):
        #  self.set_lattce_constant(self._in_lattice_constant)
        self.set_lattce_constant(3.1652000015)

        lattice_direction = [[1, 0, 0],
                             [0, 1, 0],
                             [0, 0, 1]]
        atoms = self.set_bcc_convention(in_direction=lattice_direction,
                                        in_size=(3, 3, 3))

        self.write_lmp_config_data(atoms)
        return atoms

    def set_md_lattice_potential(self, in_potential):
        self._lattice_potential = in_potential
        return

    def cal_lmp_lattice(self):
        atoms = self.gn_lattice_atoms()
        lattice_dir = "dir-lattice"

        if not os.path.isdir(lattice_dir):
            os.mkdir(lattice_dir)
        os.chdir(lattice_dir)

        self.write_lmp_config_data(atoms)
        self.gn_md_lattice("lmp_init.txt",
                           self._lattice_potential,
                           self._lattice_element)

        os.system("cp ../%s  ." % (self._lattice_potential))
        os.system("lmp_mpi -in  in.lattice")
        lattice = self.get_lmp_lattice()

        os.chdir(self.root_dir)
        return lattice

    def cal_temp_lattice(self):
        temp = 300.0
        atoms = self.gn_temp_atoms()
        lattice_dir = "dir-lattice-%4.3f" % (temp)

        if not os.path.isdir(lattice_dir):
            os.mkdir(lattice_dir)
        os.chdir(lattice_dir)

        self.write_lmp_config_data(atoms)
        self.gn_md_temp_lattice("lmp_init.txt",
                                temp,
                                self._lattice_potential,
                                self._lattice_element)

        os.system("cp ../%s  ." % (self._lattice_potential))
        os.system("lmp_mpi -in in.lattice")
        lattice = self.get_lmp_lattice()
        os.chdir(self.root_dir)
        lattice = float(lattice)
        return lattice / 5.0

    def run_lmp_lattice(self, loc_dir):
        os.chdir(loc_dir)
        os.system("lmp_linux -in in.lattice")
        os.chdir(self.root_dir)
        return

    def hcp_lattice(self):
        self.set_hcp_lattice_constant(3.1742, 5.18)
        self.set_hcp_direction()
        atoms = self.set_hcp_convention((1, 1, 1))

        ase.io.write("hcp.txt", images=atoms, format='cfg')
        self.write_lmp_config_data(atoms)

        print atoms.get_cell()
        return


usage = "usage:%prog [options] arg1 [options] arg2"
parser = OptionParser(usage=usage)
parser.add_option("-f", "--mfile",
                  action="store",
                  type="string",
                  dest="mfile",
                  default="./dummy.config.pair")

parser.add_option("-t", "--mstruct",
                  action="store",
                  type="string",
                  dest="mstruct",
                  default=None)

(options, args) = parser.parse_args()
OriDir = os.getcwd()

if __name__ == '__main__':
    if options.mstruct.lower() == 'bcc':
        Job = cal_md_lattice(in_structure='bcc',
                             lattice_constant=3.165258,
                             in_element='W',
                             in_potential="./W.set.txt")
        #  print Job.cal_temp_lattice()
        Job.cal_lmp_lattice()

    elif options.mstruct.lower() == 'hcp':
        ### ./Al-Mg.eam.fs  
        ## a  3.18421469227388
        ## c  5.18442354562984
        Job = cal_md_lattice(in_structure='hcp',
                             lattice_constant=3.18,
                             in_element='Mg',
                             in_potential="./Al-Mg.eam.fs")
        Job.hcp_lattice()

    #  Job.gn_temp_atoms()
    #           Nb                      W(W.eam.fs);   W.set.txt
    #  0K       3.307                   3.165200       3.1648492
    #  150K
    #  300K     3.3107                  3.165258
    #  450K
    #  600K
    #  750K

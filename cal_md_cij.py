#!/usr/bin/env python
# encoding: utf-8

###################################################################
#
# File Name : ./cal_md_cij.py
#
###################################################################
#
# Purpose :  calculate md cij
#
# Creation Date :
# Last Modified : Thu Mar 30 23:56:50 2017
# Created By    : Chaoming Yang
#
###################################################################

import os
try:
    from optparse import OptionParser
    import get_data
    import gn_lmp_infile
    import gn_config
    import md_pot_data

except ImportError:
    print "error during import"


class cal_md_cij(get_data.get_data,
                 gn_lmp_infile.gn_md_infile,
                 gn_config.bcc,
                 gn_config.fcc,
                 gn_config.hcp):

    def __init__(self, pot=md_pot_data.md_pot.Nb_eam):
        get_data.get_data.__init__(self)
        self.pot = pot
        self._cij_element = self.pot['element']
        self._structure = self.pot['structure']
        self._lattice_constant = self.pot['lattice']
        self._cij_potential = self.pot['file']

        gn_lmp_infile.gn_md_infile.__init__(self, self.pot)

        if self._structure == 'bcc':
            gn_config.bcc.__init__(self, self.pot)
        elif self._structure == 'fcc':
            gn_config.fcc.__init__(self, self.pot)
        elif self._structure == 'hcp':
            gn_config.hcp.__init__(self, self.pot)

        self.set_lattce_constant(self._lattice_constant)
        self.root_dir = os.getcwd()
        return

    def gn_temp_atoms(self):
        self.set_lattce_constant(self._lattice_constant)
        lattice_direction = [[1, 0, 0],
                             [0, 1, 0],
                             [0, 0, 1]]
        atoms = self.set_bcc_convention(in_direction=lattice_direction,
                                        in_size=(3, 3, 3))
        self.write_lmp_config_data(atoms)
        return atoms

    def cal_md_shear_modulus(self):
        temp = 300.0
        atoms = self.gn_temp_atoms()
        miu_dir = "dir-shearModulus-%4.3f" % (temp)

        self.mymkdir(miu_dir)
        os.chdir(miu_dir)

        self.write_lmp_config_data(atoms)

        self.gn_md_shear_lattice("lmp_init.txt",
                                 temp,
                                 self._cij_potential,
                                 self._cij_element)

        os.system("cp ../%s  ." % (self._cij_potential))
        os.system("mpirun -n 8 lmp_mpi -in in.shear")
        os.chdir(self.root_dir)
        return

    def cal_elastic_constant(self):
        self.pot = self.load_data('pot.dat')

        self.mymkdir("dir-cij")
        os.chdir("dir-cij")
        self.gn_md_cij()
        os.system("cp ../%s ." % (self.pot['file']))
        os.system("cp ../displace.mod .")
        os.system("lmp_mpi -in in.elastic > log.md_cij")
        elastic_constants = self.get_lmp_elastic_constants()
        os.chdir(self.root_dir)
        return elastic_constants


if __name__ == '__main__':
    usage = "usage:%prog [options] arg1 [options] arg2"
    parser = OptionParser(usage=usage)
    parser.add_option("-t", "--mtype",
                      action="store",
                      type="string",
                      dest="mtype",
                      help="",
                      default="read")
    (options, args) = parser.parse_args()
    drv = cal_md_cij()
    print drv.cal_elastic_constant()

    #  Job.cal_md_shear_modulus()
    #       miu 43.463638  #300K
    #       W(W.eam.fs)
    #       532.612621081359,
    #       205.015564010464,
    #       163.19875750812
    #       miu 182.26     #300K

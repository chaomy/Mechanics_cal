#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: yang37
# @Date:   2017-06-12 17:03:43
# @Last Modified by:   chaomy
# @Last Modified time: 2018-04-09 18:51:34


from optparse import OptionParser
import os
import get_data
import gn_lmp_infile
import gn_config
import md_pot_data


class cal_md_cij(get_data.get_data,
                 gn_lmp_infile.gn_md_infile,
                 gn_config.bcc):

    def __init__(self, pot=md_pot_data.md_pot.Nb_eam):
        get_data.get_data.__init__(self)
        self.pot = pot
        gn_lmp_infile.gn_md_infile.__init__(self, self.pot)
        gn_config.bcc.__init__(self, self.pot)

    def gn_temp_atoms(self):
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
        self.gn_md_shear_lattice("lmp_init.txt", temp,
                                 self.pot["file"],
                                 self.pot['element'])

        os.system("cp ../%s  ." % (self.pot["file"]))
        os.system("mpirun -n 8 lmp_mpi -in in.shear")
        os.chdir(os.pardir)

    def cal_elastic_constant(self):
        self.pot = self.load_data('pot.dat')

        self.mymkdir("dir-cij")
        os.chdir("dir-cij")
        self.gn_md_cij()
        os.system("cp ../%s ." % (self.pot['file']))
        os.system("cp ../displace.mod .")
        os.system("lmp_mpi -in in.elastic > log.md_cij")
        elastic_constants = self.get_lmp_elastic_constants()
        os.chdir(os.pardir)
        return elastic_constants


if __name__ == '__main__':
    usage = "usage:%prog [options] arg1 [options] arg2"
    parser = OptionParser(usage=usage)
    parser.add_option("-t", "--mtype",
                      action="store",
                      type="string",
                      dest="mtype")

    (options, args) = parser.parse_args()
    drv = cal_md_cij()
    print(drv.cal_elastic_constant())

    #  Job.cal_md_shear_modulus()
    #       miu 43.463638  #300K
    #       W(W.eam.fs)
    #       532.612621081359,
    #       205.015564010464,
    #       163.19875750812
    #       miu 182.26     #300K

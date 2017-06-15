#!/usr/bin/env python
# encoding: utf-8

###################################################################
#
# File Name : cal_qe_screw_dipole.py
#
###################################################################
#
# Purpose :
#
# Creation Date :
# Last Modified :
# Created By    : Chaoming Yang
#
###################################################################

try:
    import numpy as np
    import md_pot_data
    import ase
    import ase.io
    import os
    import gn_config
    import get_data
    import Intro_vasp
    import gn_pbs
    import gn_qe_inputs
    from optparse import OptionParser

except ImportError:
    print("error during import")


class qe_dislocation(get_data.get_data,
                     gn_qe_inputs.gn_qe_infile,
                     gn_pbs.gn_pbs,
                     gn_config.bcc,
                     Intro_vasp.vasp_change_box):

    def __init__(self):
        self.pot = md_pot_data.qe_pot.vca_W50Re50
        get_data.get_data.__init__(self)
        gn_qe_inputs.gn_qe_infile.__init__(self, self.pot)
        gn_pbs.gn_pbs.__init__(self)
        gn_config.bcc.__init__(self, self.pot)
        Intro_vasp.vasp_change_box.__init__(self, self.pot)
        return

    def cal_qe_dipo_screw(self,
                          in_tag="easy_easy",
                          input_s=0.0,
                          movex=0.0):

        e1 = 1. / 3. * np.array([1., 1., -2.])
        e2 = 1. / 2. * np.array([-1., 1., 0])
        e3 = np.array([0.5, 0.5, 0.5])

        n = 7
        m = 11
        atoms = self.set_bcc_convention([e1, e2, e3],
                                        (n, m, 1))  # z periodic 12
        atoms = self.cut_half_atoms(atoms)
        supercell = atoms.get_cell()
        strain = np.mat([[1.0, 0.0, 0.0],
                         [0.5, 1.0, 0.5],
                         [0.0, 0.0, 1.0]])

        supercell = strain * supercell
        atoms.set_cell(supercell)
        atoms.wrap(pbc=[1, 1, 1])
        atoms_perf = atoms.copy()
        atoms = self.intro_dipole_screw_atoms(atoms,
                                              self.pot['lattice'],
                                              movex,
                                              in_tag,
                                              input_s)
        self.write_poscar(atoms_perf)
        os.system("cp POSCAR new_perf.vasp")
        self.write_poscar(atoms)
        os.system("cp POSCAR POSCAR_init")

        self.gn_infile_dipole_screw_atoms(atoms)

        return (atoms_perf, atoms)


if __name__ == '__main__':
    usage = "usage:%prog [options] arg1 [options] arg2"
    parser = OptionParser(usage=usage)
    parser.add_option("-t", "--mtype",
                      action="store",
                      type="string",
                      dest="mtype", help="",
                      default="prp_r")
    (options, args) = parser.parse_args()

    drv = qe_dislocation()
    if options.mtype.lower() == 'dipole':
        drv.cal_qe_dipo_screw()
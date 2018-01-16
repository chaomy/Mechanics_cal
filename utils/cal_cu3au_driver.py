#!/usr/bin/env python
# encoding: utf-8
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2018-01-15 23:28:43
# @Last Modified by:   chaomy
# @Last Modified time: 2018-01-15 23:28:54


try:
    import os
    import glob
    import ase
    import ase.lattice.cubic as Cubic
    from optparse import OptionParser
    import shutil
    import Intro_vasp
    import gn_lmp_infile
    import gn_pbs
    import gn_config
    import numpy as np

except ImportError:
    print("error during import")


class cal_cu3au_driver(Intro_vasp.vasp_change_box,
                       gn_lmp_infile.gn_md_infile,
                       gn_pbs.gn_pbs,
                       gn_config.fcc):

    def __init__(self, structure='fcc'):
        self.lattice = 3.615
        gn_pbs.gn_pbs.__init__(self)
        gn_config.fcc.__init__(self)
        gn_lmp_infile.gn_md_infile.__init__(self)
        Intro_vasp.vasp_change_box.__init__(self,
                                            self.lattice)
        self.root = os.getcwd()
        self.insize = gn_lmp_infile.size
        self.insize.ydim = 40
        self.insize.zdim = 5
        return

    def prep_perf_bulk(self):
        insize = self.insize
        for i in range(1):
            insize.xdim = 80 + 10 * i
            dirname = "dir-x-%03d-z-%03d" % (insize.xdim, insize.zdim)
            self.mymkdir(dirname)
            self.write_cu3au_infile(self.lattice, insize)

            os.system("mv in.test  %s/in.cu3au" % (dirname))

            os.chdir(dirname)
            os.system("lmp_mpi -i in.cu3au")
            os.chdir(self.root)
        return

    def cal_cu3Au_dis(self):
        import atomman as am
        import atomman.lammps as lmp
        insize = self.insize

        delete_several_atoms = True
        intro_straight_dislocation = False
        intro_edge_nuclei = False
        delete_layers_atoms = False
        unitx = self.lattice * np.sqrt(2) / 2.
        unity = self.lattice * np.sqrt(2) / 2.
        unitz = self.lattice

        for i in range(1):
            xdim = 80 + 10 * i
            dirname = "dir-x-%03d-z-%03d" % (xdim, insize.zdim)

            #  atoms = ase.io.read("./%s/perf.dump" % (dirname),
            #  format='lammps-dump')
            infile = glob.glob("out/perf.dump*")[-1]
            atoms = ase.io.read(infile,
                                format='lammps-dump')

            ## introduce edge dislocation ##
            cell = atoms.get_cell()

            ## one pnt 39 to 40 ##
            xzone = [10.3 * unitx, 10.8 * unitx]
            yzone = [0.0,  4.0 * unity]
            zzone = [0.5 * cell[2, 2] - unitz, 0.5 * cell[2, 2] + 0.1]

            #  shiftAu = np.array([1.0 * unitx, 0.0, 0.0])
            ### delete some Cu atoms to create missing rows ###
            index_list = []

            for atom in atoms:
                pos = atom.position
                # three layer 7.5
                # two   layer
                if pos[1] < (3.5 * unity) and \
                    pos[0] < xzone[1] and \
                    pos[0] > xzone[0] and \
                    pos[2] < zzone[1] and \
                    pos[2] > zzone[0] and \
                        atom.symbol is 'H':
                        #  pos[0] < xzone[1] and \
                        #  pos[0] > xzone[0] and \
                        #  pos[2] < zzone[1] and \
                        #  pos[2] > zzone[0] and \
                    index_list.append(atom.index)

                elif delete_several_atoms:
                    if pos[1] < yzone[1] and \
                       pos[2] < zzone[1] and \
                       pos[2] > zzone[0] and \
                       pos[0] < xzone[1] and \
                       pos[0] > xzone[0]:
                        print atom.symbol
                        #  if atom.symbol == 'He':
                        #  index_list.append(atom.index)
                        #  atom.position += shiftAu

            print index_list
            del atoms[index_list]
            center = [0.5 * (xzone[0] + xzone[1]), yzone[1], 0.0]

            if intro_straight_dislocation:
                atoms = self.intro_single_edge_atoms(atoms,
                                                     center, 1,
                                                     [0, 1, 2])

            if intro_edge_nuclei:
                atoms = self.intro_edge_nuclei(atoms, center, -1,
                                               [0, 1, 2],
                                               zzone)

            #  cut extra half plane ##
            ratio = np.sqrt(2.) / 4.
            x_crit = ratio * self.lattice

            if delete_layers_atoms:
                for i in range(len(atoms)):
                    atom = atoms[i]
                if (atom.position[0] < (x_crit)):
                    index_list.append(atom.index)

                if index_list is not []:
                    print "delete %s atoms" % (len(index_list))
                    del atoms[index_list]

            ### output ###
            system, elements = am.convert.ase_Atoms.load(atoms)

            lmp.atom_data.dump(system, "init.txt")
            #  os.system("mv init.txt {}".format(dirname))
            #  os.system("cp in.minimize  {}".format(dirname))
        return

    def adjust(self):
        for i in range(16):
            xdim = 50 + 10 * i
            dirname = "dir-x-%03d" % (xdim)
            os.chdir(dirname)
            if not os.path.isdir("custom"):
                os.mkdir("custom")

            self.set_nnodes(1)
            self.set_ppn(12)
            self.set_wall_time(4)
            self.set_job_title("%s" % (dirname))
            self.set_main_job("mpirun lmp_linux -i in.minimize")
            self.write_pbs(od=True)

            os.chdir(self.root)
        return


if __name__ == '__main__':
    usage = "usage:%prog [options] arg1 [options] arg2"
    parser = OptionParser(usage=usage)
    parser.add_option("-t", "--mtype", action="store",
                      type="string", dest="mtype", help="",
                      default="prp_r")
    (options, args) = parser.parse_args()

    drv = cal_cu3au_driver()

    if options.mtype.lower() == 'make':
        drv.cal_cu3Au_dis()

    elif options.mtype.lower() == 'bulk':
        drv.prep_perf_bulk()

#!/usr/bin/env python
#- * - coding : utf - 8 - * -
#@Author : chaomy
#@Date : 2018 - 01 - 25 16 : 19 : 21
#@Last Modified by : chaomy
#@Last Modified time : 2018 - 01 - 25 18 : 40 : 08


from optparse import OptionParser
from gn_config import gnStructure
from copy import deepcopy
from plt_drv import plt_drv
from gsf import cal_va_gsf
import ase.lattice.cubic as cubic
import ase.lattice.orthorhombic as otho
import numpy as np
import gn_pbs
import get_data
import ase.io
import os
from numpy import sqrt
# The D03 structure is "based on FCC", but is really simple cubic
# with a basis.

latd03 = 7.46627803307887
area = 78.83577014166616


class D03Factory(cubic.SimpleCubicFactory):
    "A factory for creating Mg3Nd (D03) lattices."
    bravais_basis = [[0, 0, 0],
                     [0, 0.5, 0.5],
                     [0.5, 0, 0.5],
                     [0.5, 0.5, 0],
                     [0.25, 0.25, 0.75],
                     [0.25, 0.75, 0.25],
                     [0.75, 0.25, 0.25],
                     [0.25, 0.25, 0.25],
                     [0.75, 0.75, 0.25],
                     [0.75, 0.25, 0.75],
                     [0.25, 0.75, 0.75],
                     [0.75, 0.75, 0.75],
                     [0.5, 0.0, 0.0],
                     [0.0, 0.5, 0.0],
                     [0.0, 0.0, 0.5],
                     [0.5, 0.5, 0.5]]
    element_basis = (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1)

# x[-1, 1, 0] y[0, 0, 1] z[1, 1, 0]


class D03FactoryP110A(otho.SimpleOrthorhombicFactory):
    bravais_basis = [[0.0, 0.0, 0.0],
                     [0.5, 0.0, 0.0],
                     [0.25, 0.25, 0.0],
                     [0.75, 0.25, 0.0],
                     [0.25, 0.75, 0.0],
                     [0.75, 0.75, 0.0],
                     [0.25, 0.5, 0.5],
                     [0.75, 0.5, 0.5],
                     [0.0, 0.25, 0.5],
                     [0.0, 0.75, 0.5],
                     [0.5, 0.25, 0.5],
                     [0.5, 0.75, 0.5],
                     [0.0, 0.5, 0.0],
                     [0.5, 0.5, 0.0],
                     [0.25, 0.0, 0.5],
                     [0.75, 0.0, 0.5]]

    element_basis = (0, 0, 0, 0, 0, 0,
                     0, 0, 0, 0, 0, 0,
                     1, 1, 1, 1)


Mg3Nd = D03FactoryP110A()


class cal_d03(get_data.get_data,
              gn_pbs.gn_pbs,
              gnStructure,
              plt_drv,
              cal_va_gsf.cal_va_gsf):

    def __init__(self):
        get_data.get_data.__init__(self)
        gn_pbs.gn_pbs.__init__(self)
        gnStructure.__init__(self)
        plt_drv.__init__(self)
        cal_va_gsf.cal_va_gsf.__init__(self)
        self.disps = np.arange(0.0, 0.5, 0.05)

    def gn_displacement(self, atoms,
                        displacement_vector):
        positions = atoms.get_positions()
        atom_num = len(positions)
        displacement = deepcopy(positions)
        cut = 0.5 * np.max(positions[:, 2])
        for i in range(atom_num):
            if positions[i, 2] < cut:
                displacement[i] = [0, 0, 0]
            else:
                displacement[i] = displacement_vector
        return displacement

    def clc_lat(self):
        data = np.ndarray([20, 3])
        for i in range(20):
            mdir = "dir_{:03d}".format(i)
            os.chdir(mdir)
            (energy, vol, atoms) = self.vasp_energy_stress_vol_quick()
            data[i, 0] = atoms.get_cell()[0, 1]
            data[i, 1] = (energy)
            data[i, 2] = i
            os.chdir(os.pardir)
        np.savetxt('lat.dat', data)
        return

    def buildd03(self):
        # atoms = Mg3Nd(latticeconstant=latd03, size=(1, 1, 1),
        #               symbol=('Mg', 'Nd'))

        atoms = Mg3Nd(latticeconstant=(latd03 * sqrt(2), latd03,
                                       latd03 * sqrt(2) / 2.),
                      size=(1, 1, 18), symbol=('Mg', 'Nd'))

        # U = np.mat([[-1, 1, 0], [0, 0, 1], [0.5, 0.5, 0]])
        # Uinv = np.linalg.inv(U)
        # pos = atoms.get_scaled_positions()
        # print np.linalg.det(U)
        # ase.io.write("poscar.vasp", images=atoms, format="vasp")
        # self.write_lmp_config_data(atoms)
        return atoms

    def clc_data(self):
        disps = self.disps
        npts = len(disps)
        raw = np.ndarray([npts, 2])
        for i, disp in zip(range(npts), disps):
            mdir = 'dir_{}_{:4.3f}'.format(i, disp)
            os.chdir(mdir)
            raw[i, 0] = disp
            raw[i, 1] = np.loadtxt("out.txt")
            os.chdir(os.pardir)
        print raw
        np.savetxt("save.txt", raw)
        return

    def clc_plt(self):
        raw = np.loadtxt("save.txt")
        raw[:, 1] = (raw[:, 1] - raw[0, 1]) / area
        self.set_111plt()
        self.set_keys()
        self.ax.plot(raw[:, 0], raw[:, 1], **next(self.keysiter))
        self.fig.savefig("fig_gsf.png", **self.figsave)
        return

    def cal_stacking(self, opt="pre"):
        atoms = self.buildd03()
        perf_cells = deepcopy(atoms.get_cell())
        self.write_lmp_config_data(atoms)
        disps = self.disps
        npts = len(disps)
        for i, disp in zip(range(npts), disps):
            mdir = 'dir_{}_{:4.3f}'.format(i, disp)
            self.mymkdir(mdir)
            os.chdir(mdir)

            if opt in ["run"]:
                os.system("mpirun -n 4 lmp_mpi -i in.init")
            else:
                disp_vector = [disp, disp, 0.0]
                disp_matrix_direct = self.gn_displacement(
                    atoms.copy(), disp_vector)
                disp_matrix = deepcopy(disp_matrix_direct)

                disp_matrix[:, 0] = disp_matrix_direct[
                    :, 0] * perf_cells[0, 0]
                disp_matrix[:, 1] = disp_matrix_direct[
                    :, 1] * perf_cells[1, 1]

                local_atoms = atoms.copy()
                local_atoms.translate(disp_matrix)

                self.write_lmp_config_data(local_atoms)
                os.system("cp lmp_init.txt ../lmp.{:03d}".format(i))
                os.system("cp ../in.init .")

            os.chdir(os.pardir)
        print("area = ", latd03 * latd03 * sqrt(2))
        return


if __name__ == '__main__':
    usage = "usage:%prog [options] arg1 [options] arg2"
    parser = OptionParser(usage=usage)
    parser.add_option('-t', "--mtype", action="store",
                      type="string", dest="mtype")
    parser.add_option('-p', "--param", action="store",
                      type='string', dest="fargs")
    (options, args) = parser.parse_args()
    drv = cal_d03()
    dispatcher = {'build': drv.buildd03,
                  'gsf': drv.cal_stacking,
                  'clc': drv.clc_data,
                  'plt': drv.clc_plt}

    if options.fargs is not None:
        dispatcher[options.mtype.lower()](options.fargs)
    else:
        dispatcher[options.mtype.lower()]()

#!/usr/local/bin/python
# encoding: utf-8
import glob
import re
import numpy as np
import matplotlib.pylab as plt
from mpi4py import MPI
import gn_config
import sys


class md_crack_intro(object):

    def __init__(self):
        self.alat = 3.307898
        self.pi = 3.141592653589793
        self.burger = np.sqrt(3) / 2. * self.lattice_constant
        self.screw_coeff = self.burger / (2. * self.pi)
        self.Edge_coeff = self.burger / (2. * self.pi)
        self.P = 0.33
        return

    def Intro_Crack_cfg(self, atoms):
        alat = self.alat
        filename = sys.argv[1]
        a = alat / Ly0
        pi = self.pi
        xc, yc = 0.2, 0.5
        stress = 5.0
        k1 = stress * np.sqrt(self.pi * alat)
        Pratio = 0.333
        shearModule = 50.0
        coeff = k1 * np.sqrt(self.pi * 2) / (8 * shearModule * self.pi)
        print coeff
        pos = atoms.get_positions()
        natoms = len(pos)
        for i in range(len(atoms)):
            dx, dy = xnew[i] - xc, ynew[i] - yc
            if dx == 0:
                dx = 0.000000000001
            if dy == 0:
                dy = 0.000000000001
            r = math.sqrt(dx * dx + dy * dy)
            if ynew[i] > yc:
                if xnew[i] > xc:
                    theta = math.atan((dy / dx))
                if xnew[i] < xc:
                    theta = math.atan((-dx / dy)) + 0.5 * pi
                ux = coeff * math.sqrt(r) * ((5 - 8 * Pratio) *
                                             math.cos(0.5 * theta) - math.cos(1.5 * theta))
                uy = coeff * math.sqrt(r) * ((7 - 8 * Pratio) *
                                             math.sin(0.5 * theta) - math.sin(1.5 * theta))
            if ynew[i] < yc:
                if xnew[i] > xc:
                    theta = math.atan((-dy / dx))
                if xnew[i] < xc:
                    theta = math.atan((dx / dy)) + 0.5 * pi
                ux = coeff * math.sqrt(r) * ((5 - 8 * Pratio) *
                                             math.cos(0.5 * theta) - math.cos(1.5 * theta))
                uy = -1 * coeff * \
                    math.sqrt(r) * ((7 - 8 * Pratio) *
                                    math.sin(0.5 * theta) - math.sin(1.5 * theta))
            xnew[i] += ux
            ynew[i] += uy
        import matplotlib.pylab as plt
        from mpl_toolkits.mplot3d import Axes3D
        fig = plt.figure()
        ax = fig.gca(projection='3d')
    #    ax.set_zlim(-2.5,6.5)
        ax.scatter(xnew, ynew, znew)
        plt.show()
        AtomNumberOut = len(xnew)
        filename = "Crack.txt"
        with open(filename, mode="w") as fout:
            fout.write("#Edge_dislocation for bcc [111](110)\n")
            fout.write("\n")
            fout.write("%d atoms\n" % (AtomNumberOut))
            fout.write("1 atom types\n")
            fout.write("%f\t%f xlo xhi\n" % (0, Lx0))
            fout.write("%f\t%f ylo yhi\n" % (0, Ly0))
            fout.write("%f\t%f zlo zhi\n" % (0, Lz0))
            fout.write("xy xz yz = %8.5f %8.5f %8.5f" % (0, 0, 0))
            fout.write("\n")
            fout.write("Atoms\n")
            fout.write("\n")
            for i in range(AtomNumberOut):
                fout.write("%d  1  %12.7f %12.7f %12.7f\n" %
                           (i + 1, xnew[i], ynew[i], znew[i]))
        fout.close()
        return

    def Mpi_Intro_Crack_xyz(self,
                            k1,
                            p1, p2,
                            q1, q2,
                            u1, u2):
        from numpy.lib import scimath as SM
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        size = comm.Get_size()

        def updata_position(x, y, xc, yc):
            for i in range(len(x)):
                dx, dy = x[i] - xc, y[i] - yc
                r = np.sqrt(dx * dx + dy * dy)
                coeff = 100 * k1 * np.sqrt(2 * r / self.pi)
                if y[i] < yc:
                    theta = np.arctan2(-dy, dx)
                    ux = 1. / (u1 - u2) * \
                        (u1 * p2 * SM.sqrt(np.cos(theta) - u2 * np.sin(theta)) -
                         u2 * p1 * SM.sqrt(np.cos(theta) - u1 * np.sin(theta)))
                    uy = 1. / (u1 - u2) * \
                        (u1 * q2 * SM.sqrt(np.cos(theta) - u2 * np.sin(theta)) -
                         u2 * q1 * SM.sqrt(np.cos(theta) - u1 * np.sin(theta)))

                    ux = coeff * ux.real
                    uy = coeff * uy.real
                if y[i] >= yc:
                    theta = np.arctan2(dy, dx)
                    ux = 1. / (u1 - u2) * \
                        (u1 * p2 * SM.sqrt(np.cos(theta) - u2 * np.sin(theta)) -
                         u2 * p1 * SM.sqrt(np.cos(theta) - u1 * np.sin(theta)))

                    uy = 1. / (u1 - u2) * \
                        (u1 * q2 * SM.sqrt(np.cos(theta) - u2 * np.sin(theta)) -
                         u2 * q1 * SM.sqrt(np.cos(theta) - u1 * np.sin(theta)))
                    ux = coeff * ux.real
                    uy = -coeff * uy.real
                print ux, uy
                x[i] += ux
                y[i] += uy
            return (x, y)

        if rank == 0:
            filelist = glob.glob("./xyz/*")
            filename = filelist[-1]
            local_n = np.array([AtomNumber / size])
            print "in rank 0  x  and y size", len(x), len(y)
        else:
            local_n = np.array([0])
            x, y = None, None

        comm.Bcast(local_n, root=0)
        print "local_n of %d is " % (rank), local_n
        if rank == 0:
            if (len(x) % size == 0):
                check = True
            else:
                check = False
                body = (len(x) - np.mod(len(x), size))
                xrem, yrem = x[body:], y[body:]
                x, y = x[0:body], y[0:body]
                xrem, yrem = np.array(xrem), np.array(yrem)
            print "Now in rank 0  x  and y size", len(x), len(y)
        else:
            print "I am rank %d" % (rank)
        xc, yc = 0, 0
        x_local = np.zeros(local_n)
        y_local = np.zeros(local_n)

        if rank == 0:
            print rank,  len(x_local), len(y_local), len(x), len(y)
        if rank != 0:
            print rank,  len(x_local), len(y_local), x, y

        comm.Scatter(x, x_local, root=0)
        comm.Scatter(y, y_local, root=0)

        #print ("befor update")
        # for i in range(len(x_local)):
        #   print("[%d]  %f" %(rank, x[i]))
        (x_local, y_local) = updata_position(x_local, y_local, xc, yc)
        if rank == 0:
            if check == False:
                (xrem, yrem) = updata_position(xrem, yrem, xc, yc)
        comm.Gather(x_local, x, root=0)
        comm.Gather(y_local, y, root=0)
        if rank == 0:
            if check == False:
                x = np.concatenate([x, xrem])
                y = np.concatenate([y, yrem])
            filename = "Crack.txt"
            xlo, xhi, ylo, yhi = min(x) - 4, max(x) + 4, min(y) - 4, max(y) + 4

            with open(filename, mode="w") as fout:
                fout.write("#Edge_dislocation for bcc [111](110)\n")
                fout.write("\n")
                fout.write("%d atoms\n" % (AtomNumber))
                fout.write("1 atom types\n")
                fout.write("%f\t%f xlo xhi\n" % (xlo, xhi))
                fout.write("%f\t%f ylo yhi\n" % (ylo, yhi))
                fout.write("%f\t%f zlo zhi\n" % (zlo, zhi))
                fout.write("xy xz yz = %8.5f %8.5f %8.5f" % (0, 0, 0))
                fout.write("\n")
                fout.write("Atoms\n")
                fout.write("\n")
                for i in range(AtomNumber):
                    fout.write("%d  1  %12.7f %12.7f %12.7f\n"
                               % (Atomid[i], x[i], y[i], z[i]))
                fout.close()
        return

    def Intro_Crack_xyz(self,
                        k1,
                        p1, p2,
                        q1, q2,
                        u1, u2):
        from numpy.lib import scimath as SM
        filelist = glob.glob("./xyz/*")
        filename = filelist[-1]
        x, y, z, AtomNumber, xlo, xhi, ylo, yhi, zlo, zhi, Atomid = self.read_xyz(
            filename)
        xc, yc = 0, 0
        AtomNumber = len(x)
        for i in range(len(x)):
            dx, dy = x[i] - xc, y[i] - yc
            r = np.sqrt(dx * dx + dy * dy)
            coeff = 100 * k1 * np.sqrt(2 * r / self.pi)

            if y[i] < yc:
                theta = np.arctan2(-dy, dx)
                ux  = 1. / (u1 - u2) * \
                    (u1 * p2 * SM.sqrt(np.cos(theta) - u2 * np.sin(theta)) -
                     u2 * p1 * SM.sqrt(np.cos(theta) - u1 * np.sin(theta)))

                uy  = 1. / (u1 - u2) * \
                    (u1 * q2 * SM.sqrt(np.cos(theta) - u2 * np.sin(theta)) -
                     u2 * q1 * SM.sqrt(np.cos(theta) - u1 * np.sin(theta)))
                ux = coeff * ux.real
                uy = coeff * uy.real

            if y[i] >= yc:
                theta = np.arctan2(dy, dx)
                ux  = 1. / (u1 - u2) * \
                    (u1 * p2 * SM.sqrt(np.cos(theta) - u2 * np.sin(theta)) -
                     u2 * p1 * SM.sqrt(np.cos(theta) - u1 * np.sin(theta)))

                uy  = 1. / (u1 - u2) * \
                    (u1 * q2 * SM.sqrt(np.cos(theta) - u2 * np.sin(theta)) -
                     u2 * q1 * SM.sqrt(np.cos(theta) - u1 * np.sin(theta)))

                ux = coeff * ux.real
                uy = -coeff * uy.real
            print ux
            print uy
            x[i] += ux
            y[i] += uy
        xlo, xhi, ylo, yhi = min(x) - 4, max(x) + 4, min(y) - 4, max(y) + 4
   #     import matplotlib.pylab as plt
   #     from mpl_toolkits.mplot3d import Axes3D
   #     fig = plt.figure()
   #     ax = fig.gca(projection='3d')
   #     ax.scatter(xnew,ynew,znew)
   #     plt.show()
        filename = "Crack.txt"
        with open(filename, mode="w") as fout:
            fout.write("#Edge_dislocation for bcc [111](110)\n")
            fout.write("\n")
            fout.write("%d atoms\n" % (AtomNumber))
            fout.write("1 atom types\n")
            fout.write("%f\t%f xlo xhi\n" % (xlo, xhi))
            fout.write("%f\t%f ylo yhi\n" % (ylo, yhi))
            fout.write("%f\t%f zlo zhi\n" % (zlo, zhi))
            fout.write("xy xz yz = %8.5f %8.5f %8.5f" % (0, 0, 0))
            fout.write("\n")
            fout.write("Atoms\n")
            fout.write("\n")
            for i in range(AtomNumber):
                fout.write("%d  1  %12.7f %12.7f %12.7f\n"
                           % (Atomid[i], xnew[i], ynew[i], znew[i]))
        fout.close()
        return


class VA_ChangeBox(object):

    def __init__(self):
        self.lattice_constant = 3.3078
        self.pi = 3.141592653589793
        self.burger = np.sqrt(3) / 2. * self.lattice_constant
        self.screw_coeff = self.burger / (2. * self.pi)
        self.Edge_coeff = self.burger / (2. * self.pi)
        return

    def read_pos(self):
        with open("POSCAR", 'r') as fid:
            Raw = fid.readlines()
            fid.close()
        self.lattice_constant = float(Raw[1])
        print self.lattice_constant
        H = np.zeros([3, 3], "float")
        for i in range(3):
            for j in range(3):
                H[i, j] = float(Raw[2 + i].split()[j])
        AtomNumber = int(Raw[5])
        print AtomNumber
        Comm = str(Raw[6])
        print Comm
        AtomPosition = np.zeros([3, AtomNumber], "float")
        for j in range(AtomNumber):
            for i in range(3):
                AtomPosition[i, j] = float(Raw[7 + j].split()[i])
        print H
        print AtomPosition
        return AtomNumber, np.mat(H), Comm, AtomPosition

if __name__ == "__main__":
    M = MD_ChangeBox()
    # M.Intro_Crack_xyz(1.4828,(-0.00809419040147+9.50392012258e-19j),(-0.00142961912234+3.90049109137e-19j),(2.00783401372e-20-0.00249498017278j),(9.49992985589e-19-0.00463795644713j),(1.11022302463e-16+1.74520621177j),(1.38777878078e-16+0.572998189701j))
    # M.Intro_Screw_dipole()
    # N = VA_ChangeBox()
    # N.volume_conserving_Ortho_strain(0.5)
    # N.volume_conserving_Mono_strain(0.5)
    # M.Intro_Screw_dipole()
    M.Mpi_Intro_Crack_xyz(1.4828, (-0.00809419040147 + 9.50392012258e-19j),
                          (-0.00142961912234 + 3.90049109137e-19j),
                          (2.00783401372e-20 - 0.00249498017278j),
                          (9.49992985589e-19 - 0.00463795644713j),
                          (1.11022302463e-16 + 1.74520621177j),
                          (1.38777878078e-16 + 0.572998189701j))

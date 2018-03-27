#!/usr/bin/env python
# encoding: utf-8
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-07-05 01:09:32
# @Last Modified by:   chaomy
# @Last Modified time: 2018-03-27 15:58:47


import numpy as np
import ase.lattice
from mpi4py import MPI
from numpy.lib import scimath as SM
from numpy import arctan2
from math import pi
from math import cos
from math import sin
from math import sqrt

smsqrt = SM.sqrt


class md_crack_pre(object):

    def gn_perf_plate(self):
        atoms = ase.lattice.cubic.BodyCenteredCubic(directions=[[1, 0, 0],
                                                                [0, 1, 0],
                                                                [0, 0, 1]],
                                                    size=(120, 120, 5),
                                                    latticeconstant=self.pot[
                                                        'lattice'],
                                                    symbol=self.element,
                                                    pbc=(1, 1, 1))
        self.write_lmp_config_data(atoms)
        return atoms

    # called by mpi
    def updata_position(self, x, y, xc, yc):
        ckcoeff = self.ckcoeff
        u1, u2, p1, p2, q1, q2 = self.get_crack_coeffs()
        k1 = ckcoeff.K1
        for i in range(len(x)):
            dx, dy = x[i] - xc, y[i] - yc
            r = sqrt(dx * dx + dy * dy)
            coeff = 100 * k1 * sqrt(2 * r / pi)
            if y[i] < yc:
                theta = arctan2(-dy, dx)
                ux = 1. / (u1 - u2) * \
                    (u1 * p2 * smsqrt(cos(theta) - u2 * sin(theta)) -
                     u2 * p1 * smsqrt(cos(theta) - u1 * sin(theta)))

                uy = 1. / (u1 - u2) * \
                    (u1 * q2 * smsqrt(cos(theta) - u2 * sin(theta)) -
                     u2 * q1 * smsqrt(cos(theta) - u1 * sin(theta)))

                ux = coeff * ux.real
                uy = coeff * uy.real
            if y[i] >= yc:
                theta = arctan2(dy, dx)
                ux = 1. / (u1 - u2) * \
                    (u1 * p2 * smsqrt(cos(theta) - u2 * sin(theta)) -
                     u2 * p1 * smsqrt(cos(theta) - u1 * sin(theta)))

                uy = 1. / (u1 - u2) * \
                    (u1 * q2 * smsqrt(cos(theta) - u2 * sin(theta)) -
                     u2 * q1 * smsqrt(cos(theta) - u1 * sin(theta)))
                ux = coeff * ux.real
                uy = -coeff * uy.real
            print(ux, uy)
            x[i] += ux
            y[i] += uy
        return (x, y)

    def intro_crack_k1(self,
                       center=None,
                       atoms=None):
        u1, u2, p1, p2, q1, q2 = self.get_crack_coeffs()
        k1 = self.ckcoeff.K1
        natoms = len(atoms)
        cell = atoms.get_cell()
        pos = atoms.get_positions()

        if center is None:
            xc = 0.5 * cell[0, 0]
            yc = 0.5 * cell[1, 1]
        else:
            (xc, yc) = center

        for i in range(natoms):
            xs, ys = pos[i, 0], pos[i, 1]
            dx, dy = xs - xc, ys - yc
            r = sqrt(dx * dx + dy * dy)
            coeff = 100 * k1 * sqrt(2 * r / pi)
            if r < 140:
                if ys < yc:
                    theta = arctan2(-dy, dx)
                    ux = 1. / (u1 - u2) * \
                        (u1 * p2 * smsqrt(cos(theta) - u2 * sin(theta)) -
                         u2 * p1 * smsqrt(cos(theta) - u1 * sin(theta)))
                    uy = 1. / (u1 - u2) * \
                        (u1 * q2 * smsqrt(cos(theta) - u2 * sin(theta)) -
                         u2 * q1 * smsqrt(cos(theta) - u1 * sin(theta)))
                    ux = coeff * ux.real
                    uy = coeff * uy.real

                if ys >= yc:
                    theta = arctan2(dy, dx)
                    ux = 1. / (u1 - u2) * \
                        (u1 * p2 * smsqrt(cos(theta) - u2 * sin(theta)) -
                         u2 * p1 * smsqrt(cos(theta) - u1 * sin(theta)))

                    uy = 1. / (u1 - u2) * \
                        (u1 * q2 * smsqrt(cos(theta) - u2 * sin(theta)) -
                         u2 * q1 * smsqrt(cos(theta) - u1 * sin(theta)))
                    ux = coeff * ux.real
                    uy = -coeff * uy.real
                pos[i, 0], pos[i, 1] = xs + ux, ys + uy
        atoms.set_positions(pos)
        return atoms

    def Mpi_Intro_Crack_xyz(self, atoms):
        natoms = len(atoms)
        cell = atoms.get_cell()
        pos = atoms.get_positions()
        u1, u2, p1, p2, q1, q2 = self.get_crack_coeffs()
        k1 = self.ckcoeff.K1

        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        size = comm.Get_size()
        if rank == 0:
            local_n = np.array([natoms / size])
            print("in rank 0  x and y size", len(x), len(y))
        else:
            local_n = np.array([0])
            x, y = None, None

        comm.Bcast(local_n, root=0)
        print("local_n of %d is " % (rank), local_n)
        if rank == 0:
            if (len(x) % size == 0):
                check = True
            else:
                check = False
                body = (len(x) - np.mod(len(x), size))
                xrem, yrem = x[body:], y[body:]
                x, y = x[0:body], y[0:body]
                xrem, yrem = np.array(xrem), np.array(yrem)
            print("Now in rank 0  x  and y size", len(x), len(y))
        else:
            print("I am rank %d" % (rank))
        xc, yc = 0, 0
        x_local = np.zeros(local_n)
        y_local = np.zeros(local_n)

        if rank == 0:
            print(rank, len(x_local), len(y_local), len(x), len(y))
        if rank != 0:
            print(rank, len(x_local), len(y_local), x, y)

        comm.Scatter(x, x_local, root=0)
        comm.Scatter(y, y_local, root=0)

        # for i in range(len(x_local)):
        #   print("[%d]  %f" %(rank, x[i]))

        (x_local, y_local) = self.updata_position(x_local, y_local, xc, yc)
        if rank == 0:
            if check is False:
                (xrem, yrem) = self.updata_position(xrem, yrem, xc, yc)
        comm.Gather(x_local, x, root=0)
        comm.Gather(y_local, y, root=0)

        if rank == 0:
            if check is False:
                x = np.concatenate([x, xrem])
                y = np.concatenate([y, yrem])
            filename = "Crack.txt"
            xlo, xhi, ylo, yhi = min(x) - 4, max(x) + 4, min(y) - 4, max(y) + 4


if __name__ == "__main__":
    M.Mpi_Intro_Crack_xyz(1.4828,
                          (-0.00809419040147 + 9.50392012258e-19j),
                          (-0.00142961912234 + 3.90049109137e-19j),
                          (2.00783401372e-20 - 0.00249498017278j),
                          (9.49992985589e-19 - 0.00463795644713j),
                          (1.11022302463e-16 + 1.74520621177j),
                          (1.38777878078e-16 + 0.572998189701j))

#!/usr/bin/env python
# encoding: utf-8
#
import numpy as np
from numpy.lib import scimath as SM


class lmp_change_configs(object):

    def __init__(self):
        return

    def Intro_Crack(self, crackprams,
                    center=None,
                    atoms=None,
                    configname=None):
        k1, p1, p2, q1, q2, u1, u2 = crackprams.k1, \
            crackprams.p1, crackprams.p2, crackprams.q1, \
            crackprams.q2, crackprams.u1, crackprams.u2
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
            dx, dy = xs[i] - xc, ys[i] - yc
            r = np.sqrt(dx * dx + dy * dy)
            coeff = 100 * k1 * np.sqrt(2 * r / self.pi)
            if r < 140:
                if ys[i] < yc:
                    theta = np.arctan2(-dy, dx)
                    ux = 1. / (u1 - u2) * \
                        (u1 * p2 * SM.sqrt(np.cos(theta) - u2 * np.sin(theta)) -
                         u2 * p1 * SM.sqrt(np.cos(theta) - u1 * np.sin(theta)))
                    uy = 1. / (u1 - u2) * \
                        (u1 * q2 * SM.sqrt(np.cos(theta) - u2 * np.sin(theta)) -
                         u2 * q1 * SM.sqrt(np.cos(theta) - u1 * np.sin(theta)))
                    ux = coeff * ux.real
                    uy = coeff * uy.real

                if ys[i] >= yc:
                    theta = np.arctan2(dy, dx)
                    ux = 1. / (u1 - u2) * \
                        (u1 * p2 * SM.sqrt(np.cos(theta) - u2 * np.sin(theta)) -
                         u2 * p1 * SM.sqrt(np.cos(theta) - u1 * np.sin(theta)))

                    uy = 1. / (u1 - u2) * \
                        (u1 * q2 * SM.sqrt(np.cos(theta) - u2 * np.sin(theta)) -
                         u2 * q1 * SM.sqrt(np.cos(theta) - u1 * np.sin(theta)))
                    ux = coeff * ux.real
                    uy = -coeff * uy.real
            pos[i, 0], pos[i, 1] = xs[i] + ux, ys[i] + uy
            atoms.set_positions(pos)
        return

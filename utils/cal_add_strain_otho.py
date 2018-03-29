# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2018-02-23 15:24:11
# @Last Modified by:   chaomy
# @Last Modified time: 2018-03-04 20:00:25

# Ravindran, P., et al. "Density functional theory for calculation of
# elastic properties of orthorhombic crystals: Application to TiSi 2."
# Journal of Applied Physics 84.9 (1998): 4891-4904.


from numpy import mat, power
import copy

# def volume_conserving_ortho_strain_atoms(self, delta, atoms):
#     strain = np.mat([[1 + delta, 0, 0],
#                      [0, 1 - delta, 0],
#                      [0, 0, 1 / (1 - delta**2)]])

#     cell = atoms.get_cell()
#     pos = np.mat(atoms.get_positions())

#     cell = strain * cell
#     pos = pos * strain
#     atoms.set_positions(pos)
#     atoms.set_cell(cell)
#     return atoms


# class strain_simple(object, d, atoms, kk):
#     strainKeys = {}


class strain_otho(object):

    def strain_for_otho(self, d, atoms, kk):

        pd = 1.0 + d
        iv = 1. / power(1 - d * d, 1. / 3.)
        vd = iv * d
        pl = iv + vd
        mn = iv - vd

        strainKeys = {
            'D1': mat([[pd, 0, 0],
                       [0, 1, 0],
                       [0, 0, 1]]),

            'D2': mat([[1, 0, 0],
                       [0, pd, 0],
                       [0, 0, 1]]),

            'D3': mat([[1, 0, 0],
                       [0, 1, 0],
                       [0, 0, pd]]),

            'D4': mat([[iv, 0, 0],
                       [0, iv, vd],
                       [0, vd, iv]]),

            'D5': mat([[iv, 0, vd],
                       [0, iv, 0],
                       [vd, 0, iv]]),

            'D6': mat([[iv, vd, 0],
                       [vd, iv, 0],
                       [0, 0, iv]]),

            'D7': mat([[pl, 0, 0],
                       [0, mn, 0],
                       [0, 0, iv]]),

            'D8': mat([[pl, 0, 0],
                       [0, iv, 0],
                       [0, 0, mn]]),

            'D9': mat([[iv, 0, 0],
                       [0, pl, 0],
                       [0, 0, mn]])
        }

        strain = strainKeys[kk]

        org_cell = atoms.get_cell()
        new_cell = copy.deepcopy(org_cell)
        pos = mat(atoms.get_positions())
        pos = pos * strain

        new_cell = strain * new_cell
        atoms.set_cell(new_cell)
        atoms.set_positions(pos)
        return atoms

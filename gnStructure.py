#!/usr/bin/env python
# encoding: utf-8
#

import ase.lattice.compounds as Compounds
import ase.lattice.cubic as Cubic
import copy
import shutil
import re
import ase.io
import numpy as np


class AlloyStructure():
    def __init__(self,
                 directions,
                 size,
                 latticeconstant,
                 element,
                 structure):

        self._element = element
        self._directions = directions
        self._lattice_constant = latticeconstant
        self._structure = structure
        self._size = size
        return

    def gnL12(self):
        #### Pure Cu ###
        la_Cu = 3.49212
        la_Cu3Ag = 3.6193
        hz = 6 
        atomsCu = Cubic.FaceCenteredCubic(directions=self._directions,
                                          size=(200, 150, hz),
                                          latticeconstant=la_Cu3Ag,
                                          symbol='Cu',
                                          pbc=(1, 1, 1))

        print "Cu is", atomsCu.get_number_of_atoms()

        ase.io.write(filename="Cu.cfg",
                     images=atomsCu,
                     format='cfg')

        atoms = Compounds.AuCu3(directions=self._directions,
                                size=(100, 3, hz),
                                latticeconstant=la_Cu3Ag,
                                symbol=self._element,
                                pbc=(1, 1, 1))

        Cell = atomsCu.get_cell()
        NewCell = Cell
        atoms.set_cell(NewCell, scale_atoms=False)
        atoms.extend(atomsCu)

        ase.io.write(filename="Cu3Au.cfg",
                     images=atoms,
                     format='cfg')

        ase.io.write(filename="Cu3Au.xyz",
                     images=atoms,
                     format='extxyz')
        return

if __name__ == '__main__':
    Job = AlloyStructure(directions=[[-1, 1, 0],
                                [1, 1, 0],
                                [0, 0, -1]],
                    size=(2, 2, 2),
                    latticeconstant=3.690900,
                    element=('Au', 'Cu'),
                    structure='L12')
    Job.gnL12()

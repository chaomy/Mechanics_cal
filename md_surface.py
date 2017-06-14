#!/usr/local/bin/python
# encoding: utf-8

import os, shutil, re

class md_surface(object):
    def __init__(self,
                 in_lattice = None,
                 in_md_surface_type = None):
        self.ev_to_j = 1.60218e-19
        if in_lattice is not None:
            self._lattice =  in_lattice
        else:
            self._lattice = 3.3

        if in_md_surface_type is not None:
            self._md_surface_type = in_md_surface_type
        else:
            self._md_surface_type = '100'  # default
        self.file_name_a = 'in.surface%sA'%(self._md_surface_type)
        self.file_name_b = 'in.surface%sB'%(self._md_surface_type)
        return

    def set_md_surface_lattice(self, lattice_constant):
        self._lattice = lattice_constant
        return

    def set_md_surface_type(self, in_type):
        self._md_surface_type = in_type
        return

    def change_lattice(self, a_or_b):
        if a_or_b == 'A':
            fileName = self.file_name_a
        elif a_or_b == 'B':
            fileName = self.file_name_b

        with open(fileName, 'r') as fid:
            raw = fid.readlines()
            fid.close()
        change_lin_num = 6
        raw[change_lin_num] = \
            "variable    latparam1   equal    %f\n"\
            %(self._lattice)

        shutil.move(fileName,
                    '%s.copy'%(fileName))
        with open(fileName, 'w') as fid:
            fid.writelines(raw)
            fid.close()
        return

    def get_md_surface_energy(self, infile):
        with open(infile , 'r') as fid:
            for line in fid:
                matchA = \
                    re.search(r"Surface Area = (\d*\.\d*)", line)
                if matchA:
                    print line.split()
                    area = float(line.split()[-2])

                matchE = \
                    re.search(r"Final energy of atoms =   (-?\d*\.\d*) eV", line)
                if matchE:
                    print line.split()
                    energy = float(line.split()[-2])

                matchN = \
                    re.search("new total = \d*",line)
                if matchN :
                    print line.split()

        return (area, energy)

    def cal_surface_energy(self):
        in_file_a = self.file_name_a
        in_file_b = self.file_name_b
        if os.path.isfile(in_file_a):
            print "change Lattice A"
            self.change_lattice('A')

        if os.path.isfile(in_file_b):
            print "change Lattice B"
            self.change_lattice('B')

        os.system("lmp_linux < %s > A.log"%(in_file_a))
        os.system("lmp_linux < %s > B.log"%(in_file_b))

        area, energy_a = self.get_md_surface_energy('A.log')
        area, energy_b = self.get_md_surface_energy('B.log')

        surface_energy = 0.5 * (energy_a - energy_b) * self.ev_to_j / area

        print "surface energy is ", surface_energy
        return surface_energy

if __name__ == '__main__':
    A = md_surface('111')
    A.cal_surface_energy()

#!/usr/local/bin/python
# encoding: utf-8
#
import glob, sys, os
import numpy as np
import matplotlib.pylab as plt
import ase.io
from numpy.lib import scimath as SM
import gn_config

class lmp_change_box(object):
    def __init__(self,
                 lattice_constant = None):
        if lattice_constant is not None:
            self.lattice_constant = lattice_constant
        else:
            self.lattice_constant = 3.17
        return

    def set_intro_coeff(self):
        self.pi = 3.141592653589793
        self.burger = np.sqrt(3.) / 2. * self.lattice_constant
        self.screw_coeff = self.burger / (2. * self.pi)
        self.Edge_coeff  = self.burger / (2. * self.pi)
        self.possion = 0.33
        return

    def IntroPerf(self):
        filename = './Cu3Au.cfg'
        (elementdata, xyzdata, Cell) = self.read_cfg(filename)
        AtomNumber = len(elementdata)
        x, y, z, AtomType, AtomType2 = [], [], [], [], []
        xnew, ynew, znew = [], [], []
        count_Cu = 0
        for i in range(AtomNumber):
            x.append(float(xyzdata[i][0])*Cell[0]);
            y.append(float(xyzdata[i][1])*Cell[1]);
            z.append(float(xyzdata[i][2])*Cell[2])
            if  elementdata[i] == 'Cu':
                AtomType.append(1)
                count_Cu += 1
            elif elementdata[i] == 'Au':
                AtomType.append(2)
        xlo = 0.0; xhi = Cell[0]
        ylo = 0.0; yhi = Cell[1]
        zlo = 0.0; zhi = Cell[2]

        min_y = np.min(y)
        print "The min y new is ", np.min(y)
        unit = 3.639029845822
        Length = 2.5 * unit * np.sqrt(2)/2 - 3.2 + 0.1
        for i in range(AtomNumber):
            if y[i] - min_y <= Length and AtomType[i] == 1:
                    print "Delete", y[i], AtomType[i]
            else:
                xnew.append(x[i])
                ynew.append(y[i])
                znew.append(z[i])
                AtomType2.append(AtomType[i])

        AtomType = AtomType2
        AtomNumberOut = len(xnew)
        print "Delete ", AtomNumber - AtomNumberOut

        xlo, xhi = xlo - 5, xhi + 5
        ylo, yhi = ylo - 5, yhi + 5
        filename = 'Perf.txt'
        with open(filename,mode="w") as fout:
            fout.write("#Edge_dislocation for bcc [111](110)\n")
            fout.write("\n")
            fout.write("%d atoms\n"%(AtomNumberOut))
            fout.write("2 atom types\n")
            fout.write("%f\t%f xlo xhi\n"%(xlo,xhi))
            fout.write("%f\t%f ylo yhi\n"%(ylo,yhi))
            fout.write("%f\t%f zlo zhi\n"%(zlo,zhi))
            fout.write("xy xz yz = %8.5f %8.5f %8.5f"%(0,0,0))
            fout.write("\n")
            fout.write("Atoms\n")
            fout.write("\n")
            count = 0
            for i in range(AtomNumberOut):
                if AtomType[i] == 1:
                    count += 1
                    fout.write("%d  1  %12.7f %12.7f %12.7f\n"
                    %(count,  xnew[i] ,ynew[i] ,znew[i]))
            for i in range(AtomNumberOut):
                if AtomType[i] == 2:
                    count += 1
                    fout.write("%d  2  %12.7f %12.7f %12.7f\n"
                    %(count,  xnew[i] ,ynew[i] ,znew[i]))
        fout.close()
        return

    def Intro_Screw_cfg(self):
        filename = sys.argv[1]
        lattice_constant = 3.30789893315
        x,y,z,ccsym,energy,Sxx,Syy,Szz,Sxy,Sxz,Syz,coord1,atomid,Lx0,Ly0,Lz0 = self.read_cfg(filename)
        xnew,ynew,znew,coordN,energyN,atomidN = [],[],[],[],[],[]
        a =   float(Ly0)/ float(lattice_constant)
        dd = 0.4 /a
        xc ,yc = 0.5,0.5
        AtomNumber = len(x)
        for i in range(AtomNumber):
            if x[i] > xc :
                if  np.abs(y[i] - yc) > dd :
                    xnew.append(float(x[i])) ; ynew.append(float(y[i])); znew.append(float(z[i]))
                    coordN.append(coord1[i]) ; atomidN.append(atomid[i]); energyN.append(energy[i])
            else :
                xnew.append(float(x[i])) ; ynew.append(float(y[i])); znew.append(float(z[i]))
                coordN.append(coord1[i]) ; atomidN.append(atomid[i]); energyN.append(energy[i])
        print "atom total: " , AtomNumber
        print "atom deleted" , AtomNumber - len(xnew)
        pi = self.pi
        coeff = self.screw_coeff
        for i in range(len(xnew)):
            if ynew[i] > yc :
                if xnew[i] > xc :
                    theta = np.arctan((ynew[i] - yc)/(xnew[i] - xc))
                elif xnew[i] < xc :
                    theta = np.arctan((xc - xnew[i])/(ynew[i] - yc)) + 0.5 * pi
            elif ynew[i]< yc :
                if xnew[i] < xc :
                    theta = np.arctan((yc - ynew[i])/(xc - xnew[i])) + pi
                elif xnew[i] > xc :
                    theta = np.arctan((xnew[i] - xc)/(yc - ynew[i])) + 1.5 * pi
            dz = self.lattice_constant * coeff * theta
            print dz
            znew[i] += dz
        with open("screw.cfg",'w') as fid:
            fid.write("Number of particles = %d\n"%(len(xnew)))
            fid.write("A = 1 Angstrom (basic length-scale)\n")
            fid.write("""H0(1,1) = %f A
H0(1,2) = 0 A
H0(1,3) = 0 A
H0(2,1) = 0 A
H0(2,2) = %f A
H0(2,3) = 0 A
H0(3,1) = 0 A
H0(3,2) = 0 A
H0(3,3) = %f A
.NO_VELOCITY.
entry_count = 6
auxiliary[0] = c_coordn
auxiliary[1] = id
auxiliary[2] = energy
"""%(Lx0,Ly0,Lz0))
            for i in range(len(xnew)):
                fid.write("""92.906380
C
%6.4f %6.4f %6.4f %d %d %6.4f
"""%(xnew[i],ynew[i],znew[i],coordN[i],atomidN[i],energyN[i]))
        return

    def Intro_Screw_dipole(self):
        filelist = glob.glob('./xyz/*')
        filename = filelist[-1]
        xs , ys, zs , AtomNumber ,xlo, xhi, ylo, yhi, zlo, zhi = self.read_xyz(filename)
        xs = np.array(xs) ; ys = np.array(ys) ; zs = np.array(zs)
        xc1 , xc2 , yc = 0.25 * xhi + 1.25 * xlo , 0.75 * xhi - 0.25 * xlo , 0.5 * (yhi + ylo)
        #-------------------- cut  plane along x positive direction
        coeff  = self.screw_coeff
        xnew , ynew , znew = [] , [], []
        for i in range(AtomNumber):
            if np.abs(zs[i]) < 5.0 :
                xnew.append(float(xs[i])) ; ynew.append(float(ys[i])); znew.append(float(zs[i]))
        AtomNumberOut = len(xnew)
        print AtomNumber
        print "delete %d atoms"%(AtomNumber-AtomNumberOut)
        for i in range(AtomNumberOut):
            theta = np.arctan2((ynew[i] - yc),(xnew[i] - xc1))
            dz = coeff * theta
            znew[i] = dz + znew[i]
        for i in range(AtomNumberOut):
            theta = np.arctan2((ynew[i] - yc),(xnew[i] - xc2))
            dz = -coeff * theta
            znew[i] = dz + znew[i]
        import matplotlib.pylab as plt
        from mpl_toolkits.mplot3d import Axes3D
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        ax.set_zlim(-2.5,6.5)
        ax.scatter(xnew,ynew,znew)
        plt.show()

        filename = "screw.txt"
        with open(filename,mode="w") as fout:
            fout.write("#screw_dislocation for bcc [111](110)\n")
            fout.write("\n")
            fout.write("%d atoms\n"%(AtomNumberOut))
            fout.write("1 atom types\n")
            fout.write("%f\t%f xlo xhi\n"%(xlo,xhi))
            fout.write("%f\t%f ylo yhi\n"%(ylo,yhi))
            fout.write("%f\t%f zlo zhi\n"%(zlo,zhi))
            fout.write("xy xz yz = %8.5f %8.5f %8.5f"%(0,0,0))
            fout.write("\n")
            fout.write("Atoms\n")
            fout.write("\n")
            for i in range(AtomNumberOut):
                fout.write("%d  1  %12.7f %12.7f %12.7f\n"%(i+1, xnew[i] ,ynew[i] ,znew[i] ))
        fout.close()
        return

    def Intro_Screw(self):
        filename = './Cu3Au.cfg'
        (elementdata, xyzdata, Cell) = self.read_cfg(filename)
        AtomNumber = len(elementdata)

        x, y, z, AtomType = [], [], [], []
        count_Cu = 0
        for i in range(AtomNumber):
            x.append(float(xyzdata[i][0])*Cell[0]);
            y.append(float(xyzdata[i][1])*Cell[1]);
            z.append(float(xyzdata[i][2])*Cell[2])
            if  elementdata[i] == 'Al':
                AtomType.append(1)
                count_Cu += 1
            elif elementdata[i] == 'Ni':
                AtomType.append(2)
        xlo = 0.0; xhi = Cell[0]
        ylo = 0.0; yhi = Cell[1]
        zlo = 0.0; zhi = Cell[2]

        xs = np.array(x) ; ys = np.array(y) ; zs = np.array(z)
        print xs

        xc = 0.5 * (xlo + xhi)
        yc = 0.5 * (ylo + yhi) + 0.03

        Xc = 0.5 * (xlo + xhi)
        Yc = 0.5 * (ylo + yhi)
        #-------------------- cut  plane along x positive direction
        coeff  = self.screw_coeff

        xnew, ynew, znew = [], [], []
        for i in range(AtomNumber):
            dx, dy = xs[i] - xc, ys[i] - yc
            Dx, Dy = xs[i] - Xc, ys[i] - Yc
            r = Dx**2 + Dy**2
            if r < 10000:
                theta = np.arctan2(dy, dx)
                dz = coeff * theta
                znew.append(dz + zs[i])
                xnew.append(xs[i])
                ynew.append(ys[i])

        # ------------------- OutPut screw.txt --------------------
        xlo = xlo - 5; xhi = xhi + 5
        ylo = ylo - 5; yhi = yhi + 5
        print "Xc, Yc",  Xc, Yc

        AtomNumber = len(xnew)
        filename = "Init.txt"
        with open(filename,mode="w") as fout:
            fout.write("#screw_dislocation for bcc [111](110)\n")
            fout.write("\n")
            fout.write("%d atoms\n"%(AtomNumber))
            fout.write("2 atom types\n")
            fout.write("%f\t%f xlo xhi\n"%(xlo,xhi))
            fout.write("%f\t%f ylo yhi\n"%(ylo,yhi))
            fout.write("%f\t%f zlo zhi\n"%(zlo,zhi))
            fout.write("xy xz yz = %8.5f %8.5f %8.5f"%(0,0,0))
            fout.write("\n")
            fout.write("Atoms\n")
            fout.write("\n")
            count = 0
            for i in range(AtomNumber):
                if AtomType[i] == 1:
                    count += 1
                    fout.write("%d  1  %12.7f %12.7f %12.7f\n"
                    %(count,  xnew[i] ,ynew[i] ,znew[i]))
            for i in range(AtomNumber):
                if AtomType[i] == 2:
                    count += 1
                    fout.write("%d  2  %12.7f %12.7f %12.7f\n"
                    %(count,  xnew[i] ,ynew[i] ,znew[i]))
        fout.close()
        return

    def Intro_Edge(self):
        import math
        mtype = 'cfg'
        if mtype == 'Myrun':
            fileLists = glob.glob('xyz/perf*')
            filename = fileLists[0]
            xs, ys, zs, AtomNumber, xlo, xhi, ylo,\
                    yhi, zlo, zhi = self.read_xyz(filename)
            xlen = xhi - xlo
            ylen = yhi - ylo
            zlen = zhi - zlo
            xs = np.array(xs)
            ys = np.array(ys)
            zs = np.array(zs)

            #--------------- cut  plane along x positive direction --------
            xnew , ynew , znew, AtomType = [], [], [], []
            for i in range(AtomNumber):
                xnew.append(float(xs[i]));
                ynew.append(float(ys[i]));
                znew.append(float(zs[i]))
                AtomType.append(1)

            # ---------------- intro Edge dislocation -------------
            AtomNumberOut = len(xnew)
            print AtomNumber
            print "delete %d atoms"%(AtomNumber-AtomNumberOut)

        elif mtype == 'cfg':
            filename = './Cu3Au.cfg'
            (elementdata, xyzdata, Cell) = self.read_cfg(filename)
            AtomNumber = len(elementdata)
            xnew , ynew , znew, AtomType = [], [], [], []
            count_Cu = 0
            for i in range(AtomNumber):
                xnew.append(float(xyzdata[i][0])*Cell[0]);
                ynew.append(float(xyzdata[i][1])*Cell[1]);
                znew.append(float(xyzdata[i][2])*Cell[2])
                if  elementdata[i] == 'Cu':
                    AtomType.append(1)
                    count_Cu += 1
                elif elementdata[i] == 'Au':
                    AtomType.append(2)
            xlo = 0.0; xhi = Cell[0]
            ylo = 0.0; yhi = Cell[1]
            zlo = 0.0; zhi = Cell[2]
            AtomNumberOut = AtomNumber

        xc = 0.5 * (xlo + xhi)
        yc = 0.5 * (ylo + yhi)

        coeff  = self.Edge_coeff
        Dx, Dy, Dz = [], [], []
        for i in range(AtomNumberOut):
            dx, dy = xnew[i] - xc, ynew[i] - yc
            A = (1 - self.P) * (dx**2 + dy**2)
            if A != 0.0:
                ux = coeff * (np.arctan2(dy, dx) + (dx * dy)/(2.* A))
                uy = (1 - 2*self.P)/(4 * (1 - self.P)) * \
                        math.log(dx**2 + dy**2) + (dx**2 - dy**2)/ \
                        (4 * (1 - self.P)*(dx**2 + dy**2))
                uy *= -coeff
            elif A == 0.0:
                print "Find A == 0 "
                print xnew[i], ynew[i]
                ux = uy = 0.0

            xnew[i] += ux
            ynew[i] += uy

            Dx.append(ux)
            Dy.append(uy)
            Dz.append(0.0)

        import matplotlib.pylab as plt
        from mpl_toolkits.mplot3d import Axes3D
        fig = plt.figure()
        ax = fig.gca(projection='3d')

        #ax.set_zlim(-2.5,6.5)
        ax.quiver(xnew,ynew,znew,
                  Dx, Dy, Dz,
                  pivot = 'tail')

        ax.plot([xc,xc],
                [yc, yc],
                [zlo,zhi],
                color = 'g',
                linestyle = '--')

        ax.set_xlabel('X axis')
        ax.set_ylabel('Y axis')
        ax.set_zlabel('Z axis')

       # plt.show()
        # ------------------- OutPut Edge.txt --------------------
        xlo, xhi = xlo - 5, xhi + 5
        ylo, yhi = ylo - 5, yhi + 5

        filename = "Edge.txt"
        with open(filename,mode="w") as fout:
            fout.write("#Edge_dislocation for bcc [111](110)\n")
            fout.write("\n")
            fout.write("%d atoms\n"%(AtomNumberOut))
            fout.write("2 atom types\n")
            fout.write("%f\t%f xlo xhi\n"%(xlo,xhi))
            fout.write("%f\t%f ylo yhi\n"%(ylo,yhi))
            fout.write("%f\t%f zlo zhi\n"%(zlo,zhi))
            fout.write("xy xz yz = %8.5f %8.5f %8.5f"%(0,0,0))
            fout.write("\n")
            fout.write("Atoms\n")
            fout.write("\n")
            count = 0
            for i in range(AtomNumberOut):
                if AtomType[i] == 1:
                    count += 1
                    fout.write("%d  1  %12.7f %12.7f %12.7f\n"
                    %(count,  xnew[i] ,ynew[i] ,znew[i]))
            for i in range(AtomNumberOut):
                if AtomType[i] == 2:
                    count += 1
                    fout.write("%d  2  %12.7f %12.7f %12.7f\n"
                    %(count,  xnew[i] ,ynew[i] ,znew[i]))
        fout.close()
        return

    def Intro_Edge_cut(self):
        import math
        mtype = 'cfg'
        if mtype == 'Myrun':
            fileLists = glob.glob('xyz/perf*')
            filename = fileLists[0]
            xs, ys, zs, AtomNumber, xlo, xhi, ylo,\
                    yhi, zlo, zhi = self.read_xyz(filename)
            xlen = xhi - xlo
            ylen = yhi - ylo
            zlen = zhi - zlo
            xs = np.array(xs)
            ys = np.array(ys)
            zs = np.array(zs)
            #--------------- cut  plane along x positive direction --------
            xnew , ynew , znew, AtomType = [], [], [], []
            for i in range(AtomNumber):
                xnew.append(float(xs[i]));
                ynew.append(float(ys[i]));
                znew.append(float(zs[i]))
                AtomType.append(1)
            # ---------------- intro Edge dislocation -------------
            AtomNumberOut = len(xnew)
            print AtomNumber
            print "delete %d atoms"%(AtomNumber-AtomNumberOut)

        elif mtype == 'cfg':
            filename = './Cu3Au.cfg'
            (elementdata, xyzdata, Cell) = self.read_cfg(filename)
            AtomNumber = len(elementdata)

            x, y, z, AtomType, AtomType2 = [], [], [], [], []
            xnew, ynew, znew = [], [], []
            count_Cu = 0
            for i in range(AtomNumber):
                x.append(float(xyzdata[i][0])*Cell[0]);
                y.append(float(xyzdata[i][1])*Cell[1]);
                z.append(float(xyzdata[i][2])*Cell[2])
                if  elementdata[i] == 'Cu':
                    AtomType.append(1)
                    count_Cu += 1
                elif elementdata[i] == 'Au':
                    AtomType.append(2)
            xlo = 0.0; xhi = Cell[0]
            ylo = 0.0; yhi = Cell[1]
            zlo = 0.0; zhi = Cell[2]

            min_y = np.min(y)
            print "The min y new is ", np.min(y)
            unit = self.lattice_constant
            Length = 2.5 * unit * np.sqrt(2)/2 - 3.2 + 0.1
            for i in range(AtomNumber):
                if y[i] - min_y <= Length and AtomType[i] == 1:
                        print "Delete", y[i], AtomType[i]
                else:
                    xnew.append(x[i])
                    ynew.append(y[i])
                    znew.append(z[i])
                    AtomType2.append(AtomType[i])

            AtomType = AtomType2
            AtomNumberOut = len(xnew)
            print "Delete ", AtomNumber - AtomNumberOut

        xc = 0.5 *(xlo + xhi)
        yc = 4.5 * unit

        coeff  = self.Edge_coeff
        Dx, Dy, Dz = [], [], []
        for i in range(AtomNumberOut):
            dx, dy = xnew[i] - xc, ynew[i] - yc
            A = (1 - self.P) * (dx**2 + dy**2)
            if A != 0.0:
                ux = coeff * (np.arctan2(dy, dx) + (dx * dy)/(2.* A))
                uy = (1 - 2*self.P)/(4 * (1 - self.P)) * \
                        math.log(dx**2 + dy**2) + (dx**2 - dy**2)/ \
                        (4 * (1 - self.P)*(dx**2 + dy**2))
                uy *= -coeff
            elif A == 0.0:
                print "Find A == 0 "
                print xnew[i], ynew[i]
                ux = uy = 0.0

            xnew[i] += ux
            ynew[i] += uy

        Layer = 0.5 * np.sqrt(2) * unit
        xlo, xhi = xlo - 5  - 8 * Layer,  xhi + 5 + 8 * Layer
        ylo, yhi = ylo - 10 - 6 * Layer , yhi + 5 + 6 * Layer

        filename = "Edge.txt"
        with open(filename,mode="w") as fout:
            fout.write("#Edge_dislocation for bcc [111](110)\n")
            fout.write("\n")
            fout.write("%d atoms\n"%(AtomNumberOut))
            fout.write("2 atom types\n")
            fout.write("%f\t%f xlo xhi\n"%(xlo,xhi))
            fout.write("%f\t%f ylo yhi\n"%(ylo,yhi))
            fout.write("%f\t%f zlo zhi\n"%(zlo,zhi))
            fout.write("xy xz yz = %8.5f %8.5f %8.5f"%(0,0,0))
            fout.write("\n")
            fout.write("Atoms\n")
            fout.write("\n")
            count = 0
            for i in range(AtomNumberOut):
                if AtomType[i] == 1:
                    count += 1
                    fout.write("%d  1  %12.7f %12.7f %12.7f\n"
                    %(count,  xnew[i] ,ynew[i] ,znew[i]))
            for i in range(AtomNumberOut):
                if AtomType[i] == 2:
                    count += 1
                    fout.write("%d  2  %12.7f %12.7f %12.7f\n"
                    %(count,  xnew[i] ,ynew[i] ,znew[i]))
        fout.close()
        return

    def Intro_Edge_dipole(self):
        import math
        fileLists = glob.glob('xyz/perf*')
        filename = fileLists[0]
        xs, ys, zs, AtomNumber, xlo, xhi, ylo,\
                yhi, zlo, zhi = self.read_xyz(filename)
        xlen = xhi - xlo
        ylen = yhi - ylo
        zlen = zhi - zlo
        xs = np.array(xs)
        ys = np.array(ys)
        zs = np.array(zs)

        xc = 0.5 * (xhi + xlo);
        ycup = 0.45 * yhi + 0.25 * ylo
        ycdown = 0.25 * yhi + 0.45 * ylo

        #--------------- cut  plane along x positive direction --------
        coeff  = self.Edge_coeff
        xnew , ynew , znew = [], [], []
        for i in range(AtomNumber):
            if zs[i] >= 0:
                if xs[i] > 0.5 * xlen:
                    if  np.abs(ys[i] - yc) >= 0 :
                        xnew.append(float(xs[i]));
                        ynew.append(float(ys[i]));
                        znew.append(float(zs[i]))
                else :
                    xnew.append(float(xs[i]));
                    ynew.append(float(ys[i]));
                    znew.append(float(zs[i]))
        # ---------------- intro Edge dislocation -------------
        AtomNumberOut = len(xnew)
        print AtomNumber
        print "delete %d atoms"%(AtomNumber-AtomNumberOut)

        Dx, Dy, Dz = [], [], []
        D1x, D1y, D1z = [], [], []
        D2x, D2y, D2z = [], [], []

        for i in range(AtomNumberOut):
            dx, dy = xc- xnew[i], ycup - ynew[i]
            A = (1 - self.P) * (dx**2 + dy**2)
            if A != 0 :
                ux1 = coeff * (np.arctan2(dy, dx) + (dx * dy)/(2. * A))
                uy1 = (1 - 2*self.P)/(4 * (1 - self.P)) * \
                        math.log(dx**2 + dy**2) + (dx**2 - dy**2)/ \
                        (4 * (1 - self.P)*(dx**2 + dy**2))
                uy1 *= -coeff
            else :
                ux1 = uy1 = 0
            D1x.append(ux1)
            D1y.append(uy1)
            D1z.append(0.0)

            xnew[i] += ux1
            ynew[i] += uy1

        for i in range(AtomNumberOut):
            dx, dy = xnew[i] - xc, ynew[i] - ycdown
            A = (1 - self.P) * (dx**2 + dy**2)
            if A != 0 :
                ux2 = coeff * (np.arctan2(dy, dx) + (dx * dy)/(2. * A))
                uy2 = (1 - 2*self.P)/(4 * (1 - self.P)) * \
                        math.log(dx**2 + dy**2) + (dx**2 - dy**2)/ \
                        (4 * (1 - self.P)*(dx**2 + dy**2))
                uy2 *= -coeff
            else :
                ux2 = uy2 = 0
            xnew[i] += ux2
            ynew[i] += uy2

            D2x.append(ux2)
            D2y.append(uy2)
            D2z.append(0.0)

        import matplotlib.pylab as plt
        from mpl_toolkits.mplot3d import Axes3D
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        #ax.set_zlim(-2.5,6.5)
        #------- plot center of dislocation --------- #

        ax.plot([xc,xc],
                [ycup, ycup],
                [zlo,zhi],
                color = 'g',
                linestyle = '--')

        ax.plot([xc,xc],
                [ycdown, ycdown],
                [zlo,zhi],
                color = 'y',
                linestyle = '--')

        ax.quiver(xnew,ynew,znew,
                  D1x, D1y, D1z,
                  pivot = 'tail',
                  color = 'r')

        ax.quiver(xnew,ynew,znew,
                  D2x, D2y, D2z,
                  pivot = 'tail',
                  color = 'b')

        ax.set_xlabel('X axis')
        ax.set_ylabel('Y axis')
        ax.set_zlabel('Z axis')

        plt.show()
        xlo, xhi = xlo * 1.1 , xhi * 1.1
        ylo, yhi = ylo * 1.1 , yhi * 1.1

        # ------------------- OutPut Edge.txt --------------------
        return

    def intro_dipole_screw(self):
        atoms = ase.io.read("./lmp_init.cfg", format='cfg')
        print atoms.get_cell()
        print atoms.get_positions()
        ase.io.write("./lmp_init.xyz",
                      images=atoms,
                      format='xyz')
        M = gn_config.gnStructure()
        M.write_lmp_config_data(atoms)
        return

    def Intro_Crack(self,mtype,k1,p1,p2,q1,q2,u1,u2):
        if mtype == 'cfg':
            filename = './Init_config.cfg'
            (elementdata, xyzdata, Cell) = self.read_cfg(filename)
            AtomNumber = len(elementdata)
            xs, ys, zs = [], [], []
            xnew , ynew , znew, AtomType = [], [], [], []
            count_Cu = 0

            for i in range(AtomNumber):
                xs.append(float(xyzdata[i][0])*Cell[0]);
                ys.append(float(xyzdata[i][1])*Cell[1]);
                zs.append(float(xyzdata[i][2])*Cell[2])
                if  elementdata[i] == self.element[0]:
                    AtomType.append(1)
                    count_Cu += 1
                elif elementdata[i] == self.element[1]:
                    AtomType.append(2)

            xlo = 0.0; xhi = Cell[0]
            ylo = 0.0; yhi = Cell[1]
            zlo = 0.0; zhi = Cell[2]

            xc = 0.5 * (xlo + xhi)
            yc = 0.5 * (ylo + yhi)

        elif mtype == 'xyz':
            filelist = glob.glob("./xyz/*")
            filename = filelist[-1]
            x, y, z , AtomNumber, xlo, xhi, ylo, yhi, zlo, zhi = self.read_xyz(filename)

        for i in range(len(xs)):
            dx , dy = xs[i] - xc , ys[i] - yc
            r = np.sqrt(dx * dx + dy * dy)
            coeff = 100 * k1 * np.sqrt(2 * r / self.pi)
            if r < 140:
                if ys[i] < yc:
                    theta = np.arctan2(-dy, dx)
                    ux  = 1./(u1-u2) * (u1 * p2 * SM.sqrt(np.cos(theta) - u2 * np.sin(theta)) - u2 * p1 * SM.sqrt(np.cos(theta) - u1 * np.sin(theta)))
                    uy  = 1./(u1-u2) * (u1 * q2 * SM.sqrt(np.cos(theta) - u2 * np.sin(theta)) - u2 * q1 * SM.sqrt(np.cos(theta) - u1 * np.sin(theta)))
                    ux  = coeff * ux.real
                    uy  = coeff * uy.real
                if ys[i] >= yc:
                    theta = np.arctan2(dy, dx)
                    ux  = 1./(u1-u2) * (u1 * p2 * SM.sqrt(np.cos(theta) - u2 * np.sin(theta)) - u2 * p1 * SM.sqrt(np.cos(theta) - u1 * np.sin(theta)))
                    uy  = 1./(u1-u2) * (u1 * q2 * SM.sqrt(np.cos(theta) - u2 * np.sin(theta)) - u2 * q1 * SM.sqrt(np.cos(theta) - u1 * np.sin(theta)))
                    ux  = coeff * ux.real
                    uy  = -coeff * uy.real

                xnew.append(xs[i] + ux)
                ynew.append(ys[i] + uy)
                znew.append(zs[i])

        xlo , xhi, ylo, yhi = min(xnew) - 4 , max(xnew) + 4 , min(ynew) - 4 , max(ynew) + 4

        AtomNumberOut = len(xnew)
        filename = "Crack.txt"
        with open(filename,mode="w") as fout:
            fout.write("#Edge_dislocation for bcc [111](110)\n")
            fout.write("\n")
            fout.write("%d atoms\n"%(AtomNumberOut))
            fout.write("1 atom types\n")
            fout.write("%f\t%f xlo xhi\n"%(xlo,xhi))
            fout.write("%f\t%f ylo yhi\n"%(ylo,yhi))
            fout.write("%f\t%f zlo zhi\n"%(zlo,zhi))
            fout.write("xy xz yz = %8.5f %8.5f %8.5f"%(0,0,0))
            fout.write("\n")
            fout.write("Atoms\n")
            fout.write("\n")
            for i in range(AtomNumberOut):
                fout.write("%d  1  %12.7f %12.7f %12.7f\n"
                        %(i+1, xnew[i], ynew[i], znew[i]))
        fout.close()
        return

if __name__=="__main__":
    M = lmp_change_box()
    # M.Intro_Crack_xyz(1.4828,(-0.00809419040147+9.50392012258e-19j),(-0.00142961912234+3.90049109137e-19j),(2.00783401372e-20-0.00249498017278j),(9.49992985589e-19-0.00463795644713j),(1.11022302463e-16+1.74520621177j),(1.38777878078e-16+0.572998189701j))
    #N.Intro_Screw_dipole()
    #M.Intro_Screw()
    #M.IntroPerf()
    #M.Intro_Crack('cfg',1.4828,(-0.00809419040147+9.50392012258e-19j),(-0.00142961912234+3.90049109137e-19j),(2.00783401372e-20-0.00249498017278j),(9.49992985589e-19-0.00463795644713j),(1.11022302463e-16+1.74520621177j),(1.38777878078e-16+0.572998189701j))
    M.intro_dipole_screw()

from pylab import *
import os
import numpy as np
import pandas as pd
from fractions import Fraction
from decimal import Decimal
from itertools import cycle
import plt_drv


def fitness(v):
    return math.sqrt(v[1] * v[1] + v[2] * v[2] + v[0] * v[0])


base_path = os.getcwd()
os.chdir(base_path)

dirtemplate = 'gb-'
subdirtemplate = 'mesh-'
outtempl = 'Egb ='

delta_angle = 0.1
ang_range = 90.0

delta_incl = 0.1
incl_range = 45
# sampling the sngle space
delta_sector = 1
Nx = int(ang_range / delta_sector)


# change dirlist to erange *

dirlist = [x for x in os.listdir('.') if x.startswith(dirtemplate)]
subdirlist = []
for dir in dirlist:
    xxx = [dir + '/' +
           x for x in os.listdir(dir) if x.startswith(subdirtemplate)]
    subdirlist = subdirlist + xxx
# print dirlist
lbox_list = []
answer = []
elist = []
anglelist = []
anglelistphi = []
egblist = []
xshlist = []
zshlist = []
orientuplist = []
orientlowlist = []
gbdirlist = []


for dir in subdirlist:
    lbox_list = []
    filename = dir + '/log.lammps'
    try:
        print(dir)
        o = open(filename, 'r').readlines()
    except IOError:
        print("No file \n")
    else:
        try:
            starts = [n for n, l in enumerate(o) if l.startswith(outtempl)]
            if len(starts) > 0:
                aindex = starts[0]  # +1
        except IOError:
            print("Error")
        else:
            # relaxed energy
            if len(starts) > 0:
                lammps_filename = dir + '/gb4.in_final'
                ol = open(lammps_filename, 'r').readlines()

                startsup = [n for n, l in enumerate(
                    ol) if l.startswith('create_atoms 1 region upper')]
                aindexup = startsup[0] - 2
                startslow = [n for n, l in enumerate(
                    ol) if l.startswith('create_atoms 2 region lower')]
                aindexlow = startslow[0] - 2

                orientup = ol[aindexup]
                orientlow = ol[aindexlow]

                orientuplist.append(orientup)
                orientlowlist.append(orientlow)

                Egb = o[aindex].split()[2]
                gbdir = dir.split('/')[0]
                gbmesh = dir.split('/')[1]

                gbangle = gbdir.split('-')[1]
                # gbangle_phi=gbdir.split('-')[2]
                gbangle_phi = gbdir[
                    len(dirtemplate) + len(gbangle) + len('-'):]

                xsh = gbmesh[len(subdirtemplate):].split('-')[0]
                zsh = gbmesh[len(subdirtemplate):].split('-')[1]

#               answer.append(gbangle+' '+gbangle_phi+' '+Egb+' '+xsh+' '+zsh+'\n')
                answer.append(gbangle + ' ' + gbangle_phi + ' ' +
                              Egb + ' ' + xsh + ' ' + zsh + ' ' + orientup + '\n')

                elist.append((float(Egb), float(xsh), float(zsh)))
                anglelist.append(float(gbangle))
                anglelistphi.append(float(gbangle_phi))
                egblist.append(float(Egb))
                xshlist.append(float(xsh))
                zshlist.append(float(zsh))
#                gbdirlist.append(float(gbangle)+float(gbangle_phi))
                gbnum = int(float(gbangle) / delta_sector) + \
                    int((float(gbangle_phi) + incl_range) / delta_sector) * Nx
                gbdirlist.append(gbnum)


d = {'angle': anglelist,
     'anglephi': anglelistphi,
     'Egb': egblist,
     'xsh': xshlist,
     'zsh': zshlist,
     'orientup': orientuplist,
     'orientlow': orientlowlist,
     'gbdir': gbdirlist
     }

print(d['gbdir'])

alldataframe = pd.DataFrame(d)

# get rows indexes: group by angle minimize by fitness
idx = alldataframe.groupby(['gbdir'])['Egb'].transform(
    min) == alldataframe['Egb']

# select rows by index
gbenergy_best = alldataframe[idx].reset_index(drop=True)
gbenergy_best = gbenergy_best.sort('angle')


# plot
mplt = plt_drv.plt_drv()
mplt.set_111plt()
mplt.set_keys()
mplt.ax.plot(gbenergy_best['angle'],
             gbenergy_best['Egb'], label='Egb', **next(mplt.keysiter))
mplt.add_legends(mplt.ax)
mplt.add_x_labels(cycle(['angle']), mplt.ax)
mplt.add_y_labels(cycle(['E [eV/m^2]']), mplt.ax)
mplt.set_tick_size(mplt.ax)
mplt.fig.savefig('fig_egb.png', **mplt.figsave)


# get best energies

denergy = {'angle': gbenergy_best['angle'].tolist(),
           'anglephi': gbenergy_best['anglephi'].tolist(),
           'Egb': gbenergy_best['Egb'].tolist(),
           'xsh': gbenergy_best['xsh'].tolist(),
           'zsh': gbenergy_best['zsh'].tolist()
           }

denergydf = pd.DataFrame(denergy)

#	for k in range(aindex+10, len(o)-5):
#	    xx=o[k].split()[1]
#	    lbox_list.append(float(xx))
# aver=sum(lbox_list)/float(len(lbox_list))

# string=dir[len(dirtemplate):]+' '+str(xx)+'\n'
# answer.append(string)
# print elist
# print answer with 5 lowest energy boundaries
# for i in range(5):
#    print sorted(elist)[i]
# sorted(elist)[0] gives you the lowest energy boundary and x and z vectors

w = open('Egb.dat', 'w')
w.writelines(answer)
w.close()

np.savetxt(r'Egbbest.txt', denergydf, fmt='%f')

answer_for_buld = []
for i in range(len(gbenergy_best)):
    ang = str(gbenergy_best['angle'].tolist()[i])
    ang_phi = str(gbenergy_best['anglephi'].tolist()[i])

    xsh = gbenergy_best['xsh'].tolist()[i]
    zsh = gbenergy_best['zsh'].tolist()[i]

    u_x = gbenergy_best['orientup'].tolist()[i].split()[5] + ' ' \
        + gbenergy_best['orientup'].tolist()[i].split()[6] + ' ' \
        + gbenergy_best['orientup'].tolist()[i].split()[7]

    u_y = gbenergy_best['orientup'].tolist()[i].split()[10] + ' ' \
        + gbenergy_best['orientup'].tolist()[i].split()[11] + ' ' \
        + gbenergy_best['orientup'].tolist()[i].split()[12]

    u_z = gbenergy_best['orientup'].tolist()[i].split()[15] + ' ' \
        + gbenergy_best['orientup'].tolist()[i].split()[16] + ' ' \
        + gbenergy_best['orientup'].tolist()[i].split()[17]

    l_x = gbenergy_best['orientlow'].tolist()[i].split()[5] + ' ' \
        + gbenergy_best['orientlow'].tolist()[i].split()[6] + ' ' \
        + gbenergy_best['orientlow'].tolist()[i].split()[7]

    l_y = gbenergy_best['orientlow'].tolist()[i].split()[10] + ' ' \
        + gbenergy_best['orientlow'].tolist()[i].split()[11] + ' ' \
        + gbenergy_best['orientlow'].tolist()[i].split()[12]

    l_z = gbenergy_best['orientlow'].tolist()[i].split()[15] + ' ' \
        + gbenergy_best['orientlow'].tolist()[i].split()[16] + ' ' \
        + gbenergy_best['orientlow'].tolist()[i].split()[17]

    xuplen = fitness(list(map(int, u_x.split())))
    xlowlen = fitness(list(map(int, l_x.split())))
    x = Fraction(Decimal(xuplen / xlowlen)
                 ).limit_denominator(max_denominator=40)
    l_mult = x.numerator
    u_mult = x.denominator

    mystring = ang + ' ' + ang_phi + ' upper ' + u_x + ' ' + u_y + ' ' + u_z + ' lower ' + l_x + ' ' + l_y \
        + ' ' + l_z + ' ' + str(u_mult) + ' ' + str(l_mult) + \
        ' ' + str(xsh) + ' ' + str(zsh) + '\n'
    answer_for_buld.append(mystring)


w = open('BuildGBlist.lst_temp', 'w')
w.writelines(answer_for_buld)
w.close()

# x=Fraction(Decimal(5.0/6.0)).limit_denominator(max_denominator=40)
# print int(x.numerator)
# print int(x.denominator)

# nlist=[]
# for x in answer:
#    nlist.append(float(x.split()[0]))
#    nlist.append(float(x.split()[1]))

# narray=np.array(nlist).reshape(len(nlist)/2,2)
# print narray


# def fitfunc(x,a,b):
#    return a*x+b

#    popt, pcov =curve_fit(fitfunc,a[:,0],a[:,1])

#    print "Directory: %r " % dir
#    print aver
#    print len(lbox_list)

# a=np.array(my_list).reshape(int((maxeng-mineng)/stepeng),engaver-1,3)


#spav= np.average(a,axis=1)
#spstd= np.std(a,axis=1)


# spav[:,2]=1620-spav[:,2]

# engcol=np.array(spav[:,0]).reshape(len(spav),1)
# sputcol=np.array(spav[:,2]).reshape(len(spav),1)
# sputdevcol=np.array(spstd[:,2]).reshape(len(spav),1)

# np.savetxt('myres.dat',np.concatenate((engcol,sputcol,sputdevcol),1))

#plt.errorbar(spav[:,0],spav[:,2],yerr=spstd[:,2], fmt='-o')


# plt.xlim((mineng*0.9,maxeng*1.1))
# plt.show()


# print b
# print "STD=%r" %r np.average(a[2,:,2])
# print a.ndim


# plt.figure(1)
# plt.subplot(211)
#plt.plot(t1, f(t1), 'bo', t2, f(t2), 'k')

# plt.subplot(212)
#plt.plot(t2, np.cos(2*np.pi*t2), 'r--')
# plt.show()

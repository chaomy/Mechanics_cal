#!/usr/bin/env python 

import os
import re
import shutil
import random
import time
import sys
import subprocess as sp
#import queue_num_checker as qc


base_path = os.getcwd()
os.chdir(base_path)


argv = sys.argv
maxsubjobs = 150
checktime = 1
dirtemplate = 'gb-'  # sys.argv[1]	#'fz-'
# dirlist=[]
msublist = []


def get_subdir(dirs):
    subdirlist = []
    for dir in dirs:
        a = [dir + '/' + x for x in os.listdir(dir)]
        subdirlist.extend(a)
    return subdirlist


dirlist = [x for x in os.listdir('.') if x.startswith(dirtemplate)]
print "the list \n"
print dirlist
dirlist = get_subdir(dirlist)


for dir in dirlist:

    lammpsfile = dir + '/log.lammps'
    try:
        f = open(lammpsfile, 'r')
    except IOError:
        msublist.append(dir)
    else:
        print "There \n"
        print lammpsfile
        if os.stat(lammpsfile).st_size == 0:
            msublist.append(dir)


j = 0
# njobsrunning=qc.checkNQueue('frolov2').read().number
print msublist

# np=float(os.environ['SLURM_NPROCS'])
# pid=float(os.environ['SLURM_PROCID'])

# first=int(pid/np*len(msublist))
# last=int((pid+1)/np*len(msublist))

# j=first
while j < len(msublist):  # j <  last:
    # len(msublist):
    os.chdir(msublist[j])
#    o=open('lammps.out','w')
    p1 = sp.Popen(['lmp_mpi', '-i', 'gb4.in_final'])
    p1.wait()

    os.chdir(base_path)
    print "New j we are looping"
    print j
    print msublist[j]
    j += 1


os.chdir(base_path)

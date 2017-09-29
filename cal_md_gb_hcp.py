#!/usr/bin/env python
# encoding: utf-8
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-07-05 08:12:30
# @Last Modified by:   chaomy
# @Last Modified time: 2017-09-29 16:02:52

inhcp = """
# ---------- Initialize Simulation --------------------- 
clear 
units metal 
dimension 3 
boundary  p s p 
atom_style atomic 

# ---------- Create Atomistic Structure ---------------------
variable  alat  equal    %f  
variable  clat  equal    %f 
variable  xsize equal    10
variable  ysize equal    10
variable  zsize equal    2 
variable  xdim  equal    ${alat}*${xsize}
variable  ydim  equal    sqrt(3)*${alat}*${ysize}
variable  zdim  equal    ${clat}*${zsize}
variable  sqa   equal    sqrt(3)
variable  b1    equal    1./3.
variable  b2    equal    5./6.
variable  c     equal    ${clat}/${alat}

variable ydimmin equal   -30
variable ydimmax equal    30
variable yfreedimmin equal   -20
variable yfreedimmax equal    20

region 		whole 	block  0.0 ${xdim}  ${ydimmin} ${ydimmax} 0.000000 ${zdim} units box
create_box  4  whole 

# ---------- Create Atomistic Structure --------------------- 
#variable xdim equal xsize      #${alat}*13.92838828

region upper block INF INF -0.05 ${ydimmax} INF INF units box 
lattice   custom   ${alat} &
a1  %f  	%f 		%f     &
a2  %f      %f      %f     &
a3  %f      %f      %f     &
basis   0.0   0.0   0.0    &
basis   0.5   0.5   0.0    &
basis   0.0   ${b1} 0.5    &
basis   0.5   ${b2} 0.5   
create_atoms 1 region upper 

region lower block INF INF ${ydimmin} 0.05 INF INF units box 
lattice   custom   ${alat} &
a1  %f  	%f 		%f     &
a2  %f      %f      %f     &
a3  %f      %f      %f     &
basis   0.0   0.0   0.0    &
basis   0.5   0.5   0.0    &
basis   0.0   ${b1} 0.5    &
basis   0.5   ${b2} 0.5
create_atoms 2 region lower 

region uplayer  block INF INF ${yfreedimmax} ${ydimmax} INF INF units box 
region lowlayer block INF INF ${ydimmin} ${yfreedimmin} INF INF units box 

group upper region upper
group lower region lower

group uplayer  region uplayer
group lowlayer region lowlayer

group rigbody union uplayer lowlayer
group gbregion subtract all rigbody

set group lowlayer type 3 
set group uplayer  type 4 

#set group gbregion type 2
# ---------- Define Interatomic Potential --------------------- 
pair_style   %s 

# EAM/FS  Sun 
# pair_coeff   * *  %s   Mg  Mg  Mg  Mg 

# MEAM Kim
pair_coeff     * *  ../lib_MgNdPb.meam  Mg  ../MgNd_para.meam  Mg  Mg  Mg  Mg 

neighbor 	 2.0  bin 
neigh_modify delay 5 check yes  
# ---------- Displace atoms and delete overlapping atoms --------------------- 
displace_atoms upper move  %f  0 zmove units box

delete_atoms overlap 0.1 lower upper   # default
#delete_atoms overlap 0.1 lower upper 

# ---------- Define Settings --------------------- 
compute eng all pe/atom 
compute eatoms gbregion reduce sum c_eng 

# ---------- Run Minimization --------------------- 
fix f1 lowlayer setforce 0.0 0.0 0.0
fix f2 uplayer  aveforce 0.0 0.0 0.0

thermo  100 
thermo_style custom step pe lx ly lz press pxx pyy pzz c_eatoms 
dump	 1 all custom 100 out/rel.chkpt.*  id type mass x y z c_eng
#restart 1 out/rel.restart

min_style cg 
minimize 1e-15 1e-15 10000 10000 

variable Ecoh equal %f 
variable Ntot equal "count(gbregion)"
variable GBarea equal "lx*lz"
variable Ebulk equal "v_Ecoh*count(gbregion)"
variable Etot equal  "c_eatoms"
variable Egb equal   "(c_eatoms-v_Ecoh*count(gbregion))/lx/lz*16.021"

print "Energy per atom =  ${Ecoh}"
print "Bulk energy ${Ebulk}"
print "Etot ${Etot}"
print "AreaXY = ${GBarea}"
print "N total atoms = ${Ntot}"
print "Egb = ${Egb}" #append ${f} screen yes
print "Egb = ${Egb}" append sepfile.dat screen yes
"""
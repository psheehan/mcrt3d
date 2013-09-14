#!/usr/bin/env python

from numpy import *
import matplotlib.pyplot as plt

au = 1.496e13

f = open("grid.dat","r")
line = f.readline()
nr = int(line.split(" ")[0])
line = f.readline()
np = int(line.split(" ")[0])
line = f.readline()
nz = int(line.split(" ")[0])

rw = zeros(nr)
pw = zeros(np)
zw = zeros(nz)

for i in range(nr):
    rw[i] = float(f.readline())/au
for i in range(np):
    pw[i] = f.readline()
for i in range(nz):
    zw[i] = float(f.readline())/au

f.close()

r = zeros((nr-1,np-1,nz-1))
p = zeros((nr-1,np-1,nz-1))
z = zeros((nr-1,np-1,nz-1))
dens = zeros((nr-1,np-1,nz-1))
temp = zeros((nr-1,np-1,nz-1))

f1 = open("dens.dat","r")
f2 = open("temp.dat","r")

for i in range(nr-1):
    for j in range(np-1):
        for k in range(nz-1):
            r[i,j,k] = 0.5*(rw[i+1]+rw[i])
            p[i,j,k] = 0.5*(pw[j+1]+pw[j])
            z[i,j,k] = 0.5*(zw[k+1]+zw[k])
            dens[i,j,k] = f1.readline()
            temp[i,j,k] = f2.readline()
            
f1.close()
f2.close()

plt.imshow(temp[:,:,0],origin="lower",interpolation="nearest",vmin=temp.min(), \
        vmax=temp.max())
plt.colorbar()
plt.show()

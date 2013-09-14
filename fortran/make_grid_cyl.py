#!/usr/bin/env python

from numpy import arange, zeros, pi

au = 1.496e13

nr = 10
np = 10
nz = 10

r = arange(nr)*au/2
p = arange(np)/(np-1.)*2*pi
z = (arange(nz)-(nz-1)/2.)*au/1

dens = zeros((nr-1,np-1,nz-1)) + 1.0e-17

f = open("grid.dat","w")
f.write("{0:d}   {1:d}\n".format(nr,1))
f.write("{0:d}   {1:d}\n".format(np,1))
f.write("{0:d}   {1:d}\n".format(nz,1))
for i in range(r.size):
    f.write("{0:e}\n".format(r[i]))
for i in range(p.size):
    f.write("{0:e}\n".format(p[i]))
for i in range(z.size):
    f.write("{0:e}\n".format(z[i]))
f.close()

f = open("dens.dat","w")
for i in range(nr-1):
    for j in range(np-1):
        for k in range(nz-1):
            f.write("{0:e}\n".format(dens[i,j,k]))
f.close()

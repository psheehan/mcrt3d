#!/usr/bin/env python

from numpy import arange, zeros

au = 1.496e13

nx = 10
ny = 10
nz = 10

x = (arange(nx)-(float(nx)-1)/2)*au/1
y = (arange(ny)-(float(ny)-1)/2)*au/1
z = (arange(nz)-(float(nz)-1)/2)*au/1

dens = zeros((ny-1,nx-1,nz-1)) + 1.0e-17

f = open("grid.dat","w")
f.write("{0:d}   {1:d}\n".format(nx,1))
f.write("{0:d}   {1:d}\n".format(ny,1))
f.write("{0:d}   {1:d}\n".format(nz,1))
for i in range(x.size):
    f.write("{0:e}\n".format(x[i]))
for i in range(y.size):
    f.write("{0:e}\n".format(y[i]))
for i in range(z.size):
    f.write("{0:e}\n".format(z[i]))
f.close()

f = open("dens.dat","w")
for i in range(nx-1):
    for j in range(ny-1):
        for k in range(nz-1):
            f.write("{0:e}\n".format(dens[j,i,k]))
f.close()

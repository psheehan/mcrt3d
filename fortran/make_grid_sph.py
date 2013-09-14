#!/usr/bin/env python

from numpy import arange, zeros, pi

au = 1.496e13

nr = 10
nt = 10
np = 2

r = arange(nr)*au/2
t = arange(nt)/(nt-1.)*pi
p = arange(np)/(np-1.)*2*pi

dens = zeros((nr-1,nt-1,np-1)) + 1.0e-17

f = open("grid.dat","w")
f.write("{0:d}   {1:d}\n".format(nr,1))
f.write("{0:d}   {1:d}\n".format(nt,1))
f.write("{0:d}   {1:d}\n".format(np,1))
for i in range(r.size):
    f.write("{0:e}\n".format(r[i]))
for i in range(t.size):
    f.write("{0:e}\n".format(t[i]))
for i in range(p.size):
    f.write("{0:e}\n".format(p[i]))
f.close()

f = open("dens.dat","w")
for i in range(nr-1):
    for j in range(nt-1):
        for k in range(np-1):
            f.write("{0:e}\n".format(dens[i,j,k]))
f.close()

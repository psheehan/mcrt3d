#!/usr/bin/env python

from numpy import loadtxt, zeros, ones, arange
from hyperion.model import Model
from hyperion.model import ModelOutput
from hyperion.dust import IsotropicDust
from hyperion.util.constants import *
import matplotlib.pyplot as plt

m = Model()

data = loadtxt('../dustkappa_yso.inp', skiprows=2)
nu = c / (data[:,0].copy()[::-1] * 1.0e-4)
kabs = data[:,1].copy()[::-1]
ksca = data[:,2].copy()[::-1]
chi = kabs + ksca
albedo = ksca / chi

d = IsotropicDust(nu, albedo, chi)

nr = 10
np = 10
nz = 10

r = arange(nr)*au/2
p = arange(np)/(np-1.)*2*pi
z = (arange(nz)-(float(nz)-1)/2)*au/1

m.set_cylindrical_polar_grid(r, z, p)

dens = zeros((nr-1,np-1,nz-1)) + 1.0e-17

m.add_density_grid(dens, d)

source = m.add_spherical_source()
source.luminosity = lsun
source.radius = rsun
source.temperature = 4000.

m.set_n_photons(initial=1000000, imaging=0)
m.set_convergence(True, percentile=99., absolute=2., relative=1.02)

m.write("test_cylindrical.rtin")

#m.run("test_cylindrical.rtout", mpi=False)

n = ModelOutput('test_cylindrical.rtout')

grid = n.get_quantities()

temp = grid.quantities['temperature'][0]

for i in range(9):
    plt.imshow(temp[i,:,:],origin="lower",interpolation="nearest", \
            vmin=temp.min(),vmax=temp.max())
    plt.colorbar()
    plt.show()

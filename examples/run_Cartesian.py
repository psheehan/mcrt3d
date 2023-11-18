#!/usr/bin/env python3


from scipy.constants import au as AU
from astropy.constants import M_sun, R_sun
import matplotlib.pyplot as plt
from mcrt3d import MCRT
from mcrt3d import IsotropicDust
from time import time
import numpy
import sys

# Create a model class.

model = MCRT()

# Set up the dust.

data = numpy.loadtxt('dustkappa_yso.inp', skiprows=2)

lam = data[::-1,0].copy() * 1.0e-4
kabs = data[::-1,1].copy()
ksca = data[::-1,2].copy()

d = IsotropicDust(lam, kabs, ksca)

# Set up the grid.

nx = 10
ny = 10
nz = 10

x = (numpy.arange(nx)-(float(nx)-1)/2)*AU/1
y = (numpy.arange(ny)-(float(ny)-1)/2)*AU/1
z = (numpy.arange(nz)-(float(nz)-1)/2)*AU/1

model.set_cartesian_grid(x,y,z)

#sys.exit(0)

# Set up the density.

density = numpy.zeros((nx-1,ny-1,nz-1)) + 1.0e-16

model.grid.add_density(density, d)

# Set up the star.

model.grid.add_star(mass=M_sun.value, radius=R_sun.value, temperature=4000.)
model.grid.sources[-1].set_blackbody_spectrum(lam)

# Run the thermal simulation.

t1 = time()
model.thermal_mc(nphot=1000000, bw=True, use_mrw=False, mrw_gamma=2, \
        verbose=False)
t2 = time()
print(t2-t1)

# Run the images.

model.run_image("image", numpy.array([1.]), 256, 256, 0.1, 100000, incl=0., \
        pa=0, dpc=1.)

model.run_unstructured_image("uimage", numpy.array([1300.]), 10, 10, 2.5, \
        100000, incl=0., pa=0., dpc=1.)

# Run the spectra.

model.run_spectrum("SED", numpy.logspace(-1,4,200), 10000, incl=0, pa=0, dpc=1.)

# Plot the temperature structure.

for i in range(9):
    plt.imshow(model.grid.temperature[0][:,:,i], origin="lower",\
            interpolation="nearest", vmin=model.grid.temperature[0].min(),\
            vmax=model.grid.temperature[0].max())
    plt.colorbar()
    plt.show()

# Plot the images.

fig, ax = plt.subplots(nrows=1, ncols=1)

ax.imshow(model.images["image"].intensity[:,:,0], origin="lower", \
        interpolation="none")

plt.show()

# Plot the spectra.

fig, ax = plt.subplots(nrows=1, ncols=1)

ax.loglog(model.spectra["SED"].lam, model.spectra["SED"].intensity)

ax.set_ylim([1.0e-23,1.0e7])

plt.show()

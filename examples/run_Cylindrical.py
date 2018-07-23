#!/usr/bin/env python3


import matplotlib.pyplot as plt
from mcrt3d import MCRT, Params
from mcrt3d.grid import CylindricalGrid
from mcrt3d.dust import Dust
from mcrt3d.sources import Star
from mcrt3d.camera import Image, Spectrum
from mcrt3d.constants.astronomy import M_sun, R_sun, AU
from mcrt3d.constants.physics import c
from numpy import array, arange, pi, zeros, logspace
from time import time

# Create a model class.

model = MCRT()

# Set up the dust.

dust = Dust(filename="dustkappa_yso.inp", radmc3d=True)

# Set up the star.

star = Star(0.0,0.0,0.0,M_sun,R_sun,4000.0)
star.set_blackbody_spectrum(dust.nu)

# Set up the grid.

nr = 10
np = 10
nz = 10

r = arange(nr)*AU/2
p = arange(np)/(np-1.)*2*pi
z = (arange(nz)-(nz-1)/2.)*AU/1

density = zeros((nr-1,np-1,nz-1)) + 1.0e-16

model.set_cylindrical_grid(r,p,z)
model.grid.add_density(density, dust)
model.grid.add_source(star)

# Set the parameters for the run.

model.params.nphot = 100000
model.params.bw = True
model.params.scattering = False
model.params.verbose = False
model.params.use_mrw = True
model.params.mrw_gamma = 2

# Run the thermal simulation.

t1 = time()
model.run_thermal_mc()
t2 = time()
print(t2-t1)

# Run the images.

model.camera.nx = 256
model.camera.ny = 256
model.camera.pixel_size = AU/10
model.camera.lam = array([1300.])

image = model.run_image(incl=pi/4, pa=pi/4)

# Run the spectra.

model.camera.pixel_size = 10*AU/100
model.camera.lam = logspace(-1,4,200)

spectrum = model.run_spectrum(incl=0, pa=0)

# Plot the temperature structure.

for i in range(9):
    plt.imshow(model.grid.temperature[0][:,:,i], origin="lower",\
            interpolation="nearest", vmin=model.grid.temperature[0].min(),\
            vmax=model.grid.temperature[0].max())
    plt.colorbar()
    plt.show()

# Plot the images.

plt.imshow(image.intensity[:,:,0],origin="lower",interpolation="none")
plt.show()

# Plot the spectra.

plt.loglog(spectrum.lam, spectrum.intensity)
plt.axis([1e-1,1e4,1e-15,1e-5])
plt.show()

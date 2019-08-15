#!/usr/bin/env python3

from pdspy.constants.astronomy import AU, pc, Jy
import pdspy.modeling as modeling
import pdspy.dust as dust
import matplotlib.pyplot as plt
from numpy import array, arange, pi, zeros, logspace
import numpy

from mcrt3d import MCRT, Params
from mcrt3d.grid import CartesianGrid
from mcrt3d.dust import Dust
from mcrt3d.sources import Star
from mcrt3d.camera import Image, Spectrum
from mcrt3d.constants.astronomy import M_sun, R_sun, AU
from mcrt3d.constants.physics import c
from time import time

################################################################################
#
# Run the model with RADMC3D
#
################################################################################

m = modeling.Model()

# Read in the opacities.

data = numpy.loadtxt('../examples/dustkappa_yso.inp', skiprows=2)

lam = data[:,0].copy() * 1.0e-4
kabs = data[:,1].copy()
ksca = data[:,2].copy()

d = dust.Dust()
d.set_properties(lam, kabs, ksca)

# Set up the grid.

nx = 10
ny = 10
nz = 10

x = (numpy.arange(nx)-(float(nx)-1)/2)
y = (numpy.arange(ny)-(float(ny)-1)/2)
z = (numpy.arange(nz)-(float(nz)-1)/2)

m.grid.set_cartesian_grid(x, y, z)

# Set up the wavelength grid.

m.grid.lam = lam * 1.0e4

# Set the density.

dens = numpy.zeros((nx-1,ny-1,nz-1)) + 1.0e-16

m.grid.add_density(dens, d)

# Add a star.

source = modeling.Star(mass=0.5, luminosity=0.22658222, temperature=4000.)

m.grid.add_star(source)

# Run the thermal simulation.

t1 = time()
m.run_thermal(code="radmc3d", nphot=1e6, verbose=False)
t2 = time()
print(t2-t1)

# Run the image.

m.run_image(name="image", nphot=1e5, npix=256, pixelsize=0.1, lam="1", \
        phi=0, incl=0, code="radmc3d", dpc=1, verbose=False)

# Run the SED.

m.set_camera_wavelength(numpy.logspace(-1,4,200))

m.run_sed(name="SED", nphot=1e4, loadlambda=True, incl=0, pa=0, \
        dpc=1, code="radmc3d", camera_scatsrc_allfreq=True, \
        verbose=False)

################################################################################
#
# Run the model with MCRT3D
#
################################################################################

# Create a model class.

model = MCRT()

# Set up the dust.

dust = Dust(filename="../examples/dustkappa_yso.inp", radmc3d=True)

# Set up the star.

star = Star(0.0,0.0,0.0,M_sun,R_sun,4000.0)
star.set_blackbody_spectrum(dust.nu)

# Set up the grid.

nx = 10
ny = 10
nz = 10

x = (numpy.arange(nx)-(float(nx)-1)/2)*AU/1
y = (numpy.arange(ny)-(float(ny)-1)/2)*AU/1
z = (numpy.arange(nz)-(float(nz)-1)/2)*AU/1

density = numpy.zeros((nx-1,ny-1,nz-1)) + 1.0e-16

model.set_cartesian_grid(x,y,z)
model.grid.add_density(density, dust)
model.grid.add_source(star)

# Set the parameters for the run.

model.params.nphot = 1000000
model.params.bw = True
model.params.scattering = False
model.params.verbose = False
model.params.use_mrw = False
model.params.mrw_gamma = 2

# Run the thermal simulation.

t1 = time()
model.run_thermal_mc()
t2 = time()
print(t2-t1)

# Run the images.

model.params.nphot = 100000
model.camera.nx = 256
model.camera.ny = 256
model.camera.pixel_size = AU/10
model.camera.lam = array([1.])

image = model.run_image(incl=0, pa=0, dpc=1.)

# Run the spectra.

model.params.nphot = 10000
model.camera.pixel_size = 10*AU/100
model.camera.lam = logspace(-1,4,200)

spectrum = model.run_spectrum(incl=0, pa=0, dpc=1.)

################################################################################
#
# Make some plots.
#
################################################################################

# Plot the temperature structure.

for i in range(9):
    fig, ax = plt.subplots(nrows=1, ncols=3)

    vmin = min(m.grid.temperature[0].min(), model.grid.temperature[0].min())
    vmax = min(m.grid.temperature[0].max(), model.grid.temperature[0].max())

    diff = (m.grid.temperature[0][:,:,i] - model.grid.temperature[0][:,:,i]) / \
            m.grid.temperature[0][:,:,i] * 100

    im1 = ax[0].imshow(m.grid.temperature[0][:,:,i], origin="lower", \
            interpolation="nearest", vmin=vmin, vmax=vmax)

    im2 = ax[1].imshow(model.grid.temperature[0][:,:,i], origin="lower",\
            interpolation="nearest", vmin=vmin, vmax=vmax)

    im3 = ax[2].imshow(diff, origin="lower", interpolation="nearest")

    fig.colorbar(im1, ax=ax[1], fraction=0.046)
    fig.colorbar(im3, ax=ax[2], fraction=0.046)

    fig.subplots_adjust(left=0.1, right=0.95, wspace=0.25)
    fig.set_size_inches((9,3))

    plt.show()

# Plot the scattering function.

for i in range(9):
    fig, ax = plt.subplots(nrows=1, ncols=1)

    vmin = model.grid.scatt[0][:,:,i,0].min()
    vmax = model.grid.scatt[0][:,:,i,0].max()

    ax.imshow(model.grid.scatt[0][:,:,i,0], origin="lower", \
            interpolation="nearest", vmin=vmin, vmax=vmax)

    plt.show()

# Plot the images.

fig, ax = plt.subplots(nrows=1, ncols=3, sharex=True, sharey=True)

diff = (numpy.log10(m.images["image"].image[:,:,0,0]) - \
        numpy.log10(image.intensity[:,:,0])) / \
        numpy.log10(m.images["image"].image[:,:,0,0]) * 100

im1 = ax[0].imshow(numpy.log10(m.images["image"].image[:,:,0,0]), \
        origin="lower", interpolation="none")

im2 = ax[1].imshow(numpy.log10(image.intensity[:,:,0]), origin="lower", \
        interpolation="none")

im3 = ax[2].imshow(diff, origin="lower", interpolation="none")

fig.colorbar(im3, ax=ax[2], fraction=0.046)

fig.subplots_adjust(left=0.05, right=0.95, wspace=0.25)
fig.set_size_inches((9,3))

plt.show()

# Plot the spectra.

fig, ax = plt.subplots(nrows=1, ncols=1)

ax.loglog(m.spectra["SED"].wave, m.spectra["SED"].flux)
ax.loglog(spectrum.lam, spectrum.intensity)

plt.show()

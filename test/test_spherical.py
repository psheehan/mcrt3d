#!/usr/bin/env python3

from pdspy.constants.astronomy import AU, pc, Jy, M_sun, R_sun
import pdspy.modeling as modeling
import pdspy.dust as dust
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import numpy

from mcrt3d import MCRT
from mcrt3d import IsotropicDust
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

nr = 10
nt = 11
np = 10

r = numpy.arange(nr)/2
r[0] = 0.01
t = numpy.arange(nt)/(nt-1.)*numpy.pi
p = numpy.arange(np)/(np-1.)*2*numpy.pi

m.grid.set_spherical_grid(r, t, p)

# Set up the wavelength grid.

m.grid.lam = lam * 1.0e4

# Set the density.

dens = numpy.zeros((nr-1,nt-1,np-1)) + 1.0e-16

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

model.run_unstructured_image(numpy.array([1300.]), 10, 10, 2.5, 100000, \
        incl=0., pa=0., dpc=1.)

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

# Set up the grid.

nr = 10
nt = 11
np = 10

r = numpy.arange(nr)*AU/2
t = numpy.arange(nt)/(nt-1.)*numpy.pi
p = numpy.arange(np)/(np-1.)*2*numpy.pi

model.set_spherical_grid(r,t,p)

# Set up the dust.

data = numpy.loadtxt('../examples/dustkappa_yso.inp', skiprows=2)

lam = data[::-1,0].copy() * 1.0e-4
kabs = data[::-1,1].copy()
ksca = data[::-1,2].copy()

d = IsotropicDust(lam, kabs, ksca)

# Set up the density.

density = numpy.zeros((nr-1,nt-1,np-1)) + 1.0e-16

model.grid.add_density(density, d)

# Set up the star.

model.grid.add_star(mass=M_sun, radius=R_sun, temperature=4000.)
model.grid.sources[-1].set_blackbody_spectrum(lam)

# Run the thermal simulation.

t1 = time()
model.thermal_mc(nphot=1000000, bw=True, use_mrw=False, mrw_gamma=2, \
        verbose=False)
t2 = time()
print(t2-t1)

# Run the images.

model.run_image(numpy.array([1.]), 256, 256, 0.1, 100000, incl=0., pa=0, dpc=1.)

# Run the spectra.

model.run_spectrum(numpy.logspace(-1,4,200), 10000, incl=0, pa=0, dpc=1.)

################################################################################
#
# Make some plots.
#
################################################################################

# Plot the temperature structure.

for i in range(9):
    fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(11,3), \
            gridspec_kw=dict(left=0.05, right=0.95, wspace=0.25))

    vmin = min(m.grid.temperature[0].min(), model.grid.temperature[0].min())
    vmax = min(m.grid.temperature[0].max(), model.grid.temperature[0].max())

    diff = (m.grid.temperature[0][:,:,i] - model.grid.temperature[0][:,:,i]) / \
            m.grid.temperature[0][:,:,i] * 100

    im1 = ax[0].imshow(m.grid.temperature[0][:,:,i], origin="lower", \
            interpolation="nearest", vmin=vmin, vmax=vmax)

    im2 = ax[1].imshow(model.grid.temperature[0][:,:,i], origin="lower",\
            interpolation="nearest", vmin=vmin, vmax=vmax)

    im3 = ax[2].imshow(diff, origin="lower", interpolation="nearest")

    ax[0].set_title("RADMC-3D")
    ax[1].set_title("MCRT3D")
    ax[2].set_title("RADMC-3D - MCRT3D")

    fig.colorbar(im1, ax=ax[1], fraction=0.046)
    fig.colorbar(im3, ax=ax[2], fraction=0.046)

    plt.show()

"""
# Plot the scattering function.

for i in range(9):
    fig, ax = plt.subplots(nrows=1, ncols=1)

    with numpy.errstate(divide="ignore", invalid="ignore"):
        vmin = numpy.ma.masked_invalid(numpy.log10(\
                model.grid.scatt[0][:,:,i,0])).min()
        vmax = numpy.ma.masked_invalid(numpy.log10(\
                model.grid.scatt[0][:,:,i,0])).max()

        ax.imshow(numpy.log10(model.grid.scatt[0][:,:,i,0]), origin="lower", \
                interpolation="nearest", vmin=vmin, vmax=vmax)

    plt.show()
"""

# Plot the images.

fig, ax = plt.subplots(nrows=1, ncols=3, sharex=True, sharey=True, \
        figsize=(10,3), gridspec_kw=dict(left=0.05, right=0.95, wspace=0.25))

with numpy.errstate(divide="ignore", invalid="ignore"):
    diff = (numpy.log10(m.images["image"].image[:,:,0,0]) - \
            numpy.log10(model.images[0].intensity[:,:,0])) / \
            numpy.log10(m.images["image"].image[:,:,0,0]) * 100

    im1 = ax[0].imshow(numpy.log10(m.images["image"].image[:,:,0,0]), \
            origin="lower", interpolation="none")

    im2 = ax[1].imshow(numpy.log10(model.images[0].intensity[:,:,0]), \
            origin="lower", interpolation="none")

    im3 = ax[2].imshow(diff, origin="lower", interpolation="none")

ax[0].set_title("RADMC-3D")
ax[1].set_title("MCRT3D")
ax[2].set_title("RADMC-3D - MCRT3D")

fig.colorbar(im3, ax=ax[2], fraction=0.046)

plt.show()

# Plot the unstructured image.

triang = tri.Triangulation(model.images[1].x/AU, model.images[1].y/AU)

plt.tripcolor(triang, model.images[1].intensity[:,0], "ko-")
plt.triplot(triang, "k.-", linewidth=0.1, markersize=0.1)

plt.axes().set_aspect("equal")

plt.xlabel("x", fontsize=14)
plt.ylabel("y", fontsize=14)

plt.axes().tick_params(labelsize=14)

plt.show()

# Plot the spectra.

fig, ax = plt.subplots(nrows=1, ncols=1)

ax.loglog(m.spectra["SED"].wave, m.spectra["SED"].flux)
ax.loglog(model.spectra[0].lam, model.spectra[0].intensity)

ax.set_ylim(1.0e-23,1.0e7)

plt.show()

#!/usr/bin/env python3

from pdspy.constants.astronomy import AU, pc, Jy, M_sun, R_sun
import pdspy.modeling as modeling
import pdspy.dust as dust
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import scipy.optimize
import scipy.integrate
import numpy

from mcrt3d import MCRT
from mcrt3d import IsotropicDust
from time import time

################################################################################
#
# DEfine functions to calculate the densities.
#
################################################################################

def disk_density(r, theta, phi):
    rr = 0.5*(r[1:] + r[0:-1])
    tt = 0.5*(theta[1:] + theta[0:-1])
    pp = 0.5*(phi[1:] + phi[0:-1])

    rr, tt, pp = numpy.meshgrid(rr, tt, pp, indexing='ij')

    R = rr*numpy.sin(tt)
    z = rr*numpy.cos(tt)

    Mdisk = 0.00001
    gamma = 1.
    Rin = 0.1
    Rdisk = 150.
    h0 = 0.1
    beta = 1.1

    Sigma0 = (2-gamma)*Mdisk*M_sun / (2*numpy.pi*Rdisk**2*AU**2) * \
            numpy.exp((Rin/Rdisk)**(2-gamma))

    Sigma = Sigma0 * (R / Rdisk)**(-gamma) * numpy.exp(-(R / Rdisk)**(2-gamma))
    h = h0*(R / 1.)**beta

    rho = Sigma / (numpy.sqrt(2*numpy.pi)*h*AU) * numpy.exp(-0.5*(z / h)**2)

    rho[R < 0.1] *= 1.0e-50

    return rho

def envelope_density(r, theta, phi):
    #numpy.seterr(all='ignore')
    ##### Envelope parameters

    rin = 0.1
    rout = 1000.
    mass = 0.00001
    rcent = 150.
    cavz0 = 1.
    cavpl = 1.5
    cavrfact = 0.5
    theta_open = 60. * numpy.pi / 180.

    cava = (rout*numpy.sin(theta_open/2) / numpy.tan(theta_open/2) - \
            cavz0) * (rout*numpy.sin(theta_open/2))**-cavpl

    # Set up the coordinates.

    rr = 0.5*(r[1:] + r[0:-1])
    tt = 0.5*(theta[1:] + theta[0:-1])
    pp = 0.5*(phi[1:] + phi[0:-1])

    rr, tt, pp = numpy.meshgrid(rr, tt, pp, indexing='ij')

    mu = numpy.cos(tt)

    RR = rr * numpy.sin(tt)
    zz = rr * numpy.cos(tt)

    # Calculate mu0 at each r, theta combination.

    def func(mu0,r,mu,R_c):
        return mu0**3-mu0*(1-r/R_c)-mu*(r/R_c)

    mu0 = mu*0.
    for ir in range(rr.shape[0]):
        for it in range(rr.shape[1]):
            mu0[ir,it,0] = scipy.optimize.brenth(func, -1.0, 1.0, \
                    args=(rr[ir,it,0], mu[ir,it,0], rcent))

    ##### Make the dust density model for an Ulrich envelope.

    rho0 = 1.0

    rho = rho0 * (rr / rcent)**(-1.5) * (1 + mu/mu0)**(-0.5)* \
            (mu/mu0 + 2*mu0**2 * rcent/rr)**(-1)

    mid1 = (numpy.abs(mu) < 1.0e-10) & (rr < rcent)
    rho[mid1] = rho0 * (rr[mid1] / rcent)**(-0.5) * \
            (1. - rr[mid1] / rcent)**(-1) / 2.

    mid2 = (numpy.abs(mu) < 1.0e-10) & (rr > rcent)
    rho[mid2] = rho0 * (2.*rr[mid2]/rcent - 1)**(-0.5) * \
            (rr[mid2]/rcent - 1.)**(-1)

    rho[(rr >= rout) ^ (rr <= rin)] *= 1.0e-50

    ##### Add an outflow cavity.

    rho[numpy.abs(zz)-cavz0-cava*(RR)**cavpl > 0.0] *= cavrfact

    ##### Normalize the mass correctly.

    if theta.max() > numpy.pi/2:
        mdot = mass*M_sun/(2*numpy.pi*scipy.integrate.trapz(\
                scipy.integrate.trapz(rho*(rr*AU)**2*numpy.sin(tt), tt, \
                axis=1), rr[:,0,:]*AU, axis=0))[0]
    else:
        mdot = mass*M_sun/(4*numpy.pi*scipy.integrate.trapz(\
                scipy.integrate.trapz(rho*(rr*AU)**2*numpy.sin(tt), tt, \
                axis=1), rr[:,0,:]*AU, axis=0))[0]
    rho *= mdot

    #numpy.seterr(all='warn')

    return rho

################################################################################
#
# Run the model with RADMC3D
#
################################################################################

# Create a model class.

m = MCRT()

# Set up the grid.

nr = 100
nt = 101
np = 2

r = numpy.logspace(-1.,3.5,100)
r = numpy.concatenate(([0],r))
t = numpy.arange(nt)/(nt-1.)*numpy.pi/2
p = numpy.arange(np)/(np-1.)*2*numpy.pi

m.set_spherical_grid(r*AU,t,p)

# Set up the dust.

data = numpy.loadtxt('../examples/dustkappa_yso.inp', skiprows=2)

lam = data[::-1,0].copy() * 1.0e-4
kabs = data[::-1,1].copy()
ksca = data[::-1,2].copy()

d = IsotropicDust(lam, kabs, ksca)

# Set up the density.

density_disk = disk_density(r, t, p)

m.grid.add_density(density_disk, d)

density_envelope = envelope_density(r, t, p)

m.grid.add_density(density_envelope, d)

# Set up the star.

m.grid.add_star(mass=M_sun, radius=R_sun, temperature=4000.)
m.grid.sources[-1].set_blackbody_spectrum(lam)

# Run the thermal simulation.

t1 = time()
m.thermal_mc(nphot=1000000, bw=True, use_mrw=False, mrw_gamma=4, \
        verbose=False, nthreads=1)
t2 = time()
serial_thermal_time = t2-t1

# Run the scattering phase function calculation.

t1 = time()
m.scattering_mc(numpy.array([4.]), nphot=100000, verbose=False, nthreads=1)
t2 = time()
serial_scattering_time = t2-t1

# Run the images.

t1 = time()
m.run_image(numpy.array([4.]), 256, 256, 0.1, 100000, incl=0., pa=0, \
        dpc=140., nthreads=1)
t2 = time()
serial_image_time = t2-t1

t1 = time()
m.run_unstructured_image(numpy.array([1300.]), 25, 25, 2.0, 100000, \
        incl=0., pa=0., dpc=140., nthreads=1)
t2 = time()
serial_unimage_time = t2-t1

# Run the spectra.

t1 = time()
m.run_spectrum(numpy.logspace(-1,4,200), 10000, incl=0, pa=0, dpc=140., \
        nthreads=1)
t2 = time()
serial_sed_time = t2-t1

################################################################################
#
# Run the model with MCRT3D
#
################################################################################

# Create a model class.

model = MCRT()

# Set up the grid.

nr = 100
nt = 101
np = 2

r = numpy.logspace(-1.,3.5,100)
r = numpy.concatenate(([0],r))
t = numpy.arange(nt)/(nt-1.)*numpy.pi/2
p = numpy.arange(np)/(np-1.)*2*numpy.pi

model.set_spherical_grid(r*AU,t,p)

# Set up the dust.

data = numpy.loadtxt('../examples/dustkappa_yso.inp', skiprows=2)

lam = data[::-1,0].copy() * 1.0e-4
kabs = data[::-1,1].copy()
ksca = data[::-1,2].copy()

d = IsotropicDust(lam, kabs, ksca)

# Set up the density.

density_disk = disk_density(r, t, p)

model.grid.add_density(density_disk, d)

density_envelope = envelope_density(r, t, p)

model.grid.add_density(density_envelope, d)

# Set up the star.

model.grid.add_star(mass=M_sun, radius=R_sun, temperature=4000.)
model.grid.sources[-1].set_blackbody_spectrum(lam)

# Run the thermal simulation.

t1 = time()
model.thermal_mc(nphot=1000000, bw=True, use_mrw=False, mrw_gamma=4, \
        verbose=False, nthreads=4)
t2 = time()
parallel_thermal_time = t2-t1

# Run the scattering phase function calculation.

t1 = time()
model.scattering_mc(numpy.array([4.]), nphot=100000, verbose=False, nthreads=4)
t2 = time()
parallel_scattering_time = t2-t1

# Run the images.

t1 = time()
model.run_image(numpy.array([4.]), 256, 256, 0.1, 100000, incl=0., pa=0, \
        dpc=140., nthreads=4)
t2 = time()
parallel_image_time = t2-t1

t1 = time()
model.run_unstructured_image(numpy.array([1300.]), 25, 25, 2.0, 100000, \
        incl=0., pa=0., dpc=140., nthreads=4)
t2 = time()
parallel_unimage_time = t2-t1

# Run the spectra.

t1 = time()
model.run_spectrum(numpy.logspace(-1,4,200), 10000, incl=0, pa=0, dpc=140., \
        nthreads=4)
t2 = time()
parallel_sed_time = t2-t1

################################################################################
#
# Print the timing results.
#
################################################################################

print()
print("                          Serial      Parallel")
print("Thermal simulation        {0:6.2f}      {1:8.2f}".\
        format(serial_thermal_time, parallel_thermal_time))
print("Scattering simulation     {0:6.2f}      {1:8.2f}".\
        format(serial_scattering_time, parallel_scattering_time))
print("Imaging                   {0:6.2f}      {1:8.2f}".\
        format(serial_image_time, parallel_image_time))
print("Unstructured Imaging      {0:6.2f}      {1:8.2f}".\
        format(serial_unimage_time, parallel_unimage_time))
print("SED                       {0:6.2f}      {1:8.2f}".\
        format(serial_sed_time, parallel_sed_time))
print()

################################################################################
#
# Make some plots.
#
################################################################################

# Plot the temperature structure.

fig, ax = plt.subplots(nrows=2, ncols=3, figsize=(11,6), \
        gridspec_kw=dict(left=0.05, right=0.95, wspace=0.15, top=0.9, \
        bottom=0.05))

vmin = min(m.grid.temperature[0].min(), model.grid.temperature[0].min(), \
        m.grid.temperature[1].min(), model.grid.temperature[1].min())
vmax = max(m.grid.temperature[0].max(), model.grid.temperature[0].max(), \
        m.grid.temperature[1].max(), model.grid.temperature[1].max())

diff1 = (m.grid.temperature[0][:,:,0] - model.grid.temperature[0][:,:,0]) / \
        m.grid.temperature[0][:,:,0] * 100
diff2 = (m.grid.temperature[1][:,:,0] - model.grid.temperature[1][:,:,0]) / \
        m.grid.temperature[1][:,:,0] * 100

im1 = ax[0,0].imshow(m.grid.temperature[0][:,:,0], origin="lower", \
        interpolation="nearest", vmin=vmin, vmax=vmax)

im2 = ax[0,1].imshow(model.grid.temperature[0][:,:,0], origin="lower",\
        interpolation="nearest", vmin=vmin, vmax=vmax)

im3 = ax[0,2].imshow(diff1, origin="lower", interpolation="nearest", \
        vmin=-3, vmax=3)

im4 = ax[1,0].imshow(m.grid.temperature[1][:,:,0], origin="lower", \
        interpolation="nearest", vmin=vmin, vmax=vmax)

im5 = ax[1,1].imshow(model.grid.temperature[1][:,:,0], origin="lower",\
        interpolation="nearest", vmin=vmin, vmax=vmax)

im6 = ax[1,2].imshow(diff2, origin="lower", interpolation="nearest")

ax[0,0].set_title("Serial")
ax[0,1].set_title("Parallel")
ax[0,2].set_title("Serial - Parallel")

fig.colorbar(im1, ax=ax[0,1], fraction=0.046)
fig.colorbar(im3, ax=ax[0,2], fraction=0.046)
fig.colorbar(im4, ax=ax[1,1], fraction=0.046)
fig.colorbar(im6, ax=ax[1,2], fraction=0.046)

plt.show()

# Plot the scattering phase function.

fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(11,3), \
        gridspec_kw=dict(left=0.05, right=0.95, wspace=0.25))

with numpy.errstate(divide="ignore", invalid="ignore"):
    vmin = min(numpy.nanmin(numpy.log10(m.grid.scatt[0])[\
            numpy.isfinite(numpy.log10(m.grid.scatt[0]))]), \
            numpy.nanmin(numpy.log10(model.grid.scatt[0])[numpy.isfinite(\
            numpy.log10(model.grid.scatt[0]))]))
    vmax = max(numpy.nanmax(numpy.log10(m.grid.scatt[0])), \
            numpy.nanmax(numpy.log10(model.grid.scatt[0])))

    diff = numpy.log10(m.grid.scatt[0][:,:,0,0]) - \
            numpy.log10(model.grid.scatt[0][:,:,0,0])

    im1 = ax[0].imshow(numpy.log10(m.grid.scatt[0][:,:,0,0]), \
            origin="lower", interpolation="nearest", vmin=vmin, vmax=vmax)

    im2 = ax[1].imshow(numpy.log10(model.grid.scatt[0][:,:,0,0]), \
            origin="lower", interpolation="nearest", vmin=vmin, vmax=vmax)

    im3 = ax[2].imshow(diff, origin="lower", interpolation="nearest", \
            vmin=-0.5, vmax=0.5)

ax[0].set_title("Serial")
ax[1].set_title("Parallel")
ax[2].set_title("Serial - Parallel")

fig.colorbar(im1, ax=ax[1], fraction=0.046)
fig.colorbar(im3, ax=ax[2], fraction=0.046)

plt.show()

# Plot the images.

fig, ax = plt.subplots(nrows=1, ncols=3, sharex=True, sharey=True, \
        figsize=(10,3), gridspec_kw=dict(left=0.05, right=0.95, wspace=0.25))

with numpy.errstate(divide="ignore", invalid="ignore"):
    vmax = min(numpy.log10(numpy.nanmax(m.images[0].intensity[:,:,0])), \
            numpy.log10(numpy.nanmax(model.images[0].intensity[:,:,0])))
    vmin = vmax - 10.

    diff = (numpy.log10(m.images[0].intensity[:,:,0]) - \
            numpy.log10(model.images[0].intensity[:,:,0]))

    im1 = ax[0].imshow(numpy.log10(m.images[0].intensity[:,:,0]), \
            origin="lower", interpolation="none", vmin=vmin, vmax=vmax)

    im2 = ax[1].imshow(numpy.log10(model.images[0].intensity[:,:,0]), \
            origin="lower", interpolation="none", vmin=vmin, vmax=vmax)

    im3 = ax[2].imshow(diff, origin="lower", interpolation="none")

ax[0].set_title("Serial")
ax[1].set_title("Parallel")
ax[2].set_title("Serial - Parallel")

fig.colorbar(im3, ax=ax[2], fraction=0.046)

plt.show()

# Plot the unstructured image.

fig, ax = plt.subplots(nrows=1, ncols=2, sharex=True, sharey=True, \
        figsize=(9,4.5), gridspec_kw=dict(left=0.1, right=0.95, wspace=0.25, \
        bottom=0.15))

triang1 = tri.Triangulation(m.images[1].x/AU, m.images[1].y/AU)

ax[0].tripcolor(triang1, m.images[1].intensity[:,0], "ko-")
ax[0].triplot(triang1, "k.-", linewidth=0.1, markersize=0.1)

triang2 = tri.Triangulation(model.images[1].x/AU, model.images[1].y/AU)

ax[1].tripcolor(triang2, model.images[1].intensity[:,0], "ko-")
ax[1].triplot(triang2, "k.-", linewidth=0.1, markersize=0.1)

for axes in ax:
    axes.set_aspect("equal")

    axes.set_xlabel("x", fontsize=14)
    axes.set_ylabel("y", fontsize=14)

    axes.tick_params(labelsize=14)

ax[0].set_title("Serial")
ax[1].set_title("Parallel")

plt.show()

# Plot the spectra.

fig, ax = plt.subplots(nrows=1, ncols=1)

ax.loglog(m.spectra[0].lam, m.spectra[0].intensity)
ax.loglog(model.spectra[0].lam, model.spectra[0].intensity)

ax.set_ylim(1.0e-6,1.0e1)

plt.show()

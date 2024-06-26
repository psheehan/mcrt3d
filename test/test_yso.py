#!/usr/bin/env python3

from pdspy.constants.astronomy import AU, pc, Jy, M_sun, R_sun
from pdspy.constants.physics import m_p, G
import pdspy.modeling as modeling
import pdspy.dust as dust
import pdspy.gas as gas
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import scipy.optimize
import scipy.integrate
import numpy

from mcrt3d import MCRT
from mcrt3d import IsotropicDust
from mcrt3d import load_gas_properties_lamda
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

def disk_velocity(r, theta, phi, mstar=0.5):
    rr = 0.5*(r[1:] + r[0:-1])
    tt = 0.5*(theta[1:] + theta[0:-1])
    pp = 0.5*(phi[1:] + phi[0:-1])

    mstar *= M_sun

    rt, tt, pp = numpy.meshgrid(rr*AU, tt, pp, indexing='ij')

    rr = rt*numpy.sin(tt)
    zz = rt*numpy.cos(tt)

    v_r = numpy.zeros(rr.shape)
    v_theta = numpy.zeros(rr.shape)
    v_phi = numpy.sqrt(G*mstar*rr**2/rt**3)

    return numpy.array((v_r, v_theta, v_phi))

def disk_microturbulence(r, theta, phi):
    rr = 0.5*(r[1:] + r[0:-1])
    tt = 0.5*(theta[1:] + theta[0:-1])
    pp = 0.5*(phi[1:] + phi[0:-1])

    ##### Set up the coordinates

    rt, tt, pp = numpy.meshgrid(rr*AU, tt, pp, indexing='ij')

    rr = rt*numpy.sin(tt)
    zz = rt*numpy.cos(tt)

    ##### Make the dust density model for a protoplanetary disk.

    aturb = numpy.ones(rr.shape)*0.05*1.0e5

    return aturb

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


def envelope_velocity(r, theta, phi, mstar=0.5):
    mstar *= M_sun
    rcent = 150.*AU

    # Set up the coordinates.

    rr = 0.5*(r[1:] + r[0:-1])*AU
    tt = 0.5*(theta[1:] + theta[0:-1])
    pp = 0.5*(phi[1:] + phi[0:-1])

    rr, tt, pp = numpy.meshgrid(rr, tt, pp, indexing='ij')

    mu = numpy.cos(tt)

    # Calculate mu0 at each r, theta combination.

    def func(mu0,r,mu,R_c):
        return mu0**3-mu0*(1-r/R_c)-mu*(r/R_c)

    mu0 = mu*0.
    for ir in range(rr.shape[0]):
        for it in range(rr.shape[1]):
            mu0[ir,it,0] = scipy.optimize.brenth(func, -1.0, 1.0, \
                    args=(rr[ir,it,0], mu[ir,it,0], rcent))

    v_r = -numpy.sqrt(G*mstar/rr)*numpy.sqrt(1 + mu/mu0)
    v_theta = numpy.sqrt(G*mstar/rr) * (mu0 - mu) * \
            numpy.sqrt((mu0 + mu) / (mu0 * numpy.sin(tt)**2))
    v_phi = numpy.sqrt(G*mstar/rr) * numpy.sqrt((1 - mu0**2)/(1 - mu**2)) *\
            numpy.sqrt(1 - mu/mu0)

    return numpy.array((v_r, v_theta, v_phi))

def envelope_microturbulence(r, theta, phi):
    ##### Set up the coordinates

    rr = 0.5*(r[1:] + r[0:-1])
    tt = 0.5*(theta[1:] + theta[0:-1])
    pp = 0.5*(phi[1:] + phi[0:-1])

    rt, tt, pp = numpy.meshgrid(rr*AU, tt, pp, indexing='ij')

    rr = rt*numpy.sin(tt)
    zz = rt*numpy.cos(tt)

    ##### Make the dust density model for a protoplanetary disk.

    aturb = numpy.ones(rr.shape)*1.0e5*0.05

    return aturb

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

# Set up the gas.

g = gas.Gas()
g.set_properties_from_lambda("Data/co.dat")
abundance = 1.0e-4

# Set up the grid.

nr = 100
nt = 101
np = 2

r = numpy.logspace(-1.,3.5,100)
t = numpy.arange(nt)/(nt-1.)*numpy.pi/2
p = numpy.arange(np)/(np-1.)*2*numpy.pi

m.grid.set_spherical_grid(r, t, p)

# Set up the wavelength grid.

m.grid.lam = lam * 1.0e4

# Set the density.

density_disk = disk_density(r, t, p)

m.grid.add_density(density_disk, d)

density_envelope = envelope_density(r, t, p)

m.grid.add_density(density_envelope, d)

# Set up the gas density.

gas_density_disk = density_disk * 100 / (2.37 * m_p) * abundance
velocity_disk = disk_velocity(r, t, p)
microturbulence_disk = disk_microturbulence(r, t, p)

m.grid.add_number_density(gas_density_disk, g)
m.grid.add_velocity(velocity_disk)
m.grid.add_microturbulence(microturbulence_disk)

gas_density_envelope = density_envelope * 100 / (2.37 * m_p) * abundance
velocity_envelope = envelope_velocity(r, t, p)
microturbulence_envelope = envelope_microturbulence(r, t, p)

m.grid.add_number_density(gas_density_envelope, g)
m.grid.add_velocity(velocity_envelope)
m.grid.add_microturbulence(microturbulence_envelope)

# Add a star.

source = modeling.Star(mass=0.5, luminosity=0.22658222, temperature=4000.)

m.grid.add_star(source)

# Run the thermal simulation.

t1 = time()
m.run_thermal(code="radmc3d", nphot=1e6, verbose=False)
t2 = time()
print(t2-t1)

# Run the scattering phase function calculation.

m.set_camera_wavelength(numpy.array([4.]))

m.run_scattering(code="radmc3d", nphot=1e5, verbose=False, loadlambda=True, \
        tgas_eq_tdust=True)

# Run the image.

m.run_image(name="image", nphot=1e5, npix=256, pixelsize=0.1, lam="4", \
        phi=0, incl=0, code="radmc3d", dpc=140., tgas_eq_tdust=True, \
        verbose=False)

# Run a CO channel map.

m.run_image(name="CO2-1", nphot=0, npix=256, pixelsize=0.1, lam=None, \
        imolspec=1, iline=2, widthkms=2.5, linenlam=10, tgas_eq_tdust=True, \
        scattering_mode_max=0, incl_dust=False, incl=45, pa=0, code="radmc3d",\
        dpc=140, verbose=False, writeimage_unformatted=True)

# Run the SED.

m.set_camera_wavelength(numpy.logspace(-1,4,200))

m.run_sed(name="SED", nphot=1e4, loadlambda=True, incl=0, pa=0, \
        dpc=140., code="radmc3d", camera_scatsrc_allfreq=True, \
        tgas_eq_tdust=True, verbose=False)

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

# Set up the gas.

g = load_gas_properties_lamda("Data/co.dat")
abundance = 1.0e-4

# Set up the density.

density_disk = disk_density(r, t, p)

model.grid.add_density(density_disk, d)

density_envelope = envelope_density(r, t, p)

model.grid.add_density(density_envelope, d)

# Set up the gas density.

gas_density_disk = density_disk * 100 / (2.37 * m_p) * abundance
velocity_disk = disk_velocity(r, t, p)
microturbulence_disk = disk_microturbulence(r, t, p)

gas_density_envelope = density_envelope * 100 / (2.37 * m_p) * abundance
velocity_envelope = envelope_velocity(r, t, p)
microturbulence_envelope = envelope_microturbulence(r, t, p)

velocity = numpy.where(gas_density_disk + gas_density_envelope > 0, 
        (velocity_disk * gas_density_disk + velocity_envelope * \
        gas_density_envelope) / (gas_density_disk + gas_density_envelope), 0.)

model.grid.add_number_density(gas_density_disk, velocity, \
        microturbulence_disk, g)

model.grid.add_number_density(gas_density_envelope, velocity, \
        microturbulence_envelope, g)

# Set up the star.

model.grid.add_star(mass=M_sun, radius=R_sun, temperature=4000.)
model.grid.sources[-1].set_blackbody_spectrum(lam)

# Run the thermal simulation.

t1 = time()
model.thermal_mc(nphot=1000000, bw=True, use_mrw=False, mrw_gamma=2, \
        verbose=False)
t2 = time()
print(t2-t1)

# Run the scattering phase function calculation.

model.scattering_mc(numpy.array([4.]), nphot=100000, verbose=False)

# Run the images.

model.run_image("image", numpy.array([4.]), 256, 256, 0.1, 100000, incl=0., 
        pa=0, dpc=140.)

model.run_unstructured_image("uimage", numpy.array([1300.]), 25, 25, 2.0, \
        100000, incl=0., pa=0., dpc=140.)

model.run_image("CO2-1", m.images["CO2-1"].wave/1.0e-4, 256, 256, 0.1, 10000, 
        incl=45., pa=0, dpc=140., raytrace_dust=False, raytrace_gas=True)

# Run the spectra.

model.run_spectrum("SED", numpy.logspace(-1,4,200), 10000, incl=0, pa=0, \
        dpc=140.)

################################################################################
#
# Make some plots.
#
################################################################################

# Plot the temperature structure.

fig, ax = plt.subplots(nrows=2, ncols=3, figsize=(11,6), \
        gridspec_kw=dict(left=0.05, right=0.95, wspace=0.25))

vmin = min(m.grid.temperature[0].min(), model.grid.temperature[0].min(), \
        m.grid.temperature[1].min(), model.grid.temperature[1].min())
vmax = max(m.grid.temperature[0].max(), model.grid.temperature[0].max(), \
        m.grid.temperature[1].max(), model.grid.temperature[1].max())

diff1 = (m.grid.temperature[0][:,:,0] - model.grid.temperature[0][1:,:,0]) / \
        m.grid.temperature[0][:,:,0] * 100
diff2 = (m.grid.temperature[1][:,:,0] - model.grid.temperature[1][1:,:,0]) / \
        m.grid.temperature[1][:,:,0] * 100

im1 = ax[0,0].imshow(m.grid.temperature[0][:,:,0], origin="lower", \
        interpolation="nearest", vmin=vmin, vmax=vmax)

im2 = ax[0,1].imshow(model.grid.temperature[0][1:,:,0], origin="lower",\
        interpolation="nearest", vmin=vmin, vmax=vmax)

im3 = ax[0,2].imshow(diff1, origin="lower", interpolation="nearest")

im4 = ax[1,0].imshow(m.grid.temperature[1][:,:,0], origin="lower", \
        interpolation="nearest", vmin=vmin, vmax=vmax)

im5 = ax[1,1].imshow(model.grid.temperature[1][1:,:,0], origin="lower",\
        interpolation="nearest", vmin=vmin, vmax=vmax)

im6 = ax[1,2].imshow(diff2, origin="lower", interpolation="nearest")

ax[0,0].set_title("RADMC-3D")
ax[0,1].set_title("MCRT3D")
ax[0,2].set_title("RADMC-3D - MCRT3D")

fig.colorbar(im1, ax=ax[0,1], fraction=0.046)
fig.colorbar(im3, ax=ax[0,2], fraction=0.046)
fig.colorbar(im4, ax=ax[1,1], fraction=0.046)
fig.colorbar(im6, ax=ax[1,2], fraction=0.046)

plt.show()

# Plot the scattering phase function.

fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(11,3), \
        gridspec_kw=dict(left=0.05, right=0.95, wspace=0.25))

with numpy.errstate(divide="ignore", invalid="ignore"):
    vmin = min(numpy.nanmin(numpy.log10(m.grid.scattering_phase[0])[\
            numpy.isfinite(numpy.log10(m.grid.scattering_phase[0]))]), \
            numpy.nanmin(numpy.log10(model.grid.scatt[0]*\
            (model.grid.density[0]*d.kext[-40]+\
            model.grid.density[1]*d.kext[-40]))[numpy.isfinite(\
            numpy.log10(model.grid.scatt[0]*(model.grid.density[0]*\
            d.kext[-40] + model.grid.density[1]*d.kext[-40])))]))
    vmax = max(numpy.nanmax(numpy.log10(m.grid.scattering_phase[0])), \
            numpy.nanmax(numpy.log10(model.grid.scatt[0]*\
            (model.grid.density[0]*d.kext[-40]+\
            model.grid.density[1]*d.kext[-40]))))

    diff = numpy.log10(m.grid.scattering_phase[0][:,:,0]) - \
            numpy.log10(model.grid.scatt[0][1:,:,0,0]*\
            (model.grid.density[0][1:,:,0]*d.kext[-40] + \
            model.grid.density[1][1:,:,0]*d.kext[-40]))

    im1 = ax[0].imshow(numpy.log10(m.grid.scattering_phase[0][:,:,0]), \
            origin="lower", interpolation="nearest", vmin=vmin, vmax=vmax)

    im2 = ax[1].imshow(numpy.log10(model.grid.scatt[0][1:,:,0,0]*\
            (model.grid.density[0][1:,:,0]*d.kext[-40] + \
            model.grid.density[1][1:,:,0]*d.kext[-40])), origin="lower",\
            interpolation="nearest", vmin=vmin, vmax=vmax)

    im3 = ax[2].imshow(diff, origin="lower", interpolation="nearest", \
            vmin=-0.5, vmax=0.5)

ax[0].set_title("RADMC-3D")
ax[1].set_title("MCRT3D")
ax[2].set_title("RADMC-3D - MCRT3D")

fig.colorbar(im1, ax=ax[1], fraction=0.046)
fig.colorbar(im3, ax=ax[2], fraction=0.046)

plt.show()

# Plot the images.

fig, ax = plt.subplots(nrows=1, ncols=3, sharex=True, sharey=True, \
        figsize=(10,3), gridspec_kw=dict(left=0.05, right=0.95, wspace=0.25))

with numpy.errstate(divide="ignore", invalid="ignore"):
    vmax = min(numpy.log10(numpy.nanmax(m.images["image"].image[:,:,0,0])), \
            numpy.log10(numpy.nanmax(model.images["image"].intensity[:,:,0])))
    vmin = vmax - 10.

    diff = (numpy.log10(m.images["image"].image[:,:,0,0]) - \
            numpy.log10(model.images["image"].intensity[:,:,0]))

    im1 = ax[0].imshow(numpy.log10(m.images["image"].image[:,:,0,0]), \
            origin="lower", interpolation="none", vmin=vmin, vmax=vmax)

    im2 = ax[1].imshow(numpy.log10(model.images["image"].intensity[:,:,0]), \
            origin="lower", interpolation="none", vmin=vmin, vmax=vmax)

    im3 = ax[2].imshow(diff, origin="lower", interpolation="none")

ax[0].set_title("RADMC-3D")
ax[1].set_title("MCRT3D")
ax[2].set_title("RADMC-3D - MCRT3D")

fig.colorbar(im3, ax=ax[2], fraction=0.046)

plt.show()

# Plot the unstructured image.

triang = tri.Triangulation(model.images["uimage"].x/AU, \
        model.images["uimage"].y/AU)

plt.tripcolor(triang, model.images["uimage"].intensity[:,0], "ko-")
plt.triplot(triang, "k.-", linewidth=0.1, markersize=0.1)

plt.axes().set_aspect("equal")

plt.xlabel("x", fontsize=14)
plt.ylabel("y", fontsize=14)

plt.axes().tick_params(labelsize=14)

plt.show()

# Plot the channel maps.

fig, ax = plt.subplots(nrows=3, ncols=10, sharex=True, sharey=True, \
        figsize=(15,5), gridspec_kw=dict(hspace=0, wspace=0, left=0.05, \
        right=0.98, top=0.98, bottom=0.05))

with numpy.errstate(divide="ignore", invalid="ignore"):
    vmax = min(numpy.log10(numpy.nanmax(m.images["CO2-1"].image[:,:,:,0])), \
            numpy.log10(numpy.nanmax(model.images["CO2-1"].intensity[:,:,:])))
    vmin = vmax - 0.75

    for i in range(ax[0,:].size):
        diff = (numpy.log10(m.images["CO2-1"].image[:,:,i,0]) - \
                numpy.log10(model.images["CO2-1"].intensity[:,:,i]))

        im1 = ax[0,i].imshow(numpy.log10(m.images["CO2-1"].image[:,:,i,0]), \
                origin="lower", interpolation="none", vmin=vmin, vmax=vmax)

        im2 = ax[1,i].imshow(numpy.log10(model.images["CO2-1"].\
                intensity[:,:,i]), origin="lower", interpolation="none", \
                vmin=vmin, vmax=vmax)

        im3 = ax[2,i].imshow(diff, origin="lower", interpolation="none")

ax[0,0].set_ylabel("RADMC-3D")
ax[1,0].set_ylabel("MCRT3D")
ax[2,0].set_ylabel("RADMC-3D - MCRT3D")

plt.show()

# Plot the spectra.

fig, ax = plt.subplots(nrows=1, ncols=1)

ax.loglog(m.spectra["SED"].wave, m.spectra["SED"].flux)
ax.loglog(model.spectra["SED"].lam, model.spectra["SED"].intensity)

ax.set_ylim(1.0e-6,1.0e1)

plt.show()

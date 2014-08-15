#!/usr/bin/env python

from mcrt3d import *
from numpy import *

# Set up the grid.

def SetupParams():
    Q = Params()

    Q.set_nphot(100000)
    Q.set_bw(True)
    Q.set_scattering(False)
    Q.set_verbose(False)

    return Q

def SetupGrid():
    G = Grid()

    nr = 400
    nt = 200
    np = 2

    r = hstack([0.0,logspace(-1,log10(200.),nr-1)])*au
    t = linspace(0,pi,nt)
    p = linspace(0,2*pi,np)

    G.set_spherical_grid(r,t,p)

    dens = protoplanetary_disk(r, t, p, mass=1.0e-5) + 1.0e-21
    dust = zeros((nr-1,nt-1,np-1), dtype=int)

    G.set_physical_properties(dens,dust)

    # Set up the dust.

    dust_species = array([Dust()])

    dust_species[0].set_properties_from_file("dustkappa_yso.inp")

    G.set_dust_species(dust_species)

    # Set up the source.

    sources = array([Source()])

    sources[0].set_parameters(0.0,0.0,0.0,Msun,Rsun,4000.0)
    sources[0].set_blackbody_spectrum(dust_species[0].nu)

    G.set_sources(sources)

    return G

def protoplanetary_disk(r, theta, phi, rin=0.1, rout=200, mass=1.0e-3, \
        plrho=2.37, h=0.1, plh=58.0/45.0):
    
    ##### Disk Parameters
    
    rin *= au
    rout *= au
    mass *= Msun
    h *= au
    
    ##### Set up the coordinates
    
    rt, tt, pp = meshgrid(0.5*(r[0:r.size-1]+r[1:r.size]), \
            0.5*(theta[0:theta.size-1]+theta[1:theta.size]), \
            0.5*(phi[0:phi.size-1]+phi[1:phi.size]), indexing='ij')
    
    rr = rt*numpy.sin(tt)
    zz = rt*numpy.cos(tt)
    
    ##### Make the dust density model for a protoplanetary disk.
    
    rho0 = mass/((2*pi)**1.5*h*(1*au)**(plrho-plh))*\
        (-plrho+plh+2)/(rout**(-plrho+plh+2)-rin**(-plrho+plh+2))
    hr = h * (rr / (1*au))**(plh)
    rho0 = rho0 * (rr / (1*au))**(-plrho)
    rho1 = exp(-0.5*(zz / hr)**2)
    rho = rho0 * rho1
    rho[(rr >= rout) ^ (rr <= rin)] = 0e0
    
    return rho

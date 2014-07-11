#!/usr/bin/env python

from mcrt3d import *
from numpy import array, arange, pi, zeros

# Set up the grid.

def SetupGrid():
    G = Grid()

    nr = 10
    nt = 10
    np = 10

    r = arange(nr)*au/2
    t = arange(nt)/(nt-1.)*pi
    p = arange(np)/(np-1.)*2*pi

    G.set_spherical_grid(r,t,p)

    dens = zeros((nr-1,nt-1,np-1)) + 1.0e-17
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

def SetupImages():
    nx = 256
    ny = 256
    pixel_size = au/10

    x = (arange(nx,dtype=float)-float(nx)/2)*pixel_size
    y = (arange(ny,dtype=float)-float(ny)/2)*pixel_size

    nu = array([c_l / (1300. * 1.0e-4)])
    nnu = 1

    r = (3*4.5**2)**(1./2)*au
    incl = 0
    pa = 0

    image = Image(r, incl, pa, x, y, nx, ny, nu, pixel_size, nnu)

    return array([image])

def SetupSpectra():
    nx = 1
    ny = 1
    pixel_size = 10*au

    x = array([0.0])
    y = array([0.0])

    #nu = c_l / (logspace(-1,4,1000) * 1.0e-4)
    #nnu = 1000
    nu = c_l / array([0.1300])
    nnu = 1

    r = (3*4.5**2)**(1./2)*au
    incl = 0
    pa = 0

    spectrum = Image(r, incl, pa, x, y, nx, ny, nu, pixel_size, nnu)

    return array([spectrum])

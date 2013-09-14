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

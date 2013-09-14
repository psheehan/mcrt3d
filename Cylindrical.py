#!/usr/bin/env python

from mcrt3d import *
from numpy import array, arange, pi, zeros

# Set up the grid.

def SetupGrid():
    G = Grid()

    nr = 10
    np = 10
    nz = 10

    r = arange(nr)*au/2
    p = arange(np)/(np-1.)*2*pi
    z = (arange(nz)-(nz-1)/2.)*au/1

    G.set_cylindrical_grid(r,p,z)

    dens = zeros((nr-1,np-1,nz-1)) + 1.0e-17
    dust = zeros((nr-1,np-1,nz-1), dtype=int)

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

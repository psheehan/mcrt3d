#!/usr/bin/env python

from mcrt3d import Params
from mcrt3d.grid import CylindricalGrid
from mcrt3d.dust import Dust
from mcrt3d.sources import Star
from mcrt3d.camera import Image, Spectrum
from mcrt3d.constants.astronomy import M_sun, R_sun, AU
from mcrt3d.constants.physics import c
from numpy import array, arange, pi, zeros, logspace

# Set up the grid.

def SetupParams(model):
    model.params.nphot = 100000
    model.params.bw = True
    model.params.scattering = False
    model.params.verbose = False
    model.params.use_mrw = True
    model.params.mrw_gamma = 2

def SetupGrid(model):
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

def SetupImages():
    nx = 256
    ny = 256
    pixel_size = AU/10

    x = (arange(nx,dtype=float)-float(nx)/2)*pixel_size
    y = (arange(ny,dtype=float)-float(ny)/2)*pixel_size

    nu = array([c / (1300. * 1.0e-4)])
    nnu = 1

    r = (3*4.5**2)**(1./2)*AU
    incl = pi/4
    pa = pi/4

    image = Image(r, incl, pa, x, y, nx, ny, nu, pixel_size, nnu)

    return array([image])

def SetupSpectra():
    pixel_size = 10*AU/100

    nnu = 200
    nu = c / (logspace(-1,4,nnu) * 1.0e-4)

    r = (3*4.5**2)**(1./2)*AU
    incl = 0
    pa = 0

    spectrum = Spectrum(r, incl, pa, nu, pixel_size, nnu)

    return array([spectrum])

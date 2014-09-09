#!/usr/bin/env python

from mcrt3d import Params
from mcrt3d.grid import CartesianGrid
from mcrt3d.dust import Dust
from mcrt3d.sources import Star
from mcrt3d.camera import Image
from mcrt3d.constants.astronomy import M_sun, R_sun, AU
from mcrt3d.constants.physics import c
from numpy import array, arange, pi, zeros, logspace

# Set up the grid.

def SetupParams():
    Q = Params()

    Q.set_nphot(100000)
    Q.set_bw(True)
    Q.set_scattering(False)
    Q.set_verbose(False)
    Q.set_mrw(False)
    Q.set_mrw_gamma(4)

    return Q

def SetupGrid():
    G = SphericalGrid()

    # Set up the dust.

    dust = Dust()
    dust.set_properties_from_radmc3d("dustkappa_yso.inp")

    # Set up the star.

    star = Star()
    star.set_parameters(0.0,0.0,0.0,M_sun,R_sun,4000.0)
    star.set_blackbody_spectrum(dust.nu)

    G.add_source(star)

    # Set up the grid.

    nr = 10
    nt = 10
    np = 10

    r = arange(nr)*AU/2
    t = arange(nt)/(nt-1.)*pi
    p = arange(np)/(np-1.)*2*pi

    G.set_walls(r,t,p)

    density = zeros((nr-1,nt-1,np-1)) + 1.0e-17

    G.add_density(density, dust)

    return G

def SetupImages():
    nx = 256
    ny = 256
    pixel_size = AU/10

    x = (arange(nx,dtype=float)-float(nx)/2)*pixel_size
    y = (arange(ny,dtype=float)-float(ny)/2)*pixel_size

    nu = array([c / (1300. * 1.0e-4)])
    nnu = 1

    r = (3*4.5**2)**(1./2)*AU
    incl = 0
    pa = 0

    image = Image(r, incl, pa, x, y, nx, ny, nu, pixel_size, nnu)

    return array([image])

def SetupSpectra():
    nx = 1
    ny = 1
    pixel_size = 10*AU

    x = array([0.0])
    y = array([0.0])

    nu = c / (logspace(-1,4,1000) * 1.0e-4)
    nnu = 1000

    r = (3*4.5**2)**(1./2)*AU
    incl = pi/4
    pa = pi/4

    spectrum = Image(r, incl, pa, x, y, nx, ny, nu, pixel_size, nnu)

    return array([spectrum])

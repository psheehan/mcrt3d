#!/usr/bin/env python3

from mcrt3d import *
from mcrt3d.constants.astronomy import M_sun, R_sun, AU
from mcrt3d.constants.physics import c
from numpy import array, arange, pi, zeros, logspace

dust = Dust()
dust.set_properties_from_radmc3d("dustkappa_yso.inp")

dust.make_lookup_tables2()

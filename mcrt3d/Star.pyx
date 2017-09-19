import cython

import numpy
cimport numpy

from ..constants.physics import sigma
from .. import misc
import scipy.integrate
import numpy
import h5py

from ..mcrt3d cimport Source, Star

cdef class StarObj:
    cdef Star *obj

    cdef double x, y, z
    cdef numpy.ndarray nu, Bnu, random_nu_CPD

    property nnu:
        def __get__(self):
            return self.obj.nnu
        def __set__(self, int var):
            self.obj.nnu = var

    property mass:
        def __get__(self):
            return self.obj.mass
        def __set__(self, double var):
            self.obj.mass = var

    property radius:
        def __get__(self):
            return self.obj.radius
        def __set__(self, double var):
            self.obj.radius = var

    property temperature:
        def __get__(self):
            return self.obj.temperature
        def __set__(self, double var):
            self.obj.temperature = var

    property luminosity:
        def __get__(self):
            return self.obj.luminosity
        def __set__(self, double var):
            self.obj.luminosity = var

    def __init__(self, double x=None, double y=None, double z=None, \
            double mass=None, double radius=None, double temperature=None, \
            filename=None):

        if filename != None:
            self.read(filename)
        else:
            self.set_parameters(x, y, z, mass, radius, temperature)

        self.obj = new Star(x, y, z, mass, radius, temperature)

    def set_parameters(self, double x, double y, double z, double mass, \
            double radius, double temperature):
        self.x = x
        self.y = y
        self.z = z
        self.mass = mass
        self.radius = radius
        self.temperature = temperature
        self.luminosity = 4*numpy.pi*radius**2 * sigma*temperature**4

    def set_blackbody_spectrum(self,numpy.ndarray[double, ndim=1, mode="c"] nu):
        self.nnu = nu.size
        self.nu = nu
        self.Bnu = misc.B_nu(nu, self.temperature)
        self.luminosity = 4*numpy.pi*self.radius**2*sigma*self.temperature**4

        self.random_nu_CPD = scipy.integrate.cumtrapz(self.Bnu, x=nu, \
                initial=0) / numpy.trapz(self.Bnu, x=nu)

        cdef numpy.ndarray[double, ndim=1, mode="c"] Bnu = self.Bnu, \
                random_nu_CPD = self.random_nu_CPD

        self.obj.set_blackbody_spectrum(self.nnu, &nu[0], &Bnu[0], \
                self.luminosity, &random_nu_CPD[0])

    def read(self, filename=None, usefile=None):
        if (usefile == None):
            f = h5py.File(filename, "r")
        else:
            f = usefile

        self.mass = f['mass'].value
        self.luminosity = f['luminosity'].value
        self.temperature = f['temperature'].value
        self.radius = f['radius'].value
        self.x = f['x'].value
        self.y = f['y'].value
        self.z = f['z'].value

        if (usefile == None):
            f.close()

    def write(self, filename=None, usefile=None):
        if (usefile == None):
            f = h5py.File(filename, "w")
        else:
            f = usefile

        f['mass'] = self.mass
        f['luminosity'] = self.luminosity
        f['temperature'] = self.temperature
        f['radius'] = self.radius

        f['x'] = self.x
        f['y'] = self.y
        f['z'] = self.z

        if (usefile == None):
            f.close()

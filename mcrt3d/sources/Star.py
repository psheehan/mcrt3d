import ctypes
import numpy
import numpy.ctypeslib as npctypes
import h5py
import os
from .. import misc
from ..constants.physics import sigma

array_1d_double = npctypes.ndpointer(dtype=ctypes.c_double, ndim=1,
                flags='CONTIGUOUS')

lib = ctypes.cdll.LoadLibrary(os.path.dirname(__file__)+ \
        '/../../src/libmcrt3d.so')

lib.new_Source.restype = ctypes.c_void_p
lib.new_Source.argtypes = None

lib.set_parameters.argtypes = [ctypes.c_void_p, ctypes.c_double, \
        ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, \
        ctypes.c_double]

lib.set_blackbody_spectrum.restype = None
lib.set_blackbody_spectrum.argtypes = [ctypes.c_void_p, ctypes.c_int, \
        array_1d_double, array_1d_double, ctypes.c_double]

class Star:
    def __init__(self):
        self.obj = lib.new_Source()

    def set_parameters(self, x, y, z, mass, radius, temperature):
        self.x = x
        self.y = y
        self.z = z
        self.mass = mass
        self.radius = radius
        self.temperature = temperature
        self.luminosity = 4*numpy.pi*radius**2 * sigma*temperature**4

        lib.set_parameters(self.obj, x, y, z, mass, radius, temperature)

    def set_blackbody_spectrum(self, nu):
        self.nu = nu
        self.Bnu = misc.B_nu(nu, self.temperature)
        self.luminosity = 4*numpy.pi*self.radius**2*sigma*self.temperature**4

        lib.set_blackbody_spectrum(self.obj, self.nu.size, self.nu, \
                self.Bnu, self.luminosity)

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

import os
import ctypes
import numpy
import numpy.ctypeslib as npctypes

array_1d_double = npctypes.ndpointer(dtype=ctypes.c_double, ndim=1,
                flags='CONTIGUOUS')
array_3d_double = npctypes.ndpointer(dtype=ctypes.c_double, ndim=3,
                flags='CONTIGUOUS')

lib = ctypes.cdll.LoadLibrary(os.path.dirname(__file__)+'/../../src/libmcrt3d.so')

lib.new_Image.restype = ctypes.c_void_p
lib.new_Image.argtypes = [ctypes.c_double, ctypes.c_double, ctypes.c_double, \
        array_1d_double, array_1d_double, array_3d_double, ctypes.c_int, \
        ctypes.c_int, array_1d_double, ctypes.c_double, ctypes.c_int]

class Image:
    def __init__(self, r, incl, pa, x, y, nx, ny, nu, pixel_size, nnu):
        self.r = r
        self.incl = incl
        self.pa = pa
        self.x = x
        self.y = y
        self.intensity = numpy.zeros((nx, ny, nnu),dtype=float)
        self.nx = nx
        self.ny = ny
        self.nu = nu
        self.pixel_size = pixel_size
        self.nnu = nnu

        self.obj = lib.new_Image(self.r, self.incl, self.pa, self.x, \
                self.y, self.intensity, self.nx, self.ny, self.nu, \
                self.pixel_size, self.nnu)

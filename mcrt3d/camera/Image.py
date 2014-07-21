import os
import ctypes
import numpy

lib = ctypes.cdll.LoadLibrary(os.path.dirname(__file__)+'/../../src/libmcrt3d.so')
lib.new_Image.restype = ctypes.c_void_p

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

        self.obj = lib.new_Image(ctypes.c_double(self.r), \
                ctypes.c_double(self.incl), ctypes.c_double(self.pa), \
                self.x.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), \
                self.y.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), \
                self.intensity.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),\
                ctypes.c_int(self.nx), ctypes.c_int(self.ny), \
                self.nu.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), \
                ctypes.c_double(self.pixel_size), ctypes.c_int(self.nnu))

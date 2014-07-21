import os
import ctypes
from .dust import Dust
from .sources import Star as Source
from .grid import Grid
from .camera import Camera, Image

lib = ctypes.cdll.LoadLibrary(os.path.dirname(__file__)+'/../src/libmcrt3d.so')
lib.new_mcrt.restype = ctypes.c_void_p

class MCRT:
    def __init__(self, G):
        self.obj = lib.new_mcrt()

        lib.set_Grid(ctypes.c_void_p(self.obj),ctypes.c_void_p(G.obj))

    def thermal_mc(self, nphot, bw):
        lib.run_thermal_mc(ctypes.c_void_p(self.obj), ctypes.c_int(nphot), \
                ctypes.c_bool(bw))

    def thermal_mc_omp(self, nphot, bw, nthreads):
        lib.run_thermal_mc_omp(ctypes.c_void_p(self.obj), ctypes.c_int(nphot), \
                ctypes.c_bool(bw), ctypes.c_int(nthreads))

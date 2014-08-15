import os
import ctypes
from .dust import Dust
from .sources import Star as Source
from .grid import Grid
from .camera import Camera, Image

lib = ctypes.cdll.LoadLibrary(os.path.dirname(__file__)+'/../src/libmcrt3d.so')
lib.new_mcrt.restype = ctypes.c_void_p
lib.new_Params.restype = ctypes.c_void_p

class Params:
    def __init__(self):
        self.obj = lib.new_Params()

    def set_nphot(self, nphot):
        lib.set_nphot(ctypes.c_void_p(self.obj), ctypes.c_int(nphot))

    def set_bw(self, bw):
        lib.set_bw(ctypes.c_void_p(self.obj), ctypes.c_bool(bw))

    def set_scattering(self, scattering):
        lib.set_scattering(ctypes.c_void_p(self.obj), ctypes.c_bool(scattering))

    def set_verbose(self, verbose):
        lib.set_verbose(ctypes.c_void_p(self.obj), ctypes.c_bool(verbose))

class MCRT:
    def __init__(self, G, P, Q):
        self.obj = lib.new_mcrt()

        lib.set_Grid(ctypes.c_void_p(self.obj),ctypes.c_void_p(G.obj))
        lib.set_Params(ctypes.c_void_p(self.obj),ctypes.c_void_p(Q.obj))

    def thermal_mc(self, nphot, bw):
        lib.run_thermal_mc(ctypes.c_void_p(self.obj), ctypes.c_int(nphot), \
                ctypes.c_bool(bw))

    def thermal_mc_omp(self, nphot, bw, nthreads):
        lib.run_thermal_mc_omp(ctypes.c_void_p(self.obj), ctypes.c_int(nphot), \
                ctypes.c_bool(bw), ctypes.c_int(nthreads))

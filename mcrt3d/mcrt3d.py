import os
import ctypes
from .dust import Dust
from .sources import Star as Source
from .grid import Grid
from .camera import Camera, Image

lib = ctypes.cdll.LoadLibrary(os.path.dirname(__file__)+'/../src/libmcrt3d.so')

lib.new_mcrt.restype = ctypes.c_void_p
lib.new_mcrt.argtypes = None

lib.set_Grid.restype = None
lib.set_Grid.argtypes = [ctypes.c_void_p, ctypes.c_void_p]

lib.set_Params.restype = None
lib.set_Params.argtypes = [ctypes.c_void_p, ctypes.c_void_p]

lib.run_thermal_mc.restype = None
lib.run_thermal_mc.argtypes = [ctypes.c_void_p]

lib.new_Params.restype = ctypes.c_void_p
lib.new_Params.argtypes = None

lib.set_nphot.restype = None
lib.set_nphot.argtypes = [ctypes.c_void_p, ctypes.c_int]

lib.set_bw.restype = None
lib.set_bw.argtypes = [ctypes.c_void_p, ctypes.c_bool]

lib.set_scattering.restype = None
lib.set_scattering.argtypes = [ctypes.c_void_p, ctypes.c_bool]

lib.set_verbose.restype = None
lib.set_verbose.argtypes = [ctypes.c_void_p, ctypes.c_bool]

class Params:
    def __init__(self):
        self.obj = lib.new_Params()

    def set_nphot(self, nphot):
        lib.set_nphot(self.obj, nphot)

    def set_bw(self, bw):
        lib.set_bw(self.obj, bw)

    def set_scattering(self, scattering):
        lib.set_scattering(self.obj, scattering)

    def set_verbose(self, verbose):
        lib.set_verbose(self.obj, verbose)

class MCRT:
    def __init__(self, G, Q):
        self.obj = lib.new_mcrt()

        lib.set_Grid(self.obj, G.obj)
        lib.set_Params(self.obj, Q.obj)

    def thermal_mc(self, nphot, bw):
        lib.run_thermal_mc(self.obj)

import cython

from mcrt3d cimport *

from .grid.Grid cimport GridObj

# Define the Params and MCRT classes here because there isn't yet a better 
# place.

cdef class ParamsObj:
    cdef Params *obj

    def __init__(self):
        self.obj = new Params()
        if self.obj == NULL:
            raise MemoryError("Not enough memory.")

    def __del__(self):
        del self.obj

    property nphot:
        def __get__(self):
            return self.obj.nphot
        def __set__(self, int var):
            self.obj.set_nphot(var)

    property bw:
        def __get__(self):
            return self.obj.bw
        def __set__(self, bool var):
            self.obj.set_bw(var)

    property scattering:
        def __get__(self):
            return self.obj.scattering
        def __set__(self, bool var):
            self.obj.set_scattering(var)

    property verbose:
        def __get__(self):
            return self.obj.verbose
        def __set__(self, bool var):
            self.obj.set_verbose(var)

    property use_mrw:
        def __get__(self):
            return self.obj.use_mrw
        def __set__(self, bool var):
            self.obj.set_mrw(var)

    property mrw_gamma:
        def __get__(self):
            return self.obj.mrw_gamma
        def __set__(self, double var):
            self.obj.set_mrw_gamma(var)

cdef class MCRTObj:
    cdef MCRT *obj
    cdef GridObj grid
    cdef ParamsObj params

    def __init__(self, GridObj G, ParamsObj Q):
        self.grid = G
        self.params = Q

        self.obj = new MCRT(G.obj, Q.obj)

    def run_thermal_mc(self):
        self.obj.thermal_mc()

    def run_scattering_mc(self):
        self.obj.scattering_mc()

    def run_mc_iteration(self):
        self.obj.mc_iteration()

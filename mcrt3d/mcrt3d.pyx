import cython

import numpy
cimport numpy

from mcrt3d cimport *

from .grid.Grid cimport GridObj
from .grid.CartesianGrid cimport CartesianGridObj
from .grid.CylindricalGrid cimport CylindricalGridObj
from .grid.SphericalGrid cimport SphericalGridObj
from .camera.Image cimport ImageObj
from .camera.Spectrum cimport SpectrumObj

# Define the Params and MCRT classes here because there isn't yet a better 
# place.

cdef class ParamsObj:
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

    def __init__(self):
        self.params = ParamsObj()

    def set_cartesian_grid(self, numpy.ndarray[double, ndim=1, mode="c"] x, \
            numpy.ndarray[double, ndim=1, mode="c"] y, \
            numpy.ndarray[double, ndim=1, mode="c"] z):
        self.grid = CartesianGridObj(x, y, z)

        self.obj = new MCRT(self.grid.obj, self.params.obj)

    def set_cylindrical_grid(self, numpy.ndarray[double, ndim=1, mode="c"] r, \
            numpy.ndarray[double, ndim=1, mode="c"] p, \
            numpy.ndarray[double, ndim=1, mode="c"] z):
        self.grid = CylindricalGridObj(r, p ,z)

        self.obj = new MCRT(self.grid.obj, self.params.obj)

    def set_spherical_grid(self, numpy.ndarray[double, ndim=1, mode="c"] r, \
            numpy.ndarray[double, ndim=1, mode="c"] t, \
            numpy.ndarray[double, ndim=1, mode="c"] p):
        self.grid = SphericalGridObj(r, t, p)

        self.obj = new MCRT(self.grid.obj, self.params.obj)

    def run_thermal_mc(self):
        self.obj.thermal_mc()

    def run_image(self, ImageObj I):
        self.obj.run_image(I.obj)

    def run_spectrum(self, SpectrumObj S):
        self.obj.run_spectrum(S.obj)

    property params:
        def __get__(self):
            return self.params

    property grid:
        def __get__(self):
            return self.grid

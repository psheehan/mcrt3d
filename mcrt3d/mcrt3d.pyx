import cython

import numpy
cimport numpy

from mcrt3d cimport *

from .grid.Grid cimport GridObj
from .grid.CartesianGrid cimport CartesianGridObj
from .grid.CylindricalGrid cimport CylindricalGridObj
from .grid.SphericalGrid cimport SphericalGridObj
from .camera.Camera cimport CameraObj
from .camera.Image cimport ImageObj
from .camera.Spectrum cimport SpectrumObj

from .constants.astronomy import Jy, pc

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
    cdef CameraObj camera

    def __init__(self):
        self.params = ParamsObj()
        self.camera = CameraObj()

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

    def run_image(self, double incl=0., double pa=0., double dpc=1.):
        self.camera.x = (numpy.arange(self.camera.nx,dtype=float)-\
                float(self.camera.nx)/2)*self.camera.pixel_size
        self.camera.y = (numpy.arange(self.camera.ny,dtype=float)-\
                float(self.camera.ny)/2)*self.camera.pixel_size

        cdef double r = (4*max(self.grid.w1[-1],self.grid.w2[-1], \
                self.grid.w3[-1])**2)**0.5

        cdef ImageObj image = ImageObj(r, incl, pa, self.camera.x, \
                self.camera.y, self.camera.nx, self.camera.ny, \
                self.camera.lam, self.camera.pixel_size, self.camera.nlam)

        self.obj.run_image(image.obj)

        image.intensity *= (self.camera.pixel_size / (dpc*pc))**2 / Jy

        return image

    def run_spectrum(self, double incl=0., double pa=0., double dpc=1.):
        cdef double r = (4*max(self.grid.w1[-1],self.grid.w2[-1], \
                self.grid.w3[-1])**2)**0.5

        cdef SpectrumObj spectrum = SpectrumObj(r, incl, pa, self.camera.lam, \
                self.camera.pixel_size, self.camera.nlam)

        self.obj.run_spectrum(spectrum.obj)

        spectrum.intensity *= (self.camera.pixel_size / (dpc*pc))**2 / Jy

        return spectrum

    property params:
        def __get__(self):
            return self.params

    property camera:
        def __get__(self):
            return self.camera

    property grid:
        def __get__(self):
            return self.grid

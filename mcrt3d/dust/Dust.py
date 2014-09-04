import ctypes
import numpy
import h5py
import os
from ..constants.physics import c, sigma, k
from .. import misc

lib = ctypes.cdll.LoadLibrary(os.path.dirname(__file__)+ \
        '/../../src/libmcrt3d.so')
lib.new_Dust.restype = ctypes.c_void_p

class Dust:
    def __init__(self):
        self.obj = lib.new_Dust()

    def set_properties(self, lam, kabs, ksca):
        self.lam = lam
        self.nu = c / self.lam

        self.kabs = kabs
        self.ksca = ksca
        self.kext = self.kabs + self.ksca
        self.albedo = self.ksca / self.kext

        lib.set_optical_properties(ctypes.c_void_p(self.obj), self.lam.size, \
            self.nu.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), \
            self.lam.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), \
            self.kabs.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), \
            self.ksca.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), \
            self.kext.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), \
            self.albedo.ctypes.data_as(ctypes.POINTER(ctypes.c_double)))

        self.make_lookup_tables()

    def set_properties_from_file(self, filename=None, usefile=None):
        if (usefile == None):
            f = h5py.File(filename, "r")
        else:
            f = usefile

        self.lam = f['lam'].value
        self.nu = c / self.lam

        if ('kabs' in f):
            self.kabs = f['kabs'].value
        if ('ksca' in f):
            self.ksca = f['ksca'].value
        if (hasattr(self, 'kabs') and hasattr(self, 'ksca')):
            self.kext = self.kabs + self.ksca
            self.albedo = self.ksca / self.kext

        if (usefile == None):
            f.close()

        lib.set_optical_properties(ctypes.c_void_p(self.obj), self.lam.size, \
            self.nu.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), \
            self.lam.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), \
            self.kabs.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), \
            self.ksca.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), \
            self.kext.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), \
            self.albedo.ctypes.data_as(ctypes.POINTER(ctypes.c_double)))

        self.make_lookup_tables()

    def set_properties_from_radmc3d(self, filename):
        data = numpy.loadtxt(filename, skiprows=2)

        self.lam = data[:,0].copy() * 1.0e-4
        self.nu = c / self.lam
        self.kabs = data[:,1].copy()
        self.ksca = data[:,2].copy()
        self.kext = self.kabs + self.ksca
        self.albedo = self.ksca / self.kext

        lib.set_optical_properties(ctypes.c_void_p(self.obj), self.lam.size, \
            self.nu.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), \
            self.lam.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), \
            self.kabs.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), \
            self.ksca.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), \
            self.kext.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), \
            self.albedo.ctypes.data_as(ctypes.POINTER(ctypes.c_double)))

        self.make_lookup_tables()

    def make_lookup_tables(self):
        self.temp = numpy.logspace(-1,5,100).astype(float)
        self.ntemp = self.temp.size

        self.int_Bnukext = numpy.ones(self.temp.size)
        self.int_dBnukext = numpy.ones(self.temp.size)
        for i in range(self.temp.size-1):
            Bnu = misc.B_nu(self.nu,self.temp[i])
            dBnu = misc.dB_nu(self.nu,self.temp[i])

            self.int_Bnukext[i] = numpy.trapz(Bnu*self.kext,x=self.nu)
            self.int_dBnukext[i] = numpy.trapz(dBnu*self.kext,x=self.nu)
        self.planck_opacity = -self.int_Bnukext/(sigma*self.temp**4/numpy.pi)

        self.dplanck_opacity_dT = numpy.zeros(self.temp.size)
        self.dint_dBnukext_dT = numpy.zeros(self.temp.size)
        for i in range(self.temp.size-1):
            self.dplanck_opacity_dT[i] = (self.planck_opacity[i+1]- \
                    self.planck_opacity[i])/(self.temp[i+1]-self.temp[i])
            self.dint_dBnukext_dT[i] = (self.int_dBnukext[i+1]- \
                    self.int_dBnukext[i])/(self.temp[i+1]-self.temp[i])

        self.dkextdnu = numpy.zeros(self.nu.size)
        self.dalbedodnu = numpy.zeros(self.nu.size)
        for i in range(self.nu.size-1):
            self.dkextdnu[i] = (self.kext[i+1]-self.kext[i])/(self.nu[i+1]- \
                    self.nu[i])
            self.dalbedodnu[i] = (self.albedo[i+1]-self.albedo[i])/ \
                    (self.nu[i+1]-self.nu[i])

        self.Bnu = numpy.zeros((self.temp.size,self.nu.size))
        self.dBnu = numpy.zeros((self.temp.size,self.nu.size))
        self.dBnudT = numpy.zeros((self.temp.size,self.nu.size))
        self.ddBnudT = numpy.zeros((self.temp.size,self.nu.size))
        for i in range(self.temp.size):
            self.Bnu[i,:] = misc.B_nu(self.nu,self.temp[i])
            self.dBnu[i,:] = misc.dB_nu(self.nu,self.temp[i])

            if (i < self.temp.size-1):
                self.dBnudT[i,:] = (misc.B_nu(self.nu,self.temp[i+1])- \
                        misc.B_nu(self.nu,self.temp[i]))/ \
                        (self.temp[i+1]-self.temp[i])
                self.ddBnudT[i,:] = (misc.dB_nu(self.nu,self.temp[i+1])- \
                        misc.dB_nu(self.nu,self.temp[i]))/ \
                        (self.temp[i+1]-self.temp[i])

        lib.set_lookup_tables(ctypes.c_void_p(self.obj), self.ntemp, \
                self.temp.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), \
                self.planck_opacity.ctypes.data_as(ctypes.POINTER( \
                        ctypes.c_double)),\
                self.int_dBnukext.ctypes.data_as(ctypes.POINTER( \
                        ctypes.c_double)), \
                self.dplanck_opacity_dT.ctypes.data_as( \
                        ctypes.POINTER(ctypes.c_double)), \
                self.dint_dBnukext_dT.ctypes.data_as( \
                        ctypes.POINTER(ctypes.c_double)), \
                self.dkextdnu.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), \
                self.dalbedodnu.ctypes.data_as(ctypes.POINTER( \
                        ctypes.c_double)), \
                self.Bnu.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), \
                self.dBnu.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), \
                self.dBnudT.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), \
                self.ddBnudT.ctypes.data_as(ctypes.POINTER(ctypes.c_double)))

    def write(self, filename=None, usefile=None):
        if (usefile == None):
            f = h5py.File(filename, "w")
        else:
            f = usefile

        if hasattr(self, 'lam'):
            lam_dset = f.create_dataset("lam", (self.lam.size,), dtype='f')
            lam_dset[...] = self.lam

        if hasattr(self, 'kabs'):
            kabs_dset = f.create_dataset("kabs", (self.kabs.size,), dtype='f')
            kabs_dset[...] = self.kabs
        if hasattr(self, 'ksca'):
            ksca_dset = f.create_dataset("ksca", (self.ksca.size,), dtype='f')
            ksca_dset[...] = self.ksca
        if hasattr(self, 'g'):
            g_dset = f.create_dataset("g", (self.g.size,), dtype='f')
            g_dset[...] = self.g

        if (usefile == None):
            f.close()

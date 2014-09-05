from ..constants.physics import c, sigma, k
from ..mcrt3d import lib
from .. import misc
import numpy
import h5py

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

        lib.set_optical_properties(self.obj, self.lam.size, self.nu, self.lam, \
                self.kabs, self.ksca, self.kext, self.albedo)

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

        lib.set_optical_properties(self.obj, self.lam.size, self.nu, self.lam, \
                self.kabs, self.ksca, self.kext, self.albedo)

        self.make_lookup_tables()

    def set_properties_from_radmc3d(self, filename):
        data = numpy.loadtxt(filename, skiprows=2)

        self.lam = data[:,0].copy() * 1.0e-4
        self.nu = c / self.lam
        self.kabs = data[:,1].copy()
        self.ksca = data[:,2].copy()
        self.kext = self.kabs + self.ksca
        self.albedo = self.ksca / self.kext

        lib.set_optical_properties(self.obj, self.lam.size, self.nu, self.lam, \
                self.kabs, self.ksca, self.kext, self.albedo)

        self.make_lookup_tables()

    def make_lookup_tables(self):
        self.temp = numpy.logspace(-1,5,100).astype(float)
        self.ntemp = self.temp.size

        nu, temp = numpy.meshgrid(self.nu, self.temp)
        ksca = numpy.vstack([self.ksca for i in range(self.temp.size)])
        kext = numpy.vstack([self.kext for i in range(self.temp.size)])

        self.planck_opacity =  -numpy.pi / (sigma * self.temp**4) * \
                numpy.trapz(misc.B_nu(nu, temp)*kext, x=nu, axis=1)
        self.dplanck_opacity_dT = numpy.diff(self.planck_opacity) / \
                numpy.diff(self.temp)

        self.rosseland_extinction = -(sigma * self.temp**4 / numpy.pi) / \
                numpy.trapz(misc.B_nu(nu, temp)/kext, x=nu, axis=1)
        self.drosseland_extinction_dT = numpy.diff(self.rosseland_extinction) /\
                numpy.diff(self.temp)

        self.int_dBnukext = numpy.trapz(misc.dB_nu(nu, temp) * kext, x=nu, \
                axis=1)
        self.dint_dBnukext_dT = numpy.diff(self.int_dBnukext) / \
                numpy.diff(self.temp)

        self.dkextdnu = numpy.diff(self.kext) / numpy.diff(self.nu)
        self.dalbedodnu = numpy.diff(self.albedo) / numpy.diff(self.nu)

        self.Bnu = misc.B_nu(nu, temp)
        self.dBnu = misc.dB_nu(nu, temp)
        self.dBnudT = numpy.diff(self.Bnu, axis=0) / numpy.diff(temp, axis=0)
        self.ddBnudT = numpy.diff(self.dBnu, axis=0) / numpy.diff(temp, axis=0)

        lib.set_lookup_tables(self.obj, self.ntemp, self.temp, \
                self.planck_opacity, self.rosseland_extinction, \
                self.int_dBnukext, self.dplanck_opacity_dT, \
                self.drosseland_extinction_dT, self.dint_dBnukext_dT, \
                self.dkextdnu, self.dalbedodnu, self.Bnu, self.dBnu, \
                self.dBnudT, self.ddBnudT)

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

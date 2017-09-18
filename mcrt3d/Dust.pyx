import cython

import numpy
cimport numpy

from ..constants.physics import c, sigma, k
from .. import misc
import scipy.integrate
import h5py

from ..mcrt3d cimport Dust, IsotropicDust

cdef class DustObj:
    cdef Dust *obj

    cdef numpy.ndarray lam, nu, kabs, ksca, kext, albedo
    cdef numpy.ndarray temp, planck_opacity, dplanck_opacity_dT, \
            rosseland_extinction, drosseland_extinction_dT, random_nu_CPD, \
            random_nu_CPD_bw, drandom_nu_CPD_dT, drandom_nu_CPD_bw_dT

    def __init__(self, numpy.ndarray[double, ndim=1, mode="c"] lam=None, \
            numpy.ndarray[double, ndim=1, mode="c"] kabs=None, \
            numpy.ndarray[double, ndim=1, mode="c"] ksca=None, \
            filename=None, radmc3d=False):

        if filename != None:
            if radmc3d:
                self.set_properties_from_radmc3d(filename)
            else:
                self.set_properties_from_file(filename)
        else:
            self.set_properties(lam, kabs, ksca)

        lam = self.lam
        cdef numpy.ndarray[double, ndim=1, mode="c"] nu = self.nu
        kabs = self.kabs
        ksca = self.ksca
        cdef numpy.ndarray[double, ndim=1, mode="c"] kext = self.kext
        cdef numpy.ndarray[double, ndim=1, mode="c"] albedo = self.albedo

        self.obj = new Dust(self.obj.nlam, &nu[0], &lam[0], \
                &kabs[0], &ksca[0], &kext[0], &albedo[0])

        self.make_lookup_tables()

    property nlam:
        def __get__(self):
            return self.obj.nlam
        def __set__(self, int var):
            self.obj.nlam = var

    property ntemp:
        def __get__(self):
            return self.obj.ntemp
        def __set__(self, int var):
            self.obj.ntemp = var

    def set_properties(self, numpy.ndarray[double, ndim=1, mode="c"] lam, \
            numpy.ndarray[double, ndim=1, mode="c"] kabs, \
            numpy.ndarray[double, ndim=1, mode="c"] ksca):

        self.nlam = lam.size

        self.lam = lam
        self.nu = c / self.lam

        self.kabs = kabs
        self.ksca = ksca
        self.kext = self.kabs + self.ksca
        self.albedo = self.ksca / self.kext

    def set_properties_from_file(self, filename=None, usefile=None):
        if (usefile == None):
            f = h5py.File(filename, "r")
        else:
            f = usefile

        self.lam = f['lam'].value
        self.nu = c / self.lam

        self.nlam = self.lam.size

        if ('kabs' in f):
            self.kabs = f['kabs'].value
        if ('ksca' in f):
            self.ksca = f['ksca'].value
        if (hasattr(self, 'kabs') and hasattr(self, 'ksca')):
            self.kext = self.kabs + self.ksca
            self.albedo = self.ksca / self.kext

        if (usefile == None):
            f.close()

    def set_properties_from_radmc3d(self, filename):
        data = numpy.loadtxt(filename, skiprows=2)

        self.lam = data[:,0][::-1].copy() * 1.0e-4
        self.nu = c / self.lam
        self.kabs = data[:,1][::-1].copy()
        self.ksca = data[:,2][::-1].copy()
        self.kext = self.kabs + self.ksca
        self.albedo = self.ksca / self.kext

        self.nlam = self.lam.size

    def make_lookup_tables(self):
        self.temp = numpy.logspace(-1,5,100).astype(float)
        self.ntemp = self.temp.size

        nu, temp = numpy.meshgrid(self.nu, self.temp)
        kabs = numpy.vstack([self.kabs for i in range(self.temp.size)])
        kext = numpy.vstack([self.kext for i in range(self.temp.size)])

        # Calculate the Planck Mean Opacity and it's derivative.

        self.planck_opacity =  numpy.pi / (sigma * self.temp**4) * \
                numpy.trapz(misc.B_nu(nu, temp)*kabs, x=nu, axis=1)
        self.dplanck_opacity_dT = numpy.diff(self.planck_opacity) / \
                numpy.diff(self.temp)

        # Calculate the Rosseland Mean Extinction and it's derivative.

        self.rosseland_extinction = (sigma * self.temp**4 / numpy.pi) / \
                numpy.trapz(misc.B_nu(nu, temp)/kext, x=nu, axis=1)
        self.drosseland_extinction_dT = numpy.diff(self.rosseland_extinction) /\
                numpy.diff(self.temp)

        # Calculate the derivatives of kext and albedo.

        self.dkextdnu = numpy.diff(self.kext) / numpy.diff(self.nu)
        self.dalbedodnu = numpy.diff(self.albedo) / numpy.diff(self.nu)

        # Calculate the cumulative probability distribution that will be used
        # to generate a random nu value, both for bw and regular.

        self.random_nu_CPD = scipy.integrate.cumtrapz(kabs * \
                misc.B_nu(nu, temp), x=nu, axis=1, initial=0) / numpy.dstack( \
                [numpy.trapz(kabs * misc.B_nu(nu, temp), x=nu, axis=1) for i \
                in range(self.nu.size)])[0]
        self.drandom_nu_CPD_dT = numpy.diff(self.random_nu_CPD, axis=0) / \
                numpy.diff(temp, axis=0)

        self.random_nu_CPD_bw = scipy.integrate.cumtrapz(kabs * \
                misc.dB_nu(nu, temp), x=nu, axis=1, initial=0) / numpy.dstack( \
                [numpy.trapz(kabs * misc.dB_nu(nu, temp), x=nu, axis=1) for i \
                in range(self.nu.size)])[0]
        self.drandom_nu_CPD_bw_dT = numpy.diff(self.random_nu_CPD_bw, axis=0) /\
                numpy.diff(temp, axis=0)

        # Make sure that Cython knows they are arrays...

        cdef numpy.ndarray[double, ndim=1, mode="c"] temp_1d = self.temp, \
                rosseland_extinction = self.rosseland_extinction, \
                planck_opacity = self.plank_opacity, \
                dplanck_opacity_dT = self.dplank_opacity_dT, \
                drosseland_extinction_dT = self.drosseland_extinction_dT, \
                dkextdnu = self.dkextdnu, dalbedodnu = self.dalbedodnu
        cdef numpy.ndarray[double, ndim=2, mode="c"] \
                random_nu_CPD = self.random_nu_CPD, \
                random_nu_CPD_bw = self.random_nu_CPD_bw, \
                drandom_nu_CPD_dT = self.self.drandom_nu_CPD_dT, \
                drandom_nu_CPD_bw_dT = self.drandom_nu_CPD_bw_dT

        # Pass these arrays to the C++ code.

        self.obj.set_lookup_tables(self.obj.ntemp, &temp_1d[0], \
                &planck_opacity[0], &rosseland_extinction[0], \
                &dplanck_opacity_dT[0], &drosseland_extinction_dT[0],\
                &dkextdnu[0], &dalbedodnu[0], &random_nu_CPD[0,0], \
                &random_nu_CPD_bw[0,0], &drandom_nu_CPD_dT[0,0], \
                &drandom_nu_CPD_bw_dT[0,0])

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

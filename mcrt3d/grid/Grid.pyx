import cython

from libcpp cimport bool
from libcpp.vector cimport vector

import numpy
cimport numpy

import decimal
import h5py

from ..mcrt3d cimport Grid
from ..sources.Star cimport StarObj
from ..dust.Dust cimport DustObj

cdef class GridObj:
    def __init__(self):
        self.density = []
        self.mass = []
        self.temperature = []
        self.dust = []
        self.sources = []

        decimal.getcontext().prec = 80
        y = [i * decimal.Decimal(1) / decimal.Decimal(100) for i in range(101)]
        f = []
        for i in range(len(y)-1):
            n = 1
            while True:
                if n == 1:
                    f.append((-1)**(n+1) * y[i]**(n**2))
                else:
                    f[-1] += (-1)**(n+1) * y[i]**(n**2)
                    if (abs((f[-1] - old_f)) <= decimal.Decimal(1.0e-80) * \
                            abs(f[-1])):
                        break

                n = n+1
                old_f = f[-1]

        f.append(decimal.Decimal(1) / decimal.Decimal(2))
        f = decimal.Decimal(2) * numpy.array(f)

        """
        cdef numpy.ndarray[double, ndim=1, mode="c"] \
                dydf = numpy.diff(y) / numpy.diff(f)
        """
        dydf = numpy.diff(y) / numpy.diff(f)

        self.y = numpy.array(y, dtype=float)
        self.f = numpy.array(f, dtype=float)
        self.dydf = dydf.astype(float)
        self.ny = self.y.size

        cdef numpy.ndarray[double, ndim=1, mode="c"] yy = self.y, ff = self.f, \
                dyydff = self.dydf

        self.obj.set_mrw_tables(&yy[0], &ff[0], &dyydff[0], self.ny)

    def __del__(self):
        del self.obj

    property temperature:
        def __get__(self):
            return self.temperature

    def add_density(self, numpy.ndarray[double, ndim=3, mode="c"] density, \
            DustObj dust):

        cdef numpy.ndarray[double, ndim=3, mode="c"] \
                mass = density * self.volume, \
                temperature = numpy.ones((self.n1,self.n2,self.n3), dtype=float)

        self.density.append(density)
        self.mass.append(mass)
        self.temperature.append(temperature)
        self.dust.append(dust)

        self.obj.add_density(&density[0,0,0], &temperature[0,0,0], \
                &mass[0,0,0], dust.obj)

    def add_source(self, StarObj source):
        self.sources.append(source)

        self.obj.add_source(source.obj)

    def read(self, filename=None, usefile=None):
        if (usefile == None):
            f = h5py.File(filename, "r")
        else:
            f = usefile

        coordsystem = f['coordsystem'].value
        w1 = f['w1'].value
        w2 = f['w2'].value
        w3 = f['w3'].value

        if (coordsystem == 'cartesian'):
            self.set_cartesian_grid(w1, w2, w3)
        elif (coordsystem == 'cylindrical'):
            self.set_cylindrical_grid(w1, w2, w3)
        elif (coordsystem == 'spherical'):
            self.set_spherical_grid(w1, w2, w3)

        density = f['Density']
        for name in density:
            self.density.append(density[name].value)

        dust = f['Dust']
        for name in dust:
            d = DustObj()
            d.set_properties_from_file(usefile=dust[name])
            self.dust.append(d)

        temperature = f['Temperature']
        for name in temperature:
            self.temperature.append(temperature[name].value)

        sources = f['Sources']
        for name in sources:
            source = StarObj()
            source.read(usefile=sources[name])
            self.sources.append(source)

        if (usefile == None):
            f.close()

    def write(self, filename=None, usefile=None):
        if (usefile == None):
            f = h5py.File(filename, "w")
        else:
            f = usefile

        f['coordsystem'] = self.coordsystem
        w1_dset = f.create_dataset("w1", (self.w1.size,), dtype='f')
        w1_dset[...] = self.w1
        w2_dset = f.create_dataset("w2", (self.w2.size,), dtype='f')
        w2_dset[...] = self.w2
        w3_dset = f.create_dataset("w3", (self.w3.size,), dtype='f')
        w3_dset[...] = self.w3

        density = f.create_group("Density")
        density_dsets = []
        for i in range(len(self.density)):
            density_dsets.append(density.create_dataset( \
                    "Density{0:d}".format(i), self.density[i].shape, dtype='f'))
            density_dsets[i][...] = self.density[i]

        dust = f.create_group("Dust")
        dust_groups = []
        for i in range(len(self.dust)):
            dust_groups.append(dust.create_group("Dust{0:d}".format(i)))
            self.dust[i].write(usefile=dust_groups[i])

        temperature = f.create_group("Temperature")
        temperature_dsets = []
        for i in range(len(self.temperature)):
            temperature_dsets.append(temperature.create_dataset( \
                    "Temperature{0:d}".format(i), self.temperature[i].shape, \
                    dtype='f'))
            temperature_dsets[i][...] = self.temperature[i]

        sources = f.create_group("Sources")
        sources_groups = []
        for i in range(len(self.sources)):
            sources_groups.append(sources.create_group("Star{0:d}".format(i)))
            self.sources[i].write(usefile=sources_groups[i])

        if (usefile == None):
            f.close()

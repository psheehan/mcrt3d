import os
import ctypes
import numpy
from . import misc

au = 1.496e13
c_l = 2.99792458e10
sigma = 5.67051e-5
Msun = 1.9891e33
Lsun = 3.839e33
Rsun = 6.955e10

lib = ctypes.cdll.LoadLibrary(os.path.dirname(__file__)+'/../src/libmcrt3d.so')
lib.new_CartesianGrid.restype = ctypes.c_void_p
lib.new_CylindricalGrid.restype = ctypes.c_void_p
lib.new_SphericalGrid.restype = ctypes.c_void_p
lib.new_Dust.restype = ctypes.c_void_p
lib.new_Source.restype = ctypes.c_void_p
lib.new_Camera.restype = ctypes.c_void_p
lib.new_Image.restype = ctypes.c_void_p
lib.new_mcrt.restype = ctypes.c_void_p

class Dust:
    def __init__(self):
        self.obj = lib.new_Dust()

    def set_properties_from_file(self,filename):
        f = open(filename,"r")
        f.readline()
        self.nlam = int(f.readline())
        f.close()

        data = numpy.loadtxt(filename, skiprows=2)

        self.lam = data[:,0].copy() * 1.0e-4
        self.nu = c_l / self.lam
        self.kabs = data[:,1].copy()
        self.ksca = data[:,2].copy()
        self.kext = self.kabs + self.ksca
        self.albedo = self.ksca / self.kext

        lib.set_optical_properties(ctypes.c_void_p(self.obj), self.nlam, \
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

        self.int_Bnu = numpy.ones(self.temp.size)
        self.int_Bnukext = numpy.ones(self.temp.size)
        self.int_dBnukext = numpy.ones(self.temp.size)
        for i in range(self.temp.size-1):
            Bnu = misc.B_nu(self.nu,self.temp[i])
            dBnu = misc.dB_nu(self.nu,self.temp[i])

            self.int_Bnu[i] = numpy.trapz(Bnu,x=self.nu)
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
        self.int_Bnu_knu_nu = numpy.zeros((self.temp.size,self.nu.size))
        self.int_dBnu_knu_nu = numpy.zeros((self.temp.size,self.nu.size))
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

            for j in range(self.nu.size):
                self.int_Bnu_knu_nu[i,j] = numpy.trapz(self.Bnu[i,0:j]* \
                        self.kext[0:j],x=self.nu[0:j])
                self.int_dBnu_knu_nu[i,j] = numpy.trapz(self.dBnu[i,0:j]* \
                        self.kext[0:j],x=self.nu[0:j])

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
                self.ddBnudT.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), \
                self.int_Bnu_knu_nu.ctypes.data_as(ctypes.POINTER( \
                        ctypes.c_double)), \
                self.int_dBnu_knu_nu.ctypes.data_as(ctypes.POINTER( \
                        ctypes.c_double)))

class Source:
    def __init__(self):
        self.obj = lib.new_Source()

    def set_parameters(self, x, y, z, mass, radius, temperature):
        self.x = x
        self.y = y
        self.z = z
        self.mass = mass
        self.radius = radius
        self.temperature = temperature

        lib.set_parameters(ctypes.c_void_p(self.obj), ctypes.c_double(x), \
                ctypes.c_double(y), ctypes.c_double(z), ctypes.c_double(mass), \
                ctypes.c_double(radius), ctypes.c_double(temperature))

    def set_blackbody_spectrum(self, nu):
        self.nnu = nu.size
        self.nu = nu
        self.Bnu = misc.B_nu(nu, self.temperature)
        self.luminosity = 4*numpy.pi*self.radius**2*sigma*self.temperature**4

        lib.set_blackbody_spectrum(ctypes.c_void_p(self.obj), self.nnu, \
                self.nu.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), \
                self.Bnu.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), \
                ctypes.c_double(self.luminosity))

class Grid:
    def set_cartesian_grid(self,x,y,z):
        self.obj = lib.new_CartesianGrid()
        self.grid_type = "Cartesian"

        self.set_walls(x,y,z)

    def set_cylindrical_grid(self,r,phi,z):
        self.obj = lib.new_CylindricalGrid()
        self.grid_type = "Cylindrical"

        self.set_walls(r,phi,z)

    def set_spherical_grid(self,r,theta,phi):
        self.obj = lib.new_SphericalGrid()
        self.grid_type = "Spherical"

        self.set_walls(r,theta,phi)

    def set_walls(self,w1,w2,w3):
        self.w1 = w1
        self.w2 = w2
        self.w3 = w3

        lib.set_walls(ctypes.c_void_p(self.obj), w1.size-1, w2.size-1, \
                w3.size-1, w1.size, w2.size, w3.size, \
                w1.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), \
                w2.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), \
                w3.ctypes.data_as(ctypes.POINTER(ctypes.c_double)))

    def set_physical_properties(self, dens):
        volume = numpy.ones(dens.shape ,dtype=float)
        for i in range(volume.shape[0]):
            for j in range(volume.shape[1]):
                for k in range(volume.shape[2]):
                    if (self.grid_type == "Cartesian"):
                        volume[i,j,k] = (self.w1[i+1] - self.w1[i])* \
                            (self.w2[j+1] - self.w2[j])* \
                            (self.w3[k+1] - self.w3[k])
                    elif (self.grid_type == "Cylindrical"):
                        volume[i,j,k] = (self.w1[i+1]**2 - self.w1[i]**2)* \
                            (self.w2[j+1] - self.w2[j])* \
                            (self.w3[k+1] - self.w3[k])/2
                    elif (self.grid_type == "Spherical"):
                        volume[i,j,k] = (self.w1[i+1]**3 - self.w1[i]**3)* \
                            (self.w3[k+1] - self.w3[k])* \
                            (numpy.cos(self.w2[j]) - numpy.cos(self.w2[j+1]))/3

        mass = dens * volume
        temp = numpy.ones(dens.shape ,dtype=float)

        self.dens = numpy.array([dens])
        self.temp = numpy.array([temp])
        self.mass = numpy.array([mass])
        self.volume = volume

        lib.set_physical_properties(ctypes.c_void_p(self.obj), \
                self.dens.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), \
                self.temp.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), \
                self.mass.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), \
                volume.ctypes.data_as(ctypes.POINTER(ctypes.c_double)))

    def set_dust(self, dust):
        self.dust = dust

        lib.create_dust_array(ctypes.c_void_p(self.obj), dust.size)

        for i in range(dust.size):
            lib.set_dust(ctypes.c_void_p(self.obj), \
                    ctypes.c_void_p(dust[i].obj),i)

    def set_sources(self, sources):
        self.sources = sources

        lib.create_sources_array(ctypes.c_void_p(self.obj),sources.size)

        for i in range(sources.size):
            lib.set_sources(ctypes.c_void_p(self.obj), \
                    ctypes.c_void_p(sources[i].obj),i)

class Camera:
    def __init__(self, G):
        self.obj = lib.new_Camera()

        lib.set_camera_grid(ctypes.c_void_p(self.obj), ctypes.c_void_p(G.obj))

    def make_image(self, image):
        lib.make_image(ctypes.c_void_p(self.obj), ctypes.c_void_p(image.obj))

class Image:
    def __init__(self, r, incl, pa, x, y, nx, ny, nu, pixel_size, nnu):
        self.r = r
        self.incl = incl
        self.pa = pa
        self.x = x
        self.y = y
        self.intensity = numpy.zeros((nx, ny, nnu),dtype=float)
        self.nx = nx
        self.ny = ny
        self.nu = nu
        self.pixel_size = pixel_size
        self.nnu = nnu

        self.obj = lib.new_Image(ctypes.c_double(self.r), \
                ctypes.c_double(self.incl), ctypes.c_double(self.pa), \
                self.x.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), \
                self.y.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), \
                self.intensity.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),\
                ctypes.c_int(self.nx), ctypes.c_int(self.ny), \
                self.nu.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), \
                ctypes.c_double(self.pixel_size), ctypes.c_int(self.nnu))

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

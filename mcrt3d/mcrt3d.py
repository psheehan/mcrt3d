import os
import ctypes
import numpy.ctypeslib as npctypes

# Define a few array types to make passing numpy arrays simpler.

array_1d_double = npctypes.ndpointer(dtype=ctypes.c_double, ndim=1,
                flags='CONTIGUOUS')
array_2d_double = npctypes.ndpointer(dtype=ctypes.c_double, ndim=2,
                flags='CONTIGUOUS')
array_3d_double = npctypes.ndpointer(dtype=ctypes.c_double, ndim=3,
                flags='CONTIGUOUS')
array_4d_double = npctypes.ndpointer(dtype=ctypes.c_double, ndim=4,
                flags='CONTIGUOUS')

# Load the library.

lib = ctypes.cdll.LoadLibrary(os.path.dirname(__file__)+'/../src/libmcrt3d.so')

# Linker functions for the MCRT class.

lib.new_mcrt.restype = ctypes.c_void_p
lib.new_mcrt.argtypes = None

lib.set_Grid.restype = None
lib.set_Grid.argtypes = [ctypes.c_void_p, ctypes.c_void_p]

lib.set_Params.restype = None
lib.set_Params.argtypes = [ctypes.c_void_p, ctypes.c_void_p]

lib.run_thermal_mc.restype = None
lib.run_thermal_mc.argtypes = [ctypes.c_void_p]

# Linker functions for the Params class.

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

lib.set_mrw.restype = None
lib.set_mrw.argtypes = [ctypes.c_void_p, ctypes.c_bool]

lib.set_mrw_gamma.restype = None
lib.set_mrw_gamma.argtypes = [ctypes.c_void_p, ctypes.c_double]

# Linker functions for the Grid classes.

lib.new_CartesianGrid.restype = ctypes.c_void_p
lib.new_CartesianGrid.argtypes = None

lib.new_CylindricalGrid.restype = ctypes.c_void_p
lib.new_CylindricalGrid.argtypes = None

lib.new_SphericalGrid.restype = ctypes.c_void_p
lib.new_SphericalGrid.argtypes = None

lib.set_walls.restype = None
lib.set_walls.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_int, \
        ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_int, \
        array_1d_double, array_1d_double, array_1d_double, array_3d_double]

lib.create_dust_array.restype = None
lib.create_dust_array.argtypes = [ctypes.c_void_p, ctypes.c_int]

lib.create_physical_properties_arrays.restype = None
lib.create_physical_properties_arrays.argtypes = [ctypes.c_void_p, ctypes.c_int]

lib.set_dust.restype = None
lib.set_dust.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_int]

lib.set_physical_properties.restype = None
lib.set_physical_properties.argtypes = [ctypes.c_void_p, array_3d_double, \
        array_3d_double, array_3d_double, ctypes.c_int]

lib.create_sources_array.restype = None
lib.create_sources_array.argtypes = [ctypes.c_void_p, ctypes.c_int]

lib.set_sources.restype = None
lib.set_sources.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_int]

lib.set_mrw_tables.restype = None
lib.set_mrw_tables.argtypes = [ctypes.c_void_p, array_1d_double, \
        array_1d_double, array_1d_double, ctypes.c_int]

# Linker functions for the Dust class.

lib.new_IsotropicDust.restype = ctypes.c_void_p
lib.new_IsotropicDust.argtypes = None

lib.set_optical_properties.restype = None
lib.set_optical_properties.argtypes = [ctypes.c_void_p, ctypes.c_int, \
        array_1d_double, array_1d_double, array_1d_double, array_1d_double, \
        array_1d_double, array_1d_double]

lib.set_lookup_tables.restype = None
lib.set_lookup_tables.argtypes = [ctypes.c_void_p, ctypes.c_int, \
        array_1d_double, array_1d_double, array_1d_double, array_1d_double, \
        array_1d_double, array_1d_double, array_1d_double, array_2d_double, \
        array_2d_double, array_2d_double, array_2d_double]

# Linker functions for the Source class.

lib.new_Source.restype = ctypes.c_void_p
lib.new_Source.argtypes = None

lib.set_parameters.restype = None
lib.set_parameters.argtypes = [ctypes.c_void_p, ctypes.c_double, \
        ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, \
        ctypes.c_double]

lib.set_blackbody_spectrum.restype = None
lib.set_blackbody_spectrum.argtypes = [ctypes.c_void_p, ctypes.c_int, \
        array_1d_double, array_1d_double, ctypes.c_double, array_1d_double]

# Linker functions for the Camera class.

lib.new_Camera.restype = ctypes.c_void_p
lib.new_Camera.argtypes = None

lib.set_camera_grid.restype = None
lib.set_camera_grid.argtypes = [ctypes.c_void_p, ctypes.c_void_p]

lib.set_camera_params.restype = None
lib.set_camera_params.argtypes = [ctypes.c_void_p, ctypes.c_void_p]

lib.make_image.restype = None
lib.make_image.argtypes = [ctypes.c_void_p, ctypes.c_void_p]

# Linker functions for the Image class.

lib.new_Image.restype = ctypes.c_void_p
lib.new_Image.argtypes = [ctypes.c_double, ctypes.c_double, ctypes.c_double, \
        array_1d_double, array_1d_double, array_3d_double, ctypes.c_int, \
        ctypes.c_int, array_1d_double, ctypes.c_double, ctypes.c_int]

# Define the Params and MCRT classes here because there isn't yet a better 
# place.

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

    def set_mrw(self, use_mrw):
        lib.set_mrw(self.obj, use_mrw)

    def set_mrw_gamma(self, mrw_gamma):
        lib.set_mrw_gamma(self.obj, mrw_gamma)

class MCRT:
    def __init__(self, G, Q):
        self.obj = lib.new_mcrt()

        lib.set_Grid(self.obj, G.obj)
        lib.set_Params(self.obj, Q.obj)

    def thermal_mc(self):
        lib.run_thermal_mc(self.obj)

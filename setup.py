from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

import numpy

mcrt3d = Extension("mcrt3d.mcrt3d", \
        sources=["mcrt3d/mcrt3d.pxd","mcrt3d/mcrt3d.pyx", "src/mcrt3d.cc", "src/camera.cc", \
        "src/isotropic_dust.cc", "src/cartesian_grid.cc", \
        "src/source.cc", "src/cylindrical_grid.cc", \
        "src/misc.cc", "src/spherical_grid.cc", "src/dust.cc", "src/params.cc",\
        "src/star.cc", "src/grid.cc", "src/photon.cc"], \
        include_dirs=[numpy.get_include(),"./include"], language="c++")

dust = Extension("mcrt3d.dust.Dust", sources=["mcrt3d/Dust.pyx"], \
        include_dirs=[numpy.get_include()], language="c++")

setup(cmdclass = {'build_ext': build_ext}, ext_modules = [mcrt3d, dust])

import os, sys
import tempfile

import setuptools
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext as _build_ext

import numpy
import pybind11

def has_flag(compiler, flagname):
    """Return a boolean indicating whether a flag name is supported on
    the specified compiler.
    Borrowed from the george package.
    """
    with tempfile.NamedTemporaryFile("w", suffix=".cpp") as f:
        f.write("int main (int argc, char **argv) { return 0; }")
        try:
            compiler.compile([f.name], extra_postargs=[flagname])
        except setuptools.distutils.errors.CompileError:
            return False
    return True

# Get the compile arguments.

class build_ext(_build_ext):
    def build_extensions(self):
        # Get the include directories.

        include_dirs=[os.path.join("include"),
                pybind11.get_include(False), 
                pybind11.get_include(True)
        ]

        for ext in self.extensions:
            ext.include_dirs = include_dirs + ext.include_dirs

        # Get the proper flags.

        flags = ['-std=c++11','-stdlib=libc++','-Ofast',"-funroll-loops",\
                "-Wno-unused-function","-Wno-uninitialized",
                "-Wno-unused-local-typedefs",'-march=native',
                '-mmacosx-version-min=10.9','-fopenmp']

        opts = []
        for flag in flags:
            if has_flag(self.compiler, flag):
                opts.append(flag)

        for ext in self.extensions:
            ext.extra_compile_args = opts

        # Add link args, if appropriate.

        for ext in self.extensions:
            ext.extra_link_args += ["-mmacosx-version-min=10.9",
                    "-march=native"]
            if has_flag(self.compiler, "-fopenmp"):
                ext.extra_link_args += ["-fopenmp"]

        # Run the standard build procedure.

        _build_ext.build_extensions(self)

# Create the extension.

if __name__ == "__main__":
    r"""
    dust = Extension("mcrt3d.dust.Dust", sources=["mcrt3d/dust/Dust.pyx"], \
            include_dirs=[numpy.get_include(),"./include"], language="c++")

    star = Extension("mcrt3d.sources.Star", sources=["mcrt3d/sources/Star.pyx"], \
            include_dirs=[numpy.get_include(),"./include"], language="c++")

    grid = Extension("mcrt3d.grid.Grid", sources=["mcrt3d/grid/Grid.pyx"], \
            include_dirs=[numpy.get_include(),"./include"], language="c++")

    cartesian_grid = Extension("mcrt3d.grid.CartesianGrid", \
            sources=["mcrt3d/grid/CartesianGrid.pyx"], \
            include_dirs=[numpy.get_include(),"./include"], language="c++")

    cylindrical_grid = Extension("mcrt3d.grid.CylindricalGrid", \
            sources=["mcrt3d/grid/CylindricalGrid.pyx"], \
            include_dirs=[numpy.get_include(),"./include"], language="c++")

    spherical_grid = Extension("mcrt3d.grid.SphericalGrid", \
            sources=["mcrt3d/grid/SphericalGrid.pyx"], \
            include_dirs=[numpy.get_include(),"./include"], language="c++")

    image = Extension("mcrt3d.camera.Image", sources=["mcrt3d/camera/Image.pyx"], \
            include_dirs=[numpy.get_include(),"./include"], language="c++")

    spectrum = Extension("mcrt3d.camera.Spectrum", \
            sources=["mcrt3d/camera/Spectrum.pyx"], \
            include_dirs=[numpy.get_include(),"./include"], language="c++")

    camera = Extension("mcrt3d.camera.Camera",sources=["mcrt3d/camera/Camera.pyx"],\
            include_dirs=[numpy.get_include(),"./include"], language="c++")

    setup(cmdclass = {'build_ext': build_ext}, ext_modules = [mcrt3d, dust, star, \
            grid, cartesian_grid, cylindrical_grid, spherical_grid, image, \
            spectrum, camera])
    """

    mcrt3d = Extension("mcrt3d.mcrt3d", \
            sources=["src/mcrt3d.cc", \
            "src/camera.cc", "src/isotropic_dust.cc", "src/cartesian_grid.cc", \
            "src/source.cc", "src/cylindrical_grid.cc", "src/misc.cc", \
            "src/spherical_grid.cc", "src/dust.cc", "src/params.cc",\
            "src/star.cc", "src/grid.cc", "src/photon.cc"], \
            language="c++")

    setup(name="mcrt3d", 
            version="0.9.0", 
            author="Patrick Sheehan",
            author_email="psheehan@northwestern.edu",
            packages=["mcrt3d"], 
            ext_modules=[mcrt3d],
            description="Monte Carlo radiative transfer.",
            install_requires=["numpy","pybind11"],
            cmdclass=dict(build_ext=build_ext)
            )

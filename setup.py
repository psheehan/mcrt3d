import os, sys
import tempfile

import setuptools
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext

import numpy
import pybind11

# CUDA compilation is adapted from the source
# https://github.com/rmcgibbo/npcuda-example
# CUDA functions for compilation
def find_in_path(name, path):
    """
    Find a file in a search path
    """
    
    #adapted fom 
    #http://code.activestate.com/recipes/52224-find-a-file-given-a-search-path/
    for dir in path.split(os.pathsep):
        binpath = os.path.join(dir, name)
        if os.path.exists(binpath):
            return os.path.abspath(binpath)
    return None

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

def locate_cuda():
    """
    Locate the CUDA environment on the system
    Returns a dict with keys 'home', 'nvcc', 'include', and 'lib64'
    and values giving the absolute path to each directory.
    Starts by looking for the CUDAHOME env variable. If not found, everything
    is based on finding 'nvcc' in the PATH.
    """
    
    # first check if the CUDAHOME env variable is in use
    if 'CUDAHOME' in os.environ:
        home = os.environ['CUDAHOME']
        nvcc = os.path.join(home, 'bin', 'nvcc')
    else:
        # otherwise, search the PATH for NVCC
        nvcc = find_in_path('nvcc', os.environ['PATH'])
        if nvcc is None:
            raise EnvironmentError('The nvcc binary could not be located in '
                    'your $PATH. Either add it to your path, or set $CUDAHOME')
        home = os.path.dirname(os.path.dirname(nvcc))

    cudaconfig = {
            'home': home, 
            'nvcc': nvcc, 
            'include': os.path.join(home, 'include'), 
            'lib64': os.path.join(home, 'lib64')}
    
    for k, v in iter(cudaconfig.items()):
        if not os.path.exists(v):
            raise EnvironmentError('The CUDA %s path could not be located '
                    'in %s' % (k, v))

    return cudaconfig

def customize_compiler_for_nvcc(self):
    """
    inject deep into distutils to customize how the dispatch
    to gcc/nvcc works.
    If you subclass UnixCCompiler, it's not trivial to get your subclass
    injected in, and still have the right customizations (i.e.
    distutils.sysconfig.customize_compiler) run on it. So instead of going
    the OO route, I have this. Note, it's kindof like a wierd functional
    subclassing going on.
    """
    
    # tell the compiler it can processes .cu
    self.src_extensions.append('.cu')

    # save references to the default compiler_so and _comple methods
    default_compiler_so = self.compiler_so
    super = self._compile

    # now redefine the _compile method. This gets executed for each
    # object but distutils doesn't have the ability to change compilers
    # based on source extension: we add it.
    def _compile(obj, src, ext, cc_args, extra_postargs, pp_opts):
        if os.path.splitext(src)[1] == '.cu':
            # use the cuda for .cu files
            self.set_executable('compiler_so', CUDA['nvcc'])
            # use only a subset of the extra_postargs, which are 1-1 translated
            # from the extra_compile_args in the Extension class
            postargs = extra_postargs['nvcc']
        else:
            postargs = extra_postargs['gcc']

        super(obj, src, ext, cc_args, postargs, pp_opts)
        # reset the default compiler_so, which we might have changed for cuda
        self.compiler_so = default_compiler_so

    # inject our redefined _compile method into the class
    self._compile = _compile

# run the customize_compiler
class cuda_build_ext(build_ext):
    def build_extensions(self):
        # Get the include directories.

        include_dirs=[os.path.join("include"), 
                pybind11.get_include(False),
                pybind11.get_include(True),
        ]
        try:
            include_dirs.append(CUDA['include'])
        except:
            pass

        for ext in self.extensions:
            ext.include_dirs = include_dirs + ext.include_dirs

        # Get the proper flags.

        for ext in self.extensions:
            for flag in ext.extra_compile_args['gcc']:
                if not has_flag(self.compiler, flag):
                    ext.extra_compile_args['gcc'].remove(flag)

        # Get the proper links.

        for ext in self.extensions:
            for flag in ext.extra_link_args:
                if not has_flag(self.compiler, flag):
                    ext.extra_link_args.remove(flag)

        # Make sure that .cu files are an allowed type to compile.

        customize_compiler_for_nvcc(self.compiler)

        # Now run the standard build procedure.

        build_ext.build_extensions(self)

# Locate CUDA paths
try:
    CUDA = locate_cuda()
    cuda_found = True
except:
    cuda_found = False

# Obtain the numpy include directory. This logic works across numpy versions.
try:
    numpy_include = numpy.get_include()
except AttributeError:
    numpy_include = numpy.get_numpy_include()


if cuda_found:
    cuda = Extension("mcrt3d.mcrt3d", sources=["src/mcrt3d.cu"], language="c++",
            library_dirs=[CUDA['lib64']],
            libraries=['cudart'],
            runtime_library_dirs=[CUDA['lib64']],
            extra_compile_args={
                'gcc': [],
                'nvcc': [
                    '-O3', '-arch=sm_75', '--use_fast_math', 
                    '--ptxas-options=-v', '-c', 
                    '--compiler-options', "'-fPIC'"]},
            extra_link_args=['-lcudadevrt', '-lcudart'])

    extensions = [cuda]
else:
    cpu = Extension("mcrt3d.mcrt3d", sources=["src/mcrt3d.cc"], language="c++",
            extra_compile_args={
                'gcc': ['-std=c++11','-stdlib=libc++','-Ofast',"-funroll-loops",\
                    "-Wno-unused-function","-Wno-uninitialized",
                    "-Wno-unused-local-typedefs",'-march=native',
                    '-mmacosx-version-min=10.9','-fopenmp','-fvisibility=hidden'],
                'nvcc': []},
            extra_link_args=["-march=native",'-fopenmp',
                "-mmacosx-version-min=10.9"])

    extensions = [cpu]

setup(name="mcrt3d", 
        version="0.9.0", 
        author="Patrick Sheehan",
        author_email="psheehan@nrao.edu",
        packages=["mcrt3d"],
        ext_modules=extensions,
        description="Monte Carlo radiative transfer.",
        install_requires=["numpy","pybind11"],
        cmdclass=dict(build_ext=cuda_build_ext),
        zip_safe=False,
        )

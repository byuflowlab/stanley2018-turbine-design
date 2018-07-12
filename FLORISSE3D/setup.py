#!/usr/bin/env python
# encoding: utf-8

from numpy.distutils.core import setup, Extension

module1 = Extension('_floris', sources=['src/FLORISSE3D/floris.f90', 'src/FLORISSE3D/adStack.c', 'src/FLORISSE3D/adBuffer.f'],
                   extra_compile_args=['-O2', '-c'])

module2 = Extension('_florisDiscontinuous', sources=['src/FLORISSE3D/florisDiscontinuous.f90', 'src/FLORISSE3D/adStack.c', 'src/FLORISSE3D/adBuffer.f'],
                   extra_compile_args=['-O2', '-c'])

module3 = Extension('_shellbuckling', sources=['src/FLORISSE3D/ShellBuckling.f90'],
                   extra_compile_args=['-O2', '-c'])

module4 = Extension('_axialShear', sources=['src/FLORISSE3D/Axial_Shear.f90'],
                   extra_compile_args=['-O2', '-c'])

setup(
    name='FLORISSE3D',
    version='0.0.0',
    description='differentiable floris wake model with cosine factor',
    install_requires=['openmdao>=1.5','akima>=1.0.0'],
    package_dir={'': 'src'},
    ext_modules=[module1, module2, module3, module4],
    dependency_links=['https://github.com/andrewning/akima/tarball/master#egg=akima'],
    packages=['FLORISSE3D'],
    license='Apache License, Version 2.0',
)

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

setup(
    cmdclass = {'build_ext': build_ext},
    ext_modules = [Extension("tv_l1", 
                             sources=['tv_l1.pyx', 'tv_l1_johnson.c'])]
)

    
    #python setup.py build_ext --inplace
#!/usr/bin/env python3

import sys, platform

try:
    from setuptools import *
except ImportError:
    from distribute_setup import use_setuptools
    use_setuptools()
finally:
    from setuptools import *

from glob import glob

if sys.version_info < (2, 6):
    print >> sys.stderr, "ERROR: shmlast requires python 2.6 or greater"
    sys.exit()

import shmlast

cmdclass = {}

def main():
    setup(  name = 'shmlast',
            version = shmlast.__version__,
            description = 'reciprocal and conditional reciprocal best LAST',
            url = 'https://github.com/camillescott/shmlast',
            author = 'Camille Scott',
            author_email = 'camille.scott.w@gmail.com',
            license = 'BSD',
            test_suite = 'nose.collector',
            tests_require = ['nose'],
            packages = find_packages(),
            scripts = glob('bin/*'),
            install_requires = ['pandas',
                                'numpy',
                                ],
            zip_safe = False,
            include_package_data = True,
            cmdclass = cmdclass  )
        
if __name__ == "__main__":
    main()

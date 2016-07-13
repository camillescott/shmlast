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

if sys.version_info < (3, 3):
    print >> sys.stderr, "ERROR: shmlast requires python 3.3 or greater"
    sys.exit()

import shmlast


def main():
    setup(  name = 'shmlast',
            version = shmlast.__version__,
            description = 'reciprocal and conditional reciprocal best LAST',
            url = 'https://github.com/camillescott/shmlast',
            author = 'Camille Scott',
            author_email = 'camille.scott.w@gmail.com',
            license = 'BSD',
            packages = find_packages(),
            scripts = glob('bin/*'),
            setup_requires = ['pytest-runner'],
            tests_require = ['pytest'],
            install_requires = ['doit>=0.29.0',
                                'ficus>=0.3',
                                'matplotlib>==1.4',
                                'numpy>=1.10.0',
                                'pandas>=0.17.0',
                                'scipy>=0.17.0',
                                'screed>=0.9',
                                'seaborn>=0.6.0'],
                zip_safe = True,
                include_package_data = True )
            
if __name__ == "__main__":
    main()

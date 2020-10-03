#!/usr/bin/env python3

import sys, platform, os

try:
    from setuptools import *
except ImportError:
    from distribute_setup import use_setuptools
    use_setuptools()
finally:
    from setuptools import *

from glob import glob

if sys.version_info < (3, 5):
    print >> sys.stderr, "ERROR: shmlast requires python 3.4 or greater"
    sys.exit()

__version__ = open(os.path.join('shmlast', 'VERSION')).read().strip()

SCRIPTS=glob('bin/*')

def main():
    setup(  name = 'shmlast',
            version = __version__,
            description = 'An improved implementation of Conditional Reciprocal'\
                          ' Best Hits with LAST and Python.',
            long_description = open("README.rst").read(),
            url = 'https://github.com/camillescott/shmlast',
            author = 'Camille Scott',
            author_email = 'camille.scott.w@gmail.com',
            license = 'BSD',
            packages = find_packages(),
            scripts = SCRIPTS,
            tests_require = ['pytest',
                             'pytest-runner',
                             'codecov',
                             'pytest-benchmark'],
            install_requires = ['doit',
                                'ficus',
                                'matplotlib',
                                'numpy',
                                'pandas',
                                'scipy',
                                'screed',
                                'seaborn',
                                'filelock',
                                'ope'],
            zip_safe = False,
            include_package_data = True )
            

if __name__ == "__main__":
    main()

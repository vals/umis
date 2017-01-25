import os
from setuptools import setup, find_packages, Extension

ext = Extension('utils', ['umis/utils.pyx'])

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
        name='umis',
        version='0.5.0',
        description='Package for estimating UMI counts in Transcript Tag Counting data.',
        packages=find_packages(),
        install_requires=['click', 'pysam>=0.8.3', 'pandas', 'regex'],
        ext_modules=[ext],
        setup_requires=['cython'],
        entry_points = {
            'console_scripts': ['umis=umis.umis:umis']
        },
        url='https://github.com/vals/umis',
        author='Valentine Svensson',
        author_email='valentine@nxn.se',
        long_description=read('README.md'),
        package_data = {
        '': ['examples/*/*.json']
        }
    )

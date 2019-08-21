# encoding: utf-8
#
# @Author:    Adam Mendenhall (modified to install nonpip dependencies)
# @Date:      August 8, 2019
# @Filename:  setup.py
# @License:   BSD 3-Clause
# @Copyright: Adam Mendenhall
#

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals

from setuptools import setup, find_packages

import os
import argparse
import sys


# The NAME variable should be of the format "sdss-dewarp".
# Please check your NAME adheres to that format.
NAME = 'sdss-dewarp'
VERSION = '0.0.0'
RELEASE = 'dev' in VERSION


def run(packages, install_requires_pip, install_requires_nopip):

    setup(name=NAME,
          version=VERSION,
          license='BSD3',
          description='imagines warped images of fiducials and metrology fibers and undoes that warp',
          long_description=open('README.rst').read(),
          author='Adam Mendenhall',
          author_email='amendenhall@ucsb.edu',
          keywords='astronomy software',
          url='https://github.com/sdss/dewarp',
          include_package_data=True,
          packages=packages,
          install_requires=install_requires_pip,
          dependency_links=install_requires_nopip,
          package_dir={'': 'python'},
          scripts=['bin/dewarp'],
          classifiers=[
              'Development Status :: 4 - Beta',
              'Intended Audience :: Science/Research',
              'License :: OSI Approved :: BSD License',
              'Natural Language :: English',
              'Operating System :: OS Independent',
              'Programming Language :: Python',
              'Programming Language :: Python :: 3',
              'Topic :: Documentation :: Sphinx',
              'Topic :: Software Development :: Libraries :: Python Modules',
          ],
          )


def get_requirements(opts):
    ''' Get the proper requirements file based on the optional argument '''

    if opts.dev:
        name = 'requirements_dev.txt'
    elif opts.doc:
        name = 'requirements_doc.txt'
    else:
        name = 'requirements.txt'

    requirements_file = os.path.join(os.path.dirname(__file__), name)
    install_requires_pip = [line.strip().replace('==', '>=') for line in open(requirements_file)
                        if not line.strip().startswith('#') and line.strip() != '']

    name = 'notonpip.txt'
    requirements_file = os.path.join(os.path.dirname(__file__), name)
    install_requires_nopip = [line.strip().replace('==', '>=') for line in open(requirements_file)
                        if not line.strip().startswith('#') and line.strip() != '']
    return (install_requires_pip, install_requires_nopip)


def remove_args(parser):
    ''' Remove custom arguments from the parser '''

    arguments = []
    for action in list(parser._get_optional_actions()):
        if '--help' not in action.option_strings:
            arguments += action.option_strings

    for arg in arguments:
        if arg in sys.argv:
            sys.argv.remove(arg)


if __name__ == '__main__':

    # Custom parser to decide whether which requirements to install
    parser = argparse.ArgumentParser(prog=os.path.basename(sys.argv[0]))
    parser.add_argument('-d', '--dev', dest='dev', default=False, action='store_true',
                        help='Install all packages for development')
    parser.add_argument('-o', '--doc', dest='doc', default=False, action='store_true',
                        help='Install only core + documentation packages')

    # We use parse_known_args because we want to leave the remaining args for distutils
    args = parser.parse_known_args()[0]

    # Get the proper requirements file
    install_requires_pip, install_requires_nopip = get_requirements(args)

    # Now we remove all our custom arguments to make sure they don't interfere with distutils
    remove_args(parser)

    # Have distutils find the packages
    packages = find_packages(where='python')

    # Runs distutils
    run(packages, install_requires_pip, install_requires_nopip)

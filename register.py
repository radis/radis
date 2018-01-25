# -*- coding: utf-8 -*-
'''Automatically generate a README.rst for Pypi from my README.md, and publish
the latest version

Use
------------
    python register.py

Requirement
------------
    pandoc
    
'''

from __future__ import print_function

import os
import shutil

package_name = 'radis'

# All the rest below is standard:

if os.path.exists('README.rst'):
    os.remove('README.rst')

os.system('pandoc -s README.md --from markdown --to rst -s -o README.rst')


if os.path.exists('README.rst'):
    print('Readme generated')
    os.system('python setup.py sdist')
    os.system('python setup.py bdist_wheel --universal')
    os.system("twine upload dist/*")
    # Clean
    shutil.rmtree('dist')
    os.remove('README.rst')

print('All done')

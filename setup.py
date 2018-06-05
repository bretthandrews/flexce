import os

from setuptools import setup, find_packages
# from distutils.core import setup

NAME = 'flexce'
# do not use x.x.x-dev.  things complain.  instead use x.x.xdev
VERSION = '1.0.2'
# RELEASE = 'dev' not in VERSION

# requirements
requirements_file = os.path.join(os.path.dirname(__file__), 'requirements.txt')
install_requires = [line.strip().replace('==', '>=') for line in open(requirements_file)
                    if not line.strip().startswith('#') and line.strip() != '']

entry_points = (
    '''
        [console_scripts]
        run_flexce=flexce.run_flexce:main
    ''')

setup(
    name=NAME,
    version=VERSION,
    packages=find_packages(),
    description='Flexible Galactic Chemical Evolution Model',
    author='Brett Andrews',
    author_email='brett.h.andrews@gmail.com',
    url='https://github.com/bretthandrews/flexce',
    license='MIT',
    include_package_data=True,
    install_requires=install_requires,
    entry_points=entry_points,
)

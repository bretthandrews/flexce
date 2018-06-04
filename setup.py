from distutils.core import setup

NAME = 'flexCE'
# do not use x.x.x-dev.  things complain.  instead use x.x.xdev
VERSION = '1.0.1dev'
RELEASE = 'dev' not in VERSION

setup(name=NAME,
      version=VERSION,
      description='Flexible Galactic Chemical Evolution Model',
      author='Brett Andrews',
      author_email='brett.h.andrews@gmail.com',
      url='https://github.com/bretthandrews/flexCE',
      license='MIT',
      packages=['flexCE'])

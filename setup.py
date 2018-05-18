from setuptools import setup, find_packages

install_requirements = ['astropy', 'numpy', 'pandas', 'h5py', 'GitPython', 'cython', 'requests', 'numexpr', 'scipy', 'scikit-learn']
setup_requirements   = ['astropy', 'numpy', 'pandas', 'h5py', 'GitPython', 'cython', 'requests', 'numexpr', 'scipy', 'scikit-learn']
test_requirements    = ['pytest']

from setuptools.command.test import test as TestCommand
import sys


class PyTest(TestCommand):
    user_options = [('pytest-args=', 'a', "Arguments to pass to py.test")]

    def initialize_options(self):
        TestCommand.initialize_options(self)
        self.pytest_args = []

    def finalize_options(self):
        TestCommand.finalize_options(self)
        self.test_args = []
        self.test_suite = True

    def run_tests(self):
        #import here, cause outside the eggs aren't loaded
        import pytest
        errno = pytest.main(self.pytest_args)
        sys.exit(errno)

setup(name='cosmoutils',
      version='0.2.0',
      description='LSS Programs',
      url='http://github.com/vcalderon2009/cosmo_utils2',
      author='Victor Calderon',
      author_email='victor.calderon90@gmail.com',
      license='MIT',
      packages=find_packages(),
      install_requires=install_requirements,
      setup_requires=setup_requirements,
      tests_require=test_requirements,
      test_suite='tests',
      cmdclass = {'test': PyTest},
      zip_safe=False)

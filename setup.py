from setuptools import setup, find_packages

setup(name='atomphys',
      version='0.0.1',
      description='Atomic Physics for Python',
      url='https://gitlab.phys.ethz.ch/tiqi-projects/optical-trap/atomphys',
      author='Matt Grau',
      author_email='graum@phys.ethz.ch',
      package_data={"": ["*.json"]},
      packages=find_packages())

from setuptools import setup, find_packages

setup(name='atomphys',
      version='0.1.0',
      description='Atomic Physics for Python',
      url='https://gitlab.phys.ethz.ch/tiqi-projects/optical-trap/atomphys',
      author='Matt Grau',
      author_email='graum@phys.ethz.ch',
      package_data={'atomphys': ['data/*.json']},
      packages=find_packages(),
      requirements=[],
      extras_require={
          'dev': [
                  'pytest',
                  'pytest-cov',
                  'autopep8'
          ]
      })

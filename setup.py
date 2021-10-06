from setuptools import setup, find_packages

setup(name='atomphys',
      version='0.0.1',
      license='MIT',
      description='Atomic Physics for Python',
      keywords=['AMO', 'atomic physics'],
      author='Matt Grau',
      author_email='matt.grau@gmail.com',
      url='https://github.com/mgrau/atomphys',
      download_url='https://github.com/mgrau/atomphys/archive/refs/tags/v0.0.1-alpha.tar.gz',
      package_data={'atomphys': ['data/*.json']},
      packages=find_packages(),
      requirements=['pint'],
      extras_require={
          'dev': [
              'pytest',
              'pytest-cov',
              'autopep8',
              'mkdocs',
              'mkdocs-material',
              'mkdocstrings'
          ]
      },
      classifiers=[
          # Chose either "3 - Alpha", "4 - Beta" or "5 - Production/Stable" as the current state of your package
          'Development Status :: 3 - Alpha',
          # Define that your audience are developers
          'Intended Audience :: Science/Research',
          'Topic :: Scientific/Engineering :: Physics',
          'License :: OSI Approved :: MIT License',   # Again, pick a license
          # Specify which pyhton versions that you want to support
          'Programming Language :: Python :: 3',
      ],)

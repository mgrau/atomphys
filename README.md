# atomphys_ea

<!--intro-start-->

A Python package to help with atomic physics calculations  (fork of [atomphys](https://github.com/mgrau/atomphys))

[![Tests](https://github.com/mgrau/atomphys/actions/workflows/tests.yml/badge.svg)](https://github.com/mgrau/atomphys/actions/workflows/tests.yml)
[![Codecov](https://img.shields.io/codecov/c/github/mgrau/atomphys)](https://app.codecov.io/gh/mgrau/atomphys)
[![GitHub](https://img.shields.io/github/license/mgrau/atomphys)](LICENSE)
[![PyPI](https://img.shields.io/pypi/v/atomphys)](https://pypi.org/project/atomphys/)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/mgrau/atomphys/main?urlpath=%2Ftree%2Fexamples)

---

**Documentation**: [mgrau.github.io/atomphys/](https://mgrau.github.io/atomphys/)

**Source Code**: [github.com/mgrau/atomphys](https://github.com/mgrau/atomphys)

---

atomphys is meant to be a good starting off point for your atomic physics calculations. It can automate much of the frustrating process of searching for and compiling physical data and simple physical relations, and help you more quickly get to the good stuff.

It's designed with a natural interface and is easy to use.

## Example

```python
>>> from atomphys import Atom
>>> Rb = Atom('Rb')

>>> print(Rb('S1/2').to('P1/2').λ.to('nm'))
795 nm

>>> print(Rb('P1/2').τ.to('ns'))
27.7 ns
```

## Installation

To install atomphys, simply use pip:

```console
$ pip install atomphys
```

## Features

- Integration with [Pint](https://pint.readthedocs.io/en/stable/) for robust handling of units
- Automatically fetch energy level and transition data from the [NIST Atomic Spectra Database](https://www.nist.gov/pml/atomic-spectra-database)
- Use transition data to calculate state lifetimes, polarizabilities, transition dipole moments, cross sections, and saturation intensities

## Requirements

Python 3.6+

atomphys makes extensive use of the excellent package [Pint](https://pint.readthedocs.io/en/stable/) to handle units.

<!--intro-end-->

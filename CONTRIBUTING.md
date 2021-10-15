# Developer Guide

You can contribute to atomphys in different ways.

## Report issues

You can report bugs with the package and gaps in documentation on github in the [issue tracker](https://github.com/mgrau/atomphys/issues). Feel free to also submit feature requests, ask questions, or open general discussions.

## Building the code

To contribute to the code or documentation, you will need a suitable python development environment. If you are just getting started we recommend that you clone atomphys from github.

```console
$ git clone git@github.com:mgrau/atomphys.git
$ cd atomphys
```

We recommend create a new python virtual environment for development,

```console
$ python3 -m venv atomphys_venv
$ source atomphys_venv/bin/activate
```

and install atomphys in editable mode, with the optional dependencies for development or for building documentation.

```console
$ pip install -e .
$ pip install -r requirements_dev.txt
$ pip install -r requirements_docs.txt
```

Install the git hook scripts which performs automated checks of code style before commit

```console
$ pre-commit install
```

### Running tests

Run unit tests using pytest.

```console
$ pytest
```

You can also run the full suite of tests against every atom in the NIST ASD (warning: this is very slow!)

```console
$ pytest -m ASD
```

### Building documentation

You can build the documentatio and run it in a test server with mkdocs,

```console
$ mkdocs serve
```

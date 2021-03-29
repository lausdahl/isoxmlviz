# ISO XML VIsualiZation

This project provides a basic plotting functionality for plotting the content of ISO BUS XML files as described in ISO 11783-10
# Build the project from source

from the repo root run the following commands:

* pipenv and python 3 is required

```bash
# make the environment
pipenv install
# run the project tests in the tests folder
pipenv test
# install the build tools to build the package
pipenv build_update
# build the wheel package
pipenv build
```

# Install the tool in an insulated environment

First select a new folder which does not have a virtual environment already.
The wheel package can be taken from the `dist` folder or from the artifact produced from the job on github.

```bash
mkdir /tmp/isoxml
cd /tmp/isoxml

# install whl
pipenv install dist/isoxmlviz-1.0.0-py3-none-any.whl

pipenv run python -m "isoxmlviz" --help

usage: isoxmlviz [-h] -file FILE [-p]

optional arguments:
  -h, --help  show this help message and exit
  -file FILE  Path to a isoxml task file XML or ZIP
  -p, --pdf   Write figure to pdf

```
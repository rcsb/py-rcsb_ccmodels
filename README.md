# py-rcsb_ccmodels

Utilities for building structural models for RCSB Chemical Component definitions using
small molecule crystal structures in the CCDC.

## Introduction

This package provides workflow utilities for generating search targets for RCSB
chemical component defintions, performing substructure searches on the CCDC database for
each search target, evaluating the these search results and building chemical component
structural models, and assembling selected models into a concatenated model dictionary.

This module depends on the RCSB CCDC wrapper module `py-rcsb_utils_ccdc` which in turn
depends on the proprietary CCDC Python API and CCDC/CSD database.
The latter dependencies require licenses and separate installation that is described
with the CCDC documentation.  The 2021 version of the CCDC API is supported only in Python 3.7
and the RCSB wrapper module provides a CLI that helps to isolate the current module
from this particular version requirement.

### Installation

Download the library source software from the project repository and set the
enviroment corresponding appropriate to the CCDC and wrapper (`py-rcsb_utils_ccdc`)
in the provided script, `ccdc-api-env.sh`. This script provides an example for
a typical `pyenv` virtual configuration. To install from the source repository:

To install using pip:

```bash
# Install and configure the CCDC wrapper package
pip install rcsb.utils.ccdc
#
pip install  rcsb.ccmodels
```

Set the following environmental variables as required for you CCDC database
and Python API installation.

```bash
export CSDHOME=/Applications/CCDC/CSD_2021
# Locate your Python installation in which CCDC Python API is installed.
PYROOT=$PYENV_ROOT/versions/3.7.9
export CSD_PYTHON_ROOT_PATH=$PYROOT
#
# For Linux
export LD_LIBRARY_PATH=$PYROOT/lib:$PYROOT/lib/python3.7/site-packages/ccdc/_lib:$LD_LIBRARY_PATH
# For Macos
export DYLD_LIBRARY_PATH=$PYROOT/lib/python3.7/site-packages/ccdc/_lib
export DYLD_FRAMEWORK_PATH=$PYROOT/lib/python3.7/site-packages/ccdc/_lib
```

To install the packages from the source repository, follow these steps:

```bash
# Install the wrapper package
git clone --recurse-submodules https://github.com/rcsb/py-rcsb_utils_ccdc.git
cd py-rcsb_utils_ccdc
pip install -r requirements
pip install -e .
#

# Install the workflow utilities
git clone --recurse-submodules https://github.com/rcsb/py-rcsb_ccmodels.git
cd py-rcsb_ccmodels
pip install -r requirements
pip install -e .
# edit and set the enviroment in the following bash script ...
. ccdc-api-env.sh

```

Optionally, run test suite (currently Python versions 3.9) using
[setuptools](https://setuptools.readthedocs.io/en/latest/) or
[tox](http://tox.readthedocs.io/en/latest/example/platform.html):

```bash
  # edit and set the enviroment in the following bash script ...
  . ccdc-api-env.sh
  pip install -r requirements.txt
  python setup.py test

or simply run:

 # edit and set the enviroment in the following bash script ...
  . ccdc-api-env.sh
  tox
```

A CLI is provided to simplify access to the workflow steps.

```bash
# edit and set the enviroment in the following bash script ...
. ccdc-api-env.sh
#
python ChemCompModelExec.py --help
   -or-
cc_models_cli --help

usage: cc_models_cli [-h] [--generate] [--cc_locator CC_LOCATOR] [--bird_locator BIRD_LOCATOR]
                     [--prefix PREFIX] [--cache_path CACHE_PATH] [--use_cache]
                     [--limit_perceptions] [--search] [--update_only] [--build] [--build_align_type BUILD_ALIGN_TYPE] [--assemble] [--max_r_factor MAX_R_FACTOR]
                     [--csdhome CSDHOME] [--python_lib_path PYTHON_LIB_PATH]
                     [--python_version PYTHON_VERSION] [--num_proc NUM_PROC] [--chunk_size CHUNK_SIZE]
                     [--verbose]

optional arguments:
  -h, --help            show this help message and exit
  --generate            Generate searchable files
  --cc_locator CC_LOCATOR
                        Chemical component reference dictioanary locator
  --bird_locator BIRD_LOCATOR
                        BIRD reference dictioanary locator
  --prefix PREFIX       Prefix for identifying reference resource files (e.g. abbrev)
  --cache_path CACHE_PATH
                        Top-level cache directory path
  --use_cache           Re-use cached resource files
  --limit_perceptions   Restrict automatic OE chemical perceptions
  --search              Execute CCDC search
  --update_only         Only update current search results
  --build               Build models from CCDC search results
  --build_align_type BUILD_ALIGN_TYPE
                        Alignment criteria (default: graph-relaxed-stereo
  --assemble            Assemble models into a concatenated file
  --max_r_factor MAX_R_FACTOR
                        Maximum permissible R-value (default=10.0)
  --csdhome CSDHOME     Path to the CSD release (path to CSD_202x)
  --python_lib_path PYTHON_LIB_PATH
                        Path to Python library
  --python_version PYTHON_VERSION
                        Python library version (default: 3.7)
  --num_proc NUM_PROC   Number of processes to execute (default=2)
  --chunk_size CHUNK_SIZE
                        Number of files loaded per process
  --verbose             Verbose output
__________________________________________________
```

An example workflow script would look like the following:

```bash
#!/bin/bash
#
echo "Setup environment"
. ./ccdc-api-env.sh
#
echo "Begin search file generation"
cc_models_cli --generate --num_proc 6 --cache_path ./CACHE
#
echo "Begin search workflow"
cc_models_cli --search --num_proc 4 --cache_path ./CACHE
#
echo "Begin build workflow"
cc_models_cli --build --num_proc 4 --cache_path ./CACHE --build_align_type graph-relaxed-stereo-sdeq

echo "Begin assemble workflow"
cc_models_cli --assemble --min_r_factor 10.0  --cache_path ./CACHE
```

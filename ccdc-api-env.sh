#!/bin/bash
#
#  These settings are required to fire up the CCDC/CSD api...
#
#  Adjust these for settings for the site specific install of the
#  CSD database and the CCDC Python API
#
#  Note:  similar vscode defintions are in ./.env
#
#  Packages recommended for some of the extra ccdc examples but not required for this package.
#     pip install matplotlib jinja2 numpy biopython
#  before run, CCDC and OE_LICENSE variables should be modified as needed
case "$(uname -s)" in

   Darwin)
       echo 'Using MacOS settings'
       # Note:  the vscode defintions in ./.env
       
       # CSDS full installation folder
       export CCDC=/path_CSD_installation

       # CSD data folder
       export CSDHOME=$CCDC/CSD_2022

       # CSD Python API folder
       export PYROOT=$CCDC/Python_API_2022/miniconda
       export CSD_PYTHON_ROOT_PATH=$PYROOT
     
       # For Macos
       export DYLD_FRAMEWORK_PATH=$PYROOT/lib/python3.7/site-packages/ccdc/_lib
       export DYLD_LIBRARY_PATH=$PYROOT/lib/python3.7/site-packages/ccdc/_lib
       
       # OpenEye path and license, the license path should be modified as needed
       export OE_DIR=$PYROOT/lib/python3.7/site-packages/openeye
       export OE_LICENSE=/path_oe_license/oe_license.txt

       # Activate miniconda package installed with full CSDS package to run python 3.7 conda env
       source $PYROOT/bin/activate base
       
       ;;

   Linux)
       echo 'Using Linux settings'
       
       # CSDS full installation folder, path should be modified as needed
       export CCDC=/path_CSD_installation

       # CSD data folder
       export CSDHOME=$CCDC/CSD_2022

       # CSD Python API folder
       export PYROOT=$CCDC/Python_API_2022/miniconda
       export CSD_PYTHON_ROOT_PATH=$PYROOT
     
       # For Linux, Python library 
       export LD_LIBRARY_PATH=$PYROOT/lib/python3.7/site-packages/ccdc/_lib:$PYROOT/lib:$LD_LIBRARY_PATH

       # OpenEye path and licence
       export OE_DIR=$PYROOT/lib/python3.7/site-packages/openeye
       export OE_LICENSE=/path_oe_license/oe_license.txt

       # Activate miniconda package installed with full CSDS package to run python 3.7 conda env
       source $PYROOT/bin/activate base
       
       ;;

   *)
       echo 'Unsupported OS - no settings'
       ;;
esac

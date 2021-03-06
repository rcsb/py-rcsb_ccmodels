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
#
case "$(uname -s)" in

   Darwin)
     echo 'Using macos settings'
     #  Note:  the vscode defintions in ./.env
     export CSDHOME=/Applications/CCDC/CSD_2021
     PYROOT=$PYENV_ROOT/versions/3.7.9
     export CSD_PYTHON_ROOT_PATH=$PYROOT
     #
     # For Linux
     export LD_LIBRARY_PATH=$PYROOT/lib:$PYROOT/lib/python3.7/site-packages/ccdc/_lib:$LD_LIBRARY_PATH
     # For Macos
     export DYLD_FRAMEWORK_PATH=$PYROOT/lib/python3.7/site-packages/ccdc/_lib
     export DYLD_LIBRARY_PATH=$PYROOT/lib/python3.7/site-packages/ccdc/_lib
     #
     ;;

   Linux)
     echo 'Using Linux settings'
     export CSDHOME=/ssd/reference/cambridge-2021/CSD_2021
     export PYROOT=$PYENV_ROOT/versions/3.7.9
     export CSD_PYTHON_ROOT_PATH=$PYROOT
     #
     # For Linux
     export LD_LIBRARY_PATH=$PYROOT/lib:$PYROOT/lib/python3.7/site-packages/ccdc/_lib:$LD_LIBRARY_PATH
     # For Macos
     export DYLD_FRAMEWORK_PATH=$PYROOT/lib/python3.7/site-packages/ccdc/_lib
     export DYLD_LIBRARY_PATH=$PYROOT/lib/python3.7/site-packages/ccdc/_lib
     #
     ;;

   # Add here more strings to compare
   # See correspondence table at the bottom of this answer

   *)
     echo 'Unsupported OS - no settings'
     ;;
esac
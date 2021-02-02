##
# File: ChemCompModelExec.py
# Date: 26-Jan-2021  jdw
#
#  Execution wrapper  --  for chemical component model construction workflow
#
#  Updates:
#
##
__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Apache 2.0"

import argparse
import logging
import os
import sys

from rcsb.ccmodels.search import __version__
from rcsb.ccmodels.search.ChemCompModelBuild import ChemCompModelBuild
from rcsb.ccmodels.search.ChemCompModelGen import ChemCompModelGen
from rcsb.ccmodels.search.ChemCompModelSearch import ChemCompModelSearch

# from rcsb.utils.io.MarshalUtil import MarshalUtil

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()


def main():
    parser = argparse.ArgumentParser()
    #
    parser.add_argument("--generate", default=False, action="store_true", help="Generate searchable files")
    parser.add_argument("--cc_locator", default=None, help="Chemical component reference dictioanary locator")
    parser.add_argument("--bird_locator", default=None, help="BIRD reference dictioanary locator")
    parser.add_argument("--prefix", default=None, help="Prefix for identifying reference resource files (e.g. abbrev)")
    parser.add_argument("--cache_path", default=None, help="Top-level cache directory path")
    parser.add_argument("--use_cache", default=False, action="store_true", help="Re-use cached resource files")
    parser.add_argument("--limit_perceptions", default=False, action="store_true", help="Restrict automatic OE chemical perceptions")
    #
    parser.add_argument("--search", default=False, action="store_true", help="Execute CCDC search")
    parser.add_argument("--update_only", default=False, action="store_true", help="Only update current search results")
    #
    parser.add_argument("--build", default=False, action="store_true", help="Build models from CCDC search results")
    parser.add_argument("--assemble", default=False, action="store_true", help="Assemble models into a concatenated file")
    parser.add_argument("--max_r_factor", default=10.0, help="Maximum permissible R-value (default=10.0)")
    #
    parser.add_argument("--csdhome", default=None, help="Path to the CSD release (path to CSD_202x)")
    parser.add_argument("--python_lib_path", default=None, help="Path to Python library")
    parser.add_argument("--python_version", default=None, help="Python library version (default: 3.7)")
    #
    parser.add_argument("--num_proc", default=2, help="Number of processes to execute (default=2)")
    parser.add_argument("--chunk_size", default=10, help="Number of files loaded per process")
    parser.add_argument("--verbose", default=False, action="store_true", help="Verbose output")
    args = parser.parse_args()
    #
    try:
        pyLib = args.python_lib_path if args.python_lib_path else os.path.join(os.environ["PYENV_ROOT"], "versions", "3.7.9", "lib")
        pyRoot = os.path.dirname(pyLib)
        pyVer = args.python_version if args.python_version else "3.7"
        csdHome = args.csdhome if args.csdhome else os.environ["CSDHOME"]
        #
        doGen = args.generate
        ccLocator = args.cc_locator
        birdLocator = args.bird_locator
        prefix = args.prefix
        cachePath = args.cache_path
        numProc = int(args.num_proc)
        chunkSize = int(args.chunk_size)
        useCache = args.use_cache
        limitPerceptions = args.limit_perceptions
        #
        doSearch = args.search
        searchType = "substructure"
        updateOnly = args.update_only
        #
        doBuild = args.build
        doAssemble = args.assemble
        maxRFactor = args.max_r_factor
        verbose = args.verbose
    except Exception as e:
        logger.exception("Argument processing problem %s", str(e))
        parser.print_help(sys.stderr)
        exit(1)
    #
    try:
        os.environ["CSDHOME"] = csdHome
        os.environ["LD_LIBRARY_PATH"] = "%s:%s/python%s/site-packages/ccdc/_lib:$LD_LIBRARY_PATH" % (pyLib, pyLib, pyVer)
        os.environ["DYLD_LIBRARY_PATH"] = "%s/python%s/site-packages/ccdc/_lib" % (pyLib, pyVer)
        os.environ["DYLD_FRAMEWORK_PATH"] = "%s/python%s/site-packages/ccdc/_lib" % (pyLib, pyVer)

        logger.info("Using CSDHOME %s", os.environ["CSDHOME"])
        logger.info("Using DYLD_LIBRARY_PATH %s", os.environ["DYLD_LIBRARY_PATH"])
        logger.info("Using DYLD_FRAMEWORK_PATH %s", os.environ["DYLD_FRAMEWORK_PATH"])
        logger.info("Using source version %s", __version__)
        # -----------
        if doGen:
            ccmG = ChemCompModelGen(cachePath=cachePath, prefix=prefix)
            numMols = ccmG.buildSearchFiles(
                ccUrlTarget=ccLocator,
                birdUrlTarget=birdLocator,
                numProc=numProc,
                chunkSize=chunkSize,
                minCount=30,
                useCache=useCache,
                limitPerceptions=limitPerceptions,
            )
            logger.info("Generated %d search files using prefix %r", numMols, prefix)

        if doSearch:
            ccms = ChemCompModelSearch(cachePath=cachePath, pythonRootPath=pyRoot, csdHome=csdHome, prefix=prefix)
            rL = ccms.search(searchType, updateOnly=updateOnly, numProc=numProc, chunkSize=chunkSize)
            logger.info("Search success count %d", len(rL))

        if doBuild:
            ccmb = ChemCompModelBuild(cachePath=cachePath, prefix=prefix)
            rD = ccmb.build(alignType="relaxed-stereo", numProc=numProc, chunkSize=chunkSize, verbose=verbose)
            logger.info("Built model count %d", len(rD))

        if doAssemble:
            ccmb = ChemCompModelBuild(cachePath=cachePath, prefix=prefix)
            numAssem = ccmb.assemble(maxRFactor=maxRFactor)
            logger.info("Assembled model count %d", numAssem

    except Exception as e:
        logger.exception("Failing with %s", str(e))

    # ----------------------- - ----------------------- - ----------------------- - ----------------------- - ----------------------- -


if __name__ == "__main__":
    main()

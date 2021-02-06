##
# File:    ChemCompModelGenTests.py
# Author:  J. Westbrook
# Date:    12-Jan-2021
# Version: 0.001
#
# Update:
#
#
##
"""
Tests utilities building CCDC searchable source files for chemical component models.
"""

__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Apache 2.0"

import logging
import os
import platform
import resource
import time
import unittest

from rcsb.ccmodels.search import __version__
from rcsb.ccmodels.search.ChemCompModelGen import ChemCompModelGen

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()


class ChemCompModelGenTests(unittest.TestCase):
    skipFlag = True

    def setUp(self):
        self.__startTime = time.time()
        self.__dataPath = os.path.join(HERE, "test-data")
        self.__cachePath = os.path.join(HERE, "test-output", "CACHE")
        self.__ccUrlTarget = os.path.join(self.__dataPath, "components-abbrev.cif")
        self.__birdUrlTarget = os.path.join(self.__dataPath, "prdcc-abbrev.cif")
        #
        logger.debug("Running tests on version %s", __version__)
        logger.info("Starting %s at %s", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()))

    def tearDown(self):
        unitS = "MB" if platform.system() == "Darwin" else "GB"
        rusageMax = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        logger.info("Maximum resident memory size %.4f %s", rusageMax / 10 ** 6, unitS)
        endTime = time.time()
        logger.info("Completed %s at %s (%.4f seconds)", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - self.__startTime)

    def testMoleculeCacheFilesAbbrev(self):
        ccmG = ChemCompModelGen(cachePath=self.__cachePath, prefix="abbrev")
        numMols = ccmG.buildSearchFiles(
            ccUrlTarget=self.__ccUrlTarget,
            birdUrlTarget=self.__birdUrlTarget,
            numProc=2,
            minCount=30,
            useCache=False,
            limitPerceptions=False,
        )
        # For 156 limitPerceptions=False  and 160 for limitPerceptions=True
        self.assertGreaterEqual(numMols, 155)

    @unittest.skipIf(skipFlag, "Build search targets for the full dictionary - troubleshooting test")
    def testBuildMoleculeCacheFilesFull(self):
        ccmG = ChemCompModelGen(cachePath=self.__cachePath)
        ccmG.buildSearchFiles(
            molLimit=None,
            minCount=500,
            fpTypeList=None,
            screenTypeList=None,
            numProc=2,
            useCache=False,
            limitPerceptions=False,
        )

##
# File:    CODModelSearchTests.py
# Author:  J. Westbrook
# Date:    5-Mar-2021
# Version: 0.001
#
# Update:
##
"""
Tests of workflow for searching and COD and recovering matching data files.
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
from rcsb.ccmodels.search.CODModelSearch import CODModelSearch

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()


class CODModelSearchTests(unittest.TestCase):
    abbrevTest = False

    def setUp(self):
        self.__startTime = time.time()
        self.__dataPath = os.path.join(HERE, "test-data")
        self.__cachePath = os.path.join(HERE, "test-output", "CACHE")
        self.__ccUrlTarget = os.path.join(self.__dataPath, "components-abbrev.cif") if CODModelSearchTests.abbrevTest else None
        self.__birdUrlTarget = os.path.join(self.__dataPath, "prdcc-abbrev.cif") if CODModelSearchTests.abbrevTest else None
        #
        logger.debug("Running tests on version %s", __version__)
        logger.info("Starting %s at %s", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()))

    def tearDown(self):
        unitS = "MB" if platform.system() == "Darwin" else "GB"
        rusageMax = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        logger.info("Maximum resident memory size %.4f %s", rusageMax / 10 ** 6, unitS)
        endTime = time.time()
        logger.info("Completed %s at %s (%.4f seconds)", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - self.__startTime)

    def testSearchFilesAbbrev(self):
        csU = CODModelSearch(cachePath=self.__cachePath, ccUrlTarget=self.__ccUrlTarget, birdUrlTarget=self.__birdUrlTarget, numProc=2, prefix="abbrev")
        csU.updateDescriptors()
        numMols = csU.search(molLimit=5000)
        self.assertGreaterEqual(numMols, 3)
        csU.fetchMatchedData(useCache=False)
        #

    @unittest.skipIf(abbrevTest, "Build search targets for the full dictionary - troubleshooting test")
    def testSearchFilesFull(self):
        csU = CODModelSearch(cachePath=self.__cachePath, ccUrlTarget=self.__ccUrlTarget, birdUrlTarget=self.__birdUrlTarget, numProc=6, useCache=True)
        csU.updateDescriptors()
        numMols = csU.search(molLimit=5000)
        self.assertGreaterEqual(numMols, 500)
        numMols = csU.fetchMatchedData(useCache=True)
        self.assertGreaterEqual(numMols, 100)


def suiteSearchTests():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(CODModelSearchTests("testSearchFilesFull"))
    return suiteSelect


if __name__ == "__main__":
    #
    mySuite = suiteSearchTests()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
##
#
# File:    testChemCompModelBuild.py
# Author:  J. Westbrook
# Date:    19-Jan-2021
# Version: 0.001
#
# Updated:
#
##
"""
Test cases for to build models from CCDC search results in the chemical component model workflow.

"""
__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "john.westbrook@rcsb.org"
__license__ = "Apache 2.0"

import logging
import unittest
import time
import os
import os.path
import platform
import resource

from rcsb.ccmodels.search.ChemCompModelBuild import ChemCompModelBuild

from rcsb.ccmodels.search import __version__


HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()
logger.setLevel(logging.INFO)


class ChemCompModelBuildTests(unittest.TestCase):
    def setUp(self):
        self.__cachePath = os.path.join(HERE, "test-output", "CACHE")
        self.__prefix = "abbrev"
        self.__numProc = 2
        self.__chunkSize = 5
        #
        self.__startTime = time.time()
        logger.info("Starting %s (%s) at %s", self.id(), __version__, time.strftime("%Y %m %d %H:%M:%S", time.localtime()))

    def tearDown(self):
        unitS = "MB" if platform.system() == "Darwin" else "GB"
        rusageMax = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        logger.info("Maximum resident memory size %.4f %s", rusageMax / 10 ** 6, unitS)
        endTime = time.time()
        logger.info("Completed %s at %s (%.4f seconds)", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - self.__startTime)

    def testBuildWorkflow(self):
        """Test case: model build workflow step"""
        try:
            ccmb = ChemCompModelBuild(cachePath=self.__cachePath, prefix=self.__prefix)
            rD = ccmb.build()
            logger.info("Matched search ids %r", list(rD.keys()))
            self.assertGreaterEqual(len(rD), 4)
            qD = ccmb.fetchModelIndex()
            self.assertEqual(len(rD), len(qD))

        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()


def suiteBuildTests():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(ChemCompModelBuildTests("testBuildWorkflow"))
    return suiteSelect


if __name__ == "__main__":
    #
    mySuite = suiteBuildTests()
    unittest.TextTestRunner(verbosity=2).run(mySuite)

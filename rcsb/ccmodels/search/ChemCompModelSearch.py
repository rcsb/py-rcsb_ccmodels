##
# File:    ChemCompModelSearch.py
# Author:  J. Westbrook
# Date:    12-Jan-2021
# Version: 0.001
#
# Updated:
##
"""
CCDC search step in chemical component model workflow.


"""
__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "john.westbrook@rcsb.org"
__license__ = "Apache 2.0"

import glob
import logging
import time
import os

from rcsb.ccmodels.search.CcdcSearchExecMp import CcdcSearchExecMp
from rcsb.ccmodels.search.ChemCompModelGen import ChemCompModelGen
from rcsb.ccmodels.search import __version__

# from rcsb.utils.io.MarshalUtil import MarshalUtil

logger = logging.getLogger(__name__)


class ChemCompModelSearch(object):
    def __init__(self, cachePath, pythonRootPath, csdHome, prefix=None):
        self.__cachePath = cachePath
        self.__pythonRootPath = pythonRootPath
        self.__csdHome = csdHome
        self.__prefix = prefix
        self.__startTime = time.time()
        logger.info("Starting search (%s) at %s", __version__, time.strftime("%Y %m %d %H:%M:%S", time.localtime()))

    def getResultDirFilePath(self):
        dN = "cc-%s-result-files" % self.__prefix if self.__prefix else "cc-result-files"
        return os.path.join(self.__cachePath, dN)

    def search(self, searchType, updateOnly=False, numProc=2, chunkSize=10, fileType="sdf", timeOut=240):
        """Run the CCDC search step in chemical component model workflow.

        Args:
            searchType (str): search type: substructure|similarity
            updateOnly (bool, optional): only update search results. Defaults to False.
            numProc (int, optional): number of processors to invoke. Defaults to 2.
            chunkSize (int, optional): incremental chunk size for each subprocess. Defaults to 10.
            fileType (str, optional): molecule file format (sdf|mol2) default: sdf
            timeOut (int, optional): search time out (seconds).  defaults 240

        Returns:
            (list): searchId list for successful searches
        """
        rL = []
        try:
            resultDirPath = self.getResultDirFilePath()
            ccmg = ChemCompModelGen(self.__cachePath, prefix=self.__prefix)
            # sourceSearchDirPath = ccmg.getSearchDirFilePath()
            sTupL = ccmg.fetchPathList()
            # sTupL = self.__getSearchPathListGlob(sourceSearchDirPath, fileType=fileType)
            if updateOnly:
                rD = self.getResultIndex()
                pL = [sTup[1] for sTup in sTupL if sTup[0] not in rD and sTup[2] == fileType]
                logger.info("Updated search list length %d", len(pL))
            else:
                pL = [sTup[1] for sTup in sTupL if sTup[2] == fileType]

            if pL:
                csmp = CcdcSearchExecMp(pythonRootPath=self.__pythonRootPath, csdHome=self.__csdHome)
                resultL = csmp.runSearch(pL, resultDirPath, searchType=searchType, numProc=numProc, chunkSize=chunkSize, timeOut=timeOut)
                for result in resultL:
                    bfn, _ = os.path.splitext(os.path.basename(result))
                    rL.append(bfn)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return rL

    def getResultIndex(self):
        """Return the dictionary index of search results {identifier: idx file path}

        Args:
            resultDirPath (str): directory containing searchable mol2 files

        Returns:
            (dict): {searchId: search index file path}

        """
        rD = {}
        try:
            resultDirPath = self.getResultDirFilePath()
            logger.info("Using result directory path %r", resultDirPath)
            pthL = glob.glob(os.path.join(resultDirPath, "**", "*-index.json"), recursive=True)
            logger.info("Result path list length %d", len(pthL))
            #
            # Parse out the identifier list -
            #
            for pth in pthL:
                bfn, _ = os.path.splitext(os.path.basename(pth))
                sId = bfn.replace("-index", "")
                rD[sId] = pth
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return rD

    def __getSearchPathListGlob(self, sourceSearchDirPath, fileType="sdf"):
        """Return the path list of mol files in the input path

        Args:
            sourceSearchDirPath (str): directory containing searchable files
            fileType (str): mol2 or sdf default: sdf

        Returns:
            (list of tuples): [(searchId, molfile path)]
        """
        tupL = []
        try:
            logger.info("Using mol2 directory path %r", sourceSearchDirPath)
            fPattern = "*." + fileType
            pthL = glob.glob(os.path.join(sourceSearchDirPath, "**", fPattern), recursive=True)
            logger.info("Search list length %d", len(pthL))
            #
            # Parse out the identifier list -
            #
            for pth in pthL:
                bfn, _ = os.path.splitext(os.path.basename(pth))
                tupL.append((bfn, pth))
            #
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return tupL

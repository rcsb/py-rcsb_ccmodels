##
# File:    ChemCompModelGen.py
# Author:  J. Westbrook
# Date:    12-Jan-2021
# Version: 0.001
#
# Updated:
##
"""
Generate searchable source files for chemical component model workflow.

"""
__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "john.westbrook@rcsb.org"
__license__ = "Apache 2.0"

import logging
import time
import os

from rcsb.utils.chem.ChemCompMoleculeProvider import ChemCompMoleculeProvider
from rcsb.utils.chem.ChemCompSearchIndexProvider import ChemCompSearchIndexProvider
from rcsb.utils.chem.OeChemCompUtils import OeChemCompUtils
from rcsb.utils.chem.OeIoUtils import OeIoUtils
from rcsb.utils.chem.OeSearchMoleculeProvider import OeSearchMoleculeProvider
from rcsb.utils.io.MarshalUtil import MarshalUtil

logger = logging.getLogger(__name__)


class ChemCompModelGen(object):
    def __init__(self, cachePath, prefix=None):
        self.__cachePath = cachePath
        self.__prefix = prefix if prefix else ""

    def getSearchDirFilePath(self):
        dN = "cc-%s-search-files" % self.__prefix if self.__prefix else "cc-search-files"
        return os.path.join(self.__cachePath, dN)

    def buildSearchFiles(self, **kwargs):
        """Build cif, sdf (optional), and mol2 files for components in the chemical component search index.

        Args:
            ccUrlTarget (str): locator for source chemical component dictionary (default: full public dictionary)
            birdUrlTarget (str): locator for source BIRD dictionary (default: full public dictionary)
            limitPerceptions (bool): restrict automatic perceptions in OE molecular build operations (default: False)
            numProc (int): number of processors
            useCache (bool): use existing resource file where possible (default: True)
            molLimit (str):  limit the number to ingested chemical compont (default: None)
            quietFlag (bool): suppress output in OE library operations (default: True)

        Returns:
            (int): number molfiles generated
        """
        cachePath = self.__cachePath
        ccUrlTarget = kwargs.get("ccUrlTarget", None)
        birdUrlTarget = kwargs.get("birdUrlTarget", None)
        molLimit = kwargs.get("molLimit", None)
        quietFlag = kwargs.get("quietFlag", True)
        fpTypeList = kwargs.get("fpTypeList", [])
        screenTypeList = kwargs.get("screenTypeList", [])
        ccFileNamePrefix = "cc-%s" % self.__prefix if self.__prefix else "cc"
        oeFileNamePrefix = "oe-%s" % self.__prefix if self.__prefix else "oe"
        numProc = kwargs.get("numProc", 2)
        minCount = kwargs.get("minCount", 0)
        useCache = kwargs.get("useCache", True)
        useSdf = kwargs.get("useSdf", True)
        useMol2 = kwargs.get("useSdf", False)
        limitPerceptions = kwargs.get("limitPerceptions", False)
        logSizes = False
        #
        startTime = time.time()
        ccmP = ChemCompMoleculeProvider(
            cachePath=cachePath, useCache=useCache, ccFileNamePrefix=ccFileNamePrefix, ccUrlTarget=ccUrlTarget, birdUrlTarget=birdUrlTarget, molLimit=molLimit
        )
        ok = ccmP.testCache(minCount=minCount, logSizes=logSizes)
        logger.info("Completed chemical component provider load %r (%.4f seconds)", ok, time.time() - startTime)
        #
        startTime = time.time()
        oesmp = OeSearchMoleculeProvider(
            ccUrlTarget=ccUrlTarget,
            birdUrlTarget=birdUrlTarget,
            cachePath=cachePath,
            ccFileNamePrefix=ccFileNamePrefix,
            oeFileNamePrefix=oeFileNamePrefix,
            useCache=useCache,
            quietFlag=quietFlag,
            fpTypeList=fpTypeList,
            screenTypeList=screenTypeList,
            numProc=numProc,
            molLimit=molLimit,
            limitPerceptions=limitPerceptions,
        )
        ok = oesmp.testCache()
        logger.info("Completed OE molecule provider load %r (%.4f seconds)", ok, time.time() - startTime)
        #
        startTime = time.time()
        ccSIdxP = ChemCompSearchIndexProvider(cachePath=cachePath, useCache=useCache, ccFileNamePrefix=ccFileNamePrefix, limitPerceptions=limitPerceptions, numProc=numProc)
        ok = ccSIdxP.testCache()
        logger.info("Completed chemical component search index load %r (%.4f seconds)", ok, time.time() - startTime)
        #
        ccSIdx = ccSIdxP.getIndex() if ccSIdxP and ok else {}
        logger.info("Search index status %r index length %d", ok, len(ccSIdx))
        #
        ccIdD = {}
        mU = MarshalUtil()
        oeU = OeIoUtils(dirPath=cachePath)
        numMols = 0
        searchFileDirPath = self.getSearchDirFilePath()
        for sId in ccSIdx:
            ccId = sId.split("|")[0]
            # standard CIF definition
            if ccId not in ccIdD:
                cifPath = os.path.join(searchFileDirPath, ccId[0], ccId, ccId + ".cif")
                if not (useCache and mU.exists(cifPath)):
                    ccMol = ccmP.getMol(ccId)
                    mU.doExport(cifPath, [ccMol], fmt="mmcif")
            #
            oeMol = oesmp.getMol(sId)
            cifPath = os.path.join(searchFileDirPath, ccId[0], ccId, sId + ".cif")
            if sId != ccId and not (useCache and mU.exists(cifPath)):
                oeccU = OeChemCompUtils()
                ok = oeccU.addOeMol(sId, oeMol, missingModelXyz=True, writeIdealXyz=False)
                if ok:
                    oeccU.write(cifPath)

            if useSdf:
                molFilePath = os.path.join(searchFileDirPath, ccId[0], ccId, sId + ".sdf")
                if not (useCache and mU.exists(molFilePath)):
                    oeU.write(molFilePath, oeMol, constantMol=False, addSdTags=True)
            #
            if useMol2:
                mol2FilePath = os.path.join(searchFileDirPath, ccId[0], ccId, sId + ".mol2")
                if not (useCache and mU.exists(mol2FilePath)):
                    oeU.write(mol2FilePath, oeMol, constantMol=False, addSdTags=True)
            numMols += 1
            #
        return numMols

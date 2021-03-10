##
# File:    CODModelSearch.py
# Author:  J. Westbrook
# Date:    5-Mar-2021
# Version: 0.001
#
# Updated:
##
"""
Search the Crystallographic Open Database (COD) for matching chemical component models.

"""
__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "john.westbrook@rcsb.org"
__license__ = "Apache 2.0"

import logging
import time
import os

from rcsb.utils.chem.BatchChemSearch import BatchChemSearch
from rcsb.utils.io.MarshalUtil import MarshalUtil
from rcsb.utils.io.FileUtil import FileUtil
from rcsb.utils.multiproc.MultiProcUtil import MultiProcUtil

logger = logging.getLogger(__name__)


class CODModelSearch(object):
    def __init__(self, cachePath, **kwargs):
        self.__cachePath = cachePath
        #
        self.__useCache = kwargs.get("useCache", True)
        self.__ccUrlTarget = kwargs.get("ccUrlTarget", None)
        self.__birdUrlTarget = kwargs.get("birdUrlTarget", None)
        self.__descriptorUrlTarget = kwargs.get("descriptorUrlTarget", "http://www.crystallography.net/cod/smi/allcod.smi")
        self.__prefix = kwargs.get("prefix", None)
        self.__numProc = kwargs.get("numProc", 4)
        self.__chunkSize = kwargs.get("chunkSize", 50)
        self.__ccFileNamePrefix = "cc-%s" % self.__prefix if self.__prefix else "cc-full"
        self.__fU = FileUtil()
        # self.__ccmG = ChemCompModelGen(self.__cachePath, self.__prefix)

    def getResultIndex(self):
        mU = MarshalUtil(workPath=self.__cachePath)
        cD = mU.doImport(self.getResultFilePath(), fmt="json")
        return cD

    def getResultDetails(self, codId):
        mU = MarshalUtil(workPath=self.__cachePath)
        dD = mU.doImport(self.__getCodDetailsFilePath(codId), fmt="json")
        return dD

    def storeResultIndex(self, cD):
        mU = MarshalUtil(workPath=self.__cachePath)
        ok = mU.doExport(self.getResultFilePath(), cD, fmt="json", indent=3)
        return ok

    def getResultDirFilePath(self):
        dN = "cod-%s-result-files" % self.__prefix if self.__prefix else "cod-result-files"
        return os.path.join(self.__cachePath, dN)

    def getRawResultFilePath(self):
        dN = "cod-%s-result-files" % self.__prefix if self.__prefix else "cod-search-files"
        return os.path.join(self.__cachePath, dN, "cod-raw-result-file-index.json")

    def getResultFilePath(self):
        dN = "cod-%s-result-files" % self.__prefix if self.__prefix else "cod-search-files"
        return os.path.join(self.__cachePath, dN, "cod-result-file-index.json")

    def getDescriptorPath(self):
        fn = self.__fU.getFileName(self.__descriptorUrlTarget)
        dirPath = self.getResultDirFilePath()
        filePath = os.path.join(dirPath, fn)
        return filePath

    def updateDescriptors(self):
        self.__fetchUrl(self.__descriptorUrlTarget, filePath=self.getDescriptorPath(), useCache=False)

    def __fetchUrl(self, urlTarget, filePath, useCache=False, noRetry=False):
        ok = False
        try:
            if not (useCache and self.__fU.exists(filePath)):
                startTime = time.time()
                ok = self.__fU.get(urlTarget, filePath, noRetry=noRetry)
                endTime = time.time()
                if ok:
                    logger.debug("Fetched %s for resource file %s (status = %r) (%.4f seconds)", urlTarget, filePath, ok, endTime - startTime)
                else:
                    logger.error("Failing fetch for %s for resource file %s (status = %r) (%.4f seconds)", urlTarget, filePath, ok, endTime - startTime)
            else:
                ok = True
                logger.debug("Using cached data for %s", urlTarget)
            #
        except Exception as e:
            logger.exception("Failing for %r with %s", urlTarget, str(e))
        return ok

    def search(self, molLimit=None):
        try:
            bsw = BatchChemSearch(
                useCache=self.__useCache,
                ccUrlTarget=self.__ccUrlTarget,
                birdUrlTarget=self.__birdUrlTarget,
                ccFileNamePrefix=self.__ccFileNamePrefix,
                cachePath=self.__cachePath,
                numProc=self.__numProc,
                chunkSize=self.__chunkSize,
            )
            smiPath = self.getDescriptorPath()
            smiL = bsw.fetchDescriptorList(smiPath, swap=True)
            logger.info("Query length (%d)", len(smiL))
            #
            smiL = bsw.splitSmiles(smiL)
            retL = bsw.doQuery(smiL[:molLimit], "SMILES", matchOpts="graph-exact")
            logger.info("Result length (%d)", len(retL))
            #
            for ii, ret in enumerate(retL, 1):
                logger.debug("%5d %8s %4s (%.3f) %s: %s", ii, ret.queryId, ret.ccId, ret.fpScore, ret.queryType, ret.query)
            #
            fp = self.getRawResultFilePath()
            ok = bsw.storeMatchList(fp, retL)
            return len(retL) if ok else 0
        except Exception as e:
            logger.exception("Failing with %s", str(e))

    def __getSearchResults(self):
        """Read search results and convert to a chemical component dictionary."""
        fp = self.getRawResultFilePath()
        mU = MarshalUtil(workPath=self.__cachePath)
        rawL = mU.doImport(fp, fmt="json")
        rD = {}
        for cD in rawL:
            rD.setdefault(cD["ccId"], []).append(cD)
        return rD

    def __getCodEntryUrl(self, codId):
        # Template Examples:
        # https://molecules.crystallography.net/cod/sdf/1/00/00/1000098.sdf
        # https://molecules.crystallography.net/cod/sdf/6/00/05/6000557.sdf
        #
        baseUrl = "https://molecules.crystallography.net/cod/sdf"
        url = os.path.join(baseUrl, codId[0:1], codId[1:3], codId[3:5], codId + ".sdf")
        return url

    def __getCodDetailsUrl(self, codId):
        baseUrl = "http://www.crystallography.net/cod/optimade/structures"
        url = os.path.join(baseUrl, codId)
        return url

    def __getCodDetailsFilePath(self, codId):
        dirPath = self.getResultDirFilePath()
        fp = os.path.join(dirPath, "cod-data", codId[0:1], codId[1:3], codId[3:5], codId + ".json")
        return fp

    def __getCodEntryFilePath(self, codId):
        dirPath = self.getResultDirFilePath()
        fp = os.path.join(dirPath, "cod-data", codId[0:1], codId[1:3], codId[3:5], codId + ".sdf")
        return fp

    def fetchMatchedData(self, useCache=True):
        """Fetch COD matched entries and metadata and update the raw search index with essential COD data attrbutes.

        Args:
            useCache (bool, optional): use any cached COD data. Defaults to True.

        Returns:
            int: search result count

        """
        eCount = 0
        eSkip = 0
        rcD = {}
        cD = self.__getSearchResults()
        #
        for ccId, qDL in cD.items():
            # cifPath = self.__ccmG.getChemCompPath(ccId)
            # if not cifPath:
            #    logger.info("No CIF for %s skipping", ccId)
            #    continue
            parentId = ccId.split("|")[0]
            rqDL = []
            for qD in qDL:
                codId = qD["queryId"]
                codEntryFilePath = self.__getCodEntryFilePath(codId)
                codDetailsFilePath = self.__getCodDetailsFilePath(codId)
                ok1 = self.__fetchUrl(self.__getCodEntryUrl(codId), self.__getCodEntryFilePath(codId), useCache=useCache, noRetry=True)
                ok2 = self.__fetchUrl(self.__getCodDetailsUrl(codId), self.__getCodDetailsFilePath(codId), useCache=useCache, noRetry=True)
                tD = self.getResultDetails(codId)
                dD = tD["data"]["attributes"] if "data" in tD and "attributes" in tD["data"] else {}
                mD = tD["meta"]["implementation"] if "meta" in tD and "implementation" in tD["meta"] else {}
                if ok1 & ok2:
                    logger.info("Fetched COD entry and details for %s (%r)", codId, ok1 & ok2)
                    eCount += 1
                    qD["codEntryFilePath"] = codEntryFilePath
                    qD["codDetailsFilePath"] = codDetailsFilePath
                    # qD["cifPath"] = cifPath
                    qD["parentId"] = parentId
                    qD["chemicalName"] = dD["_cod_commonname"] if "_cod_commonname" in dD else None
                    qD["chemicalName"] = dD["_cod_chemname"] if "_cod_chemname" in dD else qD["chemicalName"]
                    qD["rValue"] = dD["_cod_Robs"] if "_cod_Robs" in dD else None
                    qD["diffrnTemp"] = dD["_cod_diffrtemp"] if "_cod_diffrtemp" in dD else None
                    qD["radiationSource"] = dD["_cod_radType"] if "_cod_radType" in dD else None
                    qD["publicationDOI"] = dD["_cod_doi"] if "_cod_doi" in dD else None
                    qD["version"] = mD["version"] if "version" in mD else None
                    qD["hasDisorder"] = "N"
                    rqDL.append(qD)
                else:
                    logger.info("Skipping entry missing data for %r at %r", codId, self.__getCodEntryUrl(codId))
                    eSkip += 1
            if rqDL:
                rcD[ccId] = rqDL
        #
        ok = self.storeResultIndex(rcD)
        logger.info("Final match result (w/sdf and metadata) (%d/%d) cod hits (%d) skipped (%d)", len(rcD), len(cD), eCount, eSkip)
        return eCount if ok else 0

    def fetchMatchedDataMp(self, numProc=6, chunkSize=5, useCache=True):
        rcD = {}
        cD = self.__getSearchResults()
        idList = list(cD.keys())
        # ---
        mpu = MultiProcUtil(verbose=True)
        mpu.setWorkingDir(self.__cachePath)
        mpu.setOptions(optionsD={"resultPath": self.__cachePath, "cD": cD, "useCache": useCache})
        mpu.set(workerObj=self, workerMethod="fetchDataWorker")

        ok, failList, resultList, _ = mpu.runMulti(dataList=idList, numProc=numProc, numResults=1, chunkSize=chunkSize)
        logger.info("Run ended with status %r success count %d failures %r", ok, len(resultList[0]), len(failList))
        for rTup in resultList[0]:
            rcD[rTup[0]] = rTup[1]
        # ---
        ok = self.storeResultIndex(rcD)
        logger.info("Final match result (w/sdf and metadata) (%d/%d)", len(rcD), len(cD))
        return True

    def fetchDataWorker(self, dataList, procName, optionsD, workingDir):
        """Worker method to fetch COD data for matched entries

        Args:
            dataList (list): list of mol2 file paths to be searched
            procName (str): processName
            optionsD (dict): dictionary of options
            workingDir (str): path to working directory (not used)

        Returns:
            (successList, resultList, []): success and result lists of mol2 paths with CCDC matches
        """
        resultPath = optionsD["resultPath"]
        cD = optionsD["cD"]
        useCache = optionsD["useCache"]
        _ = workingDir
        resultList = []
        successList = []
        startTime = time.time()
        logger.info("starting %s at %s", procName, time.strftime("%Y %m %d %H:%M:%S", time.localtime()))
        #
        eCount = 0
        eSkip = 0
        try:
            stopPath = os.path.join(resultPath, "STOP")
            logger.info("%s starting search data length %d", procName, len(dataList))
            if self.__checkStop(stopPath):
                logger.info("%s stopping", procName)
                return resultList, resultList, []
            #
            # for ccId, qDL in cD.items():
            for ccId in dataList:
                if ccId in cD:
                    qDL = cD[ccId]
                #
                parentId = ccId.split("|")[0]
                rqDL = []
                for qD in qDL:
                    codId = qD["queryId"]
                    codEntryFilePath = self.__getCodEntryFilePath(codId)
                    codDetailsFilePath = self.__getCodDetailsFilePath(codId)
                    ok1 = self.__fetchUrl(self.__getCodEntryUrl(codId), self.__getCodEntryFilePath(codId), useCache=useCache, noRetry=True)
                    ok2 = self.__fetchUrl(self.__getCodDetailsUrl(codId), self.__getCodDetailsFilePath(codId), useCache=useCache, noRetry=True)
                    tD = self.getResultDetails(codId)
                    dD = tD["data"]["attributes"] if "data" in tD and "attributes" in tD["data"] else {}
                    mD = tD["meta"]["implementation"] if "meta" in tD and "implementation" in tD["meta"] else {}
                    if ok1 & ok2:
                        logger.info("Fetched COD entry and details for %s (%r)", codId, ok1 & ok2)
                        eCount += 1
                        qD["codEntryFilePath"] = codEntryFilePath
                        qD["codDetailsFilePath"] = codDetailsFilePath
                        # qD["cifPath"] = cifPath
                        qD["parentId"] = parentId
                        qD["chemicalName"] = dD["_cod_commonname"] if "_cod_commonname" in dD else None
                        qD["chemicalName"] = dD["_cod_chemname"] if "_cod_chemname" in dD else qD["chemicalName"]
                        qD["rValue"] = dD["_cod_Robs"] if "_cod_Robs" in dD else None
                        qD["diffrnTemp"] = dD["_cod_diffrtemp"] if "_cod_diffrtemp" in dD else None
                        qD["radiationSource"] = dD["_cod_radType"] if "_cod_radType" in dD else None
                        qD["publicationDOI"] = dD["_cod_doi"] if "_cod_doi" in dD else None
                        qD["version"] = mD["version"] if "version" in mD else None
                        qD["hasDisorder"] = "N"
                        rqDL.append(qD)
                    else:
                        logger.info("Skipping entry missing data for %r at %r", codId, self.__getCodEntryUrl(codId))
                        eSkip += 1
                if rqDL:
                    resultList.append((ccId, rqDL))
                    successList.append(ccId)
        except Exception as e:
            logger.exception("Failing with %s", str(e))

        endTime = time.time()
        logger.info(
            "%s (entries %d skipped %d) (ccId result length %d) completed at %s (%.2f seconds)",
            procName,
            eCount,
            eSkip,
            len(successList),
            time.strftime("%Y %m %d %H:%M:%S", time.localtime()),
            endTime - startTime,
        )
        return successList, resultList, []

    def __checkStop(self, path):
        try:
            if os.access(path, os.F_OK):
                return True
        except Exception:
            pass
        return False

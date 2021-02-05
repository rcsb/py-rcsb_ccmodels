##
# File:    ChemCompModelAssemble.py
# Author:  J. Westbrook
# Date:     4-Feb-2021
# Version: 0.001
#
# Updated:
##
"""
Assemble model files and adjust audit records for the chemical component model workflow.

"""
__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "john.westbrook@rcsb.org"
__license__ = "Apache 2.0"

import datetime
import logging
import os
import time
from collections import defaultdict
from operator import itemgetter

from rcsb.ccmodels.search.ChemCompModelBuild import ChemCompModelBuild
from rcsb.ccmodels.search import __version__
from rcsb.utils.chemref.ChemCompModelProvider import ChemCompModelProvider
from rcsb.utils.io.MarshalUtil import MarshalUtil

logger = logging.getLogger(__name__)


class ChemCompModelAssemble(object):
    def __init__(self, cachePath, prefix=None):
        self.__cachePath = cachePath
        self.__prefix = prefix
        #
        self.__ccmb = ChemCompModelBuild(cachePath=self.__cachePath, prefix=self.__prefix)
        # self.__startTime = time.time()
        logger.info("Starting assemble (%s) at %s", __version__, time.strftime("%Y %m %d %H:%M:%S", time.localtime()))

    def assemble(self, maxRFactor=10.0):
        """Concatenate models into the input file path subject to the R value constraint.
        Relabel the models sequentially for each parent chemical component.

        Args:
            assembleModelPath (str): path for concatenated model file
            maxRFactor (float, optional): limiting R-value. Defaults to 10.0.

        Returns:
            (bool): True for success or False otherwise

        """
        dataContainerL = []
        mU = MarshalUtil(workPath=self.__cachePath)
        modelIndexD = self.__ccmb.fetchModelIndex()
        modelIndexD = self.__addPriorMatchDetails(modelIndexD)
        modelIndexD = self.__updateVariantDetails(modelIndexD)
        priorMapD = {}
        for _, mDL in modelIndexD.items():
            mDLS = sorted(mDL, key=itemgetter("priorModelId", "variantType", "rFactor"), reverse=False)

            numStd = 0
            for mD in mDLS:
                isStd = False
                if mD["variantType"].startswith("A"):
                    numStd += 1
                    isStd = True
                #
                if mD["rFactor"] > maxRFactor:
                    logger.info("Skipping model %s isStd (%r) rValue (%r)", mD["modelId"], isStd, mD["rFactor"])
                    continue
                if numStd and not isStd:
                    logger.info("Skipping model %s isStd (%r) numStd (%d)", mD["modelId"], isStd, numStd)
                    continue

                cL = mU.doImport(mD["modelPath"], fmt="mmcif")
                logger.debug("Read %d from %s", len(cL), mD["modelPath"])
                dataContainerL.extend(cL)
                if not mD["priorModelId"].startswith("Z"):
                    priorMapD[mD["modelId"]] = (mD["priorModelId"], mD["priorMatchDate"])
        #
        fn = "chem_comp_models-%s.cif" % self.__getToday()
        assembleModelPath = os.path.join(self.__ccmb.getModelDirFilePath(), fn)
        # -- relabel
        parentModelCountD = defaultdict(int)
        priorIdLD = {}
        for dataContainer in dataContainerL:
            tModelId = dataContainer.getName()
            tId = tModelId.split("_")[1]
            pId = tId.split("|")[0]
            # tCount = tModelId.split("_")[2]
            if tModelId in priorMapD:
                pCount = priorMapD[tModelId][0].split("_")[2]
                priorIdLD.setdefault(pId, []).append(int(pCount))
                self.__replaceModelId(dataContainer, tModelId, priorMapD[tModelId][0])
                self.__updateAuditDate(dataContainer, priorMapD[tModelId][1])
                parentModelCountD[pId] = sorted(priorIdLD[pId])[-1]
                logger.info("%s current model %r prior model %r count %d", pId, tModelId, priorMapD[tModelId][0], parentModelCountD[pId])
            else:
                parentModelCountD[pId] += 1
                pModelId = self.__makePublicModelId(pId, parentModelCountD[pId])
                self.__replaceModelId(dataContainer, tModelId, pModelId)

        ok = mU.doExport(assembleModelPath, dataContainerL, fmt="mmcif")
        logger.info("Assembled %d models status %r", len(dataContainerL), ok)
        self.__checkAssembledModels(assembleModelPath)
        return len(dataContainerL)

    def __getAuditDetails(self):
        """[summary]

        Returns:
            [type]: [description]

             aL.append({"audit_date": auditDate, "action_type": auditAction})
             rD.setdefault(ccId, []).append({"model_id": modelId, "db_name": dbName, "db_code": dbCode, "audit_list": aL})
        """
        ccm = ChemCompModelProvider(cachePath=self.__cachePath, useCache=False)
        rD = ccm.getAuditDetails()
        return rD

    def __addPriorMatchDetails(self, modelIndexD):
        """"""
        priorMatchLD = self.__getAuditDetails()
        for pId, mDL in modelIndexD.items():
            if pId in priorMatchLD:
                for priorMatch in priorMatchLD[pId]:
                    priorModelId = priorMatch["model_id"]
                    priorDbName = priorMatch["db_name"]
                    priorDbCode = priorMatch["db_code"]
                    priorDate = priorMatch["audit_list"][-1]["audit_date"]
                    # compare with current matches
                    numMatch = 0
                    for mD in mDL:
                        curDbName = "CSD"
                        # curDbName = mD["matchDb"]
                        curDbCode = mD["matchId"]
                        if priorDbName == curDbName and curDbCode == priorDbCode:
                            numMatch += 1
                            mD["priorModelId"] = priorModelId
                            mD["priorMatchDate"] = priorDate
                        else:
                            mD["priorModelId"] = "Znone"
                    if numMatch:
                        logger.info("%s has prior matches (%d)", pId, numMatch)
            else:
                for mD in mDL:
                    mD["priorModelId"] = "Znone"
        return modelIndexD

    def __updateVariantDetails(self, modelIndexD):
        """"""
        for _, mDL in modelIndexD.items():
            for mD in mDL:
                if not mD["variantType"]:
                    mD["variantType"] = "Anone"
        return modelIndexD

    def __makePublicModelId(self, parentId, modelNum=1):
        modelId = "M_" + parentId + "_%05d" % modelNum
        return modelId

    #  deprecated version
    def __assembleX(self, maxRFactor=10.0):
        """Concatenate models into the input file path subject to the R value constraint.
        Relabel the models sequentially for each parent chemical component.

        Args:
            assembleModelPath (str): path for concatenated model file
            maxRFactor (float, optional): limiting R-value. Defaults to 10.0.

        Returns:
            (bool): True for success or False otherwise
        """
        dataContainerL = []
        mU = MarshalUtil(workPath=self.__cachePath)
        modelIndexD = self.__ccmb.fetchModelIndex()
        for _, mDL in modelIndexD.items():
            mDLS = sorted(mDL, key=itemgetter("rFactor"), reverse=False)
            # check for non variant matches
            modelIdSelectL = []
            for mD in mDLS:
                if not mD["variantType"]:
                    modelIdSelectL.append(mD["modelId"])
            if modelIdSelectL:
                # cannonical models first -
                for mD in mDLS:
                    if mD["modelId"] in modelIdSelectL:
                        if mD["rFactor"] < maxRFactor:
                            cL = mU.doImport(mD["modelPath"], fmt="mmcif")
                            logger.debug("Read %d from %s", len(cL), mD["modelPath"])
                            dataContainerL.extend(cL)
                        else:
                            logger.info("Rejected cannonical model %s", mD["modelId"])
            else:
                # then tautomer/protomer models
                for mD in mDLS:
                    if mD["rFactor"] < maxRFactor:
                        cL = mU.doImport(mD["modelPath"], fmt="mmcif")
                        logger.debug("Read %d from %s", len(cL), mD["modelPath"])
                        dataContainerL.extend(cL)
                    else:
                        logger.info("Rejected variant model %s", mD["modelId"])
        #
        fn = "chem_comp_models-%s.cif" % self.__getToday()
        assembleModelPath = os.path.join(self.__ccmb.getModelDirFilePath(), fn)
        # -- relabel
        parentModelCountD = defaultdict(int)
        for dataContainer in dataContainerL:
            tModelId = dataContainer.getName()
            tId = tModelId.split("_")[1]
            pId = tId.split("|")[0]
            parentModelCountD[pId] += 1
            pModelId = self.__makePublicModelId(pId, parentModelCountD[pId])
            self.__replaceModelId(dataContainer, tModelId, pModelId)
        # -- relabel
        ok = mU.doExport(assembleModelPath, dataContainerL, fmt="mmcif")
        logger.info("Assembled %d models status %r", len(dataContainerL), ok)
        self.__checkAssembledModels(assembleModelPath)
        return len(dataContainerL)

    def __checkAssembledModels(self, assembleModelPath):
        catNameL = [
            "pdbx_chem_comp_model",
            "pdbx_chem_comp_model_atom",
            "pdbx_chem_comp_model_bond",
            "pdbx_chem_comp_model_descriptor",
            "pdbx_chem_comp_model_reference",
            "pdbx_chem_comp_model_feature",
            "pdbx_chem_comp_model_audit",
        ]
        mU = MarshalUtil(workPath=self.__cachePath)
        dataContainerL = mU.doImport(assembleModelPath, fmt="mmcif")
        logger.info("Read %d data containers", len(dataContainerL))
        rD = {}
        cnD = {}
        for dataContainer in dataContainerL:
            nm = dataContainer.getName()
            logger.debug("datacontainer %r", nm)
            if nm in cnD:
                logger.info("Duplicate container id %r", nm)
                cnD[nm] = True
            #
            pId = nm.split("_")[1]
            cObj = dataContainer.getObj("pdbx_chem_comp_model")
            modelId = cObj.getValue("id", 0)
            if modelId != nm:
                logger.error("modelId %r datablock %r", modelId, nm)
            #
            tD = {}
            for catName in catNameL:
                cObj = dataContainer.getObj(catName)
                nRows = cObj.getRowCount()
                tD[catName] = nRows
            cObj = dataContainer.getObj("pdbx_chem_comp_model_feature")
            skip = False
            for ii in range(cObj.getRowCount()):
                fN = cObj.getValue("feature_name", ii)
                fV = cObj.getValue("feature_value", ii)
                if fN == "heavy_atoms_only" and fV == "Y":
                    skip = True
                    break
            if not skip:
                rD.setdefault(pId, []).append(tD)
        #
        for pId, tDL in rD.items():
            for catName in catNameL:
                minV = 100000
                maxV = -1
                for tD in tDL:
                    minV = min(minV, tD[catName])
                    maxV = max(maxV, tD[catName])
                if maxV - minV > 1 and catName not in ["pdbx_chem_comp_model_feature"]:
                    logger.error("%s %s row count inconsistency %d %d", pId, catName, minV, maxV)

    def __replaceModelId(self, dataContainer, oldModelId, newModelId):
        """Update all instances of the model in the input container with the input
        replacement value.

        Args:
            dataContainerList (obj): input container list
            newModelId (str):  replacement modelId value

        """
        updateL = [
            ("pdbx_chem_comp_model", "id"),
            ("pdbx_chem_comp_model_atom", "model_id"),
            ("pdbx_chem_comp_model_bond", "model_id"),
            ("pdbx_chem_comp_model_descriptor", "model_id"),
            ("pdbx_chem_comp_model_reference", "model_id"),
            ("pdbx_chem_comp_model_feature", "model_id"),
            ("pdbx_chem_comp_model_audit", "model_id"),
        ]
        try:
            dataContainer.setName(newModelId)
            for tup in updateL:
                cObj = dataContainer.getObj(tup[0])
                cObj.replaceValue(oldModelId, newModelId, tup[1])
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return dataContainer

    def __updateAuditDate(self, dataContainer, auditDate):
        """Update the audit date for input container

        Args:
            dataContainerList (obj): input container list
            auditDate (str):  original audit date

        """
        try:
            cObj = dataContainer.getObj("pdbx_chem_comp_model_audit")
            ii = cObj.getRowCount()
            cObj.setValue(auditDate, "date", ii - 1)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return dataContainer

    def __getToday(self):
        """Return a CIF style date-timestamp value for current local time -"""
        today = datetime.datetime.today()
        # format ="%Y-%m-%d:%H:%M"
        fmt = "%Y-%m-%d"
        return str(today.strftime(fmt))

    def __makePublicModelId(self, parentId, modelNum=1):
        modelId = "M_" + parentId + "_%05d" % modelNum
        return modelId

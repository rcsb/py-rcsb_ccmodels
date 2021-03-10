##
# File:    CODModelBuild.py
# Author:  J. Westbrook
# Date:    6-Mar-2021
# Version: 0.001
#
# Updated:
##
"""
Workflow to select and build model files from matched COD entry data.

"""
__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "john.westbrook@rcsb.org"
__license__ = "Apache 2.0"

import copy
import datetime
import logging
import os
import time
from collections import namedtuple, defaultdict

from mmcif.api.DataCategory import DataCategory
from mmcif.api.PdbxContainers import DataContainer
from rcsb.ccmodels.search import __version__
from rcsb.ccmodels.search.CODModelSearch import CODModelSearch
from rcsb.utils.chem.ChemCompMoleculeProvider import ChemCompMoleculeProvider
from rcsb.utils.chem.ChemCompSearchIndexProvider import ChemCompSearchIndexProvider
from rcsb.utils.chem.OeAlignUtils import OeAlignUtils
from rcsb.utils.chem.OeChemCompUtils import OeChemCompUtils
from rcsb.utils.chem.OeDepictAlign import OeDepictSubStructureAlignMultiPage
from rcsb.utils.chem.OeDepictAlign import OeDepictMCSAlignPage
from rcsb.utils.chem.OeSearchMoleculeProvider import OeSearchMoleculeProvider
from rcsb.utils.io.decorators import timeout
from rcsb.utils.io.MarshalUtil import MarshalUtil
from rcsb.utils.multiproc.MultiProcUtil import MultiProcUtil

logger = logging.getLogger(__name__)

ComponentAtomDetails = namedtuple("ComponentAtomDetails", "atIdx atNo atName atType x y z atFormalCharge")
AlignAtomMap = namedtuple("AlignAtomMap", "refId refAtIdx refAtNo refAtName fitId fitAtIdx fitAtNo fitAtName")
AlignAtomUnMapped = namedtuple("AlignAtomUnMapped", "fitId fitAtIdx fitAtNo fitAtType fitAtName fitAtFormalCharge x y z fitNeighbors")


class CODModelBuildWorker(object):
    def __init__(self, cachePath, verbose=True):
        self.__cachePath = cachePath
        self.__verbose = verbose
        self.__ccSIdxP = None
        self.__oesmP = None
        self.__ccmP = None

    def build(self, dataList, procName, optionsD, workingDir):
        """Worker method to build chemical component models for the input search result ID list.

        Args:
            dataList (list): list of search result paths
            procName (str): processName
            optionsD (dict): dictionary of options
            workingDir (str): path to working directory (not used)

        Returns:
            (successList, resultList, []): success and result lists
        """
        _ = workingDir
        resultList = []
        successList = []
        #
        try:
            alignType = optionsD["alignType"]
            modelDirPath = optionsD["modelDirPath"]
            imageDirPath = optionsD["imageDirPath"]
            self.__ccSIdxP = optionsD["ccSIdxP"]
            self.__oesmP = optionsD["oesmP"]
            self.__ccmP = optionsD["ccmP"]
            idxIdD = optionsD["idxIdD"]
            #
            startTime = time.time()
            logger.info("%s ======== ============ ============ ", procName)
            logger.info("%s starting with length %d at %s", procName, len(dataList), time.strftime("%Y %m %d %H:%M:%S", time.localtime()))
            #
            parentD = {}
            parentModelCountD = defaultdict(int)
            for idxId in dataList:
                logger.info("%s start model build for parent Id (%s)", procName, idxId)
                #
                pairList = []
                matchDList = idxIdD[idxId]
                #
                for matchD in matchDList:
                    fitMolFilePath = matchD["codEntryFilePath"]
                    matchId = matchD["queryId"]
                    targetId = matchD["ccId"]
                    targetObj = self.__ccmP.getMol(targetId)
                    #
                    if not targetObj:
                        oeTargetObj = self.__oesmP.getMol(targetId)
                        logger.info("Starting building tautomer (%d) %r", oeTargetObj.NumAtoms(), targetId)
                        oeccU = OeChemCompUtils()
                        ok = oeccU.addOeMol(targetId, oeTargetObj, missingModelXyz=True, writeIdealXyz=False)
                        cL = oeccU.getContainerList()
                        targetObj = cL[0]
                        if not targetObj:
                            logger.error("%r has null object", targetId)
                            continue
                        logger.info("Completed building tautomer %r", targetId)
                    #
                    # matchTitle = "COD Code  " + matchId
                    # ccTitle = "Chemical Component " + targetId
                    parentId = matchD["parentId"]
                    sId = targetId
                    #
                    try:
                        nAtomsRef, refFD, nAtomsFit, fitFD, fitXyzMapD, fitAtomUnMappedL, isSkipped = self.__alignModelSubStruct(
                            targetObj,
                            fitMolFilePath,
                            alignType=alignType,
                            fitTitle=matchId,
                            refTitle=targetId,
                            onlyCloseMatches=False,
                            verbose=self.__verbose,
                            procName=procName,
                        )
                    except Exception:
                        nAtomsRef, refFD, nAtomsFit, fitFD, fitXyzMapD, fitAtomUnMappedL, isSkipped = 0, {}, 0, {}, {}, [], False
                    #
                    logger.debug(
                        ">>> %s - %s nAtomsRef %d nAtomsFit %d atommapL (%d) fitAtomUnMappedL (%d)", targetId, matchId, nAtomsRef, nAtomsFit, len(fitXyzMapD), len(fitAtomUnMappedL)
                    )
                    if nAtomsRef == 0 and nAtomsFit == 0:
                        logger.info("%s alignment fails for %s with %s", procName, targetId, matchId)
                        continue
                    if isSkipped:
                        logger.info("%s alignment skipped for %s (%d) with %s (%d)", procName, targetId, nAtomsRef, matchId, nAtomsFit)
                        continue

                    # -----
                    smilesMatch = refFD["SMILES_STEREO"] == fitFD["SMILES_STEREO"]
                    hasUnMapped = len(fitAtomUnMappedL) > 0
                    unMappedOk = self.__testUnMappedProtonation(fitAtomUnMappedL)
                    #
                    #   --- acceptance tests ---
                    acceptOk = False
                    if nAtomsRef == len(fitXyzMapD) and smilesMatch:
                        # Take all cases with matching SMILES and all atoms aligned
                        acceptOk = True
                    elif nAtomsRef > len(fitXyzMapD) and smilesMatch:
                        # Heavy atom only match
                        acceptOk = True
                    elif nAtomsRef == len(fitXyzMapD) and hasUnMapped and unMappedOk:
                        # All atoms aligned and fit includes a protonation ...
                        acceptOk = True
                    else:
                        pass

                    if len(fitXyzMapD) < 2:
                        acceptOk = False
                    #
                    if not acceptOk:
                        logger.info(
                            "%s rejecting (smilesMatch %r hasUnmapped %r (%d) unMappedOk %r nAtomsRef %d nAtomsFit %d mapped fit %d) for %s with %s",
                            procName,
                            smilesMatch,
                            hasUnMapped,
                            len(fitAtomUnMappedL),
                            unMappedOk,
                            nAtomsRef,
                            nAtomsFit,
                            len(fitXyzMapD),
                            targetId,
                            matchId,
                        )
                        continue
                    #
                    if hasUnMapped and not smilesMatch:
                        logger.info("%s SMILES for %s and %s differ with unmapped protons", procName, targetId, matchId)
                        logger.debug("Ref %-8s SMILES: %s", targetId, refFD["SMILES_STEREO"])
                        logger.debug("Fit %-8s SMILES: %s", matchId, fitFD["SMILES_STEREO"])

                    # --------- ----------------
                    #  Accept the match
                    # --------- ----------------
                    # matchId = matchD.getIdentifier()
                    # targetId = matchD.getTargetId()
                    parentD.setdefault(parentId, []).append(matchId)
                    #
                    # refImageFileName = "ref_" + targetId + "_" + matchId + ".svg"
                    # refImagePath = os.path.join(imageDirPath, sId, refImageFileName)
                    #
                    # self.__pairDepictPage(refImagePath, sId, ccTitle, refFD["OEMOL"], matchId, matchTitle, fitFD["OEMOL"], alignType=alignType)
                    # --------- ------------------
                    pairList.append((sId, refFD["OEMOL"], matchId, fitFD["OEMOL"]))
                    modelId, modelPath = self.__makeModelPath(modelDirPath, parentId, targetId, startingModelNum=parentModelCountD[parentId] + 1, maxModels=300, scanExisting=False)
                    logger.debug("targetId %r modelId %r modelPath %r", targetId, modelId, modelPath)
                    #
                    logger.info(
                        "%s accepted for %s with %s (smilesMatch %r hasUnmapped %r (%d) unMappedOk %r nAtomsRef %d nAtomsFit %d mapped fit %d) ",
                        procName,
                        targetId,
                        matchId,
                        smilesMatch,
                        hasUnMapped,
                        len(fitAtomUnMappedL),
                        unMappedOk,
                        nAtomsRef,
                        nAtomsFit,
                        len(fitXyzMapD),
                    )
                    #
                    ok, variantType = self.__writeModel(targetId, targetObj, fitFD, fitXyzMapD, fitAtomUnMappedL, matchD, modelId, modelPath)
                    if ok:
                        parentModelCountD[parentId] += 1
                        hd = matchD["hasDisorder"]
                        resultList.append(
                            {
                                "modelId": modelId,
                                "searchId": targetId,
                                "parentId": parentId,
                                "matchId": matchId,
                                "matchDB": "COD",
                                "modelPath": modelPath,
                                "rFactor": matchD["rValue"],
                                "hasDisorder": hd if hd else "N",
                                "variantType": variantType,
                            }
                        )
                #
                if pairList:
                    pdfImagePath = os.path.join(imageDirPath, sId, sId + "-all-pairs.pdf")
                    self.__depictFitList(sId, pdfImagePath, pairList, alignType=alignType)
                if resultList:
                    logger.info("%s built %d models for %s (this iteration)", procName, parentModelCountD[parentId], parentId)
                    successList.append(idxId)
                else:
                    logger.info("%s no models built for %s", procName, targetId)

            endTime = time.time()
            logger.info(
                "%s (match successes %d total models this iterations %d) completed at %s (%.2f seconds)",
                procName,
                len(successList),
                len(resultList),
                time.strftime("%Y %m %d %H:%M:%S", time.localtime()),
                endTime - startTime,
            )
            return successList, resultList, []
        except Exception as e:
            logger.exception("%s failing with %s", procName, str(e))
        return [], [], []

    def __getBuildVariant(self, sId):
        """Lookup the build type from the input chemical component search index Id.

        Args:
            sId (str): chemical component search index Id

        Returns:
            str : "tautomer_protomer" | ""
        """
        sidxD = self.__ccSIdxP.getIndexEntry(sId)
        if "protomer" in sidxD["build-type"]:
            return "tautomer_protomer"
        if "tautomer" in sidxD["build-type"]:
            return "tautomer_protomer"
        return ""

    @timeout(60)
    def __alignModelSubStruct(self, ccRefObj, molFitPath, alignType="strict", fitTitle=None, refTitle=None, onlyCloseMatches=False, verbose=False, procName="main"):
        """Align (substructure) chemical component definition search target with the candidate matching reference molecule.

        Args:
            ccRefPath (str): path to search target chemical component definition
            molFitPath (str): path to matched molfile
            alignType (str, optional):  strict|relaxed|relaxed-stereo. Defaults to "strict"
            fitTitle (str, optional): fit molecule title. Defaults to None.
            refTitle (str, optional): reference molecule title. Defaults to None.
            onlyCloseMatches (bool, optional): triage alignments for only close matches. Defaults to False
            verbose (bool, optional): enable verbose output. Defaults to False.
            procName (str, optional):  process name.  Defaults to main


        Returns:
            (tuple): (  number of atoms reference molecule,
                        reference molecule feature dictionary {"Formula": , "SMILES": , "SMILES_STEREO": , "InChI": , "InChIKey":,}
                        number of atoms in fit molecule,
                        fit molecule feature dictionary  {"Formula": , "SMILES": , "SMILES_STEREO": , "InChI": , "InChIKey":, 'xyz': }
                        atomMapD {refAtName: (fit x, fit y, fit z)}
                        fitAtomUnMappedL: [(atName, atomicNumber),(), ...],
                        isSkipped
                    )

            Note: 'xyz': [ComponentAtomDetails()]
        """
        try:
            refFD = fitFD = atomMapD = {}
            fitAtomUnMappedL = []
            isSkipped = False
            #
            logger.debug("Align target cc %s with matching model %s", ccRefObj.getName(), molFitPath)
            oesU = OeAlignUtils(workPath=self.__cachePath, verbose=verbose)
            oesU.setSearchType(sType=alignType)
            nAtomsRef = oesU.setRefObj(ccRefObj, title=refTitle)
            nAtomsFit = oesU.setFitPath(molFitPath, title=fitTitle, suppressHydrogens=False, fType="sdf", importType="3D", largestPart=True)
            if onlyCloseMatches and nAtomsFit > nAtomsRef and nAtomsFit - nAtomsRef > 2:
                isSkipped = True
                return nAtomsRef, refFD, nAtomsFit, fitFD, atomMapD, fitAtomUnMappedL, isSkipped
            (nAtomsRef, refFD, nAtomsFit, fitFD, atomMapL, fitAtomUnMappedL) = oesU.doAlignSs(unique=False)
            # -----
            # tD = {refAtName: fitAtIdx, ...}
            tD = {tup.refAtName: tup.fitAtIdx for tup in atomMapL}
            #
            # fitFD['xyz']: [ComponentAtomDetails()]
            # fitXyzD  = {atIdx1: ComponentAtomDetails(), atIdx2: ComponentAtomDetails(), ...}
            #
            fitXyzD = {tup.atIdx: tup for tup in fitFD["xyz"]}
            atomMapD = {refAtName: fitXyzD[fitAtIdx] for refAtName, fitAtIdx in tD.items() if fitAtIdx in fitXyzD}
            # -----
            return nAtomsRef, refFD, nAtomsFit, fitFD, atomMapD, fitAtomUnMappedL, isSkipped
        except Exception as e:
            logger.exception("%s failing for %r and %r with %s", procName, ccRefObj.getName(), molFitPath, str(e))
        return 0, {}, 0, {}, {}, [], False

    def __pairDepict(self, imagePath, refId, refTitle, refMol, fitId, fitTitle, fitMol, alignType="strict"):
        """Depict pairwise alignment of the input reference and fit molecules.

        Args:
            imagePath (str): path to image (format by path extension)
            refId (str): reference molecule identifier
            refTitle (str): reference molecule title
            refMol (obj): reference OE molecule object
            fitId (str): fit molecule identifier
            fitTitle (str): fit molecule title
            fitMol (obj): fit OE molecule object
            alignType (str, optional): alignment criteria (relaxed|relaxed-stereo|strict). Defaults to "strict".

        Returns:
            (list): atom mapping in all aligned figures
                    [(reference component Id, reference atom name, fit chemical component Id, fit atom name)
        """
        aML = []
        try:
            oed = OeDepictMCSAlignPage()
            oed.setSearchType(sType=alignType)

            oed.setRefMol(refMol, refId, title=refTitle)
            oed.setFitMol(fitMol, fitId, title=fitTitle)
            oed.setDisplayOptions(
                imageSizeX=2000,
                imageSizeY=1000,
                labelAtomName=True,
                labelAtomCIPStereo=True,
                labelAtomIndex=False,
                labelBondIndex=False,
                highlightStyleFit="ballAndStickInverse",
                bondDisplayWidth=0.5,
                highLightMatchColorRef="green",
                highLightNotMatchColorRef="pink",
            )
            aML = oed.alignPair(imagePath=imagePath)
            if aML:
                for (rCC, rAt, tCC, tAt) in aML:
                    logger.debug("%5s %-5s %5s %-5s", rCC, rAt, tCC, tAt)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return aML

    def __depictFitList(self, sId, pdfImagePath, pairList, alignType="exact"):
        """Depict pairwise alignments with multi-page layout in PDF format.

        Args:
            pdfImagePath (str): PDF image path
            pairList (list): [(refId, refOeMol, fitId, fitOeMol)]

        Returns:
            (list): atom mapping in all aligned figures
                    [(reference component Id, reference atom name, fit chemical component Id, fit atom name)
        """
        aML = []
        try:
            logger.debug("sId %s  pairList (%d)", sId, len(pairList))
            oed = OeDepictSubStructureAlignMultiPage()
            oed.setSearchType(sType=alignType)
            oed.setPairMolList(pairList)

            oed.setDisplayOptions(
                labelAtomName=True,
                labelAtomCIPStereo=True,
                labelAtomIndex=False,
                labelBondIndex=False,
                highlightStyleFit="ballAndStickInverse",
                pageOrientation="portrait",
                gridRows=4,
                bondDisplayWidth=0.5,
                highLightMatchColorRef="green",
                highLightNotMatchColorRef="pink",
            )
            aML = oed.alignPairListMulti(imagePath=pdfImagePath)
            logger.debug("%s atom map length (%d)", sId, len(aML))
            if aML:
                for (rCC, rAt, tCC, tAt) in aML:
                    logger.debug("%5s %-5s %5s %-5s", rCC, rAt, tCC, tAt)
            #
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return aML

    def __makeModelPath(self, modelDirPath, parentId, targetId, startingModelNum=1, maxModels=200, scanExisting=False):
        try:
            pth = None
            ii = startingModelNum
            dirPath = os.path.join(modelDirPath, parentId)
            if not os.access(dirPath, os.R_OK):
                os.makedirs(dirPath)
                modelId = self.__makeModelId(targetId, modelNum=ii)
                pth = os.path.join(dirPath, modelId + ".cif")
                return modelId, pth
            # optionally scan over any existing models before selecting a new model number.
            if scanExisting:
                while True:
                    modelId = self.__makeModelId(targetId, modelNum=ii)
                    pth = os.path.join(dirPath, modelId + ".cif")
                    if not os.access(pth, os.R_OK):
                        return modelId, pth
                    #
                    ii += 1
                    if ii > maxModels:
                        break
            else:
                modelId = self.__makeModelId(targetId, modelNum=ii)
                pth = os.path.join(dirPath, modelId + ".cif")
                return modelId, pth
        except Exception as e:
            logger.exception("Failing for %r %r %r with %s", parentId, targetId, modelDirPath, str(e))

        return None, None

    def __makeModelId(self, parentId, modelNum=1):
        modelId = "Q_" + parentId + "_%05d" % modelNum
        return modelId

    def __getToday(self):
        """Return a CIF style date-timestamp value for current local time -"""
        today = datetime.datetime.today()
        # format ="%Y-%m-%d:%H:%M"
        fmt = "%Y-%m-%d"
        return str(today.strftime(fmt))

    def __testUnMappedProtonation(self, fitAtomUnMappedL):
        # Check for the case of extra protons ...
        try:
            fitOk = True
            if fitAtomUnMappedL:
                for atU in fitAtomUnMappedL:
                    if atU.fitAtNo != 1:
                        fitOk = False
                        break
            return fitOk
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return False

    def __writeModel(self, targetId, targetObj, fitFD, fitXyzMapD, fitAtomUnMappedL, matchD, modelId, modelPath):
        """Write the chemical component model for the input chemical component Id and associated atom mapping and
        feature details --

            ComponentAtomDetails = namedtuple("ComponentAtomDetails", "index atNo name aType x y z fCharge")
            AlignAtomMap = namedtuple("AlignAtomMap", "refId refAtIdx refAtNo refAtName fitId fitAtIdx fitAtNo fitAtName")
            AlignAtomUnMapped = namedtuple("AlignAtomUnMapped", "fitId fitAtIdx fitAtNo fitAtType fitAtName fitAtFormalCharge x y z fitNeighbors")
        """
        try:
            unMappedTypeD = defaultdict(int)
            hAtomPrefix = "HEX"
            variantType = self.__getBuildVariant(targetId)
            #
            if not self.__testUnMappedProtonation(fitAtomUnMappedL):
                logger.info("Unmapped non-hydrogen atoms target %r model %r unMapped count (%d)", targetId, modelId, len(fitAtomUnMappedL))
                return False, variantType
            # Get atom partners for the unmapped atoms
            fitAtMapD = {}
            for refAtName, fAtTup in fitXyzMapD.items():
                fitAtMapD[fAtTup.atName] = refAtName
            if fitAtomUnMappedL:
                #  Check if neighbors are all mapped
                ok = True
                for fitUnTup in fitAtomUnMappedL:
                    for nAtName in fitUnTup.fitNeighbors:
                        if nAtName not in fitAtMapD:
                            ok = False
                            logger.info("Missing mapped neighbor for %r target %r model %r", nAtName, targetId, modelId)
                            break
                if not ok:
                    return False, variantType
                else:
                    logger.debug("%s match has unmapped protonation", modelId)
                    variantType = "tautomer_protomer"
            #
            #
            kList = ["xyz", "SMILES", "SMILES_STEREO", "InChI", "InChIKey"]
            for k in kList:
                if k not in fitFD:
                    logger.error("Fit feature dictionary for %s missing key %s", targetId, k)
                    return False, variantType
            # ------------
            dataContainer = DataContainer(modelId)
            #
            myContainer = targetObj
            dbName = myContainer.getName()
            if dbName.upper() != targetId.upper():
                logger.info("mismatch datablock (%r) and targetId (%r)", dbName, targetId)
            cObj = None
            if myContainer.exists("chem_comp"):
                cObj = myContainer.getObj("chem_comp")
            #
            #
            catName = "pdbx_chem_comp_model"
            if not dataContainer.exists(catName):
                dataContainer.append(DataCategory(catName, attributeNameList=["id", "comp_id"]))
            #
            parentId = targetId.split("|")[0]
            wObj = dataContainer.getObj(catName)
            wObj.setValue(modelId, "id", 0)
            wObj.setValue(parentId, "comp_id", 0)
            #
            # --------  ---------
            catName = "pdbx_chem_comp_model_atom"
            if not dataContainer.exists(catName):
                dataContainer.append(
                    DataCategory(catName, attributeNameList=["model_id", "atom_id", "type_symbol", "charge", "model_Cartn_x", "model_Cartn_y", "model_Cartn_z", "ordinal_id"])
                )
            wObj = dataContainer.getObj(catName)
            #
            if myContainer.exists("chem_comp_atom"):
                cObj = myContainer.getObj("chem_comp_atom")
            #
            #  Only write the mapped atoms in case we are missing hydrogens in the mapping
            #
            jj = 0
            for ii in range(cObj.getRowCount()):
                atName = cObj.getValue("atom_id", ii)
                atType = cObj.getValue("type_symbol", ii)
                if atName not in fitXyzMapD:
                    unMappedTypeD[atType] += 1
                    continue
                fitXyz = fitXyzMapD[atName]
                #
                # fCharge = cObj.getValue("charge", ii)
                #
                wObj.setValue(modelId, "model_id", jj)
                wObj.setValue(atName, "atom_id", jj)
                wObj.setValue(atType, "type_symbol", jj)
                #
                wObj.setValue(fitXyz.atFormalCharge, "charge", jj)
                wObj.setValue("%.4f" % fitXyz.x, "model_Cartn_x", jj)
                wObj.setValue("%.4f" % fitXyz.y, "model_Cartn_y", jj)
                wObj.setValue("%.4f" % fitXyz.z, "model_Cartn_z", jj)
                wObj.setValue(jj + 1, "ordinal_id", jj)
                jj += 1
            #
            # Add the unmapped atoms ...
            # AlignAtomUnMapped = namedtuple("AlignAtomUnMapped", "fitId fitAtIdx fitAtNo fitAtType fitAtName fitNeighbors")
            ii = wObj.getRowCount()
            for jj, uTup in enumerate(fitAtomUnMappedL):
                refAtomName = hAtomPrefix + str(jj)
                wObj.setValue(modelId, "model_id", ii)
                wObj.setValue(refAtomName, "atom_id", ii)
                wObj.setValue(uTup.fitAtType, "type_symbol", ii)
                wObj.setValue(uTup.fitAtFormalCharge, "charge", ii)
                wObj.setValue("%.4f" % uTup.x, "model_Cartn_x", ii)
                wObj.setValue("%.4f" % uTup.y, "model_Cartn_y", ii)
                wObj.setValue("%.4f" % uTup.z, "model_Cartn_z", ii)
                wObj.setValue(ii + 1, "ordinal_id", ii)
            # --------  ---------
            catName = "pdbx_chem_comp_model_bond"
            if not dataContainer.exists(catName):
                dataContainer.append(DataCategory(catName, attributeNameList=["model_id", "atom_id_1", "atom_id_2", "value_order", "ordinal_id"]))
            wObj = dataContainer.getObj(catName)
            #
            if myContainer.exists("chem_comp_bond"):
                cObj = myContainer.getObj("chem_comp_bond")
            #
            jj = 0
            for ii in range(cObj.getRowCount()):
                at1 = cObj.getValue("atom_id_1", ii)
                if at1 not in fitXyzMapD:
                    continue
                at2 = cObj.getValue("atom_id_2", ii)
                if at2 not in fitXyzMapD:
                    continue
                bType = cObj.getValue("value_order", ii)
                #
                wObj.setValue(modelId, "model_id", jj)
                wObj.setValue(at1, "atom_id_1", jj)
                wObj.setValue(at2, "atom_id_2", jj)
                wObj.setValue(bType, "value_order", jj)
                wObj.setValue(jj + 1, "ordinal_id", jj)
                jj += 1
            #
            ii = wObj.getRowCount()
            for jj, uTup in enumerate(fitAtomUnMappedL):
                at1 = hAtomPrefix + str(jj)
                for nAt in uTup.fitNeighbors:
                    at2 = fitAtMapD[nAt]
                    wObj.setValue(modelId, "model_id", ii)
                    wObj.setValue(at1, "atom_id_1", ii)
                    wObj.setValue(at2, "atom_id_2", ii)
                    wObj.setValue("SING", "value_order", ii)
                    wObj.setValue(ii + 1, "ordinal_id", ii)

            # --------  ---------
            catName = "pdbx_chem_comp_model_descriptor"
            if not dataContainer.exists(catName):
                dataContainer.append(DataCategory(catName, attributeNameList=["model_id", "type", "descriptor"]))
            wObj = dataContainer.getObj(catName)
            #
            ii = 0
            wObj.setValue(modelId, "model_id", ii)
            wObj.setValue("SMILES", "type", ii)
            wObj.setValue(fitFD["SMILES"], "descriptor", ii)
            ii += 1
            wObj.setValue(modelId, "model_id", ii)
            wObj.setValue("SMILES_CANONICAL", "type", ii)
            wObj.setValue(fitFD["SMILES_STEREO"], "descriptor", ii)
            ii += 1
            wObj.setValue(modelId, "model_id", ii)
            wObj.setValue("InChI", "type", ii)
            wObj.setValue(fitFD["InChI"], "descriptor", ii)
            ii += 1
            wObj.setValue(modelId, "model_id", ii)
            wObj.setValue("InChIKey", "type", ii)
            wObj.setValue(fitFD["InChIKey"], "descriptor", ii)
            #
            # --------  ---------
            if matchD["queryId"] is not None:
                catName = "pdbx_chem_comp_model_reference"
                if not dataContainer.exists(catName):
                    dataContainer.append(DataCategory(catName, attributeNameList=["model_id", "db_name", "db_code"]))
                wObj = dataContainer.getObj(catName)
                ii = 0
                wObj.setValue(modelId, "model_id", ii)
                wObj.setValue("COD", "db_name", ii)
                wObj.setValue(matchD["queryId"], "db_code", ii)
            #
            featureD = {}
            v = matchD["rValue"]
            vS = str(v)
            if v is not None and len(vS) > 0:
                featureD["r_factor"] = "%.3f" % float(v)
            #
            v = matchD["diffrnTemp"]
            vS = str(v)
            # remove string artifacts from temperature string ...
            if v is not None and len(vS) > 0:
                tV = vS.upper()
                try:
                    if tV.endswith("DEG.C"):
                        tV = tV.replace("AT", "")
                        tV = tV.replace("DEG.C", "")
                        tV = float(tV.strip())
                        tV = tV + 273.15
                    else:
                        tV = tV.replace("AT", "")
                        tV = tV.replace("K", "")
                        tV = float(tV.strip())
                    featureD["experiment_temperature"] = tV
                except Exception as e:
                    logger.exception("Temperature conversion fails for %s (%r) with %s", modelId, vS, tV)
            #
            v = matchD["publicationDOI"]
            vS = str(v)
            if v is not None and len(vS) > 0:
                featureD["publication_doi"] = v
            #
            v = matchD["version"]
            vS = str(v)
            if v is not None and len(vS) > 0:
                featureD["cod_version"] = v
            #
            if matchD["radiationSource"] and "neutron" in matchD["radiationSource"]:
                featureD["neutron_radiation_experiment"] = True
            if matchD["hasDisorder"] in ["Y"]:
                featureD["has_disorder"] = True
            #
            if len(unMappedTypeD) == 1 and "H" in unMappedTypeD:
                logger.info("model %r heavy_atoms_only", modelId)
                featureD["heavy_atoms_only"] = True
            else:
                featureD["all_atoms_have_sites"] = True
            # --------  ---------
            catName = "pdbx_chem_comp_model_feature"
            if not dataContainer.exists(catName):
                dataContainer.append(DataCategory(catName, attributeNameList=["model_id", "feature_name", "feature_value"]))
            wObj = dataContainer.getObj(catName)
            #
            fKeyList = ["experiment_temperature", "publication_doi", "r_factor", "csd_version"]
            ii = 0
            for fKey in fKeyList:
                if fKey in featureD:
                    wObj.setValue(modelId, "model_id", ii)
                    wObj.setValue(fKey, "feature_name", ii)
                    wObj.setValue(str(featureD[fKey]), "feature_value", ii)
                    ii += 1

            #
            boolKeyList = ["has_disorder", "neutron_radiation_experiment", "heavy_atoms_only", "all_atoms_have_sites"]
            for fKey in boolKeyList:
                if fKey in featureD:
                    if featureD[fKey]:
                        wObj.setValue(modelId, "model_id", ii)
                        wObj.setValue(fKey, "feature_name", ii)
                        wObj.setValue("Y", "feature_value", ii)
                        ii += 1
            #

            if variantType:
                wObj.setValue(modelId, "model_id", ii)
                wObj.setValue(variantType + "_match", "feature_name", ii)
                wObj.setValue("Y", "feature_value", ii)
                ii += 1

            # --------  ---------
            catName = "pdbx_chem_comp_model_audit"
            if not dataContainer.exists(catName):
                dataContainer.append(DataCategory(catName, attributeNameList=["model_id", "action_type", "date"]))
            wObj = dataContainer.getObj(catName)
            #
            ii = 0
            wObj.setValue(modelId, "model_id", ii)
            wObj.setValue("Initial release", "action_type", ii)
            wObj.setValue(self.__getToday(), "date", ii)
            # wObj.setValue('RCSB', 'processing_site',  ii)
            # wObj.setValue('JDW', 'annotator', ii)
            # wObj.setValue('?', 'details', ii)
            #
            mU = MarshalUtil(workPath=self.__cachePath)
            ok = mU.doExport(modelPath, [dataContainer], fmt="mmcif")
            return ok, variantType
        except Exception as e:
            logger.exception("Failing for %r with %s", targetId, str(e))
        return False, ""


class CODModelBuild(object):
    def __init__(self, cachePath, prefix=None, **kwargs):
        self.__cachePath = cachePath

        self.__prefix = prefix
        startTime = time.time()
        useCache = True
        self.__ccUrlTarget = kwargs.get("ccUrlTarget", None)
        self.__birdUrlTarget = kwargs.get("birdUrlTarget", None)
        ccFileNamePrefix = "cc-%s" % self.__prefix if self.__prefix else "cc"
        oeFileNamePrefix = "oe-%s" % self.__prefix if self.__prefix else "oe"
        self.__ccmP = ChemCompMoleculeProvider(
            ccUrlTarget=self.__ccUrlTarget, birdUrlTarget=self.__birdUrlTarget, cachePath=cachePath, useCache=useCache, ccFileNamePrefix=ccFileNamePrefix
        )
        ok1 = self.__ccmP.testCache()

        self.__ccSIdxP = ChemCompSearchIndexProvider(cachePath=cachePath, useCache=useCache, ccFileNamePrefix=ccFileNamePrefix)
        ok2 = self.__ccSIdxP.testCache()

        molLimit = kwargs.get("molLimit", None)
        quietFlag = kwargs.get("quietFlag", True)
        fpTypeList = kwargs.get("fpTypeList", [])
        screenTypeList = kwargs.get("screenTypeList", [])
        limitPerceptions = kwargs.get("limitPerceptions", False)
        numProc = kwargs.get("numProc", 4)
        self.__oesmP = OeSearchMoleculeProvider(
            ccUrlTarget=self.__ccUrlTarget,
            birdUrlTarget=self.__birdUrlTarget,
            cachePath=self.__cachePath,
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
        ok3 = self.__oesmP.testCache()
        self.__startTime = time.time()
        logger.info("Completed chemical component search index load %r (%.4f seconds)", ok1 & ok2 & ok3, time.time() - startTime)
        #
        logger.info("Starting model build (%s) at %s", __version__, time.strftime("%Y %m %d %H:%M:%S", time.localtime()))

    def getModelDirFilePath(self):
        dN = "cod-%s-model-files" % self.__prefix if self.__prefix else "cod-model-files"
        return os.path.join(self.__cachePath, dN)

    def getModelImageDirFilePath(self):
        dN = "cod-%s-model-image" % self.__prefix if self.__prefix else "cod-model-images"
        return os.path.join(self.__cachePath, dN)

    def __getModelIndexPath(self):
        return os.path.join(self.getModelDirFilePath(), "cod-model-index.json")

    def fetchModelIndex(self):
        mD = {}
        try:
            mU = MarshalUtil(workPath=self.__cachePath)
            fp = self.__getModelIndexPath()
            mD = mU.doImport(fp, fmt="json")
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return mD

    def storeModelIndex(self, mD):
        try:
            mU = MarshalUtil(workPath=self.__cachePath)
            fp = self.__getModelIndexPath()
            ok = mU.doExport(fp, mD, fmt="json", indent=3)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            ok = False
        return ok

    def build(self, alignType="relaxed-stereo", numProc=4, chunkSize=10, verbose=False):
        """Run the model build step in the chemical component model workflow.

        Args:
          alignType (str):  "relaxed"|"strict"| relaxed-stereo".  Default: relaxed-stereo
          numProc (int, optional): number of processes to invoke. Defaults to 4.
          chunkSize (int, optional): work chunksize. Defaults to 10.
          verbose (bool, optional): verbose logging.  Defaults to False.

        Returns:
            (dict): {searchId: [{"targetId": , "modelId": , "modelPath": ,"matchId": , "parentId": , "rFactor": , }]

        """
        retD = {}
        try:
            ccms = CODModelSearch(self.__cachePath, prefix=self.__prefix)
            modelDirPath = self.getModelDirFilePath()
            imageDirPath = self.getModelImageDirFilePath()

            #
            idxIdD = ccms.getResultIndex()
            idxIdL = list(idxIdD.keys())
            #
            logger.info("Using COD search result index length idxD (%d)", len(idxIdD))
            #
            cmbw = CODModelBuildWorker(self.__cachePath, verbose=verbose)
            mpu = MultiProcUtil(verbose=True)
            mpu.setWorkingDir(modelDirPath)
            mpu.setOptions(
                optionsD={
                    "modelDirPath": modelDirPath,
                    "imageDirPath": imageDirPath,
                    "alignType": alignType,
                    "ccSIdxP": self.__ccSIdxP,
                    "idxIdD": idxIdD,
                    "oesmP": self.__oesmP,
                    "ccmP": self.__ccmP,
                }
            )
            #
            mpu.set(workerObj=cmbw, workerMethod="build")
            ok, failList, resultList, _ = mpu.runMulti(dataList=idxIdL, numProc=numProc, numResults=1, chunkSize=chunkSize)
            logger.info("Run ended with status %r success count %d failures %r", ok, len(resultList[0]), len(failList))
            successList = copy.copy(resultList[0])
            for tD in successList:
                retD.setdefault(tD["parentId"], []).append(tD)
            #
            if retD:
                logger.info("Completed build with models for %d parent chemical definitions", len(retD))
            else:
                logger.info("No models built")
            ok = self.storeModelIndex(retD)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return retD

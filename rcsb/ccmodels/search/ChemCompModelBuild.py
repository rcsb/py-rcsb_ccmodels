##
# File:    ChemCompModelBuild.py
# Author:  J. Westbrook
# Date:    19-Jan-2021
# Version: 0.001
#
# Updated:
##
"""
Generate model files for chemical component model workflow.

"""
__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "john.westbrook@rcsb.org"
__license__ = "Apache 2.0"

import datetime
import logging
import os
import time

from mmcif.api.DataCategory import DataCategory
from mmcif.api.PdbxContainers import DataContainer
from rcsb.ccmodels.search import __version__
from rcsb.ccmodels.search.ChemCompModelSearch import ChemCompModelSearch
from rcsb.utils.chem.ChemCompSearchIndexProvider import ChemCompSearchIndexProvider
from rcsb.utils.chem.OeAlignUtils import OeAlignUtils
from rcsb.utils.chem.OeDepictAlign import OeDepictMCSAlignMultiPage
from rcsb.utils.chem.OeDepictAlign import OeDepictMCSAlignPage
from rcsb.utils.io.IndexUtils import CcdcMatchIndex
from rcsb.utils.io.MarshalUtil import MarshalUtil

logger = logging.getLogger(__name__)


class ChemCompModelBuild(object):
    def __init__(self, cachePath, prefix=None):
        self.__cachePath = cachePath

        self.__prefix = prefix
        startTime = time.time()
        useCache = True
        ccFileNamePrefix = "cc-%s" % self.__prefix if self.__prefix else "cc"
        self.__ccSIdxP = ChemCompSearchIndexProvider(cachePath=cachePath, useCache=useCache, ccFileNamePrefix=ccFileNamePrefix)
        ok = self.__ccSIdxP.testCache()
        logger.info("Completed chemical component search index load %r (%.4f seconds)", ok, time.time() - startTime)
        #
        self.__startTime = time.time()
        logger.info("Starting search (%s) at %s", __version__, time.strftime("%Y %m %d %H:%M:%S", time.localtime()))

    def getModelDirFilePath(self):
        dN = "cc-%s-model-files" % self.__prefix if self.__prefix else "cc-model-files"
        return os.path.join(self.__cachePath, dN)

    def getModelImageDirFilePath(self):
        dN = "cc-%s-model-image" % self.__prefix if self.__prefix else "cc-model-images"
        return os.path.join(self.__cachePath, dN)

    def build(self, alignType="relaxed-stereo"):
        """Run the model build step in the chemical component model workflow.

        Args:
          alignType (str):  "relaxed"|"strict"| relaxed-stereo".  Default: relaxed-stereo

        Returns:
            (dict): {searchId: [{"targetId": , "modelId": , "modelPath": ,"matchId": , "parentId": , "rFactor": , }]

        """
        retD = {}
        try:
            ccms = ChemCompModelSearch(self.__cachePath, None, None, prefix=self.__prefix)
            modelDirPath = self.getModelDirFilePath()
            imageDirPath = self.getModelImageDirFilePath()
            #
            idxPathD = ccms.getResultIndex()
            logger.info("Result index length ridxD (%d)", len(idxPathD))
            parentD = {}
            for sId, idxPath in idxPathD.items():
                logger.info("")
                logger.info("======== ============ ============ =========== ============ ===============")
                parentId = sId.split("|")[0]

                logger.info("Start model build for search Id (%s) %s", parentId, sId)
                #
                matchObjIt = CcdcMatchIndex(indexFilePath=idxPath)
                matchObjIt.sort()
                pairList = []
                startingModel = 1
                for matchObj in matchObjIt:
                    targetCifPath = matchObj.getTargetCcPath()
                    fitMolFilePath = matchObj.getMol2Path()
                    matchId = matchObj.getIdentifier()
                    targetId = matchObj.getTargetId()
                    matchTitle = "CCDC Code  " + matchId
                    ccTitle = "Chemical Component " + targetId
                    #
                    nAtomsRef, refFD, nAtomsFit, fitFD, fitXyzMapD, fitAtomExtraTupL = self.__alignModelSubStruct(
                        targetCifPath, fitMolFilePath, alignType=alignType, fitTitle=matchId, refTitle=targetId, verbose=False
                    )
                    logger.debug(
                        ">>> %s - %s nAtomsRef %d nAtomsFit %d atommapL (%d) fitAtomExtraL (%d)", targetId, matchId, nAtomsRef, nAtomsFit, len(fitXyzMapD), len(fitAtomExtraTupL)
                    )
                    # All reference atoms must match -
                    if nAtomsRef == 0 or nAtomsFit == 0 or nAtomsRef != len(fitXyzMapD):
                        logger.debug("     ----- Rejecting (incomplete/mismatch) ref %s with matchId %s", targetId, matchId)
                        continue
                    # Check for the case of extra protons ...
                    fitOk = True
                    if fitAtomExtraTupL:
                        for fitAtomExtraTup in fitAtomExtraTupL:
                            if fitAtomExtraTup[1] != 1:
                                fitOk = False
                                break
                    #
                    if not fitOk:
                        logger.debug("     ----- Rejecting (extra fit bits) ref %s with matchId %s", targetId, matchId)
                        continue
                    #
                    if refFD["SMILES_STEREO"] != fitFD["SMILES_STEREO"]:
                        logger.info("SMILES mismatch")
                        logger.info("Ref %-8s SMILES: %s", targetId, refFD["SMILES_STEREO"])
                        logger.info("Fit %-8s SMILES: %s", matchId, fitFD["SMILES_STEREO"])
                        continue
                    # --------- ----------------
                    #  Accept the match
                    # --------- ----------------
                    matchId = matchObj.getIdentifier()
                    targetId = matchObj.getTargetId()
                    parentD.setdefault(parentId, []).append(matchId)
                    #
                    refImageFileName = "ref_" + targetId + "_" + matchId + ".svg"
                    refImagePath = os.path.join(imageDirPath, sId, refImageFileName)
                    self.__pairDepictPage(refImagePath, sId, ccTitle, refFD["OEMOL"], matchId, matchTitle, fitFD["OEMOL"], alignType=alignType)
                    # --------- ------------------
                    pairList.append((sId, refFD["OEMOL"], matchId, fitFD["OEMOL"]))

                    modelId, modelPath = self.__makeModelPath(modelDirPath, parentId, startingModelNum=startingModel, maxModels=200)
                    ok = self.__writeModel(targetId, targetCifPath, fitFD, fitXyzMapD, matchObj, modelId, modelPath)
                    startingModel += 1
                    if ok:
                        retD.setdefault(sId, []).append(
                            {
                                "modelId": modelId,
                                "searchId": targetId,
                                "parentId": parentId,
                                "matchId": matchId,
                                "modelPath": modelPath,
                                "rFactor": matchObj.getRFactor(),
                                "hasDisorder": matchObj.getHasDisorder(),
                            }
                        )
                #
                if pairList:
                    pdfImagePath = os.path.join(imageDirPath, sId, sId + "-all-pairs.pdf")
                    self.__depictFitList(sId, pdfImagePath, pairList, alignType=alignType)
                if retD:
                    logger.info("Built %d models for %s", len(retD), sId)
                else:
                    logger.info("No models built for %s", sId)
            ok = self.storeModelIndex(retD)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return retD

    def __getModelIndexPath(self):
        return os.path.join(self.getModelDirFilePath(), "model-index.json")

    def fetchModelIndex(self):
        mD = {}
        try:
            mU = MarshalUtil()
            fp = self.__getModelIndexPath()
            mD = mU.doImport(fp, fmt="json")
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return mD

    def storeModelIndex(self, mD):
        try:
            mU = MarshalUtil()
            fp = self.__getModelIndexPath()
            ok = mU.doExport(fp, mD, fmt="json", indent=3)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            ok = False
        return ok

    def __getBuildVariant(self, sId):
        """Lookup the build type from the input chemical component search index Id.

        Args:
            sId (str): chemical component search index Id

        Returns:
            str : "tautomer" | ""
        """
        sidxD = self.__ccSIdxP.getIndexEntry(sId)
        if "protomer" in sidxD["build-type"]:
            return "tautomer"
        if "tautomer" in sidxD["build-type"]:
            return "tautomer"
        return ""

    def __alignModelMcss(self, ccRefPath, molFitPath, alignType="exact", fitTitle=None, refTitle=None, unique=False, minFrac=0.9, verbose=False):
        """Align (substructure) chemical component definition search target with the candidate matching reference molecule.

        Args:
            ccRefPath (str): path to search target chemical component definition
            molFitPath (str): path to matched molfile
            alignType (str, optional):  strict|relaxed|relaxed-stereo. Defaults to "strict"
            fitTitle (str, optional): fit molecule title. Defaults to None.
            refTitle (str, optional): reference molecule title. Defaults to None.
            unique (bool, optional): find unique matches. Defaults to False.
            minFrac (float, optional): minimum atom-level matching coverage fraction. Defaults to 0.9.
            verbose (bool, optional): enable verbose output. Defaults to False.

        Returns:
            (tuple): (  number of atoms reference molecule,
                        reference molecule feature dictionary {"Formula": , "SMILES": , "SMILES_STEREO": , "InChI": , "InChIKey":,}
                        number of atoms in fit molecule,
                        fit molecule feature dictionary  {"Formula": , "SMILES": , "SMILES_STEREO": , "InChI": , "InChIKey":, 'xyz': }
                        atomMapD {refAtName: (fit x, fit y, fit z)}
                        fitAtomExtraL: [(atName, atomicNumber),(), ...]
                    )

            Note: 'xyz': [(ii, atm.GetIdx(), atm.GetAtomicNum(), atm.GetName(), atm.GetType(), fitX, fitY, fitZ))
        """
        try:
            logger.debug("Align target cc %s with matching model %s", ccRefPath, molFitPath)
            oesU = OeAlignUtils(verbose=verbose)
            oesU.setSearchType(sType=alignType)
            oesU.setRefPath(ccPath=ccRefPath, title=refTitle)
            oesU.setFitPath(molFitPath, title=fitTitle, suppressHydrogens=False, fType="sdf", importType="3D")
            (nAtomsRef, refFD, nAtomsFit, fitFD, atomMapL, fitAtomExtraL) = oesU.doAlignMcss(unique=unique, minFrac=minFrac)
            # -----
            # fitFD['xyz']: [(ii, atm.GetIdx(), atm.GetAtomicNum(), atm.GetName(), atm.GetType(), x, y, z))
            #
            # tD = {refAtName: fitAtIdx, ...}
            tD = {tup[3]: tup[5] for tup in atomMapL}
            # fitXyzD  = {atIdx: (x,y,z), (x,y,z), ..}
            fitXyzD = {tup[1]: (tup[5], tup[6], tup[7]) for tup in fitFD["xyz"]}
            atomMapD = {atName: fitXyzD[fitAtIdx] for atName, fitAtIdx in tD.items() if fitAtIdx in fitXyzD}
            # -----
            return nAtomsRef, refFD, nAtomsFit, fitFD, atomMapD, fitAtomExtraL
        except Exception as e:
            logger.exception("Failing for %r and %r with %s", ccRefPath, molFitPath, str(e))
        return 0, {}, 0, {}, {}, []

    def __alignModelSubStruct(self, ccRefPath, molFitPath, alignType="strict", fitTitle=None, refTitle=None, verbose=False):
        """Align (substructure) chemical component definition search target with the candidate matching reference molecule.

        Args:
            ccRefPath (str): path to search target chemical component definition
            molFitPath (str): path to matched molfile
            alignType (str, optional):  strict|relaxed|relaxed-stereo. Defaults to "strict"
            fitTitle (str, optional): fit molecule title. Defaults to None.
            refTitle (str, optional): reference molecule title. Defaults to None.
            verbose (bool, optional): enable verbose output. Defaults to False.

        Returns:
            (tuple): (  number of atoms reference molecule,
                        reference molecule feature dictionary {"Formula": , "SMILES": , "SMILES_STEREO": , "InChI": , "InChIKey":,}
                        number of atoms in fit molecule,
                        fit molecule feature dictionary  {"Formula": , "SMILES": , "SMILES_STEREO": , "InChI": , "InChIKey":, 'xyz': }
                        atomMapD {refAtName: (fit x, fit y, fit z)}
                        fitAtomExtraL: [(atName, atomicNumber),(), ...]
                    )

            Note: 'xyz': [(ii, atm.GetIdx(), atm.GetAtomicNum(), atm.GetName(), atm.GetType(), fitX, fitY, fitZ))
        """
        try:
            logger.debug("Align target cc %s with matching model %s", ccRefPath, molFitPath)
            oesU = OeAlignUtils(verbose=verbose)
            oesU.setSearchType(sType=alignType)
            oesU.setRefPath(ccPath=ccRefPath, title=refTitle)
            oesU.setFitPath(molFitPath, title=fitTitle, suppressHydrogens=False, fType="sdf", importType="3D")
            (nAtomsRef, refFD, nAtomsFit, fitFD, atomMapL, fitAtomExtraL) = oesU.doAlignSs(unique=False)
            # -----
            # fitFD['xyz']: [(ii, atm.GetIdx(), atm.GetAtomicNum(), atm.GetName(), atm.GetType(), x, y, z))
            #
            # tD = {refAtName: fitAtIdx, ...}
            tD = {tup[3]: tup[5] for tup in atomMapL}
            # fitXyzD  = {atIdx: (x,y,z), (x,y,z), ..}
            fitXyzD = {tup[1]: (tup[5], tup[6], tup[7]) for tup in fitFD["xyz"]}
            atomMapD = {atName: fitXyzD[fitAtIdx] for atName, fitAtIdx in tD.items() if fitAtIdx in fitXyzD}
            # -----
            return nAtomsRef, refFD, nAtomsFit, fitFD, atomMapD, fitAtomExtraL
        except Exception as e:
            logger.exception("Failing for %r and %r with %s", ccRefPath, molFitPath, str(e))
        return 0, {}, 0, {}, {}, []

    def __pairDepictPage(self, imagePath, refId, refTitle, refMol, fitId, fitTitle, fitMol, alignType="strict"):
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
            oed = OeDepictMCSAlignMultiPage()
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

    def __makeModelPath(self, modelDirPath, parentId, startingModelNum=1, maxModels=200, scanExisting=False):
        pth = None
        ii = startingModelNum
        dirPath = os.path.join(modelDirPath, parentId)
        if not os.access(dirPath, os.R_OK):
            os.makedirs(dirPath)
            modelId = self.__makeModelId(parentId, modelNum=ii)
            pth = os.path.join(dirPath, modelId + ".cif")
            return modelId, pth
        # optionally scan over any existing models before selecting a new model number.
        if scanExisting:
            while True:
                modelId = self.__makeModelId(parentId, modelNum=ii)
                pth = os.path.join(dirPath, modelId + ".cif")
                if not os.access(pth, os.R_OK):
                    return modelId, pth
                #
                ii += 1
                if ii > maxModels:
                    break
        else:
            modelId = self.__makeModelId(parentId, modelNum=ii)
            pth = os.path.join(dirPath, modelId + ".cif")
            return modelId, pth

        return None, None

    def __makeModelId(self, parentId, modelNum=1):
        modelId = "M_" + parentId + "_%05d" % modelNum
        return modelId

    def __writeModel(self, targetId, targetPath, fitFD, fitXyzD, matchObj, modelId, modelPath):
        """Write the chemical component model for the input chemical component Id and associated atom mapping and
        feature details --

        """
        try:
            kList = ["xyz", "SMILES", "SMILES_STEREO", "InChI", "InChIKey"]
            for k in kList:
                if k not in fitFD:
                    logger.error("Fit feature dictionary for %s missing key %s", targetId, k)
                    return False
            #
            dataContainer = DataContainer(modelId)
            #
            mU = MarshalUtil()
            myContainerList = mU.doImport(targetPath, fmt="mmcif")
            myContainer = myContainerList[0]
            dbName = myContainer.getName()
            if dbName.upper() != targetId.upper():
                logger.error("mismatch datablock (%r) and targetId (%r)", dbName, targetId)
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
            for ii in range(cObj.getRowCount()):
                atName = cObj.getValue("atom_id", ii)
                atType = cObj.getValue("type_symbol", ii)
                fCharge = cObj.getValue("charge", ii)
                #
                fitXyz = fitXyzD[atName]
                # --------  ---------
                wObj.setValue(modelId, "model_id", ii)
                wObj.setValue(atName, "atom_id", ii)
                wObj.setValue(atType, "type_symbol", ii)
                wObj.setValue(fCharge, "charge", ii)
                wObj.setValue("%.4f" % fitXyz[0], "model_Cartn_x", ii)
                wObj.setValue("%.4f" % fitXyz[1], "model_Cartn_y", ii)
                wObj.setValue("%.4f" % fitXyz[2], "model_Cartn_z", ii)
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
            for ii in range(cObj.getRowCount()):
                at1 = cObj.getValue("atom_id_1", ii)
                at2 = cObj.getValue("atom_id_2", ii)
                bType = cObj.getValue("value_order", ii)
                wObj.setValue(modelId, "model_id", ii)
                wObj.setValue(at1, "atom_id_1", ii)
                wObj.setValue(at2, "atom_id_2", ii)
                wObj.setValue(bType, "value_order", ii)
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
            if matchObj.getIdentifier() is not None:
                catName = "pdbx_chem_comp_model_reference"
                if not dataContainer.exists(catName):
                    dataContainer.append(DataCategory(catName, attributeNameList=["model_id", "db_name", "db_code"]))
                wObj = dataContainer.getObj(catName)
                ii = 0
                wObj.setValue(modelId, "model_id", ii)
                wObj.setValue("CSD", "db_name", ii)
                wObj.setValue(matchObj.getIdentifier(), "db_code", ii)
            #
            featureD = {}
            v = matchObj.getRFactor()
            vS = str(v)
            if v is not None and len(vS) > 0:
                featureD["r_factor"] = "%.3f" % float(v)
            #
            v = matchObj.getTemperature()
            vS = str(v)
            if v is not None and len(vS) > 0:
                featureD["experiment_temperature"] = v.replace("at ", "")
            #
            v = matchObj.getCitationDOI()
            vS = str(v)
            if v is not None and len(vS) > 0:
                featureD["publication_doi"] = v
            #
            v = matchObj.getCsdVersion()
            vS = str(v)
            if v is not None and len(vS) > 0:
                featureD["csd_version"] = v
            #
            if matchObj.getRadiationSource() in ["Neutron"]:
                featureD["neutron_radiation_experiment"] = True
            if matchObj.getHasDisorder():
                featureD["has_disorder"] = True
            #
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
            boolKeyList = ["all_atoms_have_sites", "has_disorder", "neutron_radiation_experiment"]
            for fKey in boolKeyList:
                if fKey in featureD:
                    if featureD[fKey]:
                        wObj.setValue(modelId, "model_id", ii)
                        wObj.setValue(fKey, "feature_name", ii)
                        wObj.setValue("Y", "feature_value", ii)
                        ii += 1
            #
            bType = self.__getBuildVariant(targetId)
            if bType:
                wObj.setValue(modelId, "model_id", ii)
                wObj.setValue(bType + "_match", "feature_name", ii)
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
            ok = mU.doExport(modelPath, [dataContainer], fmt="mmcif")
            return ok
        except Exception as e:
            logger.exception("Failing for %r %r with %s", targetId, targetPath, str(e))
        return False

    def __getToday(self):
        """Return a CIF style date-timestamp value for current local time -"""
        today = datetime.datetime.today()
        # format ="%Y-%m-%d:%H:%M"
        fmt = "%Y-%m-%d"
        return str(today.strftime(fmt))

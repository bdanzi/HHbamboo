###------------------------------------------------------------------------------------###
###  This code needs envConfig cern.ini                                                ###
###  It creates object selection for bbmm analysis and produces root files             ###
###  Author: Brunella D'Anzi                                                           ###
###  Date: 20/08/2023                                                                  ###
###------------------------------------------------------------------------------------###

from bamboo.analysismodules import NanoAODHistoModule
from bamboo.treedecorators import NanoAODDescription
from bamboo.scalefactors import binningVariables_nano,lumiPerPeriod_default, BtagSF, get_correction, get_bTagSF_itFit, makeBtagWeightItFit

from bamboo import treedecorators as td
from bamboo import treefunctions as op
from bamboo import scalefactors

from bamboo.plots import Plot, EquidistantBinning, CutFlowReport
from bamboo import treefunctions as op
from bamboo.treeoperations import Const

from itertools import chain
from functools import partial
import sys
sys.path.append('/eos/user/b/bdanzi/LCG105/bamboo/examples/HZA_CMSDASLPC23')
from scalefactorsbbmm import ScaleFactorsbbmm
import os

from bamboo.root import loadDependency 
import logging
logger = logging.getLogger(__name__)

if os.environ["VIRTUAL_ENV"] == "":
    print("$VIRTUAL_ENV is not set. Please, activate your virtual environment")
    exit(-1)

print(os.environ["VIRTUAL_ENV"])

class category:
    def __init__(self,nMuons=0, nElectrons=0):
        self.mu  = nMuons
        self.ele = nElectrons
        self._name = self.__str__()

    def __str__(self):
        strmu  = "" if not self.mu else f"{self.mu}mu"
        strele = "" if not self.ele else f"{self.ele}ele"

        return strmu+strele

    def nMuons(self):
        return self.mu
    def nElectrons(self):
        return self.ele
    def name(self):
        return self._name

class triggerInfo:
#trigObj is a collection of triggers, trigPt is the Pt of the trigger and trigSel is the output of a op.select function
    def __init__(self,trigName,trigObj, trigPt, trigSel):
        self.tObj  = trigObj
        self.tPt   = trigPt
        self.tSel  = trigSel
        self.tName = trigName
       
    def Name(self):
        return self.tName
    def Obj(self):
        return self.tObj
    def Pt(self):
        return self.tPt
    def Sel(self):
        return self.tSel

class BaseNanoHHtobbmumu(NanoAODHistoModule):
    """ Module for HH bbmumu analysis """
    def __init__(self, args):
        super(BaseNanoHHtobbmumu, self).__init__(args)
        
    def addArgs(self,parser):
        super(BaseNanoHHtobbmumu, self).addArgs(parser)
        parser.title = """"""
        parser.add_argument("--PrintYield", action= "store_true",default= True,help="Print yield to screen (for debugging)")
        parser.add_argument("--mvaSkim", action="store_true",default= True,help="Produce MVA training skims")
        parser.add_argument("--mvaEval", action="store_true", default= False,help="Import MVA model and evaluate it on the dataframe")
        parser.add_argument("--NoSystematics", action= "store_true",default=True,help="Disable all systematic variations (default=False)")
        parser.add_argument("--Synchronization", action= "store_true",default=False,help="Produce the skims for the synchronization (without triggers, corrections of flags) if alone. If sync for specific selection, lepton, jet and channel arguments need to be used")
        parser.add_argument("--BtagReweightingOn", action= "store_true", default=False,help="Btag ratio study : Btag SF applied (without the ratio), will only do the plots for reweighting (jets and leptons args are ignored)")
        parser.add_argument("--BtagReweightingOff", action= "store_true",default=False,help="Btag ratio study : Btag Sf not applied (without the ratio), will only do the plots for reweighting (jets and leptons args are ignored)")

    #############################################################################
    #                           PrepareTree                                     #
    #############################################################################
    def prepareTree(self, tree, sample=None, sampleCfg=None, backend=None):
        
        from bamboo.treedecorators import NanoAODDescription, nanoRochesterCalc, nanoJetMETCalc, nanoJetMETCalc_METFixEE2017, nanoFatJetCalc
        from bamboo.analysisutils import configureJets, makeMultiPrimaryDatasetTriggerSelection, configureType1MET
        from bamboo.analysisutils import configureRochesterCorrection, makePileupWeight
        import ROOT
        era = sampleCfg["era"]   
        
        def isMC():
            if sampleCfg['type'] == 'data':
                return False
            elif sampleCfg['type'] in ['mc', 'signal']:
                return True
            else:
                raise RuntimeError(
                    f"The type '{sampleCfg['type']}' of {sample} dataset not understood.")
                              
        self.is_MC = isMC()
        self.era = era
        self.sample = sample
        if not self.is_MC:
            subera = sample.split("_")[1]
            print("This is not a MC sample. Subera: ",subera )
        metName = "MET"
        # nanoJetMETCalc_both = td.CalcCollectionsGroups(
        #     Jet=("pt", "mass"), changes={metName: (f"{metName}T1", f"{metName}T1Smear")},
        #     **{metName: ("pt", "phi")})
        # nanoJetMETCalc_data = td.CalcCollectionsGroups(
        #     Jet=("pt", "mass"), changes={metName: (f"{metName}T1",)},
        #     **{metName: ("pt", "phi")})
        tree,noSel,be,lumiArgs = NanoAODHistoModule.prepareTree(self, tree, sample=sample, sampleCfg=sampleCfg, description=NanoAODDescription.get('v7', year=era, isMC=self.is_MC, systVariations=[nanoRochesterCalc,td.CalcCollectionsGroups(Jet=("pt", "mass"))]), backend=backend)#td.CalcCollectionsGroups(Jet=("pt", "mass")),nanoJetMETCalc = CalcCollectionsGroups(Jet=("pt", "mass"), MET=("pt", "phi"))
        era = sampleCfg["era"]
        configureRochesterCorrection(variProxy  = tree._Muon,
                                    paramsFile = os.path.join(os.path.dirname(__file__),"RoccoR","RoccoR"+ era +"UL.txt"),
                                    isMC       = self.is_MC,
                                    backend    = be, 
                                    uName      = sample)
            
        #############################################################################
        #                           Links :                                         #  
        # JEC : https://twiki.cern.ch/twiki/bin/viewauth/CMS/JECDataMC              #
        # JER (smear) : https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution      #                          #
        #############################################################################
        #  https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution#JER_Scaling_factors_and_Uncertai
            
        if self.is_MC:
            JECTagDatabase = {"2016preVFP": "Summer19UL16APV_V7_MC", 
                              "2016postVFP": "Summer19UL16_V7_MC", 
                              "2017": "Summer19UL17_V5_MC", 
                              "2018": "Summer19UL18_V5_MC",
                              "2022": "Summer22_22Sep2023_V2_MC",
                              "2022EE": "Summer22EE_22Sep2023_V2_MC"}
            JERTagDatabase = {"2016preVFP": "Summer20UL16APV_JRV3_MC", 
                              "2016postVFP": "Summer20UL16_JRV3_MC", 
                              "2017": "Summer19UL17_JRV3_MC", 
                              "2018": "Summer19UL18_JRV2_MC",
                              "2022": "Summer22EEPrompt22_JRV1_MC",
                              "2022EE": "Summer22EEPrompt22_JRV1_MC"}
      
            if era in JECTagDatabase.keys():
                configureJets(
                    variProxy               = tree._Jet,
                    jetType                 = "AK4PFchs", #if "201" in era else "AK4PFPuppi",
                    jec                     = JECTagDatabase[era],
                    smear                   = JERTagDatabase[era],
                    jecLevels               = "default",
                    regroupTag              = "V2",
                    jesUncertaintySources   = (["Total"] if self.is_MC else None),
                    #enableSystematics     =   lambda v : not "jesTotal" in v,
                    mayWriteCache           = self.args.distributed != "worker",
                    isMC                    = self.is_MC,
                    backend                 = be,
                    uName                   = sample
                    )
                #configureType1MET(
                #    variProxy               = getattr(tree, f"_{metName}"),
                #    jec                     = JECTagDatabase[era],
                #    smear                   = JERTagDatabase[era],
                #    isT1Smear               = True,
                #    #regroupTag              = "V2",
                #    jesUncertaintySources   = (["Total"] if self.is_MC else None),
                #    mayWriteCache           = self.args.distributed != "worker",
                #    isMC                    = self.is_MC,
                #    backend                 = be,
                #    uName                   = sample
                #    )
                
            else:
                    raise RuntimeError("Could not find appropriate JEC tag for MC")
            if not self.is_MC:
                JECTagDatabase = {
                                  "2016B": "Summer19UL16APV_RunBCD_V7_DATA",
                                  "2016C": "Summer19UL16APV_RunBCD_V7_DATA",
                                  "2016D": "Summer19UL16APV_RunBCD_V7_DATA",
                                  "2016E": "Summer19UL16APV_RunEF_V7_DATA", 
                                  "2016FpreVFP": "Summer19UL16APV_RunEF_V7_DATA",                                    
                                  "2016FpostVFP": "Summer19UL16_RunFGH_V7_DATA",
                                  "2016G": "Summer19UL16_RunFGH_V7_DATA",
                                  "2016H": "Summer19UL16_RunFGH_V7_DATA",   
                                  "2017B": "Summer19UL17_RunB_V5_DATA",
                                  "2017C": "Summer19UL17_RunC_V5_DATA",
                                  "2017D": "Summer19UL17_RunD_V5_DATA",
                                  "2017E": "Summer19UL17_RunE_V5_DATA",
                                  "2017F": "Summer19UL17_RunF_V5_DATA",
                                  "2018A": "Summer19UL18_RunA_V5_DATA",
                                  "2018B": "Summer19UL18_RunB_V5_DATA",
                                  "2018C": "Summer19UL18_RunC_V5_DATA",
                                  "2018D": "Summer19UL18_RunD_V5_DATA",
                                  "2022C": "Summer22_22Sep2023_RunCD_V2_DATA",  
                                  "2022D": "Summer22_22Sep2023_RunCD_V2_DATA",
                                  "2022F": "Summer22EE_22Sep2023_RunF_V2_DATA",
                                  "2022G": "Summer22EE_22Sep2023_RunG_V2_DATA"                        
                                  }
                if subera in JECTagDatabase.keys():
                    # configureJets(variProxy, jetType, jec=None, jecLevels='default', smear=None, 
                    # useGenMatch=True, genMatchDR=0.2, genMatchDPt=3.0, jesUncertaintySources=None, regroupTag='', 
                    # uncertaintiesFallbackJetType=None, splitJER=False, addHEM2018Issue=False, enableSystematics=None, 
                    # subjets=None, mcYearForFatJets=None, isTau21DDT=False, jms=None, jmr=None, gms=None, gmr=None, 
                    # cachedir=None, mayWriteCache=False, isMC=False, backend=None, uName='')
                    configureJets(
                        variProxy               = tree._Jet,
                        jetType                 = "AK4PFchs",
                        jec                     = JECTagDatabase[subera],
                        #jecLevels               = "default",
                        #regroupTag              = "",
                        jesUncertaintySources   = (["Total"] if self.is_MC else None),
                        mayWriteCache           = self.args.distributed != "worker",
                        isMC                    = self.is_MC,
                        backend                 = be,
                        uName                   = sample  
                    )
                    #configureType1MET(
                    #   variProxy               = getattr(tree, f"_{metName}"),
                    #   jec                     = JECTagDatabase[subera],
                    #   mayWriteCache           = self.args.distributed != "worker",
                    #   isMC                    = self.is_MC,
                    #   backend                 = be,
                    #   uName                   = sample
                    #)
                    
        return tree, noSel,be,lumiArgs
 
    def definePlots(self, tree, noSel, sample=None, sampleCfg=None):
        
        from bamboo.analysisutils import configureJets, makeMultiPrimaryDatasetTriggerSelection
        from bamboo.analysisutils import configureRochesterCorrection, makePileupWeight
        
        self.plotDefaults = {  "normalized": False,
                               "legend-columns": 2,
                               "no-data": False,
                               "show-overflow": False,    
                               "show-ratio": True,
                               "y-axis-show-zero" : True,
                               "y-axis": "Events",
                               "log-y"  : "both",
                               "ratio-y-axis-range" : [0,2],
                               "save-extensions": ["pdf", "png"],
                               "ratio-y-axis" : '#frac{Data}{MC}',
                               "sort-by-yields" : True
                            }
        era = self.era  
        
        # Check era #
        if era != "2016preVFP" and era != "2016postVFP" and era != "2017" and era != "2018":
            raise RuntimeError("Unknown era {0}".format(era))
        #----- Triggers and Corrections -----#
        self.triggersPerPrimaryDataset = {}
        
        def addHLTPath(PD, HLT):
            if PD not in self.triggersPerPrimaryDataset.keys():
                self.triggersPerPrimaryDataset[PD] = []
            try:
                self.triggersPerPrimaryDataset[PD].append(
                    getattr(tree.HLT, HLT))
            except AttributeError:
                print("Couldn't find branch tree.HLT.%s, will omit it!" % HLT)
                
        #############################################################################
        #                           CutFlowReport                                   #
        #############################################################################
        self.yields = CutFlowReport("yields",printInLog=self.args.PrintYield,recursive=False)
        if self.args.PrintYield:
            self.yields.add(noSel, title="NoSelection")
        
        #############################################################################
        #                            TurnOff Syst                                   #
        #############################################################################
        if self.args.NoSystematics:
            noSel = noSel.refine('SystOff',autoSyst=False)
            if not self.args.PrintYield:
                self.yields.add(noSel, title="SystOff")
        #############################################################################
        #                            Clean DY                                       #
        #############################################################################
        #if sample == "DYJetsToLL_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8":
        #   print("0J")
        #   noSel = noSel.refine("NoSelCleaned", cut = tree.LHE.Njets == 0)
        #elif sample == "DYJetsToLL_1J_TuneCP5_13TeV-amcatnloFXFX-pythia8":
        #   print("1J")
        #   noSel = noSel.refine("NoSelCleaned", cut = tree.LHE.Njets == 1)
        #elif sample == "DYJetsToLL_1J_TuneCP5_13TeV-amcatnloFXFX-pythia8":
        #   print("2J")
        #   noSel = noSel.refine("NoSelCleaned", cut = tree.LHE.Njets == 2)
        #if self.args.PrintYield:
        #   self.yields.add(noSel, title="NoSelCleaned")

        if era == "2016postVFP" or era=="2016preVFP":
            # SingleMuon
            addHLTPath("SingleMuon","IsoMu24")
            addHLTPath("SingleMuon","IsoTkMu24")
        elif era == "2017":
            # SingleMuon #
            addHLTPath("SingleMuon","IsoMu27") 
        elif era == "2018":
            # SingleMuon #
            addHLTPath("SingleMuon","IsoMu24")
        
        #############################################################################
        #                             Pile-up                                       #
        #############################################################################
        # Get MC PU weight file #
        puWeightsFile = {"2016preVFP": "Collisions16_UltraLegacy_goldenJSON", 
                         "2016postVFP": "Collisions16_UltraLegacy_goldenJSON", 
                         "2017": "Collisions17_UltraLegacy_goldenJSON", 
                         "2018": "Collisions18_UltraLegacy_goldenJSON"}
        PUPATH = os.path.join(os.path.dirname(__file__), "POG","LUM",f'{era}_UL',"puWeights.json.gz")
        puReweight = get_correction(PUPATH, puWeightsFile[era],
                                     params={"NumTrueInteractions": lambda nTrueInt : nTrueInt},
                                     systParam="weights", systNomName="nominal", systName="pu", systVariations=("up", "down"),
                                     defineOnFirstUse=False,
                                     sel=noSel) 

        noSel = noSel.refine("NoSelwPU",weight = puReweight(tree.Pileup_nTrueInt) if self.is_MC else None)
        if self.args.PrintYield:
            self.yields.add(noSel,title="puWeight")
        #############################################################################
        #                             Gen Weight and Triggers                       #
        #############################################################################
        if self.is_MC: 
            noSel = noSel.refine('mcWeight', weight=tree.genWeight, cut=(
                op.OR(*chain.from_iterable(self.triggersPerPrimaryDataset.values()))))
        else:
            noSel = noSel.refine('trigger', cut=[makeMultiPrimaryDatasetTriggerSelection(
                sample, self.triggersPerPrimaryDataset)])
                 
        def hasAssociatedJet(lep): return lep.jet.idx != -1
        
        
        #############################################################################
        #                             Function definitions for selections           #
        #############################################################################
        # slimmedMuons after basic selection (pt > 3 && (passed('CutBasedIdLoose') || passed('SoftCutBasedId') || passed('SoftMvaId') || passed('CutBasedIdGlobalHighPt') || passed('CutBasedIdTrkHighPt')))
        def muonDef(muons):
            return op.select(muons, lambda mu : op.AND(
                    mu.pt > 20., #mu.pt > 20 in Hmm
                    op.abs(mu.eta) < 2.4,
                    mu.mediumId, # https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideMuonIdRun2#Medium_Muon
                    op.abs(mu.dxy) < 0.5, # no in Hmm
                    op.abs(mu.dz ) < 1. # no in Hmm
                    )                   
            )

        def muonDefIso(muons):
            return op.select(muons, lambda mu : op.AND(
                op.OR(op.AND(mu.fsrPhoton.idx == -1 ,op.abs(mu.pfRelIso04_all) < 0.15), # 0.25 in Hmm
                op.AND(mu.fsrPhoton.idx != -1,
                op.AND(op.deltaR(mu.p4,mu.fsrPhoton.p4)> 0.0001, op.deltaR(mu.p4,mu.fsrPhoton.p4)< 0.5), 
                mu.fsrPhoton.pt > 2., 
                op.OR(
                op.AND(op.abs(mu.fsrPhoton.eta)> 0.,op.abs(mu.fsrPhoton.eta)<1.4442),
                op.AND(op.abs(mu.fsrPhoton.eta)>1.566,op.abs(mu.fsrPhoton.eta)<2.5)
                ),
                mu.fsrPhoton.pt / mu.pt < 0.4,
                mu.fsrPhoton.relIso03 < 1.8,
                mu.fsrPhoton.dROverEt2 < 0.012,
                (op.abs(mu.pfRelIso04_all) - (mu.fsrPhoton.pt)/mu.pt) < 0.15))#mu.pfIsoId Tight, # https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideMuonIdRun2#Particle_Flow_isolation 
            )
            )
        def FSRCorrectedPt(muons):
            return op.map(muons, lambda mu: op.multiSwitch(
            (mu.fsrPhoton.idx == -1, mu.pt), # first possibility: no FSR photons
            (op.AND(
            mu.fsrPhoton.idx != -1,
            op.AND(op.deltaR(mu.p4,mu.fsrPhoton.p4)> 0.0001, op.deltaR(mu.p4,mu.fsrPhoton.p4)< 0.5), 
            mu.fsrPhoton.pt > 2., 
            op.OR(
            op.AND(op.abs(mu.fsrPhoton.eta)>0.,op.abs(mu.fsrPhoton.eta)<1.4442),
            op.AND(op.abs(mu.fsrPhoton.eta)>1.566,op.abs(mu.fsrPhoton.eta)<2.5)
            ),
            mu.fsrPhoton.pt / mu.pt < 0.4,
            mu.fsrPhoton.relIso03 < 1.8,
            mu.fsrPhoton.dROverEt2 < 0.012
            ),  
            mu.pt+(mu.fsrPhoton.pt)), # FSR photons with certain conditions, add the 
            mu.pt
        ))
            
        def elDef(electrons):
            return op.select(electrons, lambda ele: op.AND(
            ele.pt > 20., # Hmm > 20, 10 before
            op.abs(ele.eta) < 2.5,
            ele.mvaFall17V2noIso_WP90, #ele.mvaFall17V2noIso_WP80,# https://twiki.cern.ch/twiki/bin/view/CMS/MultivariateElectronIdentificationRun2#Differences_between_V1_and_V2
            op.OR( # no in Hmm
            #avoid transition region
            op.AND(op.abs(ele.eta) < 2.5,op.abs(ele.eta) > 1.57), # no in Hmm
            op.abs(ele.eta) < 1.44             # no in Hmm
            ),
            op.abs(ele.pfRelIso03_all) < 0.15, # no in Hmm
            op.abs(ele.dxy) < 0.5, # no in Hmm
            op.abs(ele.dz ) < 1. # no in Hmm

            )   
        )
            
        #############################################################################
        #                          Electron corrections                             #
        #############################################################################
        ELEPATH = os.path.join(os.path.dirname(__file__), "POG","EGM",f'{era}_UL',"electron.json.gz")
        eleId90SF = get_correction(ELEPATH, "UL-Electron-ID-SF",
                                     params={"year": era,"pt": lambda el : el.pt, "eta": lambda el : el.eta, 
                                             "WorkingPoint": "wp90noiso"},
                                     systParam="ValType", systNomName="sf", systName="eleId", 
                                     systVariations={"eleIdup": "sfup","eleIddown": "sfdown"},
                                     defineOnFirstUse=False,
                                     sel=noSel)
        eleId80SF = get_correction(ELEPATH, "UL-Electron-ID-SF",
                                     params={"year": era,"pt": lambda el : el.pt, "eta": lambda el : el.eta,
                                             "WorkingPoint": "wp80noiso"},
                                     systParam="ValType", systNomName="sf", systName="eleId", 
                                     systVariations={"eleIdup": "sfup","eleIddown": "sfdown"},
                                     defineOnFirstUse=False,
                                     sel=noSel)
        eleRecoBelowSF = get_correction(ELEPATH, "UL-Electron-ID-SF",
                                     params={"year": era,"pt": lambda el : el.pt, "eta": lambda el : el.eta,
                                             "WorkingPoint": "RecoBelow20"},
                                     systParam="ValType", systNomName="sf", systName="eleReco", 
                                     systVariations={"eleRecoup": "sfup","eleRecodown": "sfdown"},
                                     defineOnFirstUse=False,
                                     sel=noSel)
        eleRecoAboveSF = get_correction(ELEPATH, "UL-Electron-ID-SF",
                                     params={"year": era,"pt": lambda el : el.pt, "eta": lambda el : el.eta,
                                             "WorkingPoint": "RecoAbove20"},
                                     systParam="ValType", systNomName="sf", systName="eleReco", 
                                     systVariations={"eleRecoup": "sfup","eleRecodown": "sfdown"},
                                     defineOnFirstUse=False,
                                     sel=noSel)

        plots = []

        mZ = 91.1876
        mH = 125.1
        ################################################################
        # Spurious miss pT, detector noise and non-collision background#
        ################################################################
        # primary vertex filter ("Flag_goodVertices")	available in miniAOD	DONE	suggested	Work in progress, under construction	 
        # beam halo filter ("Flag_globalSuperTightHalo2016Filter")	available in miniAOD	DONE	suggested	Work in progress, under construction	Beam Halo Presentation
        # HBHE noise filter ("Flag_HBHENoiseFilter")	available in miniAOD	DONE	suggested	Work in progress, under construction	HCAL DPG Presentation
        # HBHEiso noise filter ("Flag_HBHENoiseIsoFilter")	available in miniAOD	DONE	suggested	Work in progress, under construction	same as above
        # ECAL TP filter ("Flag_EcalDeadCellTriggerPrimitiveFilter")	available in miniAOD	DONE	suggested	Work in progress, under construction	ECAL DPG Presentation
        # Bad PF Muon Filter ("Flag_BadPFMuonFilter")	available in miniAOD	DONE	suggested	Work in progress, under construction	PPD presentation, check with 2018 data
        # Bad PF Muon Dz Filter ("Flag_BadPFMuonDzFilter")	available in miniAODv2, v1*	DONE	suggested	Work in progress, under construction	New filter, not present EOY. See recipe below.
        # HF noisy hits filter ("Flag_hfNoisyHitsFilter")	available in miniAODv2	Yellow led	Yellow led	Work in progress, under construction	New filter for HF (optional). See note above.
        # Bad Charged Hadron Filter ("Flag_BadChargedCandidateFilter")	available in miniAOD	Stop	Stop	Work in progress, under construction	PPD presentation PPD update (TeV jet inefcy)
        # ee badSC noise filter ("Flag_eeBadScFilter")	available in miniAOD	DONE	suggested	Work in progress, under construction
        # ECAL bad calibration filter update ("Flag_ecalBadCalibFilter")
        # https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFiltersRun2#2018_data_EOY	
        ######################################
        #      Primary Vertex Cut            #
        ######################################
        isSpuriousBackground = noSel.refine("spuriousbkg", cut = op.AND(op.OR(era=="2016preVFP", era=="2016postVFP",era=="2018",era=="2017"),
                                             getattr(tree.Flag,"goodVertices"), 
                                             getattr(tree.Flag,"globalSuperTightHalo2016Filter"),
                                             getattr(tree.Flag,"HBHENoiseFilter"),
                                             getattr(tree.Flag,"HBHENoiseIsoFilter"),
                                             getattr(tree.Flag,"EcalDeadCellTriggerPrimitiveFilter"),
                                             getattr(tree.Flag,"BadPFMuonFilter"),
                                             getattr(tree.Flag,"BadPFMuonDzFilter"),
                                             getattr(tree.Flag,"hfNoisyHitsFilter"),
                                             getattr(tree.Flag,"eeBadScFilter"),
                                             getattr(tree.Flag,"ecalBadCalibFilter")
                                             )
                                            )
        hasGoodPrimaryVertices = isSpuriousBackground.refine("npvsGood1", cut = (tree.PV).npvsGood > 0) 
        # the one with the largest value of summed physics object pT^2 (PV.score) is taken to be the primary interaction vertex
        if self.args.PrintYield:
            self.yields.add(hasGoodPrimaryVertices, title="npvsGood1")
                
        ######################################
        #      3.2 Muons    Selection        #
        ######################################
        
        ######################################
        #      Medium Id and Eta,Pt requirement  #
        ######################################
        
        self.muonsNoIso = op.sort(muonDef(tree.Muon), lambda mu: -mu.pt)
        MuonPairsNoIso = op.combine(self.muonsNoIso, N=2, pred=lambda l1,l2 : l1.charge != l2.charge) #op.NOT(hasAssociatedJet(l1)),op.NOT(hasAssociatedJet(l2))))
        hasTwoOppositeChargeMuonsNoIso = hasGoodPrimaryVertices.refine("TwoOppositeChargeMuonsNoIso2", cut = op.rng_len(MuonPairsNoIso) > 0)
        if self.args.PrintYield:
            self.yields.add(hasTwoOppositeChargeMuonsNoIso, title="TwoOppositeChargeMuonsNoIso2")
        ########################################################################################
        #               I can have more than one pair of muons with those selections           #
        #Select the best Z asking that the invariant mass is the closest to mZ, same for Higgs #
        ########################################################################################
        bestZmm_noIso = op.rng_min_element_by(MuonPairsNoIso, lambda ll : op.abs(op.invariant_mass(ll[0].p4,ll[1].p4)-mZ))
        # Plots with only selections on muons (until medium Id step)
        if not self.args.BtagReweightingOff and not self.args.BtagReweightingOn:
            plots += self.plotPairs(bestZmm_noIso, hasTwoOppositeChargeMuonsNoIso, "TwoOppositeChargeMuonsNoIso2",self.muonsNoIso)
        
        
        #########################################################
        #      Isolation requirement                            #
        #########################################################
        #op.map(tree.Muon, lambda mu : print("Attributes in AltMuonProxy:",dir(mu.fsrPhoton)))
        self.muonsnoisocorrectedpt = FSRCorrectedPt(self.muonsNoIso)
        # Usage of FSR Corrected pt to order the reco muons
        self.muons = op.sort(muonDefIso(self.muonsNoIso), lambda mu: -(self.muonsnoisocorrectedpt[mu.idx]))
        MuonPairs_Iso = op.combine( self.muons, N=2, pred=lambda l1,l2 : l1.charge != l2.charge)
        hasTwoOppositeChargeMuonsIso = hasTwoOppositeChargeMuonsNoIso.refine("TwoOppositeChargeMuonsIso3", cut = op.rng_len(MuonPairs_Iso) == 1)
        if self.args.PrintYield:
            self.yields.add(hasTwoOppositeChargeMuonsIso, title="TwoOppositeChargeMuonsIso3")
        ########################################################################################
        # No more than one pair per event, no needed to select them                            #
        ########################################################################################
        bestZmm_Iso = op.rng_min_element_by(MuonPairs_Iso, lambda ll : op.abs(op.invariant_mass(ll[0].p4,ll[1].p4)) - mZ) # after TwoOppositeChargeMuonsIsoSelection, only one pair    
        if not self.args.BtagReweightingOff and not self.args.BtagReweightingOn:
            plots += self.plotPairs(bestZmm_Iso, hasTwoOppositeChargeMuonsIso, "TwoOppositeChargeMuonsIso3",self.muons)
        
        ######################################
        #      Trigger matching              #
        ######################################
        
        # At least one of the two muons is required to be geometrically matched to the level-3 muon
        # candidate corresponding to a single muon trigger that was fired in the event.
        # filterBit == 2 corresponds to iso for muons filterBit == 8 corresponds to isoTrk for muons 
        # https://indico.cern.ch/event/1155820/contributions/4899860/attachments/2454301/4206249/cipriani_2022June01_triggerNanoAOD.pdf
        ptThreshold = {"2016preVFP" : 26.,
                       "2016postVFP" : 26.,
                       "2017" : 29.,
                       "2018" : 26.
                       }
        ptTriggerThreshold = {"2016preVFP" : 24.,
                       "2016postVFP" : 24.,
                       "2017" : 27.,
                       "2018" : 24.
                       }
        ptL1TriggerThreshold = {"2016preVFP" : 22.,
                       "2016postVFP" : 22.,
                       "2017" : 25.,
                       "2018" : 22.
                       }
        ptL2TriggerThreshold = {"2016preVFP" : 10.,
                       "2016postVFP" : 10.,
                       "2017" : 12.,
                       "2018" : 10.
                       }
        triggerMuObj  = op.select(tree.TrigObj, lambda trgObj: op.AND( 
                        op.abs(trgObj.id) == 13,trgObj.pt > ptTriggerThreshold[era], 
                        trgObj.l1pt > ptL1TriggerThreshold[era],
                        op.OR(
                            op.AND(op.OR(trgObj.filterBits & 8,op.AND(trgObj.l2pt > ptL2TriggerThreshold[era],trgObj.filterBits & 2)),op.OR(era == "2016preVFP", era=="2016postVFP")),
                            op.AND(trgObj.l2pt > ptL2TriggerThreshold[era],trgObj.filterBits & 2)
                            )
                        )
                        )
        MuonPairs_tight = op.combine( self.muons, N=2, pred=lambda l1,l2 : op.AND(l1.charge != l2.charge,
                                                                                  op.deltaR(l1.p4,l2.p4)> 0.5,
                                                                                  op.deltaR(l1.p4,l2.p4) < 3.5,
                                                                                  op.abs(op.deltaPhi(l1.p4,l2.p4)) > 0.5,  
                                                                                  op.abs(l1.p4.Eta() - l2.p4.Eta()) < 1.5,
                                                                                  op.NOT(op.AND(op.abs(l1.p4.Eta()) > 2, op.abs(l1.p4.Eta()) < 2.5)),
                                                                                  op.NOT(op.AND(op.abs(l2.p4.Eta()) > 2, op.abs(l2.p4.Eta()) < 2.5)),    
                                                                                  op.OR(l1.tightId,l2.tightId),
                                                                                  op.OR(
                                                                                      op.AND(op.rng_any(triggerMuObj, lambda trgObj : op.deltaR(l1.p4,trgObj.p4)< 0.2), l1.pt > ptThreshold[era]),
                                                                                      op.AND(op.rng_any(triggerMuObj, lambda trgObj : op.deltaR(l2.p4,trgObj.p4)< 0.2), l2.pt > ptThreshold[era]),
                                                                                  
                                                                                  ),
                                                                                  #op.NOT(hasAssociatedJet(l1)),op.NOT(hasAssociatedJet(l2))
                                                                                  op.invariant_mass(l1.p4,l2.p4) > 40 # maybe use invariant mass with adjusted pt later?
                                                                                  ))#op.NOT(hasAssociatedJet(l1)),op.NOT(hasAssociatedJet(l2))
        hasMatchedTriggerObjectAndTightId = hasTwoOppositeChargeMuonsIso.refine("MatchedTriggerObjectAndTightId4", cut = op.rng_len(MuonPairs_tight) > 0)
        if self.args.PrintYield:
            self.yields.add(hasMatchedTriggerObjectAndTightId, title="MatchedTriggerObjectAndTightId4")
        bestZmm_tight = op.rng_min_element_by(MuonPairs_tight, lambda ll : op.abs(op.invariant_mass(ll[0].p4,ll[1].p4)-mZ))
        bestHmm_tight = op.rng_min_element_by(MuonPairs_tight, lambda ll : op.abs(op.invariant_mass(ll[0].p4,ll[1].p4)-mH))
        # Plots with only selections on muons
        if not self.args.BtagReweightingOff and not self.args.BtagReweightingOn:
            plots += self.plotPairs(bestZmm_tight, hasMatchedTriggerObjectAndTightId, "TriggerMatchingAndTightId4Zmm",self.muons)
            plots += self.plotPairs(bestHmm_tight, hasMatchedTriggerObjectAndTightId, "TriggerMatchingAndTightId4Hmm",self.muons)
        
        
        hasOneTightPairPerEvent = hasMatchedTriggerObjectAndTightId.refine("OneTightPairPerEvent5", cut = op.rng_len(MuonPairs_tight) == 1)
        if self.args.PrintYield:
            self.yields.add(hasOneTightPairPerEvent, title="OneTightPairPerEvent5")
        bestZmm_onetightpair = op.rng_min_element_by(MuonPairs_tight, lambda ll : op.abs(op.invariant_mass(ll[0].p4,ll[1].p4)-mZ))
        bestHmm_onetightpair = op.rng_min_element_by(MuonPairs_tight, lambda ll : op.abs(op.invariant_mass(ll[0].p4,ll[1].p4)-mH))
        # Plots with only selections on muons
        if not self.args.BtagReweightingOff and not self.args.BtagReweightingOn:
            plots += self.plotPairs(bestZmm_onetightpair, hasOneTightPairPerEvent, "OneTightPairPerEvent5ZmmCorrectedFSR",self.muons)
            plots += self.plotPairs(bestHmm_onetightpair, hasOneTightPairPerEvent, "OneTightPairPerEvent5HmmCorrectedFSR",self.muons)
            plots += self.plotVertices((tree.PV).npvsGood,hasOneTightPairPerEvent,"OneTightPairPerEvent5Vertices")
        #############################################################################
        #                           muon corrections                                #
        #############################################################################

        MUONPATH = os.path.join(os.path.dirname(__file__), "POG","MUO",f'{era}_UL',"muon_Z_v2.json.gz")
        muRecoSF = get_correction(MUONPATH, "NUM_TrackerMuons_DEN_genTracks", #SF are computed on Z->mumu events in the range 40 < pT < 60 GeV. The recommendation is to apply them for muons with pT in the range [10; 200] GeV. Below 10 GeV, the SF computed using Jpsi->mumu events should be used
                                     params={"pt": lambda mu : mu.pt, "abseta": lambda mu : op.abs(mu.eta)},
                                     systParam="scale_factors", systNomName="nominal", systName="muReco", 
                                     systVariations= {"muRecoup": "systup","muRecodown": "systdown"},
                                     defineOnFirstUse=False,
                                     sel=hasOneTightPairPerEvent)
       
        muIdSF = get_correction(MUONPATH, "NUM_MediumID_DEN_TrackerMuons", # computed with Dxy < 0.2 cm and Dz < 0.5, If you are using medium or loose ID in your analysis, those Dxy and Dz are not part of the selectors. In principle it does not introduce bias.
                                     params={"year": era+"_UL","pt": lambda mu : mu.pt, "abseta": lambda mu : op.abs(mu.eta)},
                                     systParam="ValType", systNomName="sf", systName="muId",
                                     systVariations={"muIdup": "systup","muIddown": "systdown"},
                                     defineOnFirstUse=False,
                                     sel=hasOneTightPairPerEvent)
        
        muIsoSF = get_correction(MUONPATH, "NUM_TightRelIso_DEN_MediumID",
                                     params={"year": era+"_UL","pt": lambda mu : mu.pt, "abseta": lambda mu : op.abs(mu.eta)},
                                     systParam="ValType", systNomName="sf", systName="muIso", 
                                     systVariations={"muIsoup": "systup","muIsodown": "systdown"},
                                     defineOnFirstUse=False,
                                     sel=hasOneTightPairPerEvent)
        if era == "2017":
            trigger_string = "NUM_IsoMu27_DEN_CutBasedIdTight_and_PFIsoTight"
        if era == "2018":
            trigger_string = "NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight"
        if era == "2016preVFP" or era =="2016postVFP":
            trigger_string = "NUM_IsoMu24_or_IsoTkMu24_DEN_CutBasedIdTight_and_PFIsoTight"
        muTriggerSF = get_correction(MUONPATH, trigger_string, #take it from https://gitlab.cern.ch/cms-muonPOG/muonefficiencies/-/blob/master/Run2/UL/2018/2018_trigger/Efficiencies_muon_generalTracks_Z_Run2018_UL_SingleMuonTriggers_schemaV2.json
                                     params={"year": era+"_UL","pt": lambda mu : mu.pt, "abseta": lambda mu : op.abs(mu.eta)},
                                     systParam="ValType", systNomName="sf", systName="muTrig", 
                                     systVariations= {"muTrigup": "systup","muTrigdown": "systdown"},
                                     defineOnFirstUse=False,
                                     sel=hasOneTightPairPerEvent)
        #############################################################################################################
        #   I have already one pair per event!                                                                      #
        #############################################################################################################
        bestZmm = op.rng_min_element_by(MuonPairs_tight, lambda ll : op.abs(op.invariant_mass(ll[0].p4,ll[1].p4)-mZ))
        bestHmm = op.rng_min_element_by(MuonPairs_tight, lambda ll : op.abs(op.invariant_mass(ll[0].p4,ll[1].p4)-mH))
        if bestHmm[0].pt > ptThreshold[era]:
            weightMuSF_1 = muRecoSF(bestHmm[0])* muTriggerSF(bestHmm[0]) * muIdSF(bestHmm[0]) * muIsoSF(bestHmm[0])
        if bestHmm[0].pt > 15. and bestHmm[0].pt <= ptThreshold[era]:
            weightMuSF_1 = muRecoSF(bestHmm[0])* muIdSF(bestHmm[0]) * muIsoSF(bestHmm[0])
        else:
            weightMuSF_1 = op.c_float(1)
        if bestHmm[1].pt > ptThreshold[era]:
            weightMuSF_2 = muRecoSF(bestHmm[1])* muTriggerSF(bestHmm[1]) * muIdSF(bestHmm[1]) * muIsoSF(bestHmm[1])
        if bestHmm[1].pt > 15. and bestHmm[1].pt <= ptThreshold[era]:
            weightMuSF_2 = muRecoSF(bestHmm[1])* muIdSF(bestHmm[1]) * muIsoSF(bestHmm[1])
        else:
            weightMuSF_2 = op.c_float(1)    
        hasMuonSFs = hasOneTightPairPerEvent.refine("MuonSFs51", cut = op.rng_len(MuonPairs_tight) == 1, weight = weightMuSF_1 * weightMuSF_2 if self.is_MC else None)
        if not self.args.BtagReweightingOff and not self.args.BtagReweightingOn:
            if self.args.PrintYield:
                self.yields.add(hasMuonSFs, title="MuonSFs51")
            # Plots with only selections on muons
            plots += self.plotPairs(bestHmm, hasMuonSFs, "MuonSFs51",self.muons)

        ######################################
        #                   FSR   Plots      #
        ######################################
        
        self.correctedPtMuontight = FSRCorrectedPt(self.muons)

        def GeoFitCorrection(bestHmm, idx,era):
            eta1 = {"2016preVFP": 411.34,
                    "2016postVFP": 411.34,
                    "2017": 582.32,
                    "2018": 650.84
                    }
            eta2 = {
                "2016preVFP": 673.40,
                "2016postVFP": 673.40,
                "2017": 974.05,
                "2018": 988.37 
            }
            eta3 = {
                "2016preVFP": 1099.0,
                "2016postVFP": 1099.0,
                "2017": 1263.4,
                "2018": 1484.6 
            }
            if  op.abs(bestHmm[idx].p4.Eta()) < 0.9: 
                pt_cor = op.c_float(op.c_float(eta1[era]) * op.c_float(bestHmm[idx].charge) * op.c_float(bestHmm[idx].dxybs) * op.c_float(bestHmm[idx].p4.Pt()) * op.c_float(bestHmm[idx].p4.Pt()) / 10000.0)
            if  op.AND(op.abs(bestHmm[idx].p4.Eta()) > 0.9, op.abs(bestHmm[idx].p4.Eta()) < 1.7):
                pt_cor = op.c_float(op.c_float(eta2[era]) * op.c_float(bestHmm[idx].charge) * op.c_float(bestHmm[idx].dxybs) * op.c_float(bestHmm[idx].p4.Pt()) * op.c_float(bestHmm[idx].p4.Pt()) / 10000.0)
            if op.abs(bestHmm[idx].p4.Eta()) > 1.7:
                pt_cor = op.c_float(op.c_float(eta3[era]) * op.c_float(bestHmm[idx].charge) * op.c_float(bestHmm[idx].dxybs) * op.c_float(bestHmm[idx].p4.Pt()) * op.c_float(bestHmm[idx].p4.Pt()) / 10000.0)
            return pt_cor
        # Pcorrected = ptRoc - coeff * d0 * pt^2 / 10000
        # Ptroc - ptcorr = coeff * d0 * pt^2/100000
        # delta P = ptRoch - ptCorr 
        # D0 = coeffic * deltaP / pt^2 
        # D0 = coeff * (Porig - Pcorr) / pOrig^2
        # 
        pt_corr1 = op.c_float(GeoFitCorrection(bestHmm,0,era))
        pt_corr2 = op.c_float(GeoFitCorrection(bestHmm,1,era))
        #lepton1 = op.construct("ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float> >",([op.c_float(op.c_float(0.)),op.c_float(0.),op.c_float(0.),op.c_float(0.)]))
        #lepton2 = op.construct("ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float> >",([op.c_float(op.c_float(0.)),op.c_float(0.),op.c_float(0.),op.c_float(0.)]))
        lepton1 = op.construct("ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float> >",([op.c_float(bestHmm[0].p4.Pt() - pt_corr1),op.c_float(bestHmm[0].p4.Eta()),op.c_float(bestHmm[0].p4.Phi()),op.c_float(bestHmm[0].p4.M())]))
        lepton2 = op.construct("ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float> >",([op.c_float(bestHmm[1].p4.Pt() - pt_corr2),op.c_float(bestHmm[1].p4.Eta()),op.c_float(bestHmm[1].p4.Phi()),op.c_float(bestHmm[1].p4.M())]))
        
        lepton1_FSR = op.construct("ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float> >",([op.c_float(op.c_float(self.correctedPtMuontight[bestHmm[0].idx]) - pt_corr1),op.c_float(bestHmm[0].p4.Eta()),op.c_float(bestHmm[0].p4.Phi()),op.c_float(bestHmm[0].p4.M())]))
        lepton2_FSR = op.construct("ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float> >",([op.c_float(op.c_float(self.correctedPtMuontight[bestHmm[1].idx]) - pt_corr2),op.c_float(bestHmm[1].p4.Eta()),op.c_float(bestHmm[1].p4.Phi()),op.c_float(bestHmm[1].p4.M())]))
        if not self.args.BtagReweightingOff and not self.args.BtagReweightingOn:
            plots += self.plotPairsCorrected(lepton1, lepton2 , hasMuonSFs, "MuonSFsGeoFit5Hmm2")
            plots += self.plotPairsCorrected(lepton1_FSR, lepton2_FSR , hasMuonSFs, "MuonSFsGeoFitFSR5Hmm2")
        
        ######################################
        #      3.8 Electron    Selection     #
        ######################################
        self.electrons = op.sort(elDef(tree.Electron), lambda el: -el.pt)
        ######################################
        #      3.4 Selection of AK4 jets     # (Jet_puId & (1 << 2))
        ######################################
        # slimmedJets, i.e. ak4 PFJets CHS with JECs applied, after basic selection (pt > 15)
        self.ak4JetsByPtNoJetPuID = op.sort(op.select(tree.Jet, lambda j : op.AND(
            j.pt > 25., #25 in Hmm
            op.abs(j.eta) < 2.5 if era == "2016preVFP" or era == "2016postVFP" else op.abs(j.eta) < 2.4, # 4.7 for Hmm jets beyond |eta| 2.5 (2.4 in 2016) cannot be b-tagged. Hence you should not include them in your calculation of b-tag event weights, and you should not use their b-tagging status/score in any MVA down the line.
            j.jetId == 6 if era == "2016preVFP" or era == "2016postVFP" else j.jetId == 6, # bit2 is tight https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID
            #op.OR((j.puId & (1 << 0)) if era == "2016preVFP" or era == "2016postVFP" else (j.puId & (1 << 2)), j.pt > 50.), # Jet PU ID bit1 is loose from JME POG 
            #op.NOT(op.rng_any(self.muons, lambda im : op.deltaR(j.p4, im.p4) < 0.4)),
            #op.NOT(op.rng_any(self.electrons, lambda im : op.deltaR(j.p4, im.p4) < 0.4)),
            #op.NOT(op.rng_any(self.muons,lambda mu: op.AND(mu.fsrPhoton.idx != -1,op.deltaR(j.p4,mu.fsrPhoton.p4)<0.4)))
            )), lambda j : -j.pt)
        self.ak4JetsByPt = op.sort(op.select(tree.Jet, lambda j : op.AND(
            j.pt > 25.,
            op.abs(j.eta) < 2.5 if era == "2016preVFP" or era == "2016postVFP" else op.abs(j.eta) < 2.4, # 4.7 for Hmm jets beyond |eta| 2.5 (2.4 in 2016) cannot be b-tagged. Hence you should not include them in your calculation of b-tag event weights, and you should not use their b-tagging status/score in any MVA down the line.
            j.jetId == 6 if era == "2016preVFP" or era == "2016postVFP" else j.jetId == 6, # bit2 is tight, 6 is tight+LepVeto https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID
            #op.OR((j.puId & (1 << 0)) if era == "2016preVFP" or era == "2016postVFP" else (j.puId & (1 << 2)), j.pt > 50.), # Jet PU ID bit1 is loose from JME POG 
            op.OR((j.puId & (1 << 0)) if era == "2016preVFP" or era == "2016postVFP" else (j.puId & (1 << 2)), j.pt > 50.), # Jet PU ID bit1 is loose from JME POG 
            #op.NOT(op.rng_any(self.muons, lambda im : op.deltaR(j.p4, im.p4) < 0.4)),
            #op.NOT(op.rng_any(self.electrons, lambda im : op.deltaR(j.p4, im.p4) < 0.4)),
            #op.NOT(op.rng_any(self.muons,lambda mu: op.AND(mu.fsrPhoton.idx != -1,op.deltaR(j.p4,mu.fsrPhoton.p4)<0.4)))
            )), lambda j : -j.pt)
        self.ak4JetsByBtagScore = op.sort(self.ak4JetsByPt,lambda j: -j.btagDeepFlavB) # Btag score ordered
        
        hasNoElectrons = hasMuonSFs.refine("VetoElectrons6", cut = op.rng_len(self.electrons) < 1)
        if not self.args.BtagReweightingOff and not self.args.BtagReweightingOn:
            if self.args.PrintYield:
                self.yields.add(hasNoElectrons, title="NoElectrons6")
            # Plots with only selections on muons
            plots += self.plotPairs(bestHmm, hasNoElectrons, "VetoElectrons6Hmm",self.muons)
            plots += self.plotPairsCorrected(lepton1, lepton2 , hasNoElectrons, "VetoElectronsGeoFit6Hmm")
            plots += self.plotPairsCorrected(lepton1_FSR, lepton2_FSR , hasNoElectrons, "VetoElectronsGeoFitFSR6Hmm")
        
        
        cleanedHighestPtJets= self.ak4JetsByPt[:op.min(op.rng_len(self.ak4JetsByPt),op.static_cast("std::size_t",op.c_int(2)))] 
        cleanedHighestBtagScoreJets = self.ak4JetsByBtagScore[:op.min(op.rng_len(self.ak4JetsByBtagScore),op.static_cast("std::size_t",op.c_int(2)))]
        if not self.args.BtagReweightingOff and not self.args.BtagReweightingOn:
            plots += self.plotJets(self.ak4JetsByPt, hasNoElectrons, "VetoElectrons6PtOrdering")
            plots += self.plotJets(self.ak4JetsByBtagScore, hasNoElectrons, "VetoElectronsBScoreOrdering")
            plots += self.plotFourObjects(self.ak4JetsByBtagScore,bestHmm,hasNoElectrons,"VetoElectronsFourObjectsBScoreOrdering")

        
        #op.map(tree.Jet, lambda jet : print("Attributes in AltJetProxy:",dir(jet)))
        JETPUIDPATH = os.path.join(os.path.dirname(__file__), "POG","JME",f'{era}_UL',"jmar.json.gz")
        JetPUid_SF_eff = get_correction(JETPUIDPATH, "PUJetID_eff", 
                                     params={"pt": lambda jet : jet.pt, "eta": lambda jet : jet.eta, "workingpoint": "L"},
                                     systParam="systematic", systNomName="nom", systName="jetpuId",
                                     systVariations={"jetpuIdup": "up","jetpuIddown": "down"},
                                     defineOnFirstUse=False,
                                     sel=hasNoElectrons)
        JetPUid_MC_eff = get_correction(JETPUIDPATH, "PUJetID_eff", 
                                     params={"pt": lambda jet : jet.pt, "eta": lambda jet : jet.eta, "workingpoint": "L"},
                                     systParam="systematic", systNomName="MCEff", systName="jetpuIdMCEff",
                                     systVariations={"jetpuIdup": "up","jetpuIddown": "down"},
                                     defineOnFirstUse=False,
                                     sel=hasNoElectrons)
        self.ak4ForPUID = op.select(self.ak4JetsByPtNoJetPuID, lambda j : j.pt<=50)
        
        if self.is_MC and not self.args.BtagReweightingOff and not self.args.BtagReweightingOn:
            # method 1.a : https://twiki.cern.ch/twiki/bin/view/CMS/BTagSFMethods#1a_Event_reweighting_using_scale
            wFail = op.extMethod("scalefactorWeightForFailingObject", returnType="double")
            # double scalefactorWeightForFailingObject(double sf, double eff) {
            #  return (1.-sf*eff)/(1.-eff);
            #  }
            #def makePuIDWeight(jets):    
            #    def JetPUIdCorrection(jet):
            #        return op.switch(op.AND(jet.pt > 12.5, jet.pt < 50., op.rng_any(tree.GenJet, lambda genjet : op.deltaR(jet.p4,genjet.p4)< 0.4)), 
            #             JetPUid_SF_eff(jet), 
            #             op.c_float(1))
            #    return op.rng_product(jets, JetPUIdCorrection)
            lambda_puid_weight = lambda jet : op.switch(op.rng_any(tree.GenJet, lambda genjet : op.deltaR(jet.p4,genjet.p4)< 0.4), # Is matching with gen jet
                                                      op.switch((jet.puId & (1 << 0)) if era == "2016preVFP" or era == "2016postVFP" else (jet.puId & (1 << 2)), # passes jet pu id cut
                                                                JetPUid_SF_eff(jet), # SF_eff
                                                                wFail(JetPUid_SF_eff(jet), JetPUid_MC_eff(jet))), # (1-SF_eff*eff) / (1-eff)
                                                      # Is PU jet : no matching with gen jets
                                                      #op.switch((j.puId & (1 << 0)) if era == "2016preVFP" or era == "2016postVFP" else (j.puId & (1 << 2)), # passed jet pu id cut
                                                      #          self.jetpuid_sf_mis(j), # SF_mis
                                                      #          wFail(self.jetpuid_sf_mis(j), self.jetpuid_mc_mis(j)))) # (1-SF_mis*mis) / (1-mis)
                                                    op.c_float(1))
            weightJetPUID_all = op.rng_product(self.ak4ForPUID, lambda j : lambda_puid_weight(j))
            #weightJetPUID_all = makePuIDWeight(self.ak4JetsByPt)
        else:
            weightJetPUID_all = op.c_float(1)
        hasJetPUid_corr = hasNoElectrons.refine("JetPUid_corr", weight = weightJetPUID_all if self.is_MC else None)
        if not self.args.BtagReweightingOff and not self.args.BtagReweightingOn:      
            plots += self.plotPairs(bestHmm, hasJetPUid_corr, "JetPUidCorrection7Hmm",self.muons)
            plots += self.plotPairsCorrected(lepton1, lepton2 , hasJetPUid_corr, "JetPUidCorrectionGeoFit7Hmm")
            plots += self.plotPairsCorrected(lepton1_FSR, lepton2_FSR , hasJetPUid_corr, "JetPUidCorrectionGeoFitFSR7Hmm")
            plots += self.plotJets(self.ak4JetsByPt, hasJetPUid_corr, "JetPUidCorrection7PtOrdering")
            plots += self.plotJets(self.ak4JetsByBtagScore, hasJetPUid_corr, "JetPUidCorrectionBScoreOrdering")
            plots += self.plotFourObjects(self.ak4JetsByBtagScore,bestHmm,hasJetPUid_corr,"JetPUidCorrectionFourObjectsBScoreOrdering")
            plots += self.plotFourObjects(self.ak4JetsByPt,bestHmm,hasJetPUid_corr,"JetPUidCorrectionFourObjectsPtOrdering")
            if self.args.PrintYield:
                self.yields.add(hasJetPUid_corr, title="hasJetPUid_corr")
        def makeBtagRatioReweighting(jsonFile, numJets, systName=None, nameHint=None):
            """ Construct a btag ratio for MC, based on the weights in a JSON file

            :param btagRatioFile: path of the JSON file with weights (binned in NumJets)
            :param numJets : expression to get the number of selected jets
            :param systName: name of the associated systematic nuisance parameter
            """
            paramVType = "Parameters::value_type::value_type"
            args = op.construct("Parameters", (op.initList("std::initializer_list<{0}>".format(paramVType), paramVType,
                (op.initList(paramVType, "float", (op.extVar("int", "BinningVariable::NumJets"), numJets)),)),))
            wFun = op.define("ILeptonScaleFactor", 'const ScaleFactor <<name>>{{"{0}"}};'.format(jsonFile), nameHint=nameHint)
            expr = wFun.get(args, op.extVar("int", "Nominal"))
            if systName:
                expr = op.systematic(expr, name=systName,
                        up  =wFun.get(args, op.extVar("int", "Up")),
                        down=wFun.get(args, op.extVar("int", "Down")))
            return op.defineOnFirstUse(expr)
        BTVPATH = os.path.join(os.path.dirname(__file__), "POG","BTV",f'{era}_UL',"btagging.json.gz")
        # btagging SF simplest method https://twiki.cern.ch/twiki/bin/viewauth/CMS/BTagShapeCalibration https://twiki.cern.ch/twiki/bin/view/CMS/JetMET
        def btvSF(flav): return get_bTagSF_itFit(BTVPATH, "deepJet", "btagDeepFlavB", flav, hasJetPUid_corr,defineOnFirstUse=False)
        if self.is_MC:
            btvWeight = makeBtagWeightItFit(self.ak4JetsByPt, btvSF)
        else:
            btvWeight = op.c_float(1)
        if self.args.BtagReweightingOn and self.args.BtagReweightingOff: 
            raise RuntimeError("Reweighting cannot be both on and off") 
        if self.args.BtagReweightingOn:
            hasDeepJetShape_corr = hasJetPUid_corr.refine("DeepJetShape_corrReweightOn", weight=btvWeight if self.is_MC else None)
            plots += self.plotJets(self.ak4JetsByPt, hasDeepJetShape_corr, "NoChannel_NoSelection_Ak4Jets_On_N")
            if self.args.PrintYield:
                self.yields.add(hasDeepJetShape_corr, title="hasDeepJetShape_corrReweightOn")
        elif self.args.BtagReweightingOff:
            hasDeepJetShape_corr = hasJetPUid_corr.refine("DeepJetShape_corrReweightOff", weight=op.c_float(1) if self.is_MC else None)
            plots += self.plotJets(self.ak4JetsByPt, hasDeepJetShape_corr, "NoChannel_NoSelection_Ak4Jets_Off_N")
            if self.args.PrintYield:
                self.yields.add(hasDeepJetShape_corr, title="hasDeepJetShape_corrReweightOff")
        else:
            if self.is_MC:
                ReweightingFileName = os.path.join(os.path.dirname(os.path.abspath(__file__)),'Reweighting_plots', f'BtagReweightingRatio_jetN_{self.sample}_{self.era}.json')
                import json
                with open(ReweightingFileName, 'r') as file:
                    data = json.load(file)
                # Check if 'binning' or 'data' is empty
                if not data.get('binning', {}).get('x') or not data.get('data'):
                    BtagRatioWeight = op.c_float(1)
                else:
                    nameHint = f'bamboo_nJetsWeight_{self.sample}'.replace('-','_')
                    if not os.path.exists(ReweightingFileName):
                        raise RuntimeError("Could not find reweighting file %s"%ReweightingFileName)
                    print('Reweighting file', ReweightingFileName)
                    BtagRatioWeight = makeBtagRatioReweighting(jsonFile=ReweightingFileName,numJets=op.rng_len(self.ak4JetsByPt),nameHint=nameHint)
            else:
                BtagRatioWeight = op.c_float(1)
            hasDeepJetShape_corr = hasJetPUid_corr.refine("DeepJetShape_corrBtagRatioWeight", cut=op.AND(op.rng_len(self.ak4JetsByPt)>=0),weight=[btvWeight,BtagRatioWeight] if self.is_MC else None)
            plots += self.plotPairs(bestHmm, hasDeepJetShape_corr, "hasDeepJetShapeCorrBtagRatioWeight8Hmm",self.muons)
            plots += self.plotPairsCorrected(lepton1, lepton2 , hasDeepJetShape_corr, "hasDeepJetShapecorrGeoFitBtagRatioWeight8Hmm")
            plots += self.plotPairsCorrected(lepton1_FSR, lepton2_FSR , hasDeepJetShape_corr, "hasDeepJetShapecorrGeoFitFSRBtagRatioWeight8Hmm")
            plots += self.plotJets(self.ak4JetsByBtagScore, hasDeepJetShape_corr, "hasDeepJetShapecorrBtagRatioWeightBScoreOrdering")
            plots += self.plotJets(self.ak4JetsByPt, hasDeepJetShape_corr, "hasDeepJetShapecorr8BtagRatioWeightPtOrdering")

        if self.is_MC and (sample.startswith("TT") or sample.startswith("tt")) and not self.args.BtagReweightingOff and not self.args.BtagReweightingOn:
            def top_pt_weight(pt):
                return op.exp(-2.02274e-01 + 1.09734e-04*pt + -1.30088e-07*pt**2 + (5.83494e+01/(pt+1.96252e+02)))

            def getTopPtWeight(tree):
                lastCopy = op.select(
                    tree.GenPart, lambda p: (op.static_cast("int", p.statusFlags) >> 13) & 1)
                tops = op.select(lastCopy, lambda p: p.pdgId == 6)
                antitops = op.select(lastCopy, lambda p: p.pdgId == -6)
                weight = op.switch(op.AND(op.rng_len(tops) >= 1, op.rng_len(antitops) >= 1),
                                   op.sqrt(top_pt_weight(
                                       tops[0].pt) * top_pt_weight(antitops[0].pt)),
                                   1.)
                return weight

            logger.info(
                "Applying Top Pt reweighting (only for TTbar samples)")

            hasDeepJetShape_corr = hasDeepJetShape_corr.refine("topPt", weight=op.systematic(
                getTopPtWeight(tree), noTopPt=op.c_float(1.)))
        else:
            hasDeepJetShape_corr = hasDeepJetShape_corr.refine("topPt", weight=op.c_float(1.))
        if self.args.PrintYield:
            self.yields.add(hasDeepJetShape_corr, "topPt reweighting")
        
        #hasDeepJetShape_corr = hasDeepJetShape_corr.refine("DeepJetShapeCut2Jets_corrBtagRatioWeight", cut=op.AND(op.rng_len(self.ak4JetsByPt)>=2,self.ak4JetsByBtagScore[0].btagDeepFlavB>=0.,self.ak4JetsByBtagScore[1].btagDeepFlavB>=0.))
        if not self.args.BtagReweightingOff and not self.args.BtagReweightingOn:   
            plots += self.plotPairs(bestHmm, hasDeepJetShape_corr, "hasDeepJetShapeTopReweightCorrBtagRatioWeight8Hmm",self.muons)
            plots += self.plotPairsCorrected(lepton1, lepton2 , hasDeepJetShape_corr, "hasDeepJetShapeTopReweightcorrGeoFitBtagRatioWeight8Hmm")
            plots += self.plotPairsCorrected(lepton1_FSR, lepton2_FSR , hasDeepJetShape_corr, "hasDeepJetShapeTopReweightcorrGeoFitFSRBtagRatioWeight8Hmm")
            plots += self.plotJets(self.ak4JetsByBtagScore, hasDeepJetShape_corr, "hasDeepJetShapecorrTopReweightBtagRatioWeightBScoreOrdering")
            plots += self.plotJets(self.ak4JetsByPt, hasDeepJetShape_corr, "hasDeepJetShapecorrTopReweight8BtagRatioWeightPtOrdering")
            plots += self.plotFourObjects(self.ak4JetsByPt,bestHmm,hasDeepJetShape_corr,"hasDeepJetShapecorrTopReweightFourObjectsPtOrdering")
            plots += self.plotFourObjectsCorrected(self.ak4JetsByPt,lepton1_FSR,lepton2_FSR,hasDeepJetShape_corr,"hasDeepJetShapecorrTopReweightPtOrderingFSRGeoFit")
            
        a = op.c_float(-99, typeName='double', cast=None)
        phi_mpi_pi  = op.extMethod("ROOT::Math::VectorUtil::Phi_mpi_pi")
        
        def outZ(dilep): # starting from a single selected pair 
            return op.NOT(op.AND((op.invariant_mass(dilep[0].p4, dilep[1].p4) > 70.),op.invariant_mass(dilep[0].p4, dilep[1].p4) < 110.))
        
        def outH(dilep): # starting from a single selected pair 
            #print("Invariant mass:",op.c_long(op.abs(op.invariant_mass(dilep[0].p4, dilep[1].p4)),cast=True) )
            return op.NOT(op.AND((op.invariant_mass(dilep[0].p4, dilep[1].p4) > 115.),op.invariant_mass(dilep[0].p4, dilep[1].p4) < 135.))
        
        def outZLeptons(l1,l2): # starting from a single selected pair 
            return op.NOT(op.AND((op.invariant_mass(l1, l2) > 70.),op.invariant_mass(l1, l2) < 110.))
        
        def outHLeptons(l1,l2): # starting from a single selected pair 
            #print("Invariant mass:",op.c_long(op.abs(op.invariant_mass(dilep[0].p4, dilep[1].p4)),cast=True) )
            return op.NOT(op.AND((op.invariant_mass(l1, l2) > 115.),op.invariant_mass(l1, l2) < 135.))
        
        def outZMass(dilepMass): # starting from a single selected pair 
            return not (dilepMass > 70. and dilepMass < 110.)
        
        def ZMass(dilepMass): # starting from a single selected pair 
            return not (dilepMass < 70. or dilepMass > 110.)
        
        def outHMass(dilepMass): # starting from a single selected pair 
            #print("Invariant mass:",op.c_long(op.abs(op.invariant_mass(dilep[0].p4, dilep[1].p4)),cast=True) )
            return not (dilepMass > 115. and dilepMass < 135.)
        
        def HMass(dilepMass): # starting from a single selected pair 
            #print("Invariant mass:",op.c_long(op.abs(op.invariant_mass(dilep[0].p4, dilep[1].p4)),cast=True) )
            return not (dilepMass < 115. or dilepMass > 135.)
        
        def range_mass_for_bjets(jets):
            return op.AND((op.invariant_mass(jets[0].p4, jets[1].p4) > 70.),op.invariant_mass(jets[0].p4, jets[1].p4) < 140.)
        
        def range_mass_for_bjets_tight(jets):
            return op.AND((op.invariant_mass(jets[0].p4, jets[1].p4) > 115.),op.invariant_mass(jets[0].p4, jets[1].p4) < 135.)
        
        # 
        # Ask to have a good Z with the refine option
        # Higgs/ Z to Muons control region
        hasZmm            = hasDeepJetShape_corr.refine("hasZmm",  cut=[op.NOT(outZLeptons(lepton1_FSR,lepton2_FSR)),
                                                         op.NOT(range_mass_for_bjets(cleanedHighestBtagScoreJets))])
        hasHmm            = hasDeepJetShape_corr.refine("hasHmm",  cut=[op.NOT(outHLeptons(lepton1_FSR,lepton2_FSR)),
                                                         op.NOT(range_mass_for_bjets(cleanedHighestBtagScoreJets))])
        hasZjets          = hasDeepJetShape_corr.refine("hasZjets", cut=[op.NOT(op.NOT(outZLeptons(lepton1_FSR,lepton2_FSR))),
                                                         range_mass_for_bjets(cleanedHighestBtagScoreJets)])
        hasHjets          = hasDeepJetShape_corr.refine("hasHjets", cut=[op.NOT(op.NOT(outHLeptons(lepton1_FSR,lepton2_FSR))),
                                                        range_mass_for_bjets(cleanedHighestBtagScoreJets)])
        isSideband2lbbZmm = hasDeepJetShape_corr.refine("isSideband2lbbZmm", cut=[op.NOT(op.NOT(outZLeptons(lepton1_FSR,lepton2_FSR))),
                                                        op.NOT(range_mass_for_bjets(cleanedHighestBtagScoreJets))])
        isSideband2lbbHmm = hasDeepJetShape_corr.refine("isSideband2lbbHmm",  cut=[op.NOT(op.NOT(outHLeptons(lepton1_FSR,lepton2_FSR))),
                                                       op.NOT(range_mass_for_bjets(cleanedHighestBtagScoreJets))])
        isSignalRegion    = hasDeepJetShape_corr.refine("isSignalRegion", cut=[op.NOT(outHLeptons(lepton1_FSR,lepton2_FSR)),
                                                        range_mass_for_bjets(cleanedHighestBtagScoreJets)])
        hasZCRDY            = hasDeepJetShape_corr.refine("hasZCRDY",  cut=[op.NOT(outZLeptons(lepton1,lepton2))])
        #hasHmmTight            = hasDeepJetShape_corr.refine("hasHmmTight",  cut=[op.NOT(outHLeptons(lepton1,lepton2)),
        #                                                 op.NOT(range_mass_for_bjets_tight(cleanedHighestBtagScoreJets))])
        #hasZjetsTight          = hasDeepJetShape_corr.refine("hasZjetsTight", cut=[op.NOT(op.NOT(outZLeptons(lepton1,lepton2))),
        #                                                 range_mass_for_bjets_tight(cleanedHighestBtagScoreJets)])
        #hasHjetsTight          = hasDeepJetShape_corr.refine("hasHjetsTight", cut=[op.NOT(op.NOT(outHLeptons(lepton1,lepton2))),
        #                                           range_mass_for_bjets_tight(cleanedHighestBtagScoreJets)])
        #isSideband2lbbZmmTight = hasDeepJetShape_corr.refine("isSideband2lbbZmmTight", cut=[op.NOT(op.NOT(outZLeptons(lepton1,lepton2))),
        #                                           op.NOT(range_mass_for_bjets_tight(cleanedHighestBtagScoreJets))])
        #isSideband2lbbHmmTight = hasDeepJetShape_corr.refine("isSideband2lbbHmmTight",  cut=[op.NOT(op.NOT(outHLeptons(lepton1,lepton2))),
        #                                           op.NOT(range_mass_for_bjets_tight(cleanedHighestBtagScoreJets))])
        #isSignalRegionTight    = hasDeepJetShape_corr.refine("isSignalRegionTight", cut=[op.NOT(outHLeptons(lepton1,lepton2)),
        #                                           range_mass_for_bjets_tight(cleanedHighestBtagScoreJets)])
        ### Save mvaVariables to be retrieved later in the postprocessor and saved in a parquet file ###
        if self.args.mvaSkim and not self.args.BtagReweightingOff and not self.args.BtagReweightingOn:
            ### Here you can specify the variables that you want to save in the .parquet file, you need to add --mvaSkim to the command line ###
            #met = tree.MET
            if self.is_MC:
                #genPartFlav
                print("This is a MC sample:",self.is_MC)
                lep1_Hcand_genPartFlav = op.c_int(bestHmm[0].genPartFlav, cast=True)
                lep2_Hcand_genPartFlav = op.c_int(bestHmm[1].genPartFlav, cast=True)
            else: 
                print("This is NOT a MC sample")
                lep1_Hcand_genPartFlav = a
                lep2_Hcand_genPartFlav = a
            mvaVariables = {
                "weight"                       : noSel.weight,
                "nmuons_cleaned"               : op.c_int(op.rng_len(self.muons), cast=True),
                "nelectrons_cleaned"           : op.c_int(op.rng_len(self.electrons), cast=True),
                # njets      
                "njets_cleaned"                : op.c_int(op.rng_len(self.ak4JetsByPt), cast=True),
                # l1
                "lep1_Hcand_genPartFlav" : lep1_Hcand_genPartFlav,
                "lep1_Hcand_px"          : bestHmm[0].p4.Px(),  
                "lep1_Hcand_py"          : bestHmm[0].p4.Py(),
                "lep1_Hcand_pt"          : bestHmm[0].p4.Pt(),
                "lep1_Hcand_ptCorrected" : lepton1.Pt(),
                "lep1_Hcand_eta"         : bestHmm[0].p4.Eta(),
                "lep1_Hcand_phi"         : bestHmm[0].p4.Phi(), 
                "lep1_Hcand_E"           : bestHmm[0].p4.E(),
                "lep1_Hcand_m"           : bestHmm[0].p4.M(),
                "lep1_Hcand_ch"          : bestHmm[0].charge,
                # l2
                "lep2_Hcand_genPartFlav" : lep2_Hcand_genPartFlav,
                "lep2_Hcand_px"          : bestHmm[1].p4.Px(),  
                "lep2_Hcand_py"          : bestHmm[1].p4.Py(),
                "lep2_Hcand_pt"          : bestHmm[1].p4.Pt(),
                "lep2_Hcand_ptCorrected" : lepton2.Pt(),
                "lep2_Hcand_eta"         : bestHmm[1].p4.Eta(),
                "lep2_Hcand_phi"         : bestHmm[1].p4.Phi(), 
                "lep2_Hcand_E"           : bestHmm[1].p4.E(),
                "lep2_Hcand_m"           : bestHmm[1].p4.M(),
                "lep2_Hcand_ch"          : bestHmm[1].charge,
                # delta between leptons
                "deltaR_l1_l2"           : op.deltaR(bestHmm[0].p4,bestHmm[1].p4),
                "deltaPhi_l1_l2"         : op.deltaPhi(bestHmm[0].p4,bestHmm[1].p4),
                "deltaEta_l1_l2"         : bestHmm[0].p4.Eta()-bestHmm[1].p4.Eta(),
                # j1 highest score        
                "j1_HighestBTagScore_pt"                  : op.switch(op.rng_len(cleanedHighestBtagScoreJets) > 0, cleanedHighestBtagScoreJets[0].p4.Pt(), a),
                "j1_HighestBTagScore_E"                   : op.switch(op.rng_len(cleanedHighestBtagScoreJets) > 0, cleanedHighestBtagScoreJets[0].p4.E(),  a),
                "j1_HighestBTagScore_m"                   : op.switch(op.rng_len(cleanedHighestBtagScoreJets) > 0, cleanedHighestBtagScoreJets[0].p4.M(),   a),
                "j1_HighestBTagScore_phi"                 : op.switch(op.rng_len(cleanedHighestBtagScoreJets) > 0, cleanedHighestBtagScoreJets[0].p4.Phi(), a),
                "j1_HighestBTagScore_eta"                 : op.switch(op.rng_len(cleanedHighestBtagScoreJets) > 0, cleanedHighestBtagScoreJets[0].p4.Eta(), a),
                "j1_HighestBTagScore_btagDeepFlavB"       : op.switch(op.rng_len(cleanedHighestBtagScoreJets) > 0, cleanedHighestBtagScoreJets[0].btagDeepFlavB, a),
                # j2  highest score   
                "j2_HighestBTagScore_pt"                  : op.switch(op.rng_len(cleanedHighestBtagScoreJets) > 0, cleanedHighestBtagScoreJets[1].p4.Pt(), a),
                "j2_HighestBTagScore_E"                   : op.switch(op.rng_len(cleanedHighestBtagScoreJets) > 0, cleanedHighestBtagScoreJets[1].p4.E(),  a),
                "j2_HighestBTagScore_m"                   : op.switch(op.rng_len(cleanedHighestBtagScoreJets) > 0, cleanedHighestBtagScoreJets[1].p4.M(),   a),
                "j2_HighestBTagScore_phi"                 : op.switch(op.rng_len(cleanedHighestBtagScoreJets) > 0, cleanedHighestBtagScoreJets[1].p4.Phi(), a),
                "j2_HighestBTagScore_eta"                 : op.switch(op.rng_len(cleanedHighestBtagScoreJets) > 0, cleanedHighestBtagScoreJets[1].p4.Eta(), a),
                "j2_HighestBTagScore_btagDeepFlavB"       : op.switch(op.rng_len(cleanedHighestBtagScoreJets) > 0, cleanedHighestBtagScoreJets[1].btagDeepFlavB, a),
                # j1 highest score        
                "j1_HighestPt_pt"                  : op.switch(op.rng_len(cleanedHighestPtJets) > 0, cleanedHighestPtJets[0].p4.Pt(),       a),
                "j1_HighestPt_E"                   : op.switch(op.rng_len(cleanedHighestPtJets) > 0, cleanedHighestPtJets[0].p4.E(),        a),
                "j1_HighestPt_m"                   : op.switch(op.rng_len(cleanedHighestPtJets) > 0, cleanedHighestPtJets[0].p4.M(),        a),
                "j1_HighestPt_phi"                 : op.switch(op.rng_len(cleanedHighestPtJets) > 0, cleanedHighestPtJets[0].p4.Phi(),      a),
                "j1_HighestPt_eta"                 : op.switch(op.rng_len(cleanedHighestPtJets) > 0, cleanedHighestPtJets[0].p4.Eta(),      a),
                "j1_HighestPt_btagDeepFlavB"       : op.switch(op.rng_len(cleanedHighestPtJets) > 0, cleanedHighestPtJets[0].btagDeepFlavB, a),
                # j2  highest score   
                "j2_HighestPt_pt"                  : op.switch(op.rng_len(cleanedHighestPtJets) > 0, cleanedHighestPtJets[1].p4.Pt(),       a),
                "j2_HighestPt_E"                   : op.switch(op.rng_len(cleanedHighestPtJets) > 0, cleanedHighestPtJets[1].p4.E(),        a),
                "j2_HighestPt_m"                   : op.switch(op.rng_len(cleanedHighestPtJets) > 0, cleanedHighestPtJets[1].p4.M(),        a),
                "j2_HighestPt_phi"                 : op.switch(op.rng_len(cleanedHighestPtJets) > 0, cleanedHighestPtJets[1].p4.Phi(),      a),
                "j2_HighestPt_eta"                 : op.switch(op.rng_len(cleanedHighestPtJets) > 0, cleanedHighestPtJets[1].p4.Eta(),      a),
                "j2_HighestPt_btagDeepFlavB"       : op.switch(op.rng_len(cleanedHighestPtJets) > 0, cleanedHighestPtJets[1].btagDeepFlavB, a),
                # delta between jets (btagscore)
                "deltaR_j1_j2_BtagScorejets"   : op.switch(op.rng_len(cleanedHighestBtagScoreJets) > 1,    op.deltaR(cleanedHighestBtagScoreJets[0].p4,cleanedHighestBtagScoreJets[1].p4),   a),
                "deltaPhi_j1_j2_BtagScorejets" : op.switch(op.rng_len(cleanedHighestBtagScoreJets) > 1,    op.deltaPhi(cleanedHighestBtagScoreJets[0].p4,cleanedHighestBtagScoreJets[1].p4), a),
                "deltaEta_j1_j2_BtagScorejets" : op.switch(op.rng_len(cleanedHighestBtagScoreJets) > 1,    cleanedHighestBtagScoreJets[0].p4.Eta() - cleanedHighestBtagScoreJets[1].p4.Eta(),a),
                # delta between jets (highest pt)
                "deltaR_j1_j2_HighestPtJets"   : op.switch(op.rng_len(cleanedHighestPtJets) > 1,    op.deltaR(cleanedHighestPtJets[0].p4,cleanedHighestPtJets[1].p4),                 a),
                "deltaPhi_j1_j2_HighestPtJets" : op.switch(op.rng_len(cleanedHighestPtJets) > 1,    op.deltaPhi(cleanedHighestPtJets[0].p4,cleanedHighestPtJets[1].p4),               a),
                "deltaEta_j1_j2_BtagScorejets" : op.switch(op.rng_len(cleanedHighestPtJets) > 1,    cleanedHighestPtJets[0].p4.Eta() - cleanedHighestPtJets[1].p4.Eta(),              a),
                # delta between lepton and (btagscore) jets
                "deltaR_l1_j1_BtagScorejets"   :  op.switch(op.rng_len(cleanedHighestBtagScoreJets) > 0,   op.deltaR(bestHmm[0].p4, cleanedHighestBtagScoreJets[0].p4),     a),
                "deltaR_l2_j2_BtagScorejets"   :  op.switch(op.rng_len(cleanedHighestBtagScoreJets) > 1,   op.deltaR(bestHmm[1].p4, cleanedHighestBtagScoreJets[1].p4),     a),
                "deltaPhi_l1_j1_BtagScorejets" :  op.switch(op.rng_len(cleanedHighestBtagScoreJets) > 0,   op.deltaPhi(bestHmm[0].p4, cleanedHighestBtagScoreJets[0].p4),   a),
                "deltaPhi_l2_j2_BtagScorejets" :  op.switch(op.rng_len(cleanedHighestBtagScoreJets) > 1,   op.deltaPhi(bestHmm[1].p4, cleanedHighestBtagScoreJets[1].p4),   a),
                "deltaEta_l1_j1_BtagScorejets" :  op.switch(op.rng_len(cleanedHighestBtagScoreJets) > 0,   bestHmm[0].p4.Eta() - cleanedHighestBtagScoreJets[0].p4.Eta(),   a),
                "deltaEta_l2_j2_BtagScorejets" :  op.switch(op.rng_len(cleanedHighestBtagScoreJets) > 1,   bestHmm[1].p4.Eta() - cleanedHighestBtagScoreJets[1].p4.Eta(),   a),
                # delta between lepton and (highest pt) jets
                "deltaR_l1_j1_HighestPtJets"   :  op.switch(op.rng_len(cleanedHighestPtJets) > 0,   op.deltaR(bestHmm[0].p4, cleanedHighestPtJets[0].p4),            a),
                "deltaR_l2_j2_HighestPtJets"   :  op.switch(op.rng_len(cleanedHighestPtJets) > 1,   op.deltaR(bestHmm[1].p4, cleanedHighestPtJets[1].p4),            a),
                "deltaPhi_l1_j1_HighestPtJets" :  op.switch(op.rng_len(cleanedHighestPtJets) > 0,   op.deltaPhi(bestHmm[0].p4, cleanedHighestPtJets[0].p4),          a),
                "deltaPhi_l2_j2_HighestPtJets" :  op.switch(op.rng_len(cleanedHighestPtJets) > 1,   op.deltaPhi(bestHmm[1].p4, cleanedHighestPtJets[1].p4),          a),
                "deltaEta_l1_j1_HighestPtJets" :  op.switch(op.rng_len(cleanedHighestPtJets) > 0,   bestHmm[0].p4.Eta() - cleanedHighestPtJets[0].p4.Eta(),          a),
                "deltaEta_l2_j2_HighestPtJets" :  op.switch(op.rng_len(cleanedHighestPtJets) > 1,   bestHmm[1].p4.Eta() - cleanedHighestPtJets[1].p4.Eta(),          a),
                # mass
                "Invariantmass_leptons_BtagScorejets" : op.switch(op.rng_len(cleanedHighestPtJets) > 1,op.invariant_mass(bestHmm[0].p4, bestHmm[1].p4, cleanedHighestPtJets[0].p4,cleanedHighestPtJets[1].p4),a),
                "Invariantmass_leptons_HighestPtJets" : op.switch(op.rng_len(cleanedHighestBtagScoreJets) > 1,op.invariant_mass(bestHmm[0].p4, bestHmm[1].p4, cleanedHighestBtagScoreJets[0].p4,cleanedHighestBtagScoreJets[1].p4),a),
                "HiggsToMMCandidateMass"              : op.invariant_mass(bestHmm[0].p4, bestHmm[1].p4),
                "HiggsToMMCandidateMassCorrected"     : (lepton1_FSR+lepton2_FSR).M(),
                #"met_pt"                       :   met.pt,
                #"met_phi"                      :   met.phi,
                #"met_deltaPhi_mu1"            :   phi_mpi_pi(met.phi - bestHmm[0].p4.Phi()),
                #"met_deltaPhi_mu2"            :   phi_mpi_pi(met.phi - bestHmm[1].p4.Phi()),

                }
            from bamboo.plots import Skim
            plots.append(Skim("ReducedMVATree", mvaVariables, isSignalRegion))
            
        if self.args.PrintYield:
            plots.append(self.yields)
            self.yields.add(hasZmm,title="ZmmCR")
            self.yields.add(hasHmm,title="HmmCR")
            self.yields.add(hasZjets,title="ZjetsCR")
            self.yields.add(hasHjets,title="HjetsCR")
            self.yields.add(isSideband2lbbHmm,title="Sideband2lbbHmm")
            self.yields.add(isSideband2lbbZmm,title="Sideband2lbbZmm")
            self.yields.add(isSignalRegion,title="SignalRegion")
            self.yields.add(isSignalRegion,title="hasZCRDY")
        
        # Plots about leptons in Hmm/Zmm and sideband control regions 
        if not self.args.BtagReweightingOff and not self.args.BtagReweightingOn:
            plots += self.plotFourObjects(self.ak4JetsByPt,bestHmm,hasZCRDY,"hasZCRDYFourObjectsPtOrdering")
            plots += self.plotFourObjectsCorrected(self.ak4JetsByPt,lepton1_FSR,lepton2_FSR,hasZCRDY,"hasZCRDYPtOrderingFSRGeoFit")
            
            plots += self.plotPairs(bestZmm, hasZmm, "CRZmmtightMuons",self.muons)
            plots += self.plotPairsCorrected(lepton1_FSR, lepton2_FSR , hasZmm, "CRZmmtightMuonsFSRGeoFit")
            plots += self.plotPairs(bestZmm, isSideband2lbbZmm,"CRSideband2lbbZmmtightMuons" ,self.muons)
            plots += self.plotPairsCorrected(lepton1_FSR, lepton2_FSR , isSideband2lbbZmm, "CRSideband2lbbZmmtightMuonsFSRGeoFit")
            
            plots += self.plotPairs(bestHmm, hasHmm, "CRHmmtightMuons",self.muons)
            plots += self.plotPairsCorrected(lepton1_FSR, lepton2_FSR , hasHmm, "CRHmmtightMuonsFSRGeoFit")
            plots += self.plotPairs(bestHmm, isSideband2lbbHmm,"CRSideband2lbbHmmtightMuons" ,self.muons)
            plots += self.plotPairsCorrected(lepton1_FSR, lepton2_FSR , isSideband2lbbHmm, "CRSideband2lbbHmmtightMuonsFSRGeoFit")
            plots += self.plotPairs(bestHmm, isSignalRegion,"SignalRegiontightMuons" ,self.muons)
            plots += self.plotPairsCorrected(lepton1_FSR, lepton2_FSR , isSignalRegion, "SignalRegiontightMuonsFSRGeoFit")
            
            
            # Plots about jets in the Zjets and Hjets control regions + sidebands
            plots += self.plotJets(self.ak4JetsByPt, hasZjets, "CRZjetsPtOrdering")
            plots += self.plotJets(self.ak4JetsByBtagScore, hasZjets, "CRZjetsBScoreOrdering")
            plots += self.plotJets(self.ak4JetsByPt, hasHjets, "CRHjetsPtOrdering")
            plots += self.plotJets(self.ak4JetsByBtagScore, hasHjets, "CRHjetsBScoreOrdering")
            plots += self.plotJets(self.ak4JetsByBtagScore, isSignalRegion, "SignalRegionBScoreOrdering")
            plots += self.plotJets(self.ak4JetsByPt, isSignalRegion, "SignalRegionPtOrdering")
            
            plots += self.plotJets(self.ak4JetsByPt, isSideband2lbbZmm, "CRSideband2lbbZmmPtOrdering")
            plots += self.plotJets(self.ak4JetsByBtagScore, isSideband2lbbZmm, "CRSideband2lbbZmmBScoreOrdering")
            plots += self.plotJets(self.ak4JetsByBtagScore, isSideband2lbbHmm, "CRSideband2lbbHmmPtOrdering")
            plots += self.plotJets(self.ak4JetsByBtagScore, isSideband2lbbHmm, "CRSideband2lbbHmmBScoreOrdering")
            
            # Plots about leptons and jets in the sidebands (jets,pair,sel,category)
            plots += self.plotFourObjects(self.ak4JetsByPt,bestZmm,isSideband2lbbZmm,"CRSideband2lbbZmmbestZmmPtOrdering")
            plots += self.plotFourObjects(self.ak4JetsByPt,bestHmm,isSignalRegion,"SignalRegionPtOrdering")
            plots += self.plotFourObjects(self.ak4JetsByPt,bestHmm,isSideband2lbbHmm,"CRSideband2lbbHmmbestHmmPtOrdering")
            plots += self.plotFourObjects(self.ak4JetsByBtagScore,bestHmm,isSideband2lbbHmm,"CRSideband2lbbHmmbestHmmBScoreOrdering")
            plots += self.plotFourObjects(self.ak4JetsByBtagScore,bestHmm,isSignalRegion,"SignalRegionBScoreOrdering")
            plots += self.plotFourObjects(self.ak4JetsByBtagScore,bestZmm,isSideband2lbbZmm,"CRSideband2lbbZmmbestZmmBScoreOrdering")
            
            plots += self.plotFourObjectsCorrected(self.ak4JetsByPt,lepton1_FSR,lepton2_FSR,isSignalRegion,"SignalRegionPtOrderingFSRGeoFit")
            plots += self.plotFourObjectsCorrected(self.ak4JetsByPt,lepton1_FSR,lepton2_FSR,isSideband2lbbHmm,"CRSideband2lbbHmmbestHmmPtOrderingFSRGeoFit")
            plots += self.plotFourObjectsCorrected(self.ak4JetsByBtagScore,lepton1_FSR,lepton2_FSR,isSideband2lbbHmm,"CRSideband2lbbHmmbestHmmBScoreOrderingFSRGeoFit")
            plots += self.plotFourObjectsCorrected(self.ak4JetsByBtagScore,lepton1_FSR,lepton2_FSR,isSignalRegion,"SignalRegionBScoreOrderingFSRGeoFit")
            
        categories = []
        categories.append(category(nMuons=2, nElectrons=0))
        
        return plots
    
    def plotVertices(self, nvertices, sel, category):
        plots = []
        plots.append(Plot.make1D(f"h_{category}_nGoodVtx", op.c_int(nvertices, cast=True),sel, EquidistantBinning(100, 0., 100.), xTitle="Number of primary vertices"))
        return plots
    def plotPairsCorrected(self,lep1, lep2, sel, category):
        plots = []
        if "SignalRegion" in category:
            plots.append(Plot.make1D(f"h_{category}_mll", op.invariant_mass(lep1,lep2), sel, EquidistantBinning(80, 40., 200.), title="Invariant Mass of two muons", xTitle= "m_{\mu\mu} (GeV)",plotopts={'no-data': True, 'show-ratio': False}))
            #if op.c_float(op.invariant_mass(lep1,lep2)) >=70. and op.c_float(op.invariant_mass(lep1,lep2)) <= 115.:
            #    plots.append(Plot.make1D(f"h_{category}_mll_narrow_Z", op.invariant_mass(lep1,lep2), sel, EquidistantBinning(50, 70., 115.), title="Invariant Mass of two muons Narrow", xTitle= "m_{\mu\mu} (GeV)",plotopts={'no-data': True, 'show-ratio': False}))
            plots.append(Plot.make1D(f"h_{category}_lep1Pt", lep1.Pt(), sel, EquidistantBinning(50, 26., 200.), title=" lepton1 pT", xTitle= "pT(l1) (GeV)",plotopts={'no-data': True, 'show-ratio': False}))
            plots.append(Plot.make1D(f"h_{category}_lep2Pt", lep2.Pt(), sel, EquidistantBinning(50, 20., 200.), title=" lepton2 pT", xTitle= "pT(l2) (GeV)",plotopts={'no-data': True, 'show-ratio': False}))
            plots.append(Plot.make1D(f"h_{category}_lep1Eta", lep1.Eta(),  sel, EquidistantBinning(10, -2.5, 2.5), title=" lepton1 Eta", xTitle= "\eta(l1)",plotopts={'no-data': True, 'show-ratio': False}))
            plots.append(Plot.make1D(f"h_{category}_lep2Eta", lep2.Eta(),  sel, EquidistantBinning(10, -2.5, 2.5), title=" lepton2 Eta", xTitle= "\eta(l2)",plotopts={'no-data': True, 'show-ratio': False}))
            plots.append(Plot.make1D(f"h_{category}_l1l2Eta", (lep1+lep2).Eta(),  sel, EquidistantBinning(10, -2.5, 2.5), title=" l1l2 Eta", xTitle= "\eta (l1l2)",plotopts={'no-data': True, 'show-ratio': False}))
            plots.append(Plot.make1D(f"h_{category}_l1l2Phi", (lep1+lep2).Phi(),  sel, EquidistantBinning(10, -2.5, 2.5), title=" l1l2 Phi", xTitle= "\phi (l1l2)",plotopts={'no-data': True, 'show-ratio': False}))
            plots.append(Plot.make1D(f"h_{category}_lep1lep2Pt", (lep1+lep2).Pt(), sel, EquidistantBinning(50, 0., 700.), title=" l1l2 pT", xTitle= "pT(l1l2) (GeV)",plotopts={'no-data': True, 'show-ratio': False}))
            plots.append(Plot.make1D(f"h_{category}_deltaRL1L2", op.deltaR(lep1,lep2),  sel, EquidistantBinning(6, 0.5, 3.5), title=" Delta R", xTitle= "\Delta R (l1,l2)" ,plotopts={'no-data': True, 'show-ratio': False}))
            plots.append(Plot.make1D(f"h_{category}_deltaPhiL1L2", op.abs(op.deltaPhi(lep1,lep2)),  sel, EquidistantBinning(5, 0.5, 2.5), title=" Delta Phi", xTitle= "\Delta \phi (l1,l2)",plotopts={'no-data': True, 'show-ratio': False} ))
            plots.append(Plot.make1D(f"h_{category}_deltaEtaL1L2", op.abs(lep1.Eta() - lep2.Eta()),  sel, EquidistantBinning(6, 0., 1.5), title=" Delta Eta", xTitle= "| \Delta \eta | (l1,l2)" ,plotopts={'no-data': True, 'show-ratio': False}))
        else:
            plots.append(Plot.make1D(f"h_{category}_mll", op.invariant_mass(lep1,lep2), sel, EquidistantBinning(80, 40., 200.), title="Invariant Mass of two muons", xTitle= "m_{\mu\mu} (GeV)"))
            if op.c_float(op.invariant_mass(lep1,lep2)) >=70. and op.c_float(op.invariant_mass(lep1,lep2)) <= 115. and "Zmm" in category and not "Zjets" in category and not "Sideband" in category:
                plots.append(Plot.make1D(f"h_{category}_mll_narrow_Z", op.invariant_mass(lep1,lep2), sel, EquidistantBinning(50, 70., 115.), title="Invariant Mass of two muons Narrow", xTitle= "m_{\mu\mu} (GeV)"))
            if op.c_float(op.invariant_mass(lep1,lep2)) >=115. and op.c_float(op.invariant_mass(lep1,lep2)) <= 135. and "Hmm" in category and not "Hjets" in category and not "Sideband" in category:
                plots.append(Plot.make1D(f"h_{category}_mll_narrow_H", op.invariant_mass(lep1,lep2), sel, EquidistantBinning(35, 115., 135.), title="Invariant Mass of two muons Narrow", xTitle= "m_{\mu\mu} (GeV)"))
            plots.append(Plot.make1D(f"h_{category}_lep1Ptnarrow", lep1.Pt(), sel, EquidistantBinning(50, 26., 200.), title=" lepton1 pT", xTitle= "pT(l1) (GeV)"))
            plots.append(Plot.make1D(f"h_{category}_lep1Pt", lep1.Pt(), sel, EquidistantBinning(50, 26., 700.), title=" lepton1 pT", xTitle= "pT(l1) (GeV)"))
            plots.append(Plot.make1D(f"h_{category}_lep2Pt", lep2.Pt(), sel, EquidistantBinning(50, 20., 380.), title=" lepton2 pT", xTitle= "pT(l2) (GeV)"))
            plots.append(Plot.make1D(f"h_{category}_lep2Ptnarrow", lep2.Pt(), sel, EquidistantBinning(50, 20., 200.), title=" lepton2 pT", xTitle= "pT(l2) (GeV)"))
            plots.append(Plot.make1D(f"h_{category}_lep1Eta", lep1.Eta(),  sel, EquidistantBinning(10, -2.5, 2.5), title=" lepton1 Eta", xTitle= "\eta (l1)"))
            plots.append(Plot.make1D(f"h_{category}_lep2Eta", lep2.Eta(),  sel, EquidistantBinning(10, -2.5, 2.5), title=" lepton2 Eta", xTitle= "\eta (l2)"))
            plots.append(Plot.make1D(f"h_{category}_l1l2Eta", (lep1+lep2).Eta(),  sel, EquidistantBinning(10, -2.5, 2.5), title=" l1l2 Eta", xTitle= "\eta (l1l2)"))
            plots.append(Plot.make1D(f"h_{category}_l1l2Phi", (lep1+lep2).Phi(),  sel, EquidistantBinning(10, -2.5, 2.5), title=" l1l2 Phi", xTitle= "\phi (l1l2)"))
            plots.append(Plot.make1D(f"h_{category}_lep1lep2Pt", (lep1+lep2).Pt(), sel, EquidistantBinning(50, 0., 700.), title=" l1l2 pT", xTitle= "pT(l1l2) (GeV)"))
            plots.append(Plot.make1D(f"h_{category}_deltaRL1L2", op.deltaR(lep1,lep2),  sel, EquidistantBinning(6, 0.5, 3.5), title=" Delta R", xTitle= "\Delta R (l1,l2)") )
            plots.append(Plot.make1D(f"h_{category}_deltaPhiL1L2", op.abs(op.deltaPhi(lep1,lep2)),  sel, EquidistantBinning(5, 0.5, 2.5), title=" Delta Phi", xTitle= "| \Delta \phi | (l1,l2)") )
            plots.append(Plot.make1D(f"h_{category}_deltaEtaL1L2", op.abs(lep1.Eta() - lep2.Eta()),  sel, EquidistantBinning(6, 0., 1.5), title=" Delta Eta", xTitle= "| \Delta \eta |(l1,l2)") )
        return plots
    
    def plotPairs(self, pair, sel, category, muons):
            plots = []
            lep1 = op.construct("ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float> >",([op.c_float(pair[0].p4.Pt()),op.c_float(pair[0].p4.Eta()),op.c_float(pair[0].p4.Phi()),op.c_float(pair[0].p4.M())]))
            lep2 = op.construct("ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float> >",([op.c_float(pair[1].p4.Pt()),op.c_float(pair[1].p4.Eta()),op.c_float(pair[1].p4.Phi()),op.c_float(pair[1].p4.M())]))
            if "SignalRegion" in category:
                plots.append(Plot.make1D(f"h_{category}_nMuons", op.c_int(op.rng_len(muons), cast=True),sel, EquidistantBinning(15, 0., 15.), xTitle="Number of cleaned muons",plotopts={'no-data': True, 'show-ratio': False}))
                plots.append(Plot.make1D(f"h_{category}_mll", op.invariant_mass(lep1,lep2), sel, EquidistantBinning(80, 40., 200.), title="Invariant Mass of two muons", xTitle= "m_{\mu\mu} (GeV)",plotopts={'no-data': True, 'show-ratio': False}))
                #if op.c_float(op.invariant_mass(lep1,lep2)) >=70. and op.c_float(op.invariant_mass(lep1,lep2)) <= 115.:
                #    plots.append(Plot.make1D(f"h_{category}_mll_narrow", op.invariant_mass(lep1,lep2), sel, EquidistantBinning(50, 70., 115.), title="Invariant Mass of two muons Narrow", xTitle= "m_{\mu\mu} (GeV)",plotopts={'no-data': True, 'show-ratio': False}))
                plots.append(Plot.make1D(f"h_{category}_lep1Pt", lep1.Pt(), sel, EquidistantBinning(50, 26., 200.), title=" lepton1 pT", xTitle= "pT(l1) (GeV)",plotopts={'no-data': True, 'show-ratio': False}))
                plots.append(Plot.make1D(f"h_{category}_lep2Pt", lep2.Pt(), sel, EquidistantBinning(50, 20., 200.), title=" lepton2 pT", xTitle= "pT(l2) (GeV)",plotopts={'no-data': True, 'show-ratio': False}))
                plots.append(Plot.make1D(f"h_{category}_lep1Eta", lep1.Eta(),  sel, EquidistantBinning(10, -2.5, 2.5), title=" lepton1 Eta", xTitle= "\eta(l1)",plotopts={'no-data': True, 'show-ratio': False}))
                plots.append(Plot.make1D(f"h_{category}_lep2Eta", lep2.Eta(),  sel, EquidistantBinning(10, -2.5, 2.5), title=" lepton2 Eta", xTitle= "\eta(l2)",plotopts={'no-data': True, 'show-ratio': False}))
                plots.append(Plot.make1D(f"h_{category}_l1l2Eta", (lep1+lep2).Eta(),  sel, EquidistantBinning(10, -2.5, 2.5), title=" l1l2 Eta", xTitle= "\eta (l1l2)",plotopts={'no-data': True, 'show-ratio': False}))
                plots.append(Plot.make1D(f"h_{category}_l1l2Phi", (lep1+lep2).Phi(),  sel, EquidistantBinning(10, -2.5, 2.5), title=" l1l2 Phi", xTitle= "\phi (l1l2)",plotopts={'no-data': True, 'show-ratio': False}))
                plots.append(Plot.make1D(f"h_{category}_lep1lep2Pt", (lep1+lep2).Pt(), sel, EquidistantBinning(50, 0., 700.), title=" l1l2 pT", xTitle= "pT(l1l2) (GeV)",plotopts={'no-data': True, 'show-ratio': False}))
                plots.append(Plot.make1D(f"h_{category}_deltaRL1L2", op.deltaR(lep1,lep2),  sel, EquidistantBinning(6, 0.5, 3.5), title=" Delta R", xTitle= "\Delta R (l1,l2)" ,plotopts={'no-data': True, 'show-ratio': False}))
                plots.append(Plot.make1D(f"h_{category}_deltaPhiL1L2", op.abs(op.deltaPhi(lep1,lep2)),  sel, EquidistantBinning(5, 0.5, 2.5), title=" Delta Phi", xTitle= "\Delta \phi (l1,l2)",plotopts={'no-data': True, 'show-ratio': False} ))
                plots.append(Plot.make1D(f"h_{category}_deltaEtaL1L2", op.abs(lep1.Eta() - lep2.Eta()),  sel, EquidistantBinning(6, 0., 1.5), title=" Delta Eta", xTitle= "| \Delta \eta | (l1,l2)" ,plotopts={'no-data': True, 'show-ratio': False}))
            else:
                plots.append(Plot.make1D(f"h_{category}_nMuons", op.c_int(op.rng_len(muons), cast=True),sel, EquidistantBinning(15, 0., 15.), xTitle="Number of cleaned muons"))
                plots.append(Plot.make1D(f"h_{category}_mll", op.invariant_mass(lep1,lep2), sel, EquidistantBinning(80, 40., 200.), title="Invariant Mass of two muons", xTitle= "m_{\mu\mu} (GeV)"))
                if op.c_float(op.invariant_mass(lep1,lep2)) >=70. and op.c_float(op.invariant_mass(lep1,lep2)) <= 115. and "Zmm" in category and not "Zjets" in category and not "Sideband" in category:
                    plots.append(Plot.make1D(f"h_{category}_mll_narrow_Z", op.invariant_mass(lep1,lep2), sel, EquidistantBinning(50, 70., 115.), title="Invariant Mass of two muons Narrow", xTitle= "m_{\mu\mu} (GeV)"))
                if op.c_float(op.invariant_mass(lep1,lep2)) >=115. and op.c_float(op.invariant_mass(lep1,lep2)) <= 135. and "Hmm" in category and not "Hjets" in category and not "Sideband" in category:
                    plots.append(Plot.make1D(f"h_{category}_mll_narrow_H", op.invariant_mass(lep1,lep2), sel, EquidistantBinning(35, 115., 135.), title="Invariant Mass of two muons Narrow", xTitle= "m_{\mu\mu} (GeV)"))
                plots.append(Plot.make1D(f"h_{category}_lep1Ptnarrow", lep1.Pt(), sel, EquidistantBinning(50, 26., 200.), title=" lepton1 pT", xTitle= "pT(l1) (GeV)"))
                plots.append(Plot.make1D(f"h_{category}_lep1Pt", lep1.Pt(), sel, EquidistantBinning(50, 26., 700.), title=" lepton1 pT", xTitle= "pT(l1) (GeV)"))
                plots.append(Plot.make1D(f"h_{category}_lep2Pt", lep2.Pt(), sel, EquidistantBinning(50, 20., 380.), title=" lepton2 pT", xTitle= "pT(l2) (GeV)"))
                plots.append(Plot.make1D(f"h_{category}_lep2Ptnarrow", lep2.Pt(), sel, EquidistantBinning(50, 20., 200.), title=" lepton2 pT", xTitle= "pT(l2) (GeV)"))
                plots.append(Plot.make1D(f"h_{category}_lep1Eta", lep1.Eta(),  sel, EquidistantBinning(10, -2.5, 2.5), title=" lepton1 Eta", xTitle= "\eta (l1)"))
                plots.append(Plot.make1D(f"h_{category}_lep2Eta", lep2.Eta(),  sel, EquidistantBinning(10, -2.5, 2.5), title=" lepton2 Eta", xTitle= "\eta (l2)"))
                plots.append(Plot.make1D(f"h_{category}_l1l2Eta", (lep1+lep2).Eta(),  sel, EquidistantBinning(10, -2.5, 2.5), title=" l1l2 Eta", xTitle= "\eta (l1l2)"))
                plots.append(Plot.make1D(f"h_{category}_l1l2Phi", (lep1+lep2).Phi(),  sel, EquidistantBinning(10, -2.5, 2.5), title=" l1l2 Phi", xTitle= "\phi (l1l2)"))
                plots.append(Plot.make1D(f"h_{category}_lep1lep2Pt", (lep1+lep2).Pt(), sel, EquidistantBinning(50, 0., 700.), title=" l1l2 pT", xTitle= "pT(l1l2) (GeV)"))
                plots.append(Plot.make1D(f"h_{category}_deltaRL1L2", op.deltaR(lep1,lep2),  sel, EquidistantBinning(6, 0.5, 3.5), title=" Delta R", xTitle= "\Delta R (l1,l2)") )
                plots.append(Plot.make1D(f"h_{category}_deltaPhiL1L2", op.abs(op.deltaPhi(lep1,lep2)),  sel, EquidistantBinning(5, 0.5, 2.5), title=" Delta Phi", xTitle= "| \Delta \phi | (l1,l2)") )
                plots.append(Plot.make1D(f"h_{category}_deltaEtaL1L2", op.abs(lep1.Eta() - lep2.Eta()),  sel, EquidistantBinning(6, 0., 1.5), title=" Delta Eta", xTitle= "| \Delta \eta |(l1,l2)") )
            return plots

    def plotJets(self, jets, sel, category):
            plots = []
            j1 = op.construct("ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float> >",([op.c_float(jets[0].pt),op.c_float(jets[0].eta),op.c_float(jets[0].phi),op.c_float(jets[0].p4.M())]))
            j2 = op.construct("ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float> >",([op.c_float(jets[1].pt),op.c_float(jets[1].eta),op.c_float(jets[1].phi),op.c_float(jets[1].p4.M())]))
            if "SignalRegion" in category:
                plots.append(Plot.make1D(f"h_{category}_nJets", op.c_int(op.rng_len(jets), cast=True),sel, EquidistantBinning(15, 0., 15.), xTitle="Number of jets",plotopts={'no-data': True, 'show-ratio': False}))
                plots.append(Plot.make1D(f"h_{category}_leadingPt", jets[0].pt,sel, EquidistantBinning(50, 30., 700.), title="pT(j1)", xTitle="pT(j1) (GeV/c)",plotopts={'no-data': True, 'show-ratio': False}))
                plots.append(Plot.make1D(f"h_{category}_leadingEta", jets[0].eta,sel, EquidistantBinning(10, -2.5, 2.5), title="eta(j1)", xTitle="\eta (j1)",plotopts={'no-data': True, 'show-ratio': False}))
                plots.append(Plot.make1D(f"h_{category}_leadingPhi", jets[0].phi,sel, EquidistantBinning(10, -2.5, 2.5), title="phi(j1)", xTitle="\phi (j1)",plotopts={'no-data': True, 'show-ratio': False}))
                plots.append(Plot.make1D(f"h_{category}_subleadingPt", jets[1].pt,sel, EquidistantBinning(50, 30., 500.), title="pT(j2)", xTitle="pT(j2) (GeV/c)",plotopts={'no-data': True, 'show-ratio': False}))
                plots.append(Plot.make1D(f"h_{category}_subleadingEta", jets[1].eta,sel, EquidistantBinning(10, -2.5, 2.5), title="eta(j2)", xTitle="\eta(j2)",plotopts={'no-data': True, 'show-ratio': False}))
                plots.append(Plot.make1D(f"h_{category}_subleadingPhi", jets[1].phi,sel, EquidistantBinning(10,-2.5, 2.5), title="phi(j2)", xTitle="\phi(j2)",plotopts={'no-data': True, 'show-ratio': False}))
                plots.append(Plot.make1D(f"h_{category}_j1j2Eta", (j1+j2).Eta() ,  sel, EquidistantBinning(10, -2.5, 2.5), title=" j1j2 Eta", xTitle= "\eta (j1j2)",plotopts={'no-data': True, 'show-ratio': False}))
                plots.append(Plot.make1D(f"h_{category}_j1j2Phi", (j1+j2).Phi() ,  sel, EquidistantBinning(10, -2.5, 2.5), title=" j1j2 Phi", xTitle= "\phi (j1j2)",plotopts={'no-data': True, 'show-ratio': False}))
                plots.append(Plot.make1D(f"h_{category}_j1j2Pt", (j1+j2).Pt() , sel, EquidistantBinning(50, 0., 700.), title=" j1j2 pT", xTitle= "pT(j1j2) (GeV)",plotopts={'no-data': True, 'show-ratio': False}))
                plots.append(Plot.make1D(f"h_{category}_deltaRJ1J2", op.deltaR(jets[0].p4,jets[1].p4) ,sel,EquidistantBinning(10, 0, 5), title="DR(j1,j2)", xTitle="\Delta R(j1,j2)",plotopts={'no-data': True, 'show-ratio': False}))
                plots.append(Plot.make1D(f"h_{category}_deltaPhiJ1J2", op.deltaPhi(jets[0].p4,jets[1].p4) ,sel,EquidistantBinning(6, -3, 3), title=" Delta Phi", xTitle= "\Delta \phi (j1,j2)" ,plotopts={'no-data': True, 'show-ratio': False}))
                plots.append(Plot.make1D(f"h_{category}_deltaEtaJ1J2", jets[0].p4.Eta() - jets[1].p4.Eta() ,  sel, EquidistantBinning(6, -3, 3), title=" Delta Eta", xTitle= "\Delta \eta (j1,j2)" ,plotopts={'no-data': True, 'show-ratio': False}))
                plots.append(Plot.make1D(f"h_{category}_mJ1J2", op.invariant_mass(jets[0].p4, jets[1].p4) ,  sel, EquidistantBinning(50, 30., 500.), title=" mJ1J2", xTitle= "m(j1,j2) (GeV)" ,plotopts={'no-data': True, 'show-ratio': False}))
                plots.append(Plot.make1D(f"h_{category}_mJ1J2_narrow", op.invariant_mass(jets[0].p4, jets[1].p4) ,  sel, EquidistantBinning(50, 30., 200.), title=" mJ1J2 narrow", xTitle= "m(j1,j2) (GeV)" ,plotopts={'no-data': True, 'show-ratio': False}))
                if "BScore" in category:
                    if jets[0].btagDeepFlavB >=0:
                        plots.append(Plot.make1D(f"h_{category}_bscoreJ1", jets[0].btagDeepFlavB,sel,EquidistantBinning(20, 0, 1.0), title=" b tagging DeepJet score jet 1", xTitle= "btag score j1",plotopts={'no-data': True, 'show-ratio': False}) )
                    if jets[1].btagDeepFlavB >=0:
                        plots.append(Plot.make1D(f"h_{category}_bscoreJ2", jets[1].btagDeepFlavB,sel,EquidistantBinning(20, 0, 1.0), title=" b tagging DeepJet score jet 2", xTitle= "btag score j2",plotopts={'no-data': True, 'show-ratio': False}) )
            else:
                plots.append(Plot.make1D(f"h_{category}_nJets", op.c_int(op.rng_len(jets), cast=True),sel, EquidistantBinning(15, 0., 15.), xTitle="Number of jets"))
                plots.append(Plot.make1D(f"h_{category}_leadingPt", jets[0].pt,sel, EquidistantBinning(50, 30., 700.), title="pT(j1)", xTitle="pT(j1) (GeV/c)"))
                plots.append(Plot.make1D(f"h_{category}_leadingEta", jets[0].eta,sel, EquidistantBinning(10, -2.5, 2.5), title="eta(j1)", xTitle=" \eta (j1)"))
                plots.append(Plot.make1D(f"h_{category}_leadingPhi", jets[0].phi,sel, EquidistantBinning(10, -2.5, 2.5), title="phi(j1)", xTitle=" \phi (j1)"))
                plots.append(Plot.make1D(f"h_{category}_leadingnMuons", jets[0].nMuons,sel, EquidistantBinning(15, 0., 15.), title="nMuons(j1)", xTitle=" nMuons (j1)"))
                if jets[0].pt < 50.:
                    plots.append(Plot.make1D(f"h_{category}_leadingpuidDisc", jets[0].puIdDisc,sel, EquidistantBinning(20, 0., 1.0), title="puidDisc(j1)", xTitle=" puidDisc (j1)"))
                plots.append(Plot.make1D(f"h_{category}_subleadingPt", jets[1].pt,sel, EquidistantBinning(50, 30., 500), title="pT(j2)", xTitle="pT(j2) (GeV/c)"))
                plots.append(Plot.make1D(f"h_{category}_subleadingEta", jets[1].eta,sel, EquidistantBinning(10, -2.5, 2.5), title="eta(j2)", xTitle="\eta (j2)"))
                plots.append(Plot.make1D(f"h_{category}_subleadingPhi", jets[1].phi,sel, EquidistantBinning(10,-2.5, 2.5), title="phi(j2)", xTitle="\phi (j2)"))
                plots.append(Plot.make1D(f"h_{category}_subleadingnMuons", jets[1].nMuons,sel, EquidistantBinning(15, 0., 15.), title="nMuons(j2)", xTitle=" nMuons (j2)"))
                if jets[1].pt < 50.:
                    plots.append(Plot.make1D(f"h_{category}_subleadingpuidDisc", jets[1].puIdDisc,sel, EquidistantBinning(20, 0., 1.0), title="puidDisc(j2)", xTitle=" puidDisc (j2)"))
                if op.c_int(op.rng_len(jets), cast=True) >= 0:    
                    plots.append(Plot.make1D(f"h_{category}_j1j2Eta", (j1+j2).Eta() ,  sel, EquidistantBinning(10, -2.5, 2.5), title=" j1j2 Eta", xTitle= "\eta (j1j2)"))
                    plots.append(Plot.make1D(f"h_{category}_j1j2Pt", (j1+j2).Pt() , sel, EquidistantBinning(50, 0., 700.), title=" j1j2 pT", xTitle= "pT(j1j2) (GeV)"))
                    plots.append(Plot.make1D(f"h_{category}_j1j2Phi", (j1+j2).Phi() ,  sel, EquidistantBinning(10, -2.5, 2.5), title=" j1j2 Phi", xTitle= "\phi (j1j2)"))
                    plots.append(Plot.make1D(f"h_{category}_deltaRJ1J2", op.deltaR(jets[0].p4,jets[1].p4) ,sel,EquidistantBinning(10, 0, 5), title="DR(j1,j2)", xTitle="\Delta R(j1,j2)"))
                    plots.append(Plot.make1D(f"h_{category}_deltaPhiJ1J2", op.deltaPhi(jets[0].p4,jets[1].p4) ,sel,EquidistantBinning(6, -3, 3), title=" Delta Phi", xTitle= "\Delta \phi (j1,j2)") )
                    plots.append(Plot.make1D(f"h_{category}_deltaEtaJ1J2", jets[0].p4.Eta() - jets[1].p4.Eta() ,  sel, EquidistantBinning(6, -3, 3), title=" Delta Eta", xTitle= "\Delta \eta (j1,j2)") )
                    plots.append(Plot.make1D(f"h_{category}_mJ1J2", op.invariant_mass(jets[0].p4, jets[1].p4) ,  sel, EquidistantBinning(50, 50., 500.), title=" mJ1J2", xTitle= "m(j1,j2) (GeV)" ))
                    plots.append(Plot.make1D(f"h_{category}_mJ1J2_narrow", op.invariant_mass(jets[0].p4, jets[1].p4) ,  sel, EquidistantBinning(50, 50., 200.), title=" mJ1J2 narrow", xTitle= "m(j1,j2) (GeV)"))
                if op.c_int(op.rng_len(jets), cast=True) == 0:    
                    plots.append(Plot.make1D(f"h_{category}_j1j2Eta0J", (j1+j2).Eta() ,  sel, EquidistantBinning(10, -2.5, 2.5), title=" j1j2 Eta", xTitle= "\eta (j1j2)"))
                    plots.append(Plot.make1D(f"h_{category}_j1j2Pt0J", (j1+j2).Pt() , sel, EquidistantBinning(50, 0., 700.), title=" j1j2 pT", xTitle= "pT(j1j2) (GeV)"))
                    plots.append(Plot.make1D(f"h_{category}_j1j2Phi0J", (j1+j2).Phi() ,  sel, EquidistantBinning(10, -2.5, 2.5), title=" j1j2 Phi", xTitle= "\phi (j1j2)"))
                    plots.append(Plot.make1D(f"h_{category}_deltaRJ1J20J", op.deltaR(jets[0].p4,jets[1].p4) ,sel,EquidistantBinning(10, 0, 5), title="DR(j1,j2)", xTitle="\Delta R(j1,j2)"))
                    plots.append(Plot.make1D(f"h_{category}_deltaPhiJ1J20J", op.deltaPhi(jets[0].p4,jets[1].p4) ,sel,EquidistantBinning(6, -3, 3), title=" Delta Phi", xTitle= "\Delta \phi (j1,j2)") )
                    plots.append(Plot.make1D(f"h_{category}_deltaEtaJ1J20J", jets[0].p4.Eta() - jets[1].p4.Eta() ,  sel, EquidistantBinning(6, -3, 3), title=" Delta Eta", xTitle= "\Delta \eta (j1,j2)") )
                    plots.append(Plot.make1D(f"h_{category}_mJ1J20J", op.invariant_mass(jets[0].p4, jets[1].p4) ,  sel, EquidistantBinning(50, 50., 500.), title=" mJ1J2", xTitle= "m(j1,j2) (GeV)" ))
                    plots.append(Plot.make1D(f"h_{category}_mJ1J2narrow0J", op.invariant_mass(jets[0].p4, jets[1].p4) ,  sel, EquidistantBinning(50, 50., 200.), title=" mJ1J2 narrow", xTitle= "m(j1,j2) (GeV)"))
                if op.c_int(op.rng_len(jets), cast=True) == 1:    
                    plots.append(Plot.make1D(f"h_{category}_j1j2Eta1J", (j1+j2).Eta() ,  sel, EquidistantBinning(10, -2.5, 2.5), title=" j1j2 Eta", xTitle= "\eta (j1j2)"))
                    plots.append(Plot.make1D(f"h_{category}_j1j2Pt1J", (j1+j2).Pt() , sel, EquidistantBinning(50, 0., 700.), title=" j1j2 pT", xTitle= "pT(j1j2) (GeV)"))
                    plots.append(Plot.make1D(f"h_{category}_j1j2Phi1J", (j1+j2).Phi() ,  sel, EquidistantBinning(10, -2.5, 2.5), title=" j1j2 Phi", xTitle= "\phi (j1j2)"))
                    plots.append(Plot.make1D(f"h_{category}_deltaRJ1J21J", op.deltaR(jets[0].p4,jets[1].p4) ,sel,EquidistantBinning(10, 0, 5), title="DR(j1,j2)", xTitle="\Delta R(j1,j2)"))
                    plots.append(Plot.make1D(f"h_{category}_deltaPhiJ1J21J", op.deltaPhi(jets[0].p4,jets[1].p4) ,sel,EquidistantBinning(6, -3, 3), title=" Delta Phi", xTitle= "\Delta \phi (j1,j2)") )
                    plots.append(Plot.make1D(f"h_{category}_deltaEtaJ1J21J", jets[0].p4.Eta() - jets[1].p4.Eta() ,  sel, EquidistantBinning(6, -3, 3), title=" Delta Eta", xTitle= "\Delta \eta (j1,j2)") )
                    plots.append(Plot.make1D(f"h_{category}_mJ1J21J", op.invariant_mass(jets[0].p4, jets[1].p4) ,  sel, EquidistantBinning(50, 50., 500.), title=" mJ1J2", xTitle= "m(j1,j2) (GeV)" ))
                    plots.append(Plot.make1D(f"h_{category}_mJ1J2narrow1J", op.invariant_mass(jets[0].p4, jets[1].p4) ,  sel, EquidistantBinning(50, 50., 200.), title=" mJ1J2 narrow", xTitle= "m(j1,j2) (GeV)"))
                if op.c_int(op.rng_len(jets), cast=True) >= 2:    
                    plots.append(Plot.make1D(f"h_{category}_j1j2Eta2J", (j1+j2).Eta() ,  sel, EquidistantBinning(10, -2.5, 2.5), title=" j1j2 Eta", xTitle= "\eta (j1j2)"))
                    plots.append(Plot.make1D(f"h_{category}_j1j2Pt2J", (j1+j2).Pt() , sel, EquidistantBinning(50, 0., 700.), title=" j1j2 pT", xTitle= "pT(j1j2) (GeV)"))
                    plots.append(Plot.make1D(f"h_{category}_j1j2Phi2J", (j1+j2).Phi() ,  sel, EquidistantBinning(10, -2.5, 2.5), title=" j1j2 Phi", xTitle= "\phi (j1j2)"))
                    plots.append(Plot.make1D(f"h_{category}_deltaRJ1J22J", op.deltaR(jets[0].p4,jets[1].p4) ,sel,EquidistantBinning(10, 0, 5), title="DR(j1,j2)", xTitle="\Delta R(j1,j2)"))
                    plots.append(Plot.make1D(f"h_{category}_deltaPhiJ1J22J", op.deltaPhi(jets[0].p4,jets[1].p4) ,sel,EquidistantBinning(6, -3, 3), title=" Delta Phi", xTitle= "\Delta \phi (j1,j2)") )
                    plots.append(Plot.make1D(f"h_{category}_deltaEtaJ1J22J", jets[0].p4.Eta() - jets[1].p4.Eta() ,  sel, EquidistantBinning(6, -3, 3), title=" Delta Eta", xTitle= "\Delta \eta (j1,j2)") )
                    plots.append(Plot.make1D(f"h_{category}_mJ1J22J", op.invariant_mass(jets[0].p4, jets[1].p4) ,  sel, EquidistantBinning(50, 50., 500.), title=" mJ1J2", xTitle= "m(j1,j2) (GeV)" ))
                    plots.append(Plot.make1D(f"h_{category}_mJ1J2narrow2J", op.invariant_mass(jets[0].p4, jets[1].p4) ,  sel, EquidistantBinning(50, 50., 200.), title=" mJ1J2 narrow", xTitle= "m(j1,j2) (GeV)"))
                if "BScore" in category:
                    if jets[0].btagDeepFlavB >=0:
                        plots.append(Plot.make1D(f"h_{category}_bscoreJ1", jets[0].btagDeepFlavB,sel,EquidistantBinning(20, 0, 1.0), title=" b tagging DeepJet score jet 1", xTitle= "btag score j1") )
                    if jets[1].btagDeepFlavB >=0:
                        plots.append(Plot.make1D(f"h_{category}_bscoreJ2", jets[1].btagDeepFlavB,sel,EquidistantBinning(20, 0, 1.0), title=" b tagging DeepJet score jet 2", xTitle= "btag score j2") )
            return plots

    def plotFourObjects(self,jets,pair,sel,category):
            plots = []
            lep1 = op.construct("ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float> >",([op.c_float(pair[0].p4.Pt()),op.c_float(pair[0].p4.Eta()),op.c_float(pair[0].p4.Phi()),op.c_float(pair[0].p4.M())]))
            lep2 = op.construct("ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float> >",([op.c_float(pair[1].p4.Pt()),op.c_float(pair[1].p4.Eta()),op.c_float(pair[1].p4.Phi()),op.c_float(pair[1].p4.M())]))
            if "SignalRegion" in category:
                plots.append(Plot.make1D(f"h_{category}_DeltaRL1J1", op.deltaR(lep1, jets[0].p4),sel, EquidistantBinning(12, 0, 6), title="DR(l1,j1)", xTitle="\Delta R(l1, j1)",plotopts={'no-data': True, 'show-ratio': False}))
                plots.append(Plot.make1D(f"h_{category}_DeltaPhiL1J1", op.deltaPhi(lep1, jets[0].p4),sel, EquidistantBinning(7, 0, 3.5), title="DeltaPhi(l1,j1)", xTitle="\Delta \phi (l1, j1)",plotopts={'no-data': True, 'show-ratio': False}))
                plots.append(Plot.make1D(f"h_{category}_deltaEtaL1J1", lep1.Eta() - jets[0].p4.Eta(),  sel, EquidistantBinning(6, -3, 3), title=" Delta Eta", xTitle= "\Delta \eta" ,plotopts={'no-data': True, 'show-ratio': False}))
                plots.append(Plot.make1D(f"h_{category}_DeltaRL2J2", op.deltaR(lep2, jets[1].p4),sel, EquidistantBinning(12, 0, 6), title="DR(l2,j2)", xTitle="\Delta R(l2, j2)",plotopts={'no-data': True, 'show-ratio': False}))
                plots.append(Plot.make1D(f"h_{category}_DeltaPhiL2J2", op.deltaPhi(lep2, jets[1].p4),sel, EquidistantBinning(7, 0, 3.5), title="DeltaPhi(l2,j2)", xTitle="\Delta \phi (l2, j2)",plotopts={'no-data': True, 'show-ratio': False}))
                plots.append(Plot.make1D(f"h_{category}_deltaEtaL2J2", lep2.Eta() - jets[1].p4.Eta(),  sel, EquidistantBinning(6, -3, 3), title=" Delta Eta", xTitle= "\Delta \eta" ,plotopts={'no-data': True, 'show-ratio': False}))
                plots.append(Plot.make1D(f"h_{category}_massLeptonsJets", op.invariant_mass(lep1, lep2, jets[0].p4,jets[1].p4) ,sel, EquidistantBinning(50, 150., 700.), title="m(j1,j2,l1,l2)", xTitle="m(j1,j2,l1,l2)",plotopts={'no-data': True, 'show-ratio': False}))
            else:
                if op.c_int(op.rng_len(jets), cast=True) >= 0:    
                    plots.append(Plot.make1D(f"h_{category}_massLeptonsJets", op.invariant_mass(lep1, lep2, jets[0].p4,jets[1].p4) ,sel, EquidistantBinning(50, 150., 700.), title="m(j1,j2,l1,l2)", xTitle="m(j1,j2,l1,l2)"))
                    plots.append(Plot.make1D(f"h_{category}_mll", op.invariant_mass(lep1,lep2), sel, EquidistantBinning(80, 40., 200.), title="Invariant Mass of two muons", xTitle= "m_{\mu\mu} (GeV)"))
                    plots.append(Plot.make1D(f"h_{category}_lep1lep2Pt", (lep1+lep2).Pt(), sel, EquidistantBinning(50, 0., 700.), title=" l1l2 pT", xTitle= "pT(l1l2) (GeV)"))
                    plots.append(Plot.make1D(f"h_{category}_DeltaRL1J1", op.deltaR(lep1, jets[0].p4),sel, EquidistantBinning(12, 0, 6), title="DR(l1,j1)", xTitle="\Delta R(l1, j1)"))
                    plots.append(Plot.make1D(f"h_{category}_DeltaPhiL1J1", op.deltaPhi(lep1, jets[0].p4),sel, EquidistantBinning(7, 0, 3.5), title="DeltaPhi(l1,j1)", xTitle="\Delta \phi (l1, j1)"))
                    plots.append(Plot.make1D(f"h_{category}_deltaEtaL1J1", lep1.Eta() - jets[0].p4.Eta(),  sel, EquidistantBinning(6, -3, 3), title=" Delta Eta", xTitle= "\Delta \eta (l1,j1)") )
                    plots.append(Plot.make1D(f"h_{category}_DeltaPhiL2J2", op.deltaPhi(lep2, jets[1].p4),sel, EquidistantBinning(7, 0, 3.5), title="DeltaPhi(l2,j2)", xTitle="\Delta \phi (l2, j2)"))
                    plots.append(Plot.make1D(f"h_{category}_DeltaRL2J2", op.deltaR(lep2, jets[1].p4),sel, EquidistantBinning(12, 0, 6), title="DR(l2,j2)", xTitle="\Delta R(l2, j2)"))
                    plots.append(Plot.make1D(f"h_{category}_deltaEtaL2J2", lep2.Eta() - jets[1].p4.Eta(),  sel, EquidistantBinning(6, -3, 3), title=" Delta Eta", xTitle= "\Delta \eta (l2,j2)") )
                if op.c_int(op.rng_len(jets), cast=True) == 0:    
                    plots.append(Plot.make1D(f"h_{category}_massLeptonsJets0J", op.invariant_mass(lep1, lep2, jets[0].p4,jets[1].p4) ,sel, EquidistantBinning(50, 150., 700.), title="m(j1,j2,l1,l2)", xTitle="m(j1,j2,l1,l2)"))
                    plots.append(Plot.make1D(f"h_{category}_mll0J", op.invariant_mass(lep1,lep2), sel, EquidistantBinning(80, 40., 200.), title="Invariant Mass of two muons", xTitle= "m_{\mu\mu} (GeV)"))
                    plots.append(Plot.make1D(f"h_{category}_lep1lep2Pt0J", (lep1+lep2).Pt(), sel, EquidistantBinning(50, 0., 700.), title=" l1l2 pT", xTitle= "pT(l1l2) (GeV)"))
                    plots.append(Plot.make1D(f"h_{category}_DeltaRL1J10J", op.deltaR(lep1, jets[0].p4),sel, EquidistantBinning(12, 0, 6), title="DR(l1,j1)", xTitle="\Delta R(l1, j1)"))
                    plots.append(Plot.make1D(f"h_{category}_DeltaPhiL1J10J", op.deltaPhi(lep1, jets[0].p4),sel, EquidistantBinning(7, 0, 3.5), title="DeltaPhi(l1,j1)", xTitle="\Delta \phi (l1, j1)"))
                    plots.append(Plot.make1D(f"h_{category}_deltaEtaL1J10J", lep1.Eta() - jets[0].p4.Eta(),  sel, EquidistantBinning(6, -3, 3), title=" Delta Eta", xTitle= "\Delta \eta (l1,j1)") )
                    plots.append(Plot.make1D(f"h_{category}_DeltaPhiL2J20J", op.deltaPhi(lep2, jets[1].p4),sel, EquidistantBinning(7, 0, 3.5), title="DeltaPhi(l2,j2)", xTitle="\Delta \phi (l2, j2)"))
                    plots.append(Plot.make1D(f"h_{category}_DeltaRL2J20J", op.deltaR(lep2, jets[1].p4),sel, EquidistantBinning(12, 0, 6), title="DR(l2,j2)", xTitle="\Delta R(l2, j2)"))
                    plots.append(Plot.make1D(f"h_{category}_deltaEtaL2J20J", lep2.Eta() - jets[1].p4.Eta(),  sel, EquidistantBinning(6, -3, 3), title=" Delta Eta", xTitle= "\Delta \eta (l2,j2)") )
                if op.c_int(op.rng_len(jets), cast=True) == 1:    
                    plots.append(Plot.make1D(f"h_{category}_massLeptonsJets1J", op.invariant_mass(lep1, lep2, jets[0].p4,jets[1].p4) ,sel, EquidistantBinning(50, 150., 700.), title="m(j1,j2,l1,l2)", xTitle="m(j1,j2,l1,l2)"))    
                    plots.append(Plot.make1D(f"h_{category}_mll1J", op.invariant_mass(lep1,lep2), sel, EquidistantBinning(80, 40., 200.), title="Invariant Mass of two muons", xTitle= "m_{\mu\mu} (GeV)"))
                    plots.append(Plot.make1D(f"h_{category}_lep1lep2Pt1J", (lep1+lep2).Pt(), sel, EquidistantBinning(50, 0., 700.), title=" l1l2 pT", xTitle= "pT(l1l2) (GeV)"))
                    plots.append(Plot.make1D(f"h_{category}_DeltaRL1J11J", op.deltaR(lep1, jets[0].p4),sel, EquidistantBinning(12, 0, 6), title="DR(l1,j1)", xTitle="\Delta R(l1, j1)"))
                    plots.append(Plot.make1D(f"h_{category}_DeltaPhiL1J11J", op.deltaPhi(lep1, jets[0].p4),sel, EquidistantBinning(7, 0, 3.5), title="DeltaPhi(l1,j1)", xTitle="\Delta \phi (l1, j1)"))
                    plots.append(Plot.make1D(f"h_{category}_deltaEtaL1J11J", lep1.Eta() - jets[0].p4.Eta(),  sel, EquidistantBinning(6, -3, 3), title=" Delta Eta", xTitle= "\Delta \eta (l1,j1)") )
                    plots.append(Plot.make1D(f"h_{category}_DeltaPhiL2J21J", op.deltaPhi(lep2, jets[1].p4),sel, EquidistantBinning(7, 0, 3.5), title="DeltaPhi(l2,j2)", xTitle="\Delta \phi (l2, j2)"))
                    plots.append(Plot.make1D(f"h_{category}_DeltaRL2J21J", op.deltaR(lep2, jets[1].p4),sel, EquidistantBinning(12, 0, 6), title="DR(l2,j2)", xTitle="\Delta R(l2, j2)"))
                    plots.append(Plot.make1D(f"h_{category}_deltaEtaL2J21J", lep2.Eta() - jets[1].p4.Eta(),  sel, EquidistantBinning(6, -3, 3), title=" Delta Eta", xTitle= "\Delta \eta (l2,j2)") )
                if op.c_int(op.rng_len(jets), cast=True) >= 2:    
                    plots.append(Plot.make1D(f"h_{category}_massLeptonsJets2J", op.invariant_mass(lep1, lep2, jets[0].p4,jets[1].p4) ,sel, EquidistantBinning(50, 150., 700.), title="m(j1,j2,l1,l2)", xTitle="m(j1,j2,l1,l2)"))
                    plots.append(Plot.make1D(f"h_{category}_mll2J", op.invariant_mass(lep1,lep2), sel, EquidistantBinning(80, 40., 200.), title="Invariant Mass of two muons", xTitle= "m_{\mu\mu} (GeV)"))
                    plots.append(Plot.make1D(f"h_{category}_lep1lep2Pt2J", (lep1+lep2).Pt(), sel, EquidistantBinning(50, 0., 700.), title=" l1l2 pT", xTitle= "pT(l1l2) (GeV)"))
                    plots.append(Plot.make1D(f"h_{category}_DeltaRL1J12J", op.deltaR(lep1, jets[0].p4),sel, EquidistantBinning(12, 0, 6), title="DR(l1,j1)", xTitle="\Delta R(l1, j1)"))
                    plots.append(Plot.make1D(f"h_{category}_DeltaPhiL1J12J", op.deltaPhi(lep1, jets[0].p4),sel, EquidistantBinning(7, 0, 3.5), title="DeltaPhi(l1,j1)", xTitle="\Delta \phi (l1, j1)"))
                    plots.append(Plot.make1D(f"h_{category}_deltaEtaL1J12J", lep1.Eta() - jets[0].p4.Eta(),  sel, EquidistantBinning(6, -3, 3), title=" Delta Eta", xTitle= "\Delta \eta (l1,j1)") )
                    plots.append(Plot.make1D(f"h_{category}_DeltaPhiL2J22J", op.deltaPhi(lep2, jets[1].p4),sel, EquidistantBinning(7, 0, 3.5), title="DeltaPhi(l2,j2)", xTitle="\Delta \phi (l2, j2)"))
                    plots.append(Plot.make1D(f"h_{category}_DeltaRL2J22J", op.deltaR(lep2, jets[1].p4),sel, EquidistantBinning(12, 0, 6), title="DR(l2,j2)", xTitle="\Delta R(l2, j2)"))
                    plots.append(Plot.make1D(f"h_{category}_deltaEtaL2J22J", lep2.Eta() - jets[1].p4.Eta(),  sel, EquidistantBinning(6, -3, 3), title=" Delta Eta", xTitle= "\Delta \eta (l2,j2)") )
            return plots
    def plotFourObjectsCorrected(self,jets,lep1,lep2,sel,category):
            plots = []
            if "SignalRegion" in category:
                plots.append(Plot.make1D(f"h_{category}_DeltaRL1J1", op.deltaR(lep1, jets[0].p4),sel, EquidistantBinning(12, 0, 6), title="DR(l1,j1)", xTitle="\Delta R(l1, j1)",plotopts={'no-data': True, 'show-ratio': False}))
                plots.append(Plot.make1D(f"h_{category}_DeltaPhiL1J1", op.deltaPhi(lep1, jets[0].p4),sel, EquidistantBinning(7, 0, 3.5), title="DeltaPhi(l1,j1)", xTitle="\Delta \phi (l1, j1)",plotopts={'no-data': True, 'show-ratio': False}))
                plots.append(Plot.make1D(f"h_{category}_deltaEtaL1J1", lep1.Eta() - jets[0].p4.Eta(),  sel, EquidistantBinning(6, -3, 3), title=" Delta Eta", xTitle= "\Delta \eta" ,plotopts={'no-data': True, 'show-ratio': False}))
                plots.append(Plot.make1D(f"h_{category}_DeltaRL2J2", op.deltaR(lep2, jets[1].p4),sel, EquidistantBinning(12, 0, 6), title="DR(l2,j2)", xTitle="\Delta R(l2, j2)",plotopts={'no-data': True, 'show-ratio': False}))
                plots.append(Plot.make1D(f"h_{category}_DeltaPhiL2J2", op.deltaPhi(lep2, jets[1].p4),sel, EquidistantBinning(7, 0, 3.5), title="DeltaPhi(l2,j2)", xTitle="\Delta \phi (l2, j2)",plotopts={'no-data': True, 'show-ratio': False}))
                plots.append(Plot.make1D(f"h_{category}_deltaEtaL2J2", lep2.Eta() - jets[1].p4.Eta(),  sel, EquidistantBinning(6, -3, 3), title=" Delta Eta", xTitle= "\Delta \eta" ,plotopts={'no-data': True, 'show-ratio': False}))
                plots.append(Plot.make1D(f"h_{category}_massLeptonsJets", op.invariant_mass(lep1, lep2, jets[0].p4,jets[1].p4) ,sel, EquidistantBinning(50, 150., 700.), title="m(j1,j2,l1,l2)", xTitle="m(j1,j2,l1,l2)",plotopts={'no-data': True, 'show-ratio': False}))
            else:
                if op.c_int(op.rng_len(jets), cast=True) >= 0:    
                    plots.append(Plot.make1D(f"h_{category}_massLeptonsJets", op.invariant_mass(lep1, lep2, jets[0].p4,jets[1].p4) ,sel, EquidistantBinning(50, 150., 700.), title="m(j1,j2,l1,l2)", xTitle="m(j1,j2,l1,l2)"))
                    plots.append(Plot.make1D(f"h_{category}_mll", op.invariant_mass(lep1,lep2), sel, EquidistantBinning(80, 40., 200.), title="Invariant Mass of two muons", xTitle= "m_{\mu\mu} (GeV)"))
                    plots.append(Plot.make1D(f"h_{category}_lep1lep2Pt", (lep1+lep2).Pt(), sel, EquidistantBinning(50, 0., 700.), title=" l1l2 pT", xTitle= "pT(l1l2) (GeV)"))
                    plots.append(Plot.make1D(f"h_{category}_DeltaRL1J1", op.deltaR(lep1, jets[0].p4),sel, EquidistantBinning(12, 0, 6), title="DR(l1,j1)", xTitle="\Delta R(l1, j1)"))
                    plots.append(Plot.make1D(f"h_{category}_DeltaPhiL1J1", op.deltaPhi(lep1, jets[0].p4),sel, EquidistantBinning(7, 0, 3.5), title="DeltaPhi(l1,j1)", xTitle="\Delta \phi (l1, j1)"))
                    plots.append(Plot.make1D(f"h_{category}_deltaEtaL1J1", lep1.Eta() - jets[0].p4.Eta(),  sel, EquidistantBinning(6, -3, 3), title=" Delta Eta", xTitle= "\Delta \eta (l1,j1)") )
                    plots.append(Plot.make1D(f"h_{category}_DeltaPhiL2J2", op.deltaPhi(lep2, jets[1].p4),sel, EquidistantBinning(7, 0, 3.5), title="DeltaPhi(l2,j2)", xTitle="\Delta \phi (l2, j2)"))
                    plots.append(Plot.make1D(f"h_{category}_DeltaRL2J2", op.deltaR(lep2, jets[1].p4),sel, EquidistantBinning(12, 0, 6), title="DR(l2,j2)", xTitle="\Delta R(l2, j2)"))
                    plots.append(Plot.make1D(f"h_{category}_deltaEtaL2J2", lep2.Eta() - jets[1].p4.Eta(),  sel, EquidistantBinning(6, -3, 3), title=" Delta Eta", xTitle= "\Delta \eta (l2,j2)") )
                if op.c_int(op.rng_len(jets), cast=True) == 0:    
                    plots.append(Plot.make1D(f"h_{category}_massLeptonsJets0J", op.invariant_mass(lep1, lep2, jets[0].p4,jets[1].p4) ,sel, EquidistantBinning(50, 150., 700.), title="m(j1,j2,l1,l2)", xTitle="m(j1,j2,l1,l2)"))
                    plots.append(Plot.make1D(f"h_{category}_mll0J", op.invariant_mass(lep1,lep2), sel, EquidistantBinning(80, 40., 200.), title="Invariant Mass of two muons", xTitle= "m_{\mu\mu} (GeV)"))
                    plots.append(Plot.make1D(f"h_{category}_lep1lep2Pt0J", (lep1+lep2).Pt(), sel, EquidistantBinning(50, 0., 700.), title=" l1l2 pT", xTitle= "pT(l1l2) (GeV)"))
                    plots.append(Plot.make1D(f"h_{category}_DeltaRL1J10J", op.deltaR(lep1, jets[0].p4),sel, EquidistantBinning(12, 0, 6), title="DR(l1,j1)", xTitle="\Delta R(l1, j1)"))
                    plots.append(Plot.make1D(f"h_{category}_DeltaPhiL1J10J", op.deltaPhi(lep1, jets[0].p4),sel, EquidistantBinning(7, 0, 3.5), title="DeltaPhi(l1,j1)", xTitle="\Delta \phi (l1, j1)"))
                    plots.append(Plot.make1D(f"h_{category}_deltaEtaL1J10J", lep1.Eta() - jets[0].p4.Eta(),  sel, EquidistantBinning(6, -3, 3), title=" Delta Eta", xTitle= "\Delta \eta (l1,j1)") )
                    plots.append(Plot.make1D(f"h_{category}_DeltaPhiL2J20J", op.deltaPhi(lep2, jets[1].p4),sel, EquidistantBinning(7, 0, 3.5), title="DeltaPhi(l2,j2)", xTitle="\Delta \phi (l2, j2)"))
                    plots.append(Plot.make1D(f"h_{category}_DeltaRL2J20J", op.deltaR(lep2, jets[1].p4),sel, EquidistantBinning(12, 0, 6), title="DR(l2,j2)", xTitle="\Delta R(l2, j2)"))
                    plots.append(Plot.make1D(f"h_{category}_deltaEtaL2J20J", lep2.Eta() - jets[1].p4.Eta(),  sel, EquidistantBinning(6, -3, 3), title=" Delta Eta", xTitle= "\Delta \eta (l2,j2)") )
                if op.c_int(op.rng_len(jets), cast=True) == 1:    
                    plots.append(Plot.make1D(f"h_{category}_massLeptonsJets1J", op.invariant_mass(lep1, lep2, jets[0].p4,jets[1].p4) ,sel, EquidistantBinning(50, 150., 700.), title="m(j1,j2,l1,l2)", xTitle="m(j1,j2,l1,l2)"))    
                    plots.append(Plot.make1D(f"h_{category}_mll1J", op.invariant_mass(lep1,lep2), sel, EquidistantBinning(80, 40., 200.), title="Invariant Mass of two muons", xTitle= "m_{\mu\mu} (GeV)"))
                    plots.append(Plot.make1D(f"h_{category}_lep1lep2Pt1J", (lep1+lep2).Pt(), sel, EquidistantBinning(50, 0., 700.), title=" l1l2 pT", xTitle= "pT(l1l2) (GeV)"))
                    plots.append(Plot.make1D(f"h_{category}_DeltaRL1J11J", op.deltaR(lep1, jets[0].p4),sel, EquidistantBinning(12, 0, 6), title="DR(l1,j1)", xTitle="\Delta R(l1, j1)"))
                    plots.append(Plot.make1D(f"h_{category}_DeltaPhiL1J11J", op.deltaPhi(lep1, jets[0].p4),sel, EquidistantBinning(7, 0, 3.5), title="DeltaPhi(l1,j1)", xTitle="\Delta \phi (l1, j1)"))
                    plots.append(Plot.make1D(f"h_{category}_deltaEtaL1J11J", lep1.Eta() - jets[0].p4.Eta(),  sel, EquidistantBinning(6, -3, 3), title=" Delta Eta", xTitle= "\Delta \eta (l1,j1)") )
                    plots.append(Plot.make1D(f"h_{category}_DeltaPhiL2J21J", op.deltaPhi(lep2, jets[1].p4),sel, EquidistantBinning(7, 0, 3.5), title="DeltaPhi(l2,j2)", xTitle="\Delta \phi (l2, j2)"))
                    plots.append(Plot.make1D(f"h_{category}_DeltaRL2J21J", op.deltaR(lep2, jets[1].p4),sel, EquidistantBinning(12, 0, 6), title="DR(l2,j2)", xTitle="\Delta R(l2, j2)"))
                    plots.append(Plot.make1D(f"h_{category}_deltaEtaL2J21J", lep2.Eta() - jets[1].p4.Eta(),  sel, EquidistantBinning(6, -3, 3), title=" Delta Eta", xTitle= "\Delta \eta (l2,j2)") )
                if op.c_int(op.rng_len(jets), cast=True) >= 2:    
                    plots.append(Plot.make1D(f"h_{category}_massLeptonsJets2J", op.invariant_mass(lep1, lep2, jets[0].p4,jets[1].p4) ,sel, EquidistantBinning(50, 150., 700.), title="m(j1,j2,l1,l2)", xTitle="m(j1,j2,l1,l2)"))
                    plots.append(Plot.make1D(f"h_{category}_mll2J", op.invariant_mass(lep1,lep2), sel, EquidistantBinning(80, 40., 200.), title="Invariant Mass of two muons", xTitle= "m_{\mu\mu} (GeV)"))
                    plots.append(Plot.make1D(f"h_{category}_lep1lep2Pt2J", (lep1+lep2).Pt(), sel, EquidistantBinning(50, 0., 700.), title=" l1l2 pT", xTitle= "pT(l1l2) (GeV)"))
                    plots.append(Plot.make1D(f"h_{category}_DeltaRL1J12J", op.deltaR(lep1, jets[0].p4),sel, EquidistantBinning(12, 0, 6), title="DR(l1,j1)", xTitle="\Delta R(l1, j1)"))
                    plots.append(Plot.make1D(f"h_{category}_DeltaPhiL1J12J", op.deltaPhi(lep1, jets[0].p4),sel, EquidistantBinning(7, 0, 3.5), title="DeltaPhi(l1,j1)", xTitle="\Delta \phi (l1, j1)"))
                    plots.append(Plot.make1D(f"h_{category}_deltaEtaL1J12J", lep1.Eta() - jets[0].p4.Eta(),  sel, EquidistantBinning(6, -3, 3), title=" Delta Eta", xTitle= "\Delta \eta (l1,j1)") )
                    plots.append(Plot.make1D(f"h_{category}_DeltaPhiL2J22J", op.deltaPhi(lep2, jets[1].p4),sel, EquidistantBinning(7, 0, 3.5), title="DeltaPhi(l2,j2)", xTitle="\Delta \phi (l2, j2)"))
                    plots.append(Plot.make1D(f"h_{category}_DeltaRL2J22J", op.deltaR(lep2, jets[1].p4),sel, EquidistantBinning(12, 0, 6), title="DR(l2,j2)", xTitle="\Delta R(l2, j2)"))
                    plots.append(Plot.make1D(f"h_{category}_deltaEtaL2J22J", lep2.Eta() - jets[1].p4.Eta(),  sel, EquidistantBinning(6, -3, 3), title=" Delta Eta", xTitle= "\Delta \eta (l2,j2)") )
            return plots
    def postProcess(self, taskList, config=None, workdir=None, resultsdir=None):
        super(BaseNanoHHtobbmumu, self).postProcess(taskList, config=config, workdir=workdir, resultsdir=resultsdir)
        if not self.plotList:
            self.plotList = self.getPlotList(resultsdir=resultsdir, config=config)
        from bamboo.plots import Plot, DerivedPlot, CutFlowReport
        plotList_cutflowreport = [ap for ap in self.plotList if isinstance(ap, CutFlowReport)]
        plotList_plotIt = [ap for ap in self.plotList
                           if (isinstance(ap, Plot) or isinstance(ap, DerivedPlot))
                           and len(ap.binnings) == 1]
        eraMode, eras = self.args.eras
        if eras is None:
            eras = list(config["eras"].keys())
        if plotList_cutflowreport:
            from bamboo.analysisutils import printCutFlowReports
            printCutFlowReports(
                config, plotList_cutflowreport, workdir=workdir, resultsdir=resultsdir,
                readCounters=self.readCounters, eras=(eraMode, eras), verbose=self.args.verbose)
        if plotList_plotIt:
            from bamboo.analysisutils import writePlotIt, runPlotIt
            import os
            cfgName = os.path.join(workdir, "plots.yml")
            writePlotIt(
                config, plotList_plotIt, cfgName, eras=eras, workdir=workdir, resultsdir=resultsdir,
                readCounters=self.readCounters, plotDefaults=self.plotDefaults,
                vetoFileAttributes=self.__class__.CustomSampleAttributes)
            runPlotIt(
                cfgName, workdir=workdir, plotIt=self.args.plotIt, eras=(eraMode, eras),
                verbose=self.args.verbose)
        
        from bamboo.plots import Skim
        skims = [ap for ap in self.plotList if isinstance(ap, Skim)]

        from bamboo.analysisutils import loadPlotIt
        p_config, samples, _, systematics, legend = loadPlotIt(config, [], eras=self.args.eras[1], workdir=workdir, resultsdir=resultsdir, readCounters=self.readCounters, vetoFileAttributes=self.__class__.CustomSampleAttributes)
        
        from bamboo.root import gbl
        import math
        
        ########################################################################################################
        ######## save a TTree with invariant mass of gg and bb for the different categories#####################
        ########################################################################################################
        import os.path
        from bamboo.plots import Skim
        skims = [ap for ap in self.plotList if isinstance(ap, Skim)]
        if self.args.mvaEval and skims:
            from bamboo.analysisutils import loadPlotIt
            p_config, samples, _, systematics, legend = loadPlotIt(config, [], eras=self.args.eras[1], workdir=workdir, resultsdir=resultsdir, readCounters=self.readCounters, vetoFileAttributes=self.__class__.CustomSampleAttributes)
            try:
                from bamboo.root import gbl
                import pandas as pd
                for skim in skims:
                    frames = []
                    for smp in samples:
                        print(smp.name)
                        for cb in (smp.files if hasattr(smp, "files") else [smp]):  # could be a helper in plotit
                            tree = cb.tFile.Get(skim.treeName)
                            if not tree:
                               print( f"KEY TTree {skim.treeName} does not exist, we are gonna skip this {smp}\n")
                            else:
                               cols = gbl.ROOT.RDataFrame(cb.tFile.Get(skim.treeName)).AsNumpy()
                               cols["weight"] *= cb.scale
                               cols["process"] = [smp.name]*len(cols["weight"])
                               frames.append(pd.DataFrame(cols))
                    df = pd.concat(frames)
                    df["process"] = pd.Categorical(df["process"], categories=pd.unique(df["process"]), ordered=False)
                    pqoutname = os.path.join(resultsdir, f"{skim.name}.parquet")
                    df.to_parquet(pqoutname)
                    logger.info(f"Dataframe for skim {skim.name} saved to {pqoutname}")
            except ImportError as ex:
                logger.error("Could not import pandas, no dataframes will be saved")

        ########################################################################################################
        ######## save a parquet file with the variables needed for MVA training ################################
        ########################################################################################################
        
        #mvaSkim
        if self.args.mvaSkim and skims:
            from bamboo.analysisutils import loadPlotIt
            try:
                from bamboo.root import gbl
                import pandas as pd
                for skim in skims:
                    frames = []
                    for smp in samples:
                        print("Parquet file creation")
                        print(smp.name)
                        for cb in (smp.files if hasattr(smp, "files") else [smp]):  # could be a helper in plotit
                            tree = cb.tFile.Get(skim.treeName)
                            if not tree:
                               print( f"KEY TTree {skim.treeName} does not exist, we are gonna skip this {smp}\n")
                            else:
                               cols = gbl.ROOT.RDataFrame(cb.tFile.Get(skim.treeName)).AsNumpy()
                               cols["weight"] *= cb.scale
                               cols["process"] = [smp.name]*len(cols["weight"])
                               frames.append(pd.DataFrame(cols))
                    df = pd.concat(frames)
                    df["process"] = pd.Categorical(df["process"], categories=pd.unique(df["process"]), ordered=False)
                    pqoutname = os.path.join(resultsdir, f"{skim.name}.parquet")
                    df.to_parquet(pqoutname)
                    logger.info(f"Dataframe for skim {skim.name} saved to {pqoutname}")
            except ImportError as ex:
                logger.error("Could not import pandas, no dataframes will be saved")




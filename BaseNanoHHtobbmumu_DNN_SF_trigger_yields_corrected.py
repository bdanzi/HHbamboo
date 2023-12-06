###------------------------------------------------------------------------------------###
###  This code needs envConfig cern.ini                                                ###
###  It creates object selection for bbmm analysis and produces root files             ###
###  Author: Brunella D'Anzi                                                           ###
###  Date: 20/08/2023                                                                  ###
###------------------------------------------------------------------------------------###

from bamboo.analysismodules import NanoAODHistoModule
from bamboo.treedecorators import NanoAODDescription
from bamboo.scalefactors import binningVariables_nano,lumiPerPeriod_default, BtagSF, get_scalefactor

from bamboo import treedecorators as td
from bamboo import treefunctions as op
from bamboo import scalefactors

from bamboo.plots import Plot, EquidistantBinning, CutFlowReport
from bamboo import treefunctions as op
from bamboo.treeoperations import Const

from itertools import chain
from functools import partial

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
    #-------------------------------------------------------------------------------------------#
    #                                       prepareTree                                         #
    #-------------------------------------------------------------------------------------------#
        
    def prepareTree(self, tree, sample=None, sampleCfg=None, backend=None):
        
        from bamboo.treedecorators import NanoAODDescription, nanoRochesterCalc, nanoJetMETCalc, nanoJetMETCalc_METFixEE2017, nanoFatJetCalc
        from bamboo.analysisutils import configureJets, makeMultiPrimaryDatasetTriggerSelection
        from bamboo.analysisutils import configureRochesterCorrection, makePileupWeight
        
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
        if not self.is_MC:
            subera = sample.split("_")[1]
            print("This is not a MC sample. Subera: ",subera )
        
        tree,noSel,be,lumiArgs = NanoAODHistoModule.prepareTree(self, tree, sample=sample, sampleCfg=sampleCfg, description=NanoAODDescription.get('v7', year=era, isMC=self.is_MC, systVariations=[nanoRochesterCalc,td.CalcCollectionsGroups(Jet=("pt", "mass"))]), backend=backend)
        
        if era == "2016":
            configureRochesterCorrection(variProxy  = tree._Muon,
                                         paramsFile = os.path.join(os.path.dirname(__file__),"RoccoR","RoccoR"+ subera +"UL.txt"),
                                         isMC       = self.is_MC,
                                         backend    = be, 
                                         uName      = sample)
        elif era == "2017":
            configureRochesterCorrection(variProxy  = tree._Muon,
                                         paramsFile = os.path.join(os.path.dirname(__file__),"RoccoR","RoccoR"+ era +"UL.txt"),
                                         isMC       = self.is_MC,
                                         backend    = be, 
                                         uName      = sample)
        elif era == "2018":
            configureRochesterCorrection(variProxy  = tree._Muon,
                                         paramsFile = os.path.join(os.path.dirname(__file__),"RoccoR","RoccoR"+ era +"UL.txt"),
                                         isMC       = self.is_MC,
                                         backend    = be, 
                                         uName      = sample)
        # Links : 
            # JEC : https://twiki.cern.ch/twiki/bin/viewauth/CMS/JECDataMC
            # JER (smear) : https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution
            
        if self.is_MC:
            JECTagDatabase = {"2016APV": "Summer19UL16APV_V7_MC", 
                              "2016": "Summer19UL16_V7_MC", 
                              "2017": "Summer19UL17_V5_MC", 
                              "2018": "Summer19UL18_V5_MC"}
            JERTagDatabase = {"2016APV": "Summer20UL16APV_JRV3_MC", 
                              "2016": "Summer20UL16_JRV3_MC", 
                              "2017": "Summer19UL17_JRV3_MC", 
                              "2018": "Summer19UL18_JRV2_MC"}
            if era in JECTagDatabase.keys():
                configureJets(
                    variProxy               = tree._Jet,
                    jetType                 = "AK4PFchs",
                    jec                     = JECTagDatabase[era],
                    smear                   = JERTagDatabase[era],
                    jecLevels               = "default",
                    regroupTag              = "V2",
                    jesUncertaintySources   = "All",
                    mayWriteCache           = self.args.distributed != "worker",
                    isMC                    = self.is_MC,
                    backend                 = be,
                    uName                   = sample
                    )
            else:
                    raise RuntimeError("Could not find appropriate JEC tag for MC")
        if not self.is_MC:
            JECTagDatabase = {
                              "2016APV": "Summer19UL16_RunBCDEF_V7_DATA", 
                              "2016postAPV": "Summer19UL16_RunFGH_V7_DATA", 
                              "2017B": "Summer19UL17_RunB_V5_DATA",
                              "2017C": "Summer19UL17_RunC_V5_DATA",
                              "2017D": "Summer19UL17_RunD_V5_DATA",
                              "2017E": "Summer19UL17_RunE_V5_DATA",
                              "2018A": "Summer19UL18_RunA_V5_DATA",
                              "2018B": "Summer19UL18_RunB_V5_DATA",
                              "2018C": "Summer19UL18_RunC_V5_DATA",
                              "2018D": "Summer19UL18_RunD_V5_DATA"
                              }
            for era in JECTagDatabase.keys():
                if era in sampleCfg['db'] or era in sampleCfg['db'][0]:
                    configureJets(
                        variProxy               = tree._Jet,
                        jetType                 = "AK4PFchs",
                        jec                     = JECTagDatabase[era],
                        jecLevels               = "default",
                        regroupTag              = "V2",
                        jesUncertaintySources   = "All",
                        mayWriteCache           = self.args.distributed != "worker",
                        isMC                    = self.is_MC,
                        backend                 = be,
                        uName                   = sample
                        )
                    
        return tree, noSel,be,lumiArgs
 
    def definePlots(self, tree, noSel, sample=None, sampleCfg=None):
        
        from bamboo.analysisutils import configureJets, makeMultiPrimaryDatasetTriggerSelection
        from bamboo.analysisutils import configureRochesterCorrection, makePileupWeight
        
        self.plotDefaults = {  "normalized": False,
                               "legend-columns": 2,
                               "no-data": False,
                               "show-overflow": True,    
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
        if era != "2016" and era != "2017" and era != "2018":
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
                
        #----- CutFlow report -----#
        self.yields = CutFlowReport("yields",printInLog=self.args.PrintYield,recursive=self.args.PrintYield)
        if self.args.PrintYield:
            self.yields.add(noSel, title="NoSelection")
            
        #----- Turn off systs -----#
        if self.args.NoSystematics:
            noSel = noSel.refine('SystOff',autoSyst=False)
            if self.args.PrintYield:
                self.yields.add(noSel, title="SystOff")
        
        #----- Theory uncertainties -----#
        # PS weights #
        if self.is_MC:
            self.psISRSyst = op.switch(op.rng_len(tree.PSWeight) == 4,
                                       op.systematic(op.c_float(1.), name="psISR", up=tree.PSWeight[2], down=tree.PSWeight[0]),
                                       op.systematic(op.c_float(1.), name="psISR", up=op.c_float(1.), down=op.c_float(1.)))
            self.psFSRSyst = op.switch(op.rng_len(tree.PSWeight) == 4,
                                       op.systematic(op.c_float(1.), name="psFSR", up=tree.PSWeight[3], down=tree.PSWeight[1]),
                                       op.systematic(op.c_float(1.), name="psFSR", up=op.c_float(1.), down=op.c_float(1.)))
            noSel = noSel.refine("PSweights", weight = [self.psISRSyst, self.psFSRSyst])
            if self.args.PrintYield:
                self.yields.add(noSel, title="PSweights")

        
        if era == "2016":
            # SingleMuon
            addHLTPath("SingleMuon","IsoMu22")
            addHLTPath("SingleMuon","IsoTkMu22")
            addHLTPath("SingleMuon","IsoMu22_eta2p1")
            addHLTPath("SingleMuon","IsoTkMu22_eta2p1")
            addHLTPath("SingleMuon","IsoMu24")
            addHLTPath("SingleMuon","IsoTkMu24")
            # SingleElectron
            addHLTPath("SingleElectron","Ele27_WPTight_Gsf")
            addHLTPath("SingleElectron","Ele25_eta2p1_WPTight_Gsf")
            addHLTPath("SingleElectron","Ele27_eta2p1_WPLoose_Gsf")
            # DoubleMuon
            addHLTPath("DoubleMuon","Mu17_TrkIsoVVL_Mu8_TrkIsoVVL")
            addHLTPath("DoubleMuon","Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ")
            addHLTPath("DoubleMuon","Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL")
            addHLTPath("DoubleMuon","Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ")
            # DoubleEGamma
            addHLTPath("DoubleEGamma","Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ")
            # MuonEG
            addHLTPath("MuonEG","Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL")
            addHLTPath("MuonEG","Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ")
            addHLTPath("MuonEG","Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL")
            addHLTPath("MuonEG","Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ")
           # addHLTPath("MuonEG","Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL") #additional one 
           # addHLTPath("MuonEG","Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL") #additional one 
           # addHLTPath("MuonEG","Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL") #additional one 
           # addHLTPath("MuonEG","Mu8_DiEle12_CaloIdL_TrackIdL") #additional one 
           # addHLTPath("MuonEG","DiMu9_Ele9_CaloIdL_TrackIdL")  #additional one 
            
        elif era == "2017":
            # SingleMuon #
            addHLTPath("SingleMuon","IsoMu24")
            addHLTPath("SingleMuon","IsoMu27") 
            # SingleElectron #
            addHLTPath("SingleElectron","Ele35_WPTight_Gsf")
            addHLTPath("SingleElectron","Ele32_WPTight_Gsf")
            #addHLTPath("SingleElectron","Ele38_WPTight_Gsf") #additional one
            #addHLTPath("SingleElectron","Ele40_WPTight_Gsf") #additional one
            # DoubleMuon #
            addHLTPath("DoubleMuon","Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8")
            addHLTPath("DoubleMuon","Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8")
            #addHLTPath("DoubleMuon","TripleMu_12_10_5") #additional one
            #addHLTPath("DoubleMuon","TripleMu_10_5_5_D2") #additional one
            # DoubleEGamma #
            addHLTPath("DoubleEGamma","Ele23_Ele12_CaloIdL_TrackIdL_IsoVL")
            #addHLTPath("DoubleEGamma","DoubleEle33_CaloIdL_GsfTrkIdVL") #additional one
            #addHLTPath("Ele16_Ele12_Ele8_CaloIdL_TrackIdL") #additional one
            # MuonEG #
            addHLTPath("MuonEG","Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ")
            addHLTPath("MuonEG","Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ")
            addHLTPath("MuonEG","Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL")
            addHLTPath("MuonEG","Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ")
            #addHLTPath("MuonEG","DiMu9_Ele9_CaloIdL_TrackIdL_DZ") #additional one
            #addHLTPath("MuonEG","Mu8_DiEle12_CaloIdL_TrackIdL") #additional one
            #addHLTPath("MuonEG","Mu8_DiEle12_CaloIdL_TrackIdL_DZ") #additional one
        elif era == "2018":
            # SingleMuon #
            addHLTPath("SingleMuon","IsoMu24")
            addHLTPath("SingleMuon","IsoMu27") 
            # DoubleMuon #
            addHLTPath("DoubleMuon","Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8")
            # EGamma composed by SingleElectron + DoubleEGamma #
            addHLTPath("EGamma","Ele23_Ele12_CaloIdL_TrackIdL_IsoVL")
            addHLTPath("EGamma","Ele32_WPTight_Gsf")
            #addHLTPath("EGamma","DoubleEle25_CaloIdL_MW") #additional one
            # MuonEG #
            addHLTPath("MuonEG","Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ")
            addHLTPath("MuonEG","Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ")
            addHLTPath("MuonEG","Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL")
            #addHLTPath("MuonEG","DiMu9_Ele9_CaloIdL_TrackIdL_DZ") #additional one
        
        #############################################################################
        #                            Pre-firing rates                               #
        #############################################################################
        if era in ["2016","2017","2018"] and self.is_MC and hasattr(tree,'L1PreFiringWeight_Nom'):
            self.L1Prefiring = op.systematic(tree.L1PreFiringWeight_Nom,
                                             name = "L1PreFiring",
                                             up   = tree.L1PreFiringWeight_Up,
                                             down = tree.L1PreFiringWeight_Dn)
            
            noSel = noSel.refine("L1PreFiringRate", weight = self.L1Prefiring)
            if self.args.PrintYield:
                self.yields.add(noSel,title="Prefiring")
        #############################################################################
        #                             Pile-up                                       #
        #############################################################################
        # Get MC PU weight file #
        puWeightsFile = None
        if self.is_MC:
            if 'related-sample' in sampleCfg.keys():
                puWeightsFile = os.path.join(os.path.dirname(__file__), "PileUpWeights",f'{sampleCfg["related-sample"]}_{era}.json')
                nameHint = f'puweightFromFile{sampleCfg["related-sample"]}'.replace('-','_')
            else:
                puWeightsFile = os.path.join(os.path.dirname(__file__), "PileUpWeights",f'puweights_{era}.json')
                nameHint = f'puweightFromFile{sample}'.replace('-','_')
            if not os.path.exists(puWeightsFile):
                raise RuntimeError("Could not find pileup file %s"%puWeightsFile)
            self.PUWeight = makePileupWeight(puWeightsFile, tree.Pileup_nTrueInt, systName="pileup",nameHint=nameHint)
            noSel = noSel.refine("puWeight", weight = self.PUWeight)
            if self.args.PrintYield:
                self.yields.add(noSel,title="puWeight")
        # Gen Weight and Triggers
        if self.is_MC: 
            noSel = noSel.refine('mcWeight', weight=tree.genWeight, cut=(
                op.OR(*chain.from_iterable(self.triggersPerPrimaryDataset.values()))))
        else:
            noSel = noSel.refine('trigger', cut=[makeMultiPrimaryDatasetTriggerSelection(
                sample, self.triggersPerPrimaryDataset)])
        if self.args.PrintYield:
                self.yields.add(noSel,title="TriggerFired")
                
        def hasAssociatedJet(lep): return lep.jet.idx != -1
        
        def muonDef(muons):
            return op.select(muons, lambda mu : op.AND(
                    mu.pt >= 5.,
                    op.abs(mu.eta) <= 2.4,
                    op.abs(mu.dxy) <= 0.05,
                    op.abs(mu.dz) <= 0.1,
                    mu.miniPFRelIso_all <= 0.4, # mini PF relative isolation, total (with scaled rho*EA PU corrections)
                    mu.sip3d <= 8,
                    mu.looseId
                    )                   
            )
        
        def muonConePt(muons):
            return op.map(muons, lambda lep: op.multiSwitch(
            (op.AND(op.abs(lep.pdgId) != 11, op.abs(lep.pdgId) != 13), lep.pt),
            (op.AND(op.abs(lep.pdgId) == 13, lep.mediumId, lep.mvaTTH > 0.50), lep.pt),
            0.9*lep.pt*(1.+lep.jetRelIso)
        ))


        def elDef(electrons):
            return op.select(electrons, lambda ele: op.AND(
                    ele.pt >= 7.,
                    op.abs(ele.eta) <= 2.5,
                    op.abs(ele.dxy) <= 0.05,
                    op.abs(ele.dz) <= 0.1,
                    ele.miniPFRelIso_all <= 0.4, # mini PF relative isolation, total (with scaled rho*EA PU cor
                    ele.sip3d <= 8,
                    ele.mvaFall17V2noIso_WPL, 
                    ele.lostHits <=1 # number of missing inner hits
                        )   
                    )

        def elConePt(electrons):
            return op.map(electrons, lambda lep: op.multiSwitch(
                (op.AND(op.abs(lep.pdgId) != 11, op.abs(lep.pdgId) != 13), lep.pt),
                (op.AND(op.abs(lep.pdgId) == 11, lep.mvaTTH > 0.30), lep.pt),
                0.9*lep.pt*(1.+lep.jetRelIso)
            ))

        def lepton_associatedJetLessThanMediumBtag(lep, era): 
            return op.OR(op.AND(era=="2016",op.OR(op.NOT(hasAssociatedJet(lep)), lep.jet.btagDeepFlavB <= 0.3093)),
                         op.AND(era=="2017",op.OR(op.NOT(hasAssociatedJet(lep)), lep.jet.btagDeepFlavB <= 0.3033)),
                         op.AND(era=="2018",op.OR(op.NOT(hasAssociatedJet(lep)), lep.jet.btagDeepFlavB <= 0.2770)))  


        def lepton_associatedJetLessThanTightBtag(lep,era): 
            return op.OR(op.AND(era=="2016",op.OR(op.NOT(hasAssociatedJet(lep)), lep.jet.btagDeepFlavB <= 0.7221)),
                         op.AND(era=="2017",op.OR(op.NOT(hasAssociatedJet(lep)), lep.jet.btagDeepFlavB <= 0.7489)),
                         op.AND(era=="2018",op.OR(op.NOT(hasAssociatedJet(lep)), lep.jet.btagDeepFlavB <= 0.7264)))
        
        def muon_x(mu): 
            return op.min(op.max(0., (0.9*mu.pt*(1+mu.jetRelIso))-20.)/(45.-20.), 1.)


        def muon_btagInterpolation(mu,era):
            if era=="2016":
                return muon_x(mu)*0.0614 + (1-muon_x(mu))*0.3093
            elif era=="2017":
                return muon_x(mu)*0.0521 + (1-muon_x(mu))*0.3033
            elif era=="2018":
                return muon_x(mu)*0.0494 + (1-muon_x(mu))*0.2770
        
        def muon_deepJetInterpIfMvaFailed(mu,era): 
            return op.OR(op.NOT(hasAssociatedJet(mu)), mu.jet.btagDeepFlavB < muon_btagInterpolation(mu,era))
        
        def muonFakeSel(muons,era):
            return op.select(muons, lambda mu: op.AND(
            muonConePt(muons)[mu.idx] >= 10.,
            op.OR(lepton_associatedJetLessThanMediumBtag(mu,era), op.AND(mu.jetRelIso < 0.8, muon_deepJetInterpIfMvaFailed(mu,era))))
            )

        def muonTightSel(muons,era): 
            return op.select(muons, lambda mu: op.AND(
                    muonConePt(muons)[mu.idx] >= 10.,
                    lepton_associatedJetLessThanMediumBtag(mu,era),
                    mu.mvaTTH >= 0.50,
                    mu.mediumId
            ))

        def elFakeSel(electrons,era):
            return op.select(electrons, lambda el: op.AND(
                elConePt(electrons)[el.idx] >= 10,
                op.OR(
                    op.AND(op.abs(el.eta+el.deltaEtaSC) <= 1.479, el.sieie <= 0.011),
                    op.AND(op.abs(el.eta+el.deltaEtaSC) > 1.479, el.sieie <= 0.030)
                ),
                el.hoe <= 0.10,
                el.eInvMinusPInv >= -0.04,
                op.OR(el.mvaTTH >= 0.30, op.AND(el.jetRelIso < 0.7, el.mvaFall17V2noIso_WP90)), # Lepton MVA id from ttH
                op.switch(
                    el.mvaTTH < 0.30,
                    lepton_associatedJetLessThanTightBtag(el,era),
                    lepton_associatedJetLessThanMediumBtag(el,era)),
                el.lostHits == 0,  # number of missing inner hits
                el.convVeto # Passes conversion veto
                ))


        def elTightSel(electrons,era): 
            return op.select(electrons, lambda el: op.AND(
            elConePt(electrons)[el.idx] >= 10.,
            op.OR(
                op.AND(op.abs(el.eta+el.deltaEtaSC) <= 1.479, el.sieie <= 0.011),
                op.AND(op.abs(el.eta+el.deltaEtaSC) > 1.479, el.sieie <= 0.030)
            ),
            el.hoe <= 0.10,
            el.eInvMinusPInv >= -0.04,
            el.convVeto,
            el.mvaTTH >= 0.30,
            el.lostHits == 0,
            lepton_associatedJetLessThanMediumBtag(el,era),
            ))

        def cleanElectrons(electrons, muons):
            return op.select(electrons, lambda el: op.NOT(
                op.rng_any( # test if any item in a range passes a selection
                    muons, lambda mu: op.deltaR(el.p4, mu.p4) <= 0.3))
            )
        
        def cleanMuons(electrons, muons):
            return op.select(muons, lambda mu: op.NOT(
                op.rng_any( # test if any item in a range passes a selection
                    electrons, lambda el: op.deltaR(el.p4, mu.p4) <= 0.3))
            )


        def ak4jetDef(jets,era):
            return op.select(jets, lambda jet: op.AND(
                jet.jetId & 1 if era == "2016" else jet.jetId & 2,  # Jet ID flags bit1 is loose, bit2 is tight, bit3 is tightLepVeto
                jet.pt >= 25.,
                op.abs(jet.eta) <= 2.4,
                op.OR(((jet.puId >> 2) & 1), jet.pt > 50.) # Jet PU ID bit1 is loose # no puId in Run3 so far
            ))

        def cleaningWithRespectToLeadingLeptons(jets,electrons, muons, DR):
            return op.select(jets, lambda j: 
               op.multiSwitch(
                # If Only electrons check
                (op.AND(op.rng_len(electrons) >= 2, op.rng_len(muons) == 0), op.AND(op.deltaR(j.p4, electrons[0].p4) >= DR, op.deltaR(j.p4, electrons[1].p4) >= DR)),
                # Elif only muons
                (op.AND(op.rng_len(electrons) == 0, op.rng_len(muons) >= 2), op.AND(op.deltaR(j.p4, muons[0].p4) >= DR, op.deltaR(j.p4, muons[1].p4) >= DR)),
                # Elif one electron + one muon
                (op.AND(op.rng_len(electrons) == 1, op.rng_len(muons) == 1), op.AND(op.deltaR(j.p4, electrons[0].p4) >= DR, op.deltaR(j.p4, muons[0].p4) >= DR)),
                # Elif at least one electron + at least one muon
                (op.AND(op.rng_len(electrons) >= 1, op.rng_len(muons) >= 1), op.switch(
                # Elif Electron is the leading lepton
                elConePt(electrons)[0] > muonConePt(muons)[0],
                op.switch(
                    # if one electron
                    op.rng_len(electrons) == 1,
                    op.AND(op.deltaR(j.p4, electrons[0].p4) >= DR, op.deltaR(j.p4, muons[0].p4) >= DR),
                    # elif more than one electron     
                    op.switch(
                        elConePt(electrons)[1] > muonConePt(muons)[0],
                        op.AND(op.deltaR(j.p4, electrons[0].p4) >= DR, op.deltaR(j.p4, electrons[1].p4) >= DR),
                                op.AND(op.deltaR(j.p4, electrons[0].p4) >= DR, op.deltaR(j.p4, muons[0].p4) >= DR)
                            )
                        ),
                # Muon is the leading lepton
                op.switch(
                     # if one muon
                    op.rng_len(muons) == 1,
                    op.AND(op.deltaR(j.p4, muons[0].p4) >= DR, op.deltaR(j.p4, electrons[0].p4) >= DR),
                    op.switch(
                        muonConePt(muons)[1] > elConePt(electrons)[0],
                        op.AND(op.deltaR(j.p4, muons[0].p4) >= DR, op.deltaR(j.p4, muons[1].p4) >= DR),
                                op.AND(op.deltaR(j.p4, muons[0].p4) >= DR, op.deltaR(j.p4, electrons[0].p4) >= DR)))
                            )
                        ), 
                op.c_bool(True)
                )
            )

        plots = []

        # lepton cone-pt definitions
        self.muon_conept = muonConePt(tree.Muon)
        self.electron_conept = elConePt(tree.Electron)
        
        # lepton definitions sorted by their cone-pt
        self.muons = op.sort(muonDef(tree.Muon), lambda mu: -self.muon_conept[mu.idx])
        self.electrons = op.sort(elDef(tree.Electron), lambda el: -self.electron_conept[el.idx])
        
        # cleaning electrons wrt muons
        self.cleanedElectrons = cleanElectrons(self.electrons,self.muons)
        cleanedElectrons = self.cleanedElectrons
        # cleaning muons wrt electrons
        self.cleanedMuons = cleanMuons(self.electrons,self.muons)
        cleanedMuons = self.cleanedMuons
        
        # Fakeable leptons
        self.fakeMuons = muonFakeSel(self.muons,era)
        self.fakeElectrons = elFakeSel(self.electrons,era)
        
        # tight leptons
        self.tightMuons = muonTightSel(self.cleanedMuons,era)
        self.tightElectrons = elTightSel(self.cleanedElectrons,era)

        # Ak4 jets sorted by their pt
        self.ak4JetsByPt = op.sort(tree.Jet, lambda jet: -jet.pt)
        # Preselection #
        self.lambda_ak4JetsPreSel = lambda j : op.AND(j.jetId & 1 if era == "2016" else j.jetId & 2, # Jet ID flags bit1 is loose, bit2 is tight, bit3 is tightLepVeto
                                                      j.pt >= 25.,
                                                      op.abs(j.eta) <= 2.4)
        # AK4 Jets sorted by their pt and selected
        ak4JetsPreSel = ak4jetDef(self.ak4JetsByPt,era)
        
        def returnLambdaCleaningWithRespectToLeadingLeptons(DR):
                return lambda j : op.multiSwitch(
                      (op.AND(op.rng_len(self.fakeElectrons) >= 2,op.rng_len(self.fakeMuons) == 0), 
                          # Only electrons 
                          op.AND(op.deltaR(j.p4, self.fakeElectrons[0].p4)>=DR, op.deltaR(j.p4, self.fakeElectrons[1].p4)>=DR)),
                      (op.AND(op.rng_len(self.fakeElectrons) == 0,op.rng_len(self.fakeMuons) >= 2), 
                          # Only muons  
                          op.AND(op.deltaR(j.p4, self.fakeMuons[0].p4)>=DR, op.deltaR(j.p4, self.fakeMuons[1].p4)>=DR)),
                      (op.AND(op.rng_len(self.fakeElectrons) == 1,op.rng_len(self.fakeMuons) == 1),
                          # One electron + one muon
                          op.AND(op.deltaR(j.p4, self.fakeElectrons[0].p4)>=DR, op.deltaR(j.p4, self.fakeMuons[0].p4)>=DR)),
                      (op.AND(op.rng_len(self.fakeElectrons) >= 1,op.rng_len(self.fakeMuons) >= 1),
                          # At least one electron + at least one muon
                       op.switch(self.electron_conept[self.fakeElectrons[0].idx] > self.muon_conept[self.fakeMuons[0].idx],
                                 # Electron is leading #
                                 op.switch(op.rng_len(self.fakeElectrons) == 1,
                                           op.AND(op.deltaR(j.p4, self.fakeElectrons[0].p4)>=DR, op.deltaR(j.p4, self.fakeMuons[0].p4)>=DR),
                                           op.switch(self.electron_conept[self.fakeElectrons[1].idx] > self.muon_conept[self.fakeMuons[0].idx],
                                                     op.AND(op.deltaR(j.p4, self.fakeElectrons[0].p4)>=DR, op.deltaR(j.p4, self.fakeElectrons[1].p4)>=DR),
                                                     op.AND(op.deltaR(j.p4, self.fakeElectrons[0].p4)>=DR, op.deltaR(j.p4, self.fakeMuons[0].p4)>=DR))),
                                 # Muon is leading #
                                 op.switch(op.rng_len(self.fakeMuons) == 1,
                                           op.AND(op.deltaR(j.p4, self.fakeMuons[0].p4)>=DR, op.deltaR(j.p4, self.fakeElectrons[0].p4)>=DR),
                                           op.switch(self.muon_conept[self.fakeMuons[1].idx] > self.electron_conept[self.fakeElectrons[0].idx],
                                                     op.AND(op.deltaR(j.p4, self.fakeMuons[0].p4)>=DR, op.deltaR(j.p4, self.fakeMuons[1].p4)>=DR),
                                                     op.AND(op.deltaR(j.p4, self.fakeMuons[0].p4)>=DR, op.deltaR(j.p4, self.fakeElectrons[0].p4)>=DR))))),
                       op.c_bool(True))
        self.lambda_cleanAk4Jets = returnLambdaCleaningWithRespectToLeadingLeptons(0.4)
        # remove jets within cone of DR<0.4 of leading lept
        self.ak4Jets = cleaningWithRespectToLeadingLeptons(ak4JetsPreSel, self.fakeElectrons, self.fakeMuons, 0.4) # Pt ordered
        ak4Jets = self.ak4Jets
        ############     Btagging     #############
        # The pfDeepFlavour (DeepJet) algorithm is used
        if era == "2016": 
            self.lambda_ak4BtagLoose =   lambda jet    : jet.btagDeepFlavB > 0.0614
            self.lambda_ak4Btag =   lambda jet    : jet.btagDeepFlavB > 0.3093
            self.lambda_ak4NoBtag = lambda jet    : jet.btagDeepFlavB <= 0.3093
        elif era =="2017":
            self.lambda_ak4BtagLoose =   lambda jet    : jet.btagDeepFlavB > 0.0521
            self.lambda_ak4Btag =   lambda jet    : jet.btagDeepFlavB > 0.3033
            self.lambda_ak4NoBtag = lambda jet    : jet.btagDeepFlavB <= 0.3033
        elif era == "2018":
            self.lambda_ak4BtagLoose =   lambda jet    : jet.btagDeepFlavB > 0.0494
            self.lambda_ak4Btag =   lambda jet    : jet.btagDeepFlavB > 0.2770
            self.lambda_ak4NoBtag = lambda jet    : jet.btagDeepFlavB <= 0.2770

        self.ak4JetsByBtagScore = op.sort(self.ak4Jets,lambda j: -j.btagDeepFlavB) # Btag score ordered
        # Ak4 Jets selected with loose requirement on bscore and ordered by it
        self.ak4BJets = op.select(self.ak4Jets, self.lambda_ak4Btag)
        self.ak4BJetsLoose = op.select(self.ak4Jets, self.lambda_ak4BtagLoose)
        self.ak4LightJetsByPt = op.select(self.ak4Jets, self.lambda_ak4NoBtag)
        self.ak4LightJetsByBtagScore = op.sort(self.ak4LightJetsByPt, lambda jet : -jet.btagDeepFlavB)
        
        cleanedHighestBtagScoreJets = self.ak4JetsByBtagScore[:op.min(op.rng_len(self.ak4JetsByBtagScore),op.static_cast("std::size_t",op.c_int(2)))]         
        cleanedHighestPtJets = self.ak4Jets[:op.min(op.rng_len(self.ak4Jets),op.static_cast("std::size_t",op.c_int(2)))] 
        
        #############################################################################
        #                                 PU Jets                                   #
        #############################################################################
        # Correction to be computed later on all the ak4 jets (regular + VBF) on which the cut is applied 
        # To avoid double counting : OR(regular jet, VBF jet)
        self.ak4ForPUID = op.select(self.ak4JetsByPt, lambda j : op.AND(j.pt<=50.,
                                                                        op.OR(self.lambda_ak4JetsPreSel(j),
                                                                        self.lambda_cleanAk4Jets(j))))
        
        #############################################################################
        #                             Scalefactors                                  #
        #############################################################################
        
        mZ = 91.1876
        mH = 125.1
        
        MuMuLooseSel = op.combine(self.cleanedMuons, N=2, pred=lambda l1,l2 : op.AND(l1.charge != l2.charge, op.OR(op.AND(l1.pt > 25., l2.pt> 15.),op.AND(l1.pt> 15., l2.pt > 25.))))
        ElElLooseSel = op.combine(self.cleanedElectrons, N=2, pred=lambda l1,l2 : op.AND(l1.charge != l2.charge, op.OR(op.AND(l1.pt > 25., l2.pt> 15.),op.AND(l1.pt> 15., l2.pt > 25.))))
        ElMuLooseSel = op.combine((self.cleanedElectrons, self.cleanedMuons), N=2, pred=lambda l1,l2 : op.AND(l1.charge != l2.charge, op.OR(op.AND(l1.pt > 25., l2.pt> 15.),op.AND(l1.pt> 15., l2.pt > 25.))))
        
        # OS tight dilepton collections
        ElElTightSel = op.combine(self.tightElectrons, N=2, pred=lambda lep1, lep2 : lep1.charge != lep2.charge)
        MuMuTightSel = op.combine(self.tightMuons, N=2, pred=lambda lep1, lep2 : lep1.charge != lep2.charge)
        ElMuTightSel = op.combine((self.tightElectrons, self.tightMuons), N=2, pred= lambda el, mu : el.charge != mu.charge)
        
        def lowMllCut(dileptons): return op.NOT(op.rng_any(
            dileptons, lambda dilep: op.invariant_mass(dilep[0].p4, dilep[1].p4) < 12.))
        
        # low Mll cut : reject events with dilepton mass below 12 GeV
        mllCut = op.AND(lowMllCut(ElElLooseSel), lowMllCut(MuMuLooseSel), lowMllCut(ElMuLooseSel))
        #mllCut = op.AND(lowMllCut(MuMuLooseSel))
        leptonMultiplicityCut_mumu = op.AND(op.rng_len(ElElTightSel) == 0, op.rng_len(MuMuTightSel) == 1, op.rng_len(ElMuTightSel) == 0)
        
        #--> Select the best Z asking that the invariant mass is the closest to mZ, same for Higgs
        bestZmm = op.rng_min_element_by(MuMuTightSel, lambda ll : op.abs(op.invariant_mass(ll[0].p4,ll[1].p4)-mZ))
        bestHmm = op.rng_min_element_by(MuMuTightSel, lambda ll : op.abs(op.invariant_mass(ll[0].p4,ll[1].p4)-mH))
        
        a = op.c_float(-99, typeName='double', cast=None)
        
        def outZ(dilep): # starting from a single selected pair 
            return op.NOT(op.abs(op.invariant_mass(dilep[0].p4, dilep[1].p4) - mZ) <= 10.)
        
        def outH(dilep): # starting from a single selected pair 
            #print("Invariant mass:",op.c_long(op.abs(op.invariant_mass(dilep[0].p4, dilep[1].p4)),cast=True) )
            return op.NOT(op.abs(op.invariant_mass(dilep[0].p4, dilep[1].p4) - mH) <= 10.)
        
        def out_of_range_mass(jets, min,max):
            return op.OR(op.switch(op.rng_len(jets) > 1, op.invariant_mass(jets[0].p4,jets[1].p4), a) < min, op.switch(op.rng_len(jets) > 1, op.invariant_mass(jets[0].p4,jets[1].p4), a) > max )
        
        min_mass_bjets = 70.
        max_mass_bjets = 140.
        # Ask to have a good Z with the refine option
        # Higgs/ Z to Muons control region
        hasZmm      = noSel.refine("hasZmm",  cut=[ mllCut, leptonMultiplicityCut_mumu,
                                                   op.rng_len(cleanedHighestBtagScoreJets) > 1,
                                                   op.NOT(outZ(bestZmm)),
                                                   out_of_range_mass(cleanedHighestBtagScoreJets,min_mass_bjets,max_mass_bjets)])
        
        hasHmm      = noSel.refine("hasHmm",  cut=[ mllCut, leptonMultiplicityCut_mumu,
                                                   op.rng_len(cleanedHighestBtagScoreJets) > 1,
                                                   op.NOT(outH(bestHmm)),
                                                   out_of_range_mass(cleanedHighestBtagScoreJets,min_mass_bjets,max_mass_bjets) ])
        hasZjets    = noSel.refine("hasZjets",  cut=[ mllCut, leptonMultiplicityCut_mumu,
                                                   op.rng_len(cleanedHighestBtagScoreJets) > 1,
                                                   outZ(bestZmm),
                                                   op.NOT(out_of_range_mass(cleanedHighestBtagScoreJets,min_mass_bjets,max_mass_bjets)) ])
        
        hasHjets    = noSel.refine("hasHjets",  cut=[ mllCut, leptonMultiplicityCut_mumu,
                                                   op.rng_len(cleanedHighestBtagScoreJets) > 1,
                                                   op.NOT(outH(bestHmm)),
                                                   out_of_range_mass(cleanedHighestBtagScoreJets,min_mass_bjets,max_mass_bjets) ])
        isSideband2lbbZmm = noSel.refine("isSideband2lbbZmm", cut=[ mllCut, leptonMultiplicityCut_mumu,
                                                   op.rng_len(cleanedHighestBtagScoreJets) > 1,
                                                   outZ(bestZmm),
                                                   out_of_range_mass(cleanedHighestBtagScoreJets,min_mass_bjets,max_mass_bjets)])
        isSideband2lbbHmm = noSel.refine("isSideband2lbbHmm",  cut=[ mllCut,leptonMultiplicityCut_mumu,
                                                   op.rng_len(cleanedHighestBtagScoreJets) > 1,
                                                   outH(bestHmm),
                                                   out_of_range_mass(cleanedHighestBtagScoreJets,min_mass_bjets,max_mass_bjets) ])
        isSignalRegion = noSel.refine("isSignalRegion", 
                                                   cut=[
                                                   mllCut, leptonMultiplicityCut_mumu,
                                                   op.rng_len(cleanedHighestBtagScoreJets) > 1,
                                                   op.NOT(outH(bestHmm)),
                                                   op.NOT(out_of_range_mass(cleanedHighestBtagScoreJets,min_mass_bjets,max_mass_bjets)) ])
        
        plots = []

        ### Save mvaVariables to be retrieved later in the postprocessor and saved in a parquet file ###
        if self.args.mvaSkim:
            ### Here you can specify the variables that you want to save in the .parquet file, you need to add --mvaSkim to the command line ###
            if self.is_MC:
                #genPartFlav
                print("This is a MC sample:",self.is_MC)
                lep1_Hcand_genPartFlav = op.c_int(bestHmm[0].genPartFlav, cast=True)
                lep2_Hcand_genPartFlav = op.c_int(bestHmm[1].genPartFlav, cast=True)
            else: 
                print("This is NOT a MC sample:")
                lep1_Hcand_genPartFlav = a
                lep2_Hcand_genPartFlav = a
            mvaVariables = {
                "weight"                       : noSel.weight,
                "nmuons_cleaned"               : op.c_int(op.rng_len(cleanedMuons), cast=True),
                "nelectrons_cleaned"           : op.c_int(op.rng_len(cleanedElectrons), cast=True),
                "nmuonPairs_tight"             : op.c_int(op.rng_len(MuMuTightSel), cast=True),
                "nelectrons_tight"             : op.c_int(op.rng_len(ElElTightSel), cast=True),
                # njets      
                "njets_cleaned"                : op.c_int(op.rng_len(ak4Jets), cast=True),
                # l1
                "lep1_Hcand_genPartFlav" : lep1_Hcand_genPartFlav,
                "lep1_Hcand_px"          : bestHmm[0].p4.Px(),  
                "lep1_Hcand_py"          : bestHmm[0].p4.Py(),
                "lep1_Hcand_pt"          : bestHmm[0].p4.Pt(),
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
            if self.is_MC:
                self.yields.add(isSignalRegion,title="SignalRegion")
        
        # Plots about leptons in Hmm/Zmm and sideband control regions 
        plots += self.plotPairs(bestZmm, hasZmm, "CRZmm_cleanedMuons",self.cleanedMuons)
        plots += self.plotPairs(bestHmm, hasHmm, "CRHmm_cleanedMuons",self.cleanedMuons)
        plots += self.plotPairs(bestZmm, isSideband2lbbZmm,"Sideband2lbbZmm_cleanedMuons" ,self.cleanedMuons)
        plots += self.plotPairs(bestHmm, isSideband2lbbHmm,"Sideband2lbbHmm_cleanedMuons" ,self.cleanedMuons)
        
        plots += self.plotPairs(bestZmm, hasZmm, "CRZmm_tightMuons",self.tightMuons)
        plots += self.plotPairs(bestHmm, hasHmm, "CRHmm_tightMuons",self.tightMuons)
        plots += self.plotPairs(bestZmm, isSideband2lbbZmm,"Sideband2lbbZmm_tightMuons" ,self.tightMuons)
        plots += self.plotPairs(bestHmm, isSideband2lbbHmm,"Sideband2lbbHmm_tightMuons" ,self.tightMuons)
        
        # Plots about jets in the Zjets and Hjets control regions + sidebands
        plots += self.plotJets(cleanedHighestPtJets, hasZjets, "Z_jets_CR_PtOrdering")
        plots += self.plotJets(cleanedHighestBtagScoreJets, hasZjets, "Z_jets_CR_BScoreOrdering")
        plots += self.plotJets(cleanedHighestBtagScoreJets, hasHjets, "H_jets_CR_PtOrdering")
        plots += self.plotJets(cleanedHighestBtagScoreJets, hasHjets, "H_jets_CR_BScoreOrdering")
        
        plots += self.plotJets(cleanedHighestPtJets, isSideband2lbbZmm, "Sideband2lbbZmm_CR_PtOrdering")
        plots += self.plotJets(cleanedHighestBtagScoreJets, isSideband2lbbZmm, "Sideband2lbbZmm_CR_BScoreOrdering")
        plots += self.plotJets(cleanedHighestBtagScoreJets, isSideband2lbbHmm, "Sideband2lbbHmm_CR_PtOrdering")
        plots += self.plotJets(cleanedHighestBtagScoreJets, isSideband2lbbHmm, "Sideband2lbbHmm_CR_BScoreOrdering")
    
        # Plots about leptons and jets in the sidebands (jets,pair,sel,category)
        plots += self.plotFourObjects(cleanedHighestPtJets,bestZmm,isSideband2lbbZmm,"Sideband2lbbZmm_bestZmm_PtOrdering")
        plots += self.plotFourObjects(cleanedHighestBtagScoreJets,bestZmm,isSideband2lbbZmm,"Sideband2lbbZmm_bestZmm_BSCoreOrdering")
        plots += self.plotFourObjects(cleanedHighestPtJets,bestHmm,isSideband2lbbHmm,"Sideband2lbbZmm_bestHmm_PtOrdering")
        plots += self.plotFourObjects(cleanedHighestBtagScoreJets,bestHmm,isSideband2lbbHmm,"Sideband2lbbZmm_bestHmm_BSCoreOrdering")
        
        if self.is_MC:
            plots += self.plotPairs(bestHmm, isSignalRegion,"SignalRegion_tightMuons" ,self.tightMuons)
            plots += self.plotPairs(bestHmm, isSignalRegion,"SignalRegion_cleanedMuons" ,self.cleanedMuons)
            plots += self.plotJets(cleanedHighestBtagScoreJets, isSignalRegion, "SignalRegion")
            plots += self.plotFourObjects(cleanedHighestPtJets,bestHmm,isSignalRegion,"SignalRegion")
            
        categories = []
        categories.append(category(nMuons=2, nElectrons=0))
        
        return plots

    def plotPairs(self, pair, sel, category, muons):
            plots = []
            if "SignalRegion" in category:
                plots.append(Plot.make1D(f"h_{category}_nMuons", op.rng_len(muons),sel, EquidistantBinning(15, 0., 15.), xTitle="Number of cleaned muons",plotopts={'no-data': True, 'show-ratio': False}))
                plots.append(Plot.make1D(f"h_{category}_mll", op.invariant_mass(pair[0].p4, pair[1].p4), sel, EquidistantBinning(50, 20., 200.), title="Invariant Mass of two muons", xTitle= "m_{\mu\mu} (GeV)",plotopts={'no-data': True, 'show-ratio': False}))
                plots.append(Plot.make1D(f"h_{category}_lep1_pT", pair[0].p4.Pt(), sel, EquidistantBinning(50, 20., 200.), title=" lepton1 pT", xTitle= "pT (GeV)",plotopts={'no-data': True, 'show-ratio': False}))
                plots.append(Plot.make1D(f"h_{category}_lep2_pT", pair[1].p4.Pt(), sel, EquidistantBinning(50, 10., 200.), title=" lepton2 pT", xTitle= "pT (GeV)",plotopts={'no-data': True, 'show-ratio': False}))
                plots.append(Plot.make1D(f"h_{category}_lep1_eta", pair[0].p4.Eta(),  sel, EquidistantBinning(5, -2.5, 2.5), title=" lepton1 Eta", xTitle= "\eta(l1)",plotopts={'no-data': True, 'show-ratio': False}))
                plots.append(Plot.make1D(f"h_{category}_lep2_eta", pair[1].p4.Eta(),  sel, EquidistantBinning(5, -2.5, 2.5), title=" lepton2 Eta", xTitle= "\eta(l2)",plotopts={'no-data': True, 'show-ratio': False}))
                plots.append(Plot.make1D(f"h_{category}_deltaR_l1_l2", op.deltaR(pair[0].p4,pair[1].p4),  sel, EquidistantBinning(6, 0, 6), title=" Delta R", xTitle= "\Delta R" ,plotopts={'no-data': True, 'show-ratio': False}))
                plots.append(Plot.make1D(f"h_{category}_deltaPhi_l1_l2", op.deltaPhi(pair[0].p4,pair[1].p4),  sel, EquidistantBinning(6, -3, 3), title=" Delta Phi", xTitle= "\Delta \phi",plotopts={'no-data': True, 'show-ratio': False} ))
                plots.append(Plot.make1D(f"h_{category}_deltaEta_l1_l2", pair[0].p4.Eta() - pair[1].p4.Eta(),  sel, EquidistantBinning(6, -3, 3), title=" Delta Eta", xTitle= "\Delta \eta" ,plotopts={'no-data': True, 'show-ratio': False}))
            else:
                plots.append(Plot.make1D(f"h_{category}_nMuons", op.rng_len(muons),sel, EquidistantBinning(15, 0., 15.), xTitle="Number of cleaned muons"))
                plots.append(Plot.make1D(f"h_{category}_mll", op.invariant_mass(pair[0].p4, pair[1].p4), sel, EquidistantBinning(50, 20., 200.), title="Invariant Mass of two muons", xTitle= "m_{\mu\mu} (GeV)"))
                plots.append(Plot.make1D(f"h_{category}_lep1_pT", pair[0].p4.Pt(), sel, EquidistantBinning(50, 20., 200.), title=" lepton1 pT", xTitle= "pT (GeV)"))
                plots.append(Plot.make1D(f"h_{category}_lep2_pT", pair[1].p4.Pt(), sel, EquidistantBinning(50, 10., 200.), title=" lepton2 pT", xTitle= "pT (GeV)"))
                plots.append(Plot.make1D(f"h_{category}_lep1_eta", pair[0].p4.Eta(),  sel, EquidistantBinning(5, -2.5, 2.5), title=" lepton1 Eta", xTitle= "\eta(l1)"))
                plots.append(Plot.make1D(f"h_{category}_lep2_eta", pair[1].p4.Eta(),  sel, EquidistantBinning(5, -2.5, 2.5), title=" lepton2 Eta", xTitle= "\eta(l2)"))
                plots.append(Plot.make1D(f"h_{category}_deltaR_l1_l2", op.deltaR(pair[0].p4,pair[1].p4),  sel, EquidistantBinning(6, 0, 6), title=" Delta R", xTitle= "\Delta R") )
                plots.append(Plot.make1D(f"h_{category}_deltaPhi_l1_l2", op.deltaPhi(pair[0].p4,pair[1].p4),  sel, EquidistantBinning(6, -3, 3), title=" Delta Phi", xTitle= "\Delta \phi") )
                plots.append(Plot.make1D(f"h_{category}_deltaEta_l1_l2", pair[0].p4.Eta() - pair[1].p4.Eta(),  sel, EquidistantBinning(6, -3, 3), title=" Delta Eta", xTitle= "\Delta \eta") )
            return plots

    def plotJets(self, jets, sel, category):
            plots = []
            if "SignalRegion" in category:
                plots.append(Plot.make1D(f"h_{category}_nJets", op.rng_len(jets),sel, EquidistantBinning(15, 0., 15.), xTitle="Number of b-jets",plotopts={'no-data': True, 'show-ratio': False}))
                plots.append(Plot.make1D(f"h_{category}_leading_pt", jets[0].pt,sel, EquidistantBinning(50, 0, 500), title="pT(j1)", xTitle="pT(j1) (GeV/c)",plotopts={'no-data': True, 'show-ratio': False}))
                plots.append(Plot.make1D(f"h_{category}_leading_eta", jets[0].eta,sel, EquidistantBinning(25, -2.5, 2.5), title="eta(j1)", xTitle="bjet \eta",plotopts={'no-data': True, 'show-ratio': False}))
                plots.append(Plot.make1D(f"h_{category}_leading_phi", jets[0].phi,sel, EquidistantBinning(25, -2.5, 2.5), title="phi(j1)", xTitle="bjet \phi",plotopts={'no-data': True, 'show-ratio': False}))
                plots.append(Plot.make1D(f"h_{category}_subleading_pt", jets[1].pt,sel, EquidistantBinning(50, 0, 500), title="pT(j2)", xTitle="pT(j2) (GeV/c)",plotopts={'no-data': True, 'show-ratio': False}))
                plots.append(Plot.make1D(f"h_{category}_subleading_eta", jets[1].eta,sel, EquidistantBinning(25, -2.5, 2.5), title="eta(j2)", xTitle="\eta(j2)",plotopts={'no-data': True, 'show-ratio': False}))
                plots.append(Plot.make1D(f"h_{category}_subleading_phi", jets[1].phi,sel, EquidistantBinning(25,-2.5, 2.5), title="phi(j2)", xTitle="\phi(j2)",plotopts={'no-data': True, 'show-ratio': False}))
                plots.append(Plot.make1D(f"h_{category}_deltaR_j1_j2", op.deltaR(jets[0].p4,jets[1].p4),sel,EquidistantBinning(10, 0, 5), title="DR(j1,j2)", xTitle="\Delta R(j1,j2)",plotopts={'no-data': True, 'show-ratio': False}))
                plots.append(Plot.make1D(f"h_{category}_deltaPhi_j1_j2", op.deltaPhi(jets[0].p4,jets[1].p4),sel,EquidistantBinning(7, -3.3, 3.3), title=" Delta Phi", xTitle= "\Delta \phi" ,plotopts={'no-data': True, 'show-ratio': False}))
                plots.append(Plot.make1D(f"h_{category}_deltaEta_j1_j2", jets[0].p4.Eta() - jets[1].p4.Eta(),  sel, EquidistantBinning(6, -3, 3), title=" Delta Eta", xTitle= "\Delta \eta" ,plotopts={'no-data': True, 'show-ratio': False}))
            else:
                plots.append(Plot.make1D(f"h_{category}_nJets", op.rng_len(jets),sel, EquidistantBinning(15, 0., 15.), xTitle="Number of b-jets"))
                plots.append(Plot.make1D(f"h_{category}_leading_pt", jets[0].pt,sel, EquidistantBinning(50, 0, 500), title="pT(j1)", xTitle="pT(j1) (GeV/c)"))
                plots.append(Plot.make1D(f"h_{category}_leading_eta", jets[0].eta,sel, EquidistantBinning(25, -2.5, 2.5), title="eta(j1)", xTitle="bjet \eta"))
                plots.append(Plot.make1D(f"h_{category}_leading_phi", jets[0].phi,sel, EquidistantBinning(25, -2.5, 2.5), title="phi(j1)", xTitle="bjet \phi"))
                plots.append(Plot.make1D(f"h_{category}_subleading_pt", jets[1].pt,sel, EquidistantBinning(50, 0, 500), title="pT(j2)", xTitle="pT(j2) (GeV/c)"))
                plots.append(Plot.make1D(f"h_{category}_subleading_eta", jets[1].eta,sel, EquidistantBinning(25, -2.5, 2.5), title="eta(j2)", xTitle="\eta(j2)"))
                plots.append(Plot.make1D(f"h_{category}_subleading_phi", jets[1].phi,sel, EquidistantBinning(25,-2.5, 2.5), title="phi(j2)", xTitle="\phi(j2)"))
                plots.append(Plot.make1D(f"h_{category}_deltaR_j1_j2", op.deltaR(jets[0].p4,jets[1].p4),sel,EquidistantBinning(10, 0, 5), title="DR(j1,j2)", xTitle="\Delta R(j1,j2)"))
                plots.append(Plot.make1D(f"h_{category}_deltaPhi_j1_j2", op.deltaPhi(jets[0].p4,jets[1].p4),sel,EquidistantBinning(7, -3.3, 3.3), title=" Delta Phi", xTitle= "\Delta \phi") )
                plots.append(Plot.make1D(f"h_{category}_deltaEta_j1_j2", jets[0].p4.Eta() - jets[1].p4.Eta(),  sel, EquidistantBinning(6, -3, 3), title=" Delta Eta", xTitle= "\Delta \eta") )
                if "BScoreOrdering" in category:
                    plots.append(Plot.make1D(f"h_{category}_bscore_j1", jets[0].btagDeepFlavB,sel,EquidistantBinning(10, 0, 1.0), title=" b tagging DeepJet score jet 1", xTitle= "btag score j1") )
                    plots.append(Plot.make1D(f"h_{category}_bscore_j2", jets[1].btagDeepFlavB,sel,EquidistantBinning(10, 0, 1.0), title=" b tagging DeepJet score jet 2", xTitle= "btag score j2") )
            return plots

    def plotFourObjects(self,jets,pair,sel,category):
            plots = []
            if "SignalRegion" in category:
                plots.append(Plot.make1D(f"h_{category}_DeltaR_l1_j1", op.deltaR(pair[0].p4, jets[0].p4),sel, EquidistantBinning(10, 0, 7), title="DR(l1,j1)", xTitle="\Delta R(l1, ak4bjet1)",plotopts={'no-data': True, 'show-ratio': False}))
                plots.append(Plot.make1D(f"h_{category}_DeltaR_l2_j2", op.deltaR(pair[1].p4, jets[1].p4),sel, EquidistantBinning(10, 0, 7), title="DR(l2,j2)", xTitle="\Delta R(l2, ak4bjet2)",plotopts={'no-data': True, 'show-ratio': False}))
                plots.append(Plot.make1D(f"h_{category}_DeltaPhi_l1_j1", op.deltaPhi(pair[0].p4, jets[0].p4),sel, EquidistantBinning(10, 0, 7), title="DeltaPhi(l1,j1)", xTitle="\Delta \phi (l1, ak4bjet1)",plotopts={'no-data': True, 'show-ratio': False}))
                plots.append(Plot.make1D(f"h_{category}_DeltaPhi_l2_j2", op.deltaPhi(pair[1].p4, jets[1].p4),sel, EquidistantBinning(10, 0, 7), title="DeltaPhi(l2,j2)", xTitle="\Delta \phi (l2, ak4bjet2)",plotopts={'no-data': True, 'show-ratio': False}))
                plots.append(Plot.make1D(f"h_{category}_deltaEta_l1_j1", pair[0].p4.Eta() - jets[0].p4.Eta(),  sel, EquidistantBinning(6, -3, 3), title=" Delta Eta", xTitle= "\Delta \eta" ,plotopts={'no-data': True, 'show-ratio': False}))
                plots.append(Plot.make1D(f"h_{category}_deltaEta_l2_j2", pair[1].p4.Eta() - jets[1].p4.Eta(),  sel, EquidistantBinning(6, -3, 3), title=" Delta Eta", xTitle= "\Delta \eta" ,plotopts={'no-data': True, 'show-ratio': False}))
                plots.append(Plot.make1D(f"h_{category}_mass_leptonsANDjets", op.invariant_mass(pair[0].p4, pair[1].p4, jets[0].p4,jets[1].p4),sel, EquidistantBinning(1000, 0, 500), title="m(j1,j2,l1,l2)", xTitle="m(j1,j2,l1,l2)",plotopts={'no-data': True, 'show-ratio': False}))
            else:
                plots.append(Plot.make1D(f"h_{category}_DeltaR_l1_j1", op.deltaR(pair[0].p4, jets[0].p4),sel, EquidistantBinning(10, 0, 7), title="DR(l1,j1)", xTitle="\Delta R(l1, ak4bjet1)"))
                plots.append(Plot.make1D(f"h_{category}_DeltaR_l2_j2", op.deltaR(pair[1].p4, jets[1].p4),sel, EquidistantBinning(10, 0, 7), title="DR(l2,j2)", xTitle="\Delta R(l2, ak4bjet2)"))
                plots.append(Plot.make1D(f"h_{category}_DeltaPhi_l1_j1", op.deltaPhi(pair[0].p4, jets[0].p4),sel, EquidistantBinning(10, 0, 7), title="DeltaPhi(l1,j1)", xTitle="\Delta \phi (l1, ak4bjet1)"))
                plots.append(Plot.make1D(f"h_{category}_DeltaPhi_l2_j2", op.deltaPhi(pair[1].p4, jets[1].p4),sel, EquidistantBinning(10, 0, 7), title="DeltaPhi(l2,j2)", xTitle="\Delta \phi (l2, ak4bjet2)"))
                plots.append(Plot.make1D(f"h_{category}_deltaEta_l1_j1", pair[0].p4.Eta() - jets[0].p4.Eta(),  sel, EquidistantBinning(6, -3, 3), title=" Delta Eta", xTitle= "\Delta \eta") )
                plots.append(Plot.make1D(f"h_{category}_deltaEta_l2_j2", pair[1].p4.Eta() - jets[1].p4.Eta(),  sel, EquidistantBinning(6, -3, 3), title=" Delta Eta", xTitle= "\Delta \eta") )
                plots.append(Plot.make1D(f"h_{category}_mass_leptonsANDjets", op.invariant_mass(pair[0].p4, pair[1].p4, jets[0].p4,jets[1].p4),sel, EquidistantBinning(1000, 0, 500), title="m(j1,j2,l1,l2)", xTitle="m(j1,j2,l1,l2)"))
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




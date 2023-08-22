from bamboo.analysismodules import NanoAODHistoModule
from bamboo.treedecorators import NanoAODDescription
from bamboo.scalefactors import binningVariables_nano,lumiPerPeriod_default

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

if os.environ["VIRTUAL_ENV"] == "":
    print("$VIRTUAL_ENV is not set. Please, activate your virtual environment")
    exit(-1)

print(os.environ["VIRTUAL_ENV"])
#print("{0}/../ClassicSVfit/interface".format(os.environ["VIRTUAL_ENV"]))

#loadDependency( headers="SVfitBambooProducer.h",
#                includePath="{0}/../ClassicSVfit/interface".format(os.environ["VIRTUAL_ENV"]),
#                dynamicPath="{0}/../build-ClassicSVfit-FastMTT/src".format(os.environ["VIRTUAL_ENV"]), 
#                libraries="ClassicSVfit")
                
#svFitBambooProducer = op.define("SVfitBambooProducer", 'auto <<name>> = SVfitBambooProducer();')


class category:
    def __init__(self,nMuons=0, nElectrons=0, nTaus=0):
        #FIXME: add checks on total number of leptons

        self.mu  = nMuons
        self.ele = nElectrons
        self.tau = nTaus  
        self._name = self.__str__()

    def __str__(self):
        strmu  = "" if not self.mu else f"{self.mu}mu"
        strele = "" if not self.ele else f"{self.ele}ele"
        strtau = "" if not self.tau else f"{self.tau}tau"  

        return strmu+strele+strtau

    def nMuons(self):
        return self.mu
    def nElectrons(self):
        return self.ele
    def nTaus(self):
        return self.tau
    def name(self):
        return self._name


class ZhToLLTauTau(NanoAODHistoModule):
    """ Module for Zh lltautau analysis"""
    def __init__(self, args):
        super(ZhToLLTauTau, self).__init__(args)

    def isMC(self, sampleName):
        return sampleName.split("_")[0] not in ("SingleMuon", "DoubleMuon", "SingleEGamma")
    def prepareTree(self, tree, sample=None, sampleCfg=None, backend=None):
        from bamboo.treedecorators import nanoRochesterCalc
        from bamboo.analysisutils import configureJets
        from bamboo.analysisutils import configureRochesterCorrection
        tree,noSel,be,lumiArgs = NanoAODHistoModule.prepareTree(self, tree, sample=sample, sampleCfg=sampleCfg, description=NanoAODDescription.get('v7', year='2018', isMC=True, systVariations=[nanoRochesterCalc]), backend=backend)
        era = sampleCfg["era"]
        if era == "2018":
            configureRochesterCorrection(tree._Muon, "./RoccoR2018UL.txt", isMC=self.isMC(sample), backend=be)

        return tree,noSel,be,lumiArgs
 
    def definePlots(self,t, noSel, sample=None, sampleCfg=None):

        plots = []

        if self.isMC(sample):
            noSel = noSel.refine("mcWeight", weight=t.genWeight)

       
        yields = CutFlowReport("yields", printInLog=True)
        era = sampleCfg["era"]

        if era == "2018":
            singleMuonTriggers= [ t.HLT.IsoMu24 ] 


        isoMuFilterMask = 0xA

        #Trigger
        hasTriggerFired = noSel.refine("passSingleMuonHLT", cut=op.OR(*(singleMuonTriggers))) 

        plots = []

        triggerObj = op.select(t.TrigObj, lambda trgObj: op.AND( trgObj.id == 13,
                                                                 (trgObj.filterBits & isoMuFilterMask) )) 
        muons = op.sort(op.select(t.Muon, lambda mu : op.AND(
            mu.pt > 5.,
            op.abs(mu.p4.Eta()) < 2.4 ,
            op.abs(mu.pfRelIso04_all ) < 0.40,
            op.abs(mu.dz ) < 1.,
            op.abs(mu.dxy) < 0.5,
            op.OR(mu.isGlobal,op.AND(mu.isTracker,mu.nStations>0)),
            op.sqrt(mu.dxy**2 + mu.dz**2)/op.sqrt(mu.dxyErr**2+mu.dzErr**2) < 4 ## SIP3D
            )), lambda mu : -mu.pt)
        electrons = op.sort(op.select(t.Electron, lambda el : op.AND(
            el.pt > 7.,
            op.abs(el.p4.Eta()) < 2.5,
            op.abs(el.dxy ) < 0.5,
            op.abs(el.dz ) < 1.,
            op.sqrt(el.dxy**2 + el.dz**2)/op.sqrt(el.dxyErr**2+el.dzErr**2) < 4, ## SIP3D
            op.abs(el.pfRelIso03_all) < 0.40
            )), lambda el : -el.pt)

        taus = op.select(t.Tau, lambda tau : op.AND(
            tau.p4.Pt() > 20., 
            op.abs(tau.p4.Eta()) < 2.4 ,
            op.abs(tau.dz) < 0.2
        ))

        #Select only taus not in the proximity of electrons/muons 
        cleanedTaus = op.select(taus, lambda it : op.AND(
                                                  op.NOT(op.rng_any(electrons, lambda ie : op.deltaR(it.p4, ie.p4) < 0.3 )),
                                                  op.NOT(op.rng_any(muons, lambda im : op.deltaR(it.p4, im.p4) < 0.3 ))
                                ))
 
         
        mZ = 91.1876
        mH = 125.1
        
        pairsMuMu = op.combine(muons, N=2, pred=lambda l1,l2 : op.AND(l1.charge != l2.charge, l1.pt > 20., l2.pt> 10.))

        #--> Select the best Z asking that the invariant mass is the closest to mZ
        bestZ = op.rng_min_element_by(pairsMuMu, lambda ll : op.abs(op.invariant_mass(ll[0].p4,ll[1].p4)-mZ))
        bestH = op.rng_min_element_by(pairsMuMu, lambda ll : op.abs(op.invariant_mass(ll[0].p4,ll[1].p4)-mH))
        #Ask to have a good Z with the refine option
        hasZmm      = hasTriggerFired.refine("hasZmm",  cut=[op.rng_len(pairsMuMu) > 0] )
        hasHmm      = hasTriggerFired.refine("hasHmm",  cut=[op.rng_len(pairsMuMu) > 0] )
        
        
        hasZmmTight = hasZmm.refine("hasZmmTight",  cut=[ bestZ[0].pfRelIso04_all < 0.15,
                                                          bestZ[0].tightId ])
       
     
        plots += self.plotPairs(bestZ, hasZmm, "ZosMuMu")
        plots += self.plotPairs(bestH, hasHmm, "HosMuMu")
        plots += self.plotPairs(bestZ, hasZmmTight, "osMuMu_Tight")
    

        categories = []
        categories.append(category(nMuons=1, nElectrons=0, nTaus=1))  #--> In the same way, add the other categories
       

        plots.append(yields)
        yields.add(noSel, title="NoSel")
        yields.add(hasTriggerFired, title="TriggerFired")
        yields.add(hasZmm, title="Zloose")
        #-->add other yields


        for cat in categories:
            #muons not used for the Z candidate
            otherLeptons = op.select(muons, partial(lambda l,oz=None : op.AND(l.idx != oz[0].idx, l.idx != oz[1].idx), oz=bestZ))
            
            #Muon category
            if cat.nMuons() > 0:
                oslep3lep4 = op.combine((otherLeptons, cleanedTaus), pred=lambda mu,tau : mu.charge != tau.charge)
            #--> Add the other categories
            
            bestH = op.rng_max_element_by(oslep3lep4, lambda ll : ll[0].p4.Pt()+ll[1].p4.Pt())
            hasSeconPair = hasZmmTight.refine(f"hasSeconPairCategory_{cat}", cut=[op.rng_len(oslep3lep4) > 0])
            yields.add(hasSeconPair, title=f"With a higgs pair cadidate {cat}")

        return plots

    def plotPairs(self, pair, sel, category):
        
            #-->Feel free to add whatever plot you are interested in
            plots = []
            plots.append(Plot.make1D(f"h_{category}_mll", op.invariant_mass(pair[0].p4, pair[1].p4), sel, EquidistantBinning(100, 20., 200.), title="Invariant Mass of two muons", xTitle= "m_{\mu\mu} (GeV)"))
            plots.append(Plot.make1D(f"h_{category}_lep1_pT", pair[0].p4.Pt(), sel, EquidistantBinning(100, 20., 200.), title=" lepton1 pT", xTitle= "pT (GeV)"))
            plots.append(Plot.make1D(f"h_{category}_lep2_pT", pair[1].p4.Pt(), sel, EquidistantBinning(100, 10., 200.), title=" lepton2 pT", xTitle= "pT (GeV)"))
            plots.append(Plot.make1D(f"h_{category}_lep1_eta", pair[0].p4.Eta(),  sel, EquidistantBinning(100, -3, 3), title=" lepton1 Eta", xTitle= "eta"))
            plots.append(Plot.make1D(f"h_{category}_lep2_eta", pair[1].p4.Eta(),  sel, EquidistantBinning(100, -3, 3), title=" lepton2 Eta", xTitle= "eta"))
            plots.append(Plot.make1D(f"h_{category}_deltaR", op.deltaR(pair[0].p4,pair[1].p4),  sel, EquidistantBinning(20, 0, 10), title=" Delta R", xTitle= "deltaR") )
            plots.append(Plot.make1D(f"h_{category}_deltaPhi", op.deltaPhi(pair[0].p4,pair[1].p4),  sel, EquidistantBinning(50, -3.3, 3.3), title=" Delta Phi", xTitle= "deltaPhi") )

            return plots


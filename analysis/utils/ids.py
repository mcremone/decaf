import numpy as np
from coffea.util import save
import awkward as ak
import os

######
## Electron
## Electron_cutBased Int_t cut-based ID Fall17 V2
## (0:fail, 1:veto, 2:loose, 3:medium, 4:tight)
## https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2
######


def isLooseElectron(e):
    
    mask = (
        (e.pt > 10)
        & (abs(e.eta+e.deltaEtaSC) < 1.4442)
        & (abs(e.dxy) < 0.05)
        & (abs(e.dz) < 0.1)
        & (e.cutBased >= 2)
    ) | (
        (e.pt > 10)
        & (abs(e.eta+e.deltaEtaSC) > 1.5660)
        & (abs(e.eta+e.deltaEtaSC) < 2.5)
        & (abs(e.dxy) < 0.1)
        & (abs(e.dz) < 0.2)
        & (e.cutBased >= 2)
    )
    
    return mask


def isTightElectron(e):
    
    mask = (
        (e.pt > 40)
        & (abs(e.eta+e.deltaEtaSC) < 1.4442)
        & (abs(e.dxy) < 0.05)
        & (abs(e.dz) < 0.1)
        & (e.cutBased == 4)
    ) | (
        (e.pt > 40)
        & (abs(e.eta+e.deltaEtaSC) > 1.5660)
        & (abs(e.eta+e.deltaEtaSC) < 2.5)
        & (abs(e.dxy) < 0.1)
        & (abs(e.dz) < 0.2)
        & (e.cutBased == 4)
    )
    
    return mask


#######
## Muon
## Muon ID WPs:
## https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2#Muon_selectors_Since_9_4_X
## Muon isolation WPs:
## https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonSelection#Muon_Isolation
#######


def isLooseMuon(mu):
    
    mask = (mu.pt > 15) & (abs(mu.eta) < 2.4) & mu.looseId & (mu.pfRelIso04_all < 0.25)
    
    return mask


def isTightMuon(mu):
    
    mask = (mu.pt > 30) & (abs(mu.eta) < 2.4) & mu.tightId & (mu.pfRelIso04_all < 0.15)
    
    return mask


def isSoftMuon(mu):
    
    mask = (mu.pt > 5) & (abs(mu.eta) < 2.4) & mu.tightId & (mu.pfRelIso04_all > 0.15)
    
    return mask


######
## Tau
## https://twiki.cern.ch/twiki/bin/viewauth/CMS/TauIDRecommendationForRun2
## The decayModeFindingNewDMs: recommended for use with DeepTauv2p1, where decay
## modes 5 and 6 should be explicitly rejected.
## This should already be applied in NanoAOD.
##
## Tau_idDeepTau2017v2p1VSe ID working points (bitmask):
## 1 = VVVLoose, 2 = VVLoose, 4 = VLoose, 8 = Loose,
## 16 = Medium, 32 = Tight, 64 = VTight, 128 = VVTight
##
## Tau_idDeepTau2017v2p1VSjet ID working points (bitmask):
## 1 = VVVLoose, 2 = VVLoose, 4 = VLoose, 8 = Loose,
## 16 = Medium, 32 = Tight, 64 = VTight, 128 = VVTight
##
## Tau_idDeepTau2017v2p1VSmu ID working points (bitmask):
## 1 = VLoose, 2 = Loose, 4 = Medium, 8 = Tight
######


def isLooseTau(tau):
    
    try:
        decayModeDMs=tau.decayModeFindingNewDMs
    except:
        decayModeDMs=~np.isnan(ak.ones_like(tau.pt))

    mask = (
        (tau.pt > 20)
        & (abs(tau.eta) < 2.3)
        #& ~(tau.decayMode == 5)
        #& ~(tau.decayMode == 6)
        & decayModeDMs
        #& ((tau.idDeepTau2017v2p1VSe & 16) == 16)
        & ((tau.idDeepTau2017v2p1VSjet & 4) == 4)
        #& ((tau.idDeepTau2017v2p1VSmu & 2) == 2)
    )
    
    return mask


######
## Photon
## https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2
## Photon_cutBased Int_t cut-based ID bitmap, Fall17V2,
## (0:fail, 1:loose, 2:medium, 3:tight)
## Note: Photon IDs are integers, not bit masks
######


def isLoosePhoton(pho):
    
    mask = (
        (pho.pt > 20)
        & ~((abs(pho.eta) > 1.4442) & (abs(pho.eta) < 1.5660))
        & (abs(pho.eta) < 2.5)
        & (pho.cutBased >= 1)
    )
    
    return mask&(pho.electronVeto)


def isTightPhoton(pho):
    
    mask = (pho.pt > 230) & (pho.cutBased == 3)
    
    return mask&(pho.isScEtaEB)&(pho.electronVeto) #tight photons are barrel only


######
## Fatjet
## https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVUL
## Tight working point including lepton veto (TightLepVeto)
######


def isGoodAK15(fj):
    
    mask = (fj.pt > 160) & (abs(fj.eta) < 2.4) & ((fj.jetId & 2) == 2)
    
    return mask


######
## Jet
## https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVUL
## Tight working point including lepton veto (TightLepVeto)
##
## For Jet ID flags, bit1 is Loose (always false in 2017 since it does not
## exist), bit2 is Tight, bit3 is TightLepVeto. The POG recommendation is to
## use Tight Jet ID as the standard Jet ID.
######
## PileupJetID
## https://twiki.cern.ch/twiki/bin/view/CMS/PileupJetIDUL
## Using Loose Pileup ID
##
## Note: There is a bug in 2016 UL in which bit values for Loose and Tight Jet
## Pileup IDs are accidentally flipped relative to 2017 UL and 2018 UL.
##
## For 2016 UL,
## Jet_puId = (passtightID*4 + passmediumID*2 + passlooseID*1).
##
## For 2017 UL and 2018 UL,
## Jet_puId = (passlooseID*4 + passmediumID*2 + passtightID*1).
######


def isGoodAK4(j, year):

    puId_value = 4
    if '2016' in year:
        puId_value =1
    
    mask = (
        (j.pt > 30) 
        & (abs(j.eta) < 2.4) 
        & ((j.jetId & 2) == 2)
        & ((j.pt >= 50) | ((j.puId & puId_value) == puId_value))
        )

    return mask

def isSoftAK4(j, year):

    def puId_cut_low_pt(jet_pt):
        puId = (0.85-0.7)*(jet_pt-30)/(30-8) + 0.85
        return puId

    mask = (
        (j.pt > 15) 
        & (abs(j.eta) < 2.4) 
        & ((j.jetId & 2) == 2)
        & (j.puIdDisc > puId_cut_low_pt(j.pt))
        )

    return (mask&~(j.pt > 30))|isGoodAK4(j, year)


######
## HEM
######


def isHEMJet(j):

    mask = (j.pt > 30) & (j.eta > -3.0) & (j.eta < -1.3) & (j.phi > -1.57) & (j.phi < -0.87)
    return mask


ids = {}
ids["isLooseElectron"] = isLooseElectron
ids["isTightElectron"] = isTightElectron
ids["isLooseMuon"] = isLooseMuon
ids["isTightMuon"] = isTightMuon
ids["isSoftMuon"] = isSoftMuon
ids["isLooseTau"] = isLooseTau
ids["isLoosePhoton"] = isLoosePhoton
ids["isTightPhoton"] = isTightPhoton
ids["isGoodAK4"] = isGoodAK4
ids["isSoftAK4"] = isSoftAK4
ids["isGoodAK15"] = isGoodAK15
ids["isHEMJet"] = isHEMJet
path = "decaf/analysis/data/" if "srv" in os.getcwd() else "data/"   ### to make it run with coffea4bees
save(ids, f"{path}/ids.coffea")

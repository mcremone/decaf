#!/usr/bin/env python
import logging
import numpy as np
import awkward as ak
import json
import copy
from collections import defaultdict
from coffea import processor
import hist
from coffea.analysis_tools import Weights, PackedSelection
from coffea.lumi_tools import LumiMask
from coffea.util import load, save
from optparse import OptionParser
from uproot_methods import TVector2Array

class AnalysisProcessor(processor.ProcessorABC):

    lumis = { 
        #Values from https://twiki.cern.ch/twiki/bin/view/CMS/LumiRecommendationsRun2                                                      
        '2016postVFP': 36.31,
        '2016preVFP': 36.31,
        '2017': 41.48,
        '2018': 59.83
    }

    lumiMasks = {
        '2016postVFP': LumiMask("data/jsons/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt"),
        '2016preVFP': LumiMask("data/jsons/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt"),
        '2017': LumiMask("data/jsons/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt"),
        '2018"': LumiMask("data/jsons/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt"),
    }
    
    met_filters = {
        # https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFiltersRun2
        '2016postVFP': [
                'goodVertices',
                'globalSuperTightHalo2016Filter',
                'HBHENoiseFilter',
                'HBHENoiseIsoFilter',
                'EcalDeadCellTriggerPrimitiveFilter',
                'BadPFMuonFilter',
                'BadPFMuonDzFilter',
                'eeBadScFilter'
                ],

        '2016preVFP': [
                'goodVertices',
                'globalSuperTightHalo2016Filter',
                'HBHENoiseFilter',
                'HBHENoiseIsoFilter',
                'EcalDeadCellTriggerPrimitiveFilter',
                'BadPFMuonFilter',
                'BadPFMuonDzFilter',
                'eeBadScFilter'
                ],
        
        '2017': [
                'goodVertices', 
                'globalSuperTightHalo2016Filter', 
                'HBHENoiseFilter', 
                'HBHENoiseIsoFilter', 
                'EcalDeadCellTriggerPrimitiveFilter', 
                'BadPFMuonFilter', 
                'BadPFMuonDzFilter', 
                'eeBadScFilter', 
                'ecalBadCalibFilter'
                ],

        '2018': [
                'goodVertices', 
                'globalSuperTightHalo2016Filter', 
                'HBHENoiseFilter', 
                'HBHENoiseIsoFilter', 
                'EcalDeadCellTriggerPrimitiveFilter', 
                'BadPFMuonFilter', 
                'BadPFMuonDzFilter', 
                'eeBadScFilter', 
                'ecalBadCalibFilter'
                ]
    }
            
    def __init__(self, year, xsec, corrections, ids, common):

        self._year = year
        self._lumi = 1000.*float(AnalysisProcessor.lumis[year])
        self._xsec = xsec
        self._systematics = True
        self._skipJER = False

        self._samples = {
            'sr':('ZJets','WJets','DY','TT','ST','WW','WZ','ZZ','QCD','HToBB','HTobb','MET','mhs'),
            'wmcr':('WJets','DY','TT','ST','WW','WZ','ZZ','QCD','HToBB','HTobb','MET'),
            'tmcr':('WJets','DY','TT','ST','WW','WZ','ZZ','QCD','HToBB','HTobb','MET'),
            'wecr':('WJets','DY','TT','ST','WW','WZ','ZZ','QCD','HToBB','HTobb','SingleElectron','EGamma'),
            'tecr':('WJets','DY','TT','ST','WW','WZ','ZZ','QCD','HToBB','HTobb','SingleElectron','EGamma'),
            'qcdcr':('ZJets','WJets','TT','ST','WW','WZ','ZZ','QCD','MET'),
        }
        
        self._ZHbbvsQCDwp = {
            '2016': 0.53,
            '2017': 0.61,
            '2018': 0.65
        }

        self._met_triggers = {
            '2016postVFP': [
                'PFMETNoMu90_PFMHTNoMu90_IDTight',
                'PFMETNoMu100_PFMHTNoMu100_IDTight',
                'PFMETNoMu110_PFMHTNoMu110_IDTight',
                'PFMETNoMu120_PFMHTNoMu120_IDTight'
            ],
            '2016preVFP': [
                'PFMETNoMu90_PFMHTNoMu90_IDTight',
                'PFMETNoMu100_PFMHTNoMu100_IDTight',
                'PFMETNoMu110_PFMHTNoMu110_IDTight',
                'PFMETNoMu120_PFMHTNoMu120_IDTight'
            ],
            '2017': [
                'PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60',
                'PFMETNoMu120_PFMHTNoMu120_IDTight'
            ],
            '2018': [
                'PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60',
                'PFMETNoMu120_PFMHTNoMu120_IDTight'
            ]
        }

        self._singleelectron_triggers = { #2017 and 2018 from monojet, applying dedicated trigger weights
            '2016postVFP': [
                'Ele27_WPTight_Gsf',
                'Ele105_CaloIdVT_GsfTrkIdT'
            ],
            '2016preVFP': [
                'Ele27_WPTight_Gsf',
                'Ele105_CaloIdVT_GsfTrkIdT'
            ],
            '2017': [
                'Ele35_WPTight_Gsf',
                'Ele115_CaloIdVT_GsfTrkIdT',
                'Photon200'
            ],
            '2018': [
                'Ele32_WPTight_Gsf',
                'Ele115_CaloIdVT_GsfTrkIdT',
                'Photon200'
            ]
        }

        self._corrections = corrections
        self._ids = ids
        self._common = common

        ptbins=[30.0, 
                60.0, 
                90.0, 
                120.0, 
                150.0, 
                180.0, 
                210.0, 
                250.0, 
                280.0, 
                310.0, 
                340.0, 
                370.0, 
                400.0, 
                430.0, 
                470.0, 
                510.0, 
                550.0, 
                590.0, 
                640.0, 
                690.0, 
                740.0, 
                790.0, 
                840.0, 
                900.0, 
                960.0, 
                1020.0, 
                1090.0, 
                1160.0, 
                1250.0]

        self.make_output = lambda: {
            'sumw': 0.,
            'template': hist.Hist(
                hist.axis.StrCategory([], name='region', growth=True),
                hist.axis.StrCategory([], name='systematic', growth=True),
                hist.axis.Variable([250,310,370,470,590,840,1020,1250,3000], name='recoil', label=r'$U$ [GeV]'),
                hist.axis.Variable([40,50,60,70,80,90,100,110,120,130,150,160,180,200,220,240,300], name='fjmass', label=r'AK15 Jet $m_{sd}$'),
                hist.axis.Variable([0, self._ZHbbvsQCDwp[self._year], 1], name='ZHbbvsQCD', label='ZHbbvsQCD', flow=False),
                hist.storage.Weight(),
            ),
            'ZHbbvsQCD': hist.Hist(
                hist.axis.StrCategory([], name='region', growth=True),
                hist.axis.Regular(15,0,1, name='ZHbbvsQCD', label='ZHbbvsQCD'),
                hist.storage.Weight(),
            ),
            'mindphirecoil': hist.Hist(
                hist.axis.StrCategory([], name='region', growth=True),
                hist.axis.Regular(30,0,3.5, name='mindphirecoil', label='Min dPhi(Recoil,AK4s)'),
                hist.axis.Variable([0, self._ZHbbvsQCDwp[self._year], 1], name='ZHbbvsQCD', label='ZHbbvsQCD', flow=False),
                hist.storage.Weight(),
            ),
            'minDphirecoil': hist.Hist(
                hist.axis.StrCategory([], name='region', growth=True),
                hist.axis.Regular(30,0,3.5, name='minDphirecoil', label='Min dPhi(Recoil,AK15s)'),
                hist.axis.Variable([0, self._ZHbbvsQCDwp[self._year], 1], name='ZHbbvsQCD', label='ZHbbvsQCD', flow=False),
                hist.storage.Weight(),
            ),
            'CaloMinusPfOverRecoil': hist.Hist(
                hist.axis.StrCategory([], name='region', growth=True),
                hist.axis.Regular(35,0,1, name='CaloMinusPfOverRecoil', label='Calo - Pf / Recoil'),
                hist.axis.Variable([0, self._ZHbbvsQCDwp[self._year], 1], name='ZHbbvsQCD', label='ZHbbvsQCD', flow=False),
                hist.storage.Weight(),
            ),
            'met': hist.Hist(
                hist.axis.StrCategory([], name='region', growth=True),
                hist.axis.Regular(30,0,600, name='met', label='MET'),
                hist.axis.Variable([0, self._ZHbbvsQCDwp[self._year], 1], name='ZHbbvsQCD', label='ZHbbvsQCD', flow=False),
                hist.storage.Weight(),
            ),
            'metphi': hist.Hist(
                hist.axis.StrCategory([], name='region', growth=True),
                hist.axis.Regular(35,-3.5,3.5, name='metphi', label='MET phi'),
                hist.axis.Variable([0, self._ZHbbvsQCDwp[self._year], 1], name='ZHbbvsQCD', label='ZHbbvsQCD', flow=False),
                hist.storage.Weight(),
            ),
            'mindphimet': hist.Hist(
                hist.axis.StrCategory([], name='region', growth=True),
                hist.axis.Regular(30,0,3.5, name='mindphimet', label='Min dPhi(MET,AK4s)'),
                hist.axis.Variable([0, self._ZHbbvsQCDwp[self._year], 1], name='ZHbbvsQCD', label='ZHbbvsQCD', flow=False),
                hist.storage.Weight(),
            ),
            'minDphimet': hist.Hist(
                hist.axis.StrCategory([], name='region', growth=True),
                hist.axis.Regular(30,0,3.5, name='minDphimet', label='Min dPhi(MET,AK15s)'),
                hist.axis.Variable([0, self._ZHbbvsQCDwp[self._year], 1], name='ZHbbvsQCD', label='ZHbbvsQCD', flow=False),
                hist.storage.Weight(),
            ),
            'j1pt': hist.Hist(
                hist.axis.StrCategory([], name='region', growth=True),
                hist.axis.Variable(ptbins, name='j1pt', label='AK4 Leading Jet Pt'),
                hist.axis.Variable([0, self._ZHbbvsQCDwp[self._year], 1], name='ZHbbvsQCD', label='ZHbbvsQCD', flow=False),
                hist.storage.Weight(),
            ),
            'j1eta': hist.Hist(
                hist.axis.StrCategory([], name='region', growth=True),
                hist.axis.Regular(35,-3.5,3.5, name='j1eta', label='AK4 Leading Jet Eta'),
                hist.axis.Variable([0, self._ZHbbvsQCDwp[self._year], 1], name='ZHbbvsQCD', label='ZHbbvsQCD', flow=False),
                hist.storage.Weight(),
            ),
            'j1phi': hist.Hist(
                hist.axis.StrCategory([], name='region', growth=True),
                hist.axis.Regular(35,-3.5,3.5, name='j1phi', label='AK4 Leading Jet Phi'),
                hist.axis.Variable([0, self._ZHbbvsQCDwp[self._year], 1], name='ZHbbvsQCD', label='ZHbbvsQCD', flow=False),
                hist.storage.Weight(),
            ),
            'fj1pt': hist.Hist(
                hist.axis.StrCategory([], name='region', growth=True),
                hist.axis.Variable(ptbins, name='fj1pt', label='AK15 Leading SoftDrop Jet Pt'),
                hist.axis.Variable([0, self._ZHbbvsQCDwp[self._year], 1], name='ZHbbvsQCD', label='ZHbbvsQCD', flow=False),
                hist.storage.Weight(),
            ),
            'fj1eta': hist.Hist(
                hist.axis.StrCategory([], name='region', growth=True),
                hist.axis.Regular(35,-3.5,3.5, name='fj1eta', label='AK15 Leading SoftDrop Jet Eta'),
                hist.axis.Variable([0, self._ZHbbvsQCDwp[self._year], 1], name='ZHbbvsQCD', label='ZHbbvsQCD', flow=False),
                hist.storage.Weight(),
            ),
            'fj1phi': hist.Hist(
                hist.axis.StrCategory([], name='region', growth=True),
                hist.axis.Regular(35,-3.5,3.5, name='fj1phi', label='AK15 Leading SoftDrop Jet Phi'),
                hist.axis.Variable([0, self._ZHbbvsQCDwp[self._year], 1], name='ZHbbvsQCD', label='ZHbbvsQCD', flow=False),
                hist.storage.Weight(),
            ),
            'njets': hist.Hist(
                hist.axis.StrCategory([], name='region', growth=True),
                hist.axis.IntCategory([0, 1, 2, 3, 4, 5, 6], name='njets', label='AK4 Number of Jets'),
                hist.axis.Variable([0, self._ZHbbvsQCDwp[self._year], 1], name='ZHbbvsQCD', label='ZHbbvsQCD', flow=False),
                hist.storage.Weight(),
            ),
            'ndflvL': hist.Hist(
                hist.axis.StrCategory([], name='region', growth=True),
                hist.axis.IntCategory([0, 1, 2, 3, 4, 5, 6], name='ndflvL', label='AK4 Number of deepFlavor Loose Jets'),
                hist.axis.Variable([0, self._ZHbbvsQCDwp[self._year], 1], name='ZHbbvsQCD', label='ZHbbvsQCD', flow=False),
                hist.storage.Weight(),
            ),
            'nfjclean': hist.Hist(
                hist.axis.StrCategory([], name='region', growth=True),
                hist.axis.IntCategory([0, 1, 2, 3, 4], name='nfjclean', label='AK15 Number of Cleaned Jets'),
                hist.axis.Variable([0, self._ZHbbvsQCDwp[self._year], 1], name='ZHbbvsQCD', label='ZHbbvsQCD', flow=False),
                hist.storage.Weight(),
            ),
            'mT': hist.Hist(
                hist.axis.StrCategory([], name='region', growth=True),
                hist.axis.Regular(20,0,600, name='mT', label='Transverse Mass'),
                hist.axis.Variable([0, self._ZHbbvsQCDwp[self._year], 1], name='ZHbbvsQCD', label='ZHbbvsQCD', flow=False),
                hist.storage.Weight(),
            ),
            'l1pt': hist.Hist(
                hist.axis.StrCategory([], name='region', growth=True),
                hist.axis.Variable(ptbins, name='l1pt', label='Leading Lepton Pt'),
                hist.axis.Variable([0, self._ZHbbvsQCDwp[self._year], 1], name='ZHbbvsQCD', label='ZHbbvsQCD', flow=False),
                hist.storage.Weight(),
            ),
            'l1eta': hist.Hist(
                hist.axis.StrCategory([], name='region', growth=True),
                hist.axis.Regular(48,-2.4,2.4, name='l1eta', label='Leading Lepton Eta'),
                hist.axis.Variable([0, self._ZHbbvsQCDwp[self._year], 1], name='ZHbbvsQCD', label='ZHbbvsQCD', flow=False),
                hist.storage.Weight(),
            ),
            'l1phi': hist.Hist(
                hist.axis.StrCategory([], name='region', growth=True),
                hist.axis.Regular(64,-3.2,3.2, name='l1phi', label='Leading Lepton Phi'),
                hist.axis.Variable([0, self._ZHbbvsQCDwp[self._year], 1], name='ZHbbvsQCD', label='ZHbbvsQCD', flow=False),
                hist.storage.Weight(),
            ),
    }

    def process(self, events):
        isData = not hasattr(events, "genWeight")
        if isData:
            # Nominal JEC are already applied in data
            return self.process_shift(events, None)

        import cachetools
        jec_cache = cachetools.Cache(np.inf)
    
        nojer = "NOJER" if self._skipJER else ""
        thekey = f"{self._year}mc{nojer}"

        def add_jec_variables(jets, event_rho):
            jets["pt_raw"] = (1 - jets.rawFactor)*jets.pt
            jets["mass_raw"] = (1 - jets.rawFactor)*jets.mass
            jets["pt_gen"] = ak.values_astype(ak.fill_none(jets.matched_gen.pt, 0), np.float32)
            jets["event_rho"] = ak.broadcast_arrays(event_rho, jets.pt)[0]
            return jets
        
        jets = jet_factory[thekey].build(add_jec_variables(events.Jet, events.fixedGridRhoFastjetAll), jec_cache)
        subjets = jet_factory[thekey].build(add_jec_variables(events.AK15PFPuppiSubjet, events.fixedGridRhoFastjetAll), jec_cache)
        met = met_factory.build(events.MET, jets, {})

        shifts = [({"Jet": jets, "AK15PFPuppiSubjet": subjets, "MET": met}, None)]
        if self._systematics:
            shifts.extend([
                ({"Jet": jets.JES_jes.up, "AK15PFPuppiSubjet": subjets.JES_jes.up, "MET": met.JES_jes.up}, "JESUp"),
                ({"Jet": jets.JES_jes.down, "AK15PFPuppiSubjet": subjets.JES_jes.down, "MET": met.JES_jes.down}, "JESDown"),
                ({"Jet": jets, "AK15PFPuppiSubjet": subjets, "MET": met.MET_UnclusteredEnergy.up}, "UESUp"),
                ({"Jet": jets, "AK15PFPuppiSubjet": subjets, "MET": met.MET_UnclusteredEnergy.down}, "UESDown"),
            ])
            if not self._skipJER:
                shifts.extend([
                    ({"Jet": jets.JER.up, "AK15PFPuppiSubjet": subjets.JER.up, "MET": met.JER.up}, "JERUp"),
                    ({"Jet": jets.JER.down, "AK15PFPuppiSubjet": subjets.JER.down, "MET": met.JER.down}, "JERDown"),
                ])
        return processor.accumulate(self.process_shift(update(events, collections), name) for collections, name in shifts)

    def process_shift(self, events, shift_name):

        dataset = events.metadata['dataset']

        selected_regions = []
        for region, samples in self._samples.items():
            for sample in samples:
                if sample not in dataset: continue
                selected_regions.append(region)

        isData = not hasattr(events, "genWeight")
        selection = PackedSelection()
        weights = Weights(len(events), storeIndividual=True)
        output = self.make_output()
        if shift_name is None and not isData:
            output['sumw'] = ak.sum(events.genWeight)

        ###
        #Getting corrections, ids from .coffea files
        ###

        get_met_trig_weight      = self._corrections['get_met_trig_weight'][self._year]
        get_ele_loose_id_sf      = self._corrections['get_ele_loose_id_sf']
        get_ele_tight_id_sf      = self._corrections['get_ele_tight_id_sf']
        get_ele_trig_weight      = self._corrections['get_ele_trig_weight'][self._year]
        get_ele_reco_sf_below20  = self._corrections['get_ele_reco_sf_below20'][self._year]
        get_ele_reco_err_below20 = self._corrections['get_ele_reco_err_below20'][self._year]
        get_ele_reco_sf_above20  = self._corrections['get_ele_reco_sf_above20'][self._year]
        get_ele_reco_err_above20 = self._corrections['get_ele_reco_err_above20'][self._year]
        get_muon_loose_id_sf     = self._corrections['get_muon_loose_id_sf']
        get_muon_tight_id_sf     = self._corrections['get_muon_tight_id_sf']
        get_muon_loose_iso_sf    = self._corrections['get_muon_loose_iso_sf']
        get_muon_tight_iso_sf    = self._corrections['get_muon_tight_iso_sf']
        get_mu_rochester_sf      = self._corrections['get_mu_rochester_sf'][self._year]
        get_met_xy_correction    = self._corrections['get_met_xy_correction']
        get_pu_weight            = self._corrections['get_pu_weight']    
        get_nlo_ewk_weight       = self._corrections['get_nlo_ewk_weight']    
        get_nnlo_nlo_weight      = self._corrections['get_nnlo_nlo_weight'][self._year]
        get_msd_corr             = self._corrections['get_msd_corr']
        get_deepflav_weight      = self._corrections['get_btag_weight']['deepflav'][self._year]
        jet_factory              = self._corrections['jet_factory']
        fatjet_factory           = self._corrections['fatjet_factory']
        met_factory              = self._corrections['met_factory']
        
        isLooseElectron = self._ids['isLooseElectron'] 
        isTightElectron = self._ids['isTightElectron'] 
        isLooseMuon     = self._ids['isLooseMuon']     
        isTightMuon     = self._ids['isTightMuon']     
        isLooseTau      = self._ids['isLooseTau']      
        isLoosePhoton   = self._ids['isLoosePhoton']   
        isTightPhoton   = self._ids['isTightPhoton']   
        isGoodJet       = self._ids['isGoodJet']       
        isGoodFatJet    = self._ids['isGoodFatJet']    
        isHEMJet        = self._ids['isHEMJet']        
        
        match = self._common['match']
        deepflavWPs = self._common['btagWPs']['deepflav'][self._year]
        deepcsvWPs = self._common['btagWPs']['deepcsv'][self._year]

        ###
        #Initialize global quantities (MET ecc.)
        ###

        npv = events.PV.npvsGood
        run = events.run
        calomet = events.CaloMET
        met = events.MET
        corrected_pt, corrected_phi = get_met_xy_correction(self._year, npv, run, met.pt, met.phi, isData)
        met['pt'] = corrected_pt
        met['phi'] = corrected_phi
        

        ###
        #Initialize physics objects
        ###

        mu = events.Muon
        if isData:
            _k = get_mu_rochester_sf.kScaleDT(mu.charge, mu.pt, mu.eta, mu.phi)
        else:
            hasgen = (mu.matched_gen.pt.fillna(-1) > 0)
            u = awkward.JaggedArray.fromoffsets(mu.pt.offsets, np.random.rand(*_pt.flatten().shape))
            kspread = get_mu_rochester_sf.kSpreadMC(mu.charge[hasgen], mu.pt[hasgen], mu.eta[hasgen], mu.phi[hasgen],
                                           mu.matched_gen.pt[hasgen])
            ksmear = get_mu_rochester_sf.kSmearMC(mu.charge[~hasgen], mu.pt[~hasgen], mu.eta[~hasgen], mu.phi[~hasgen],
                                         mu.nTrackerLayers[~hasgen], u[~hasgen])
            k = mu.pt.ones_like()
            k[hasgen] = kspread
            k[~hasgen] = ksmear
        mask = _pt < 200
        rochester_pt = mu.pt.ones_like()
        rochester_pt[~mask] = mu.pt[~mask]
        rochester_pt[mask] = (k * mu.pt)[mask]
        mu['pt'] = rochester_pt
        mu['isloose'] = isLooseMuon(mu.pt,mu.eta,mu.pfRelIso04_all,mu.looseId,self._year)
        mu['istight'] = isTightMuon(mu.pt,mu.eta,mu.pfRelIso04_all,mu.tightId,self._year)
        mu['T'] = TVector2Array.from_polar(mu.pt, mu.phi)
        mu_loose=mu[mu.isloose]
        mu_tight=mu[mu.istight]
        mu_ntot = mu.counts
        mu_nloose = mu_loose.counts
        mu_ntight = mu_tight.counts
        leading_mu = mu[mu.pt.argmax()]
        leading_mu = leading_mu[leading_mu.istight]

        e = events.Electron
        e['isclean'] = ~match(e,mu_loose,0.3) 
        e['isloose'] = isLooseElectron(e.pt,e.eta+e.deltaEtaSC,e.dxy,e.dz,e.cutBased,self._year)
        e['istight'] = isTightElectron(e.pt,e.eta+e.deltaEtaSC,e.dxy,e.dz,e.cutBased,self._year)
        e['T'] = TVector2Array.from_polar(e.pt, e.phi)
        e_clean = e[e.isclean]
        e_loose = e_clean[e_clean.isloose]
        e_tight = e_clean[e_clean.istight]
        e_ntot = e.counts
        e_nloose = e_loose.counts
        e_ntight = e_tight.counts
        leading_e = e[e.pt.argmax()]
        leading_e = leading_e[leading_e.isclean]
        leading_e = leading_e[leading_e.istight]

        tau = events.Tau
        tau['isclean']=~match(tau,mu_loose,0.4)&~match(tau,e_loose,0.4)
        tau['isloose']=isLooseTau(tau.pt,tau.eta,tau.idDecayMode,tau.idDecayModeNewDMs,tau.idDeepTau2017v2p1VSe,tau.idDeepTau2017v2p1VSjet,tau.idDeepTau2017v2p1VSmu,self._year)
        tau_clean=tau[tau.isclean]
        tau_loose=tau_clean[tau_clean.isloose]
        tau_ntot=tau.counts
        tau_nloose=tau_loose.counts

        pho = events.Photon
        pho['isclean']=~match(pho,mu_loose,0.5)&~match(pho,e_loose,0.5)&~match(pho,tau_loose,0.5)
        _id = 'cutBasedBitmap'
        if self._year=='2016': 
            _id = 'cutBased'
        pho['isloose']=isLoosePhoton(pho.pt,pho.eta,pho[_id],self._year)&(pho.electronVeto) #added electron veto flag
        pho['T'] = TVector2Array.from_polar(pho.pt, pho.phi)
        pho_clean=pho[pho.isclean]
        pho_loose=pho_clean[pho_clean.isloose]
        pho_ntot=pho.counts
        pho_nloose=pho_loose.counts

        fj = events.AK15PFPuppiJet
        fj['sd'] = fj.subjets.sum()
        fj['isclean'] =~match(fj.sd,pho_loose,1.5)&~match(fj.sd,mu_loose,1.5)&~match(fj.sd,e_loose,1.5)&~match(fj.sd,tau_loose,1.5)
        fj['isgood'] = isGoodFatJet(fj.sd.pt, fj.sd.eta, fj.jetId)
        fj['msd_corr'] = get_msd_corr(fj)
        probQCD=fj.probQCDbb+fj.probQCDcc+fj.probQCDb+fj.probQCDc+fj.probQCDothers
        probZHbb=fj.probZbb+fj.probHbb
        fj['ZHbbvsQCD'] = probZHbb/(probZHbb+probQCD)
        fj_good = fj[fj.isgood]
        fj_clean = fj_good[fj_good.isclean]
        fj_ntot = fj.counts
        fj_ngood = fj_good.counts
        fj_nclean = fj_clean.counts

        j = events.Jet
        j['isgood'] = isGoodJet(j.pt, j.eta, j.jetId, j.puId, j.neHEF, j.chHEF)
        j['isHEM'] = isHEMJet(j.pt, j.eta, j.phi)
        j['isclean'] = ~match(j,e_loose,0.4)&~match(j,mu_loose,0.4)&~match(j,pho_loose,0.4)&~match(j,tau_loose,0.4)
        j['isiso'] = ~match(j,fj_clean[fj_clean.pt.argmax()],1.5)
        j['isdcsvL'] = (j.btagDeepB>deepcsvWPs['loose'])
        j['isdflvL'] = (j.btagDeepFlavB>deepflavWPs['loose'])
        j_good = j[j.isgood]
        j_clean = j_good[j_good.isclean]
        j_iso = j_clean[j_clean.isiso]
        j_dcsvL = j_iso[j_iso.isdcsvL]
        j_dflvL = j_iso[j_iso.isdflvL]
        j_HEM = j[j.isHEM]
        j_ntot=j.counts
        j_ngood=j_good.counts
        j_nclean=j_clean.counts
        j_niso=j_iso.counts
        j_ndcsvL=j_dcsvL.counts
        j_ndflvL=j_dflvL.counts
        j_nHEM = j_HEM.counts
        leading_j = j[j.pt.argmax()]
        leading_j = leading_j[leading_j.isgood]
        leading_j = leading_j[leading_j.isclean]

        ###
        # Calculate recoil and transverse mass
        ###

        met['T']  = TVector2Array.from_polar(met.pt, met.phi)

        u = {
            'sr'    : met.T,
            'wecr'  : met.T+leading_e.T.sum(),
            'tecr'  : met.T+leading_e.T.sum(),
            'wmcr'  : met.T+leading_mu.T.sum(),
            'tmcr'  : met.T+leading_mu.T.sum(),
            'qcdcr' : met.T,
        }

        mT = {
            'wecr'  : np.sqrt(2*leading_e.pt.sum()*met.pt*(1-np.cos(met.T.delta_phi(leading_e.T.sum())))),
            'tecr'  : np.sqrt(2*leading_e.pt.sum()*met.pt*(1-np.cos(met.T.delta_phi(leading_e.T.sum())))),
            'wmcr'  : np.sqrt(2*leading_mu.pt.sum()*met.pt*(1-np.cos(met.T.delta_phi(leading_mu.T.sum())))),
            'tmcr'  : np.sqrt(2*leading_mu.pt.sum()*met.pt*(1-np.cos(met.T.delta_phi(leading_mu.T.sum())))) 
        }

        ###
        #Calculating weights
        ###
        if not isData:
            
            gen = events.GenPart

            gen['isb'] = (abs(gen.pdgId)==5)&gen.hasFlags(['fromHardProcess', 'isLastCopy'])
            gen['isc'] = (abs(gen.pdgId)==4)&gen.hasFlags(['fromHardProcess', 'isLastCopy'])
            gen['isTop'] = (abs(gen.pdgId)==6)&gen.hasFlags(['fromHardProcess', 'isLastCopy'])
            genTops = gen[gen.isTop]
            nlo = np.ones(events.size)
            if('TTJets' in dataset): 
                nlo = np.sqrt(get_ttbar_weight(genTops[:,0].pt.sum()) * get_ttbar_weight(genTops[:,1].pt.sum()))
                
            gen['isW'] = (abs(gen.pdgId)==24)&gen.hasFlags(['fromHardProcess', 'isLastCopy'])
            gen['isZ'] = (abs(gen.pdgId)==23)&gen.hasFlags(['fromHardProcess', 'isLastCopy'])
            
            genWs = gen[gen.isW] 
            genZs = gen[gen.isZ]
            genDYs = gen[gen.isZ&(gen.mass>30)]
            
            nnlo_nlo = {}
            nlo_qcd = np.ones(events.size)
            nlo_ewk = np.ones(events.size)
            if('WJets' in dataset): 
                nlo_qcd = get_nlo_qcd_weight['w'](genWs.pt.max())
                nlo_ewk = get_nlo_ewk_weight['w'](genWs.pt.max())
                for systematic in get_nnlo_nlo_weight['w']:
                    nnlo_nlo[systematic] = get_nnlo_nlo_weight['w'][systematic](genWs.pt.max())*((genWs.counts>0)&(genWs.pt.max()>=100)) + \
                                           (~((genWs.counts>0)&(genWs.pt.max()>=100))).astype(np.int)
            elif('DY' in dataset): 
                nlo_qcd = get_nlo_qcd_weight['dy'](genDYs.pt.max())
                nlo_ewk = get_nlo_ewk_weight['dy'](genDYs.pt.max())
                for systematic in get_nnlo_nlo_weight['dy']:
                    nnlo_nlo[systematic] = get_nnlo_nlo_weight['dy'][systematic](genDYs.pt.max())*((genDYs.counts>0)&(genDYs.pt.max()>=100)) + \
                                           (~((genDYs.counts>0)&(genDYs.pt.max()>=100))).astype(np.int)
            elif('ZJets' in dataset): 
                nlo_qcd = get_nlo_qcd_weight['z'](genZs.pt.max())
                nlo_ewk = get_nlo_ewk_weight['z'](genZs.pt.max())
                for systematic in get_nnlo_nlo_weight['z']:
                    nnlo_nlo[systematic] = get_nnlo_nlo_weight['z'][systematic](genZs.pt.max())*((genZs.counts>0)&(genZs.pt.max()>=100)) + \
                                           (~((genZs.counts>0)&(genZs.pt.max()>=100))).astype(np.int)

            ###
            # Calculate PU weight and systematic variations
            ###

            pu = get_pu_weight(events.Pileup.nTrueInt)

            ###
            # Trigger efficiency weight
            ###

            trig = {
                'sr':   get_met_trig_weight(met.pt),
                'wmcr': get_met_trig_weight(u['wmcr'].mag),
                'tmcr': get_met_trig_weight(u['tmcr'].mag),
                'wecr': get_ele_trig_weight(leading_e.eta.sum()+leading_e.deltaEtaSC.sum(), leading_e.pt.sum()),
                'tecr': get_ele_trig_weight(leading_e.eta.sum()+leading_e.deltaEtaSC.sum(), leading_e.pt.sum()),
                'qcdcr':   get_met_trig_weight(met.pt),
            }

            ### 
            # Calculating electron and muon ID weights
            ###

            mueta = abs(leading_mu.eta.sum())
            if self._year=='2016':
                mueta=leading_mu.eta.sum()
            ids ={
                'sr':  np.ones(events.size),
                'wmcr': get_mu_tight_id_sf(mueta,leading_mu.pt.sum()),
                'tmcr': get_mu_tight_id_sf(mueta,leading_mu.pt.sum()),
                'wecr': get_ele_tight_id_sf(leading_e.eta.sum()+leading_e.deltaEtaSC.sum(),leading_e.pt.sum()),
                'tecr': get_ele_tight_id_sf(leading_e.eta.sum()+leading_e.deltaEtaSC.sum(),leading_e.pt.sum()),
                'qcdcr':  np.ones(events.size),
            }

            ###
            # Reconstruction weights for electrons
            ###

            def ele_reco_sf(pt, eta):#2017 has separate weights for low/high pT (threshold at 20 GeV)
                return get_ele_reco_sf(eta, pt)*(pt>20).astype(np.int) + get_ele_reco_lowet_sf(eta, pt)*(~(pt>20)).astype(np.int)

            if self._year == '2017':
                sf = ele_reco_sf
            else:
                sf = get_ele_reco_sf

            reco = {
                'sr': np.ones(events.size),
                'wmcr': np.ones(events.size),
                'tmcr': np.ones(events.size),
                'wecr': sf(leading_e.eta.sum()+leading_e.deltaEtaSC.sum(),leading_e.pt.sum()),
                'tecr': sf(leading_e.eta.sum()+leading_e.deltaEtaSC.sum(),leading_e.pt.sum()),
                'qcdcr': np.ones(events.size),
            }

            ###
            # Isolation weights for muons
            ###

            isolation = {
                'sr'  : np.ones(events.size),
                'wmcr': get_mu_tight_iso_sf(mueta,leading_mu.pt.sum()),
                'tmcr': get_mu_tight_iso_sf(mueta,leading_mu.pt.sum()),
                'wecr': np.ones(events.size),
                'tecr': np.ones(events.size),
                'qcdcr'  : np.ones(events.size),
            }

            ###
            # AK4 b-tagging weights
            ###

            btagSF, \
            btagSFbc_correlatedUp, \
            btagSFbc_correlatedDown, \
            btagSFbc_uncorrelatedUp, \
            btagSFbc_uncorrelatedDown, \
            btagSFlight_correlatedUp, \
            btagSFlight_correlatedDown, \
            btagSFlight_uncorrelatedUp, \
            btagSFlight_uncorrelatedDown  = get_deepflav_weight['loose'](j_iso.pt,j_iso.eta,j_iso.hadronFlavour,j_iso.isdflvL)

            if 'L1PreFiringWeight' in events.columns: 
                weights.add('prefiring', events.L1PreFiringWeight.Nom, events.L1PreFiringWeight.Up, events.L1PreFiringWeight.Dn)
            weights.add('genw',events.genWeight)
            weights.add('nlo_qcd',nlo_qcd)
            weights.add('nlo_ewk',nlo_ewk)
            if 'cen' in nnlo_nlo:
                #weights.add('nnlo_nlo',nnlo_nlo['cen'])
                weights.add('qcd1',np.ones(events.size), nnlo_nlo['qcd1up']/nnlo_nlo['cen'], nnlo_nlo['qcd1do']/nnlo_nlo['cen'])
                weights.add('qcd2',np.ones(events.size), nnlo_nlo['qcd2up']/nnlo_nlo['cen'], nnlo_nlo['qcd2do']/nnlo_nlo['cen'])
                weights.add('qcd3',np.ones(events.size), nnlo_nlo['qcd3up']/nnlo_nlo['cen'], nnlo_nlo['qcd3do']/nnlo_nlo['cen'])
                weights.add('ew1',np.ones(events.size), nnlo_nlo['ew1up']/nnlo_nlo['cen'], nnlo_nlo['ew1do']/nnlo_nlo['cen'])
                weights.add('ew2G',np.ones(events.size), nnlo_nlo['ew2Gup']/nnlo_nlo['cen'], nnlo_nlo['ew2Gdo']/nnlo_nlo['cen'])
                weights.add('ew3G',np.ones(events.size), nnlo_nlo['ew3Gup']/nnlo_nlo['cen'], nnlo_nlo['ew3Gdo']/nnlo_nlo['cen'])
                weights.add('ew2W',np.ones(events.size), nnlo_nlo['ew2Wup']/nnlo_nlo['cen'], nnlo_nlo['ew2Wdo']/nnlo_nlo['cen'])
                weights.add('ew3W',np.ones(events.size), nnlo_nlo['ew3Wup']/nnlo_nlo['cen'], nnlo_nlo['ew3Wdo']/nnlo_nlo['cen'])
                weights.add('ew2Z',np.ones(events.size), nnlo_nlo['ew2Zup']/nnlo_nlo['cen'], nnlo_nlo['ew2Zdo']/nnlo_nlo['cen'])
                weights.add('ew3Z',np.ones(events.size), nnlo_nlo['ew3Zup']/nnlo_nlo['cen'], nnlo_nlo['ew3Zdo']/nnlo_nlo['cen'])
                weights.add('mix',np.ones(events.size), nnlo_nlo['mixup']/nnlo_nlo['cen'], nnlo_nlo['mixdo']/nnlo_nlo['cen'])
                weights.add('muF',np.ones(events.size), nnlo_nlo['muFup']/nnlo_nlo['cen'], nnlo_nlo['muFdo']/nnlo_nlo['cen'])
                weights.add('muR',np.ones(events.size), nnlo_nlo['muRup']/nnlo_nlo['cen'], nnlo_nlo['muRdo']/nnlo_nlo['cen'])
            weights.add('pileup',pu)
            weights.add('trig', trig[region])
            weights.add('ids', ids[region])
            weights.add('reco', reco[region])
            weights.add('isolation', isolation[region])
            weights.add('btagSF',btagSF)
            weights.add('btagSFbc_correlated',np.ones(events.size), btagSFbc_correlatedUp/btagSF, btagSFbc_correlatedDown/btagSF)
            weights.add('btagSFbc_uncorrelated',np.ones(events.size), btagSFbc_uncorrelatedUp/btagSF, btagSFbc_uncorrelatedDown/btagSF)
            weights.add('btagSFlight_correlated',np.ones(events.size), btagSFlight_correlatedUp/btagSF, btagSFlight_correlatedDown/btagSF)
            weights.add('btagSFlight_uncorrelated',np.ones(events.size), btagSFlight_uncorrelatedUp/btagSF, btagSFlight_uncorrelatedDown/btagSF)

        
        ###
        # Selections
        ###

        lumimask = np.ones(events.size, dtype=np.bool)
        if Data:
            lumimask = lumiMasks[self._year](events.run, events.luminosityBlock)
        selection.add('lumimask', lumimask)

        met_filters =  np.ones(events.size, dtype=np.bool)
        if isData: met_filters = met_filters & events.Flag['eeBadScFilter']#this filter is recommended for data only
        for flag in AnalysisProcessor.met_filter_flags[self._year]:
            met_filters = met_filters & events.Flag[flag]
        selection.add('met_filters',met_filters)

        triggers = np.zeros(events.size, dtype=np.bool)
        for path in self._met_triggers[self._year]:
            if path not in events.HLT.columns: continue
            triggers = triggers | events.HLT[path]
        selection.add('met_triggers', triggers)

        triggers = np.zeros(events.size, dtype=np.bool)
        for path in self._singleelectron_triggers[self._year]:
            if path not in events.HLT.columns: continue
            triggers = triggers | events.HLT[path]
        selection.add('singleelectron_triggers', triggers)

        noHEMj = np.ones(events.size, dtype=np.bool)
        if self._year=='2018': noHEMj = (j_nHEM==0)

        noHEMmet = np.ones(events.size, dtype=np.bool)
        if self._year=='2018': noHEMmet = (met.pt>470)|(met.phi>-0.62)|(met.phi<-1.62)

        leading_fj = fj[fj.sd.pt.argmax()]
        leading_fj = leading_fj[leading_fj.isgood]
        leading_fj = leading_fj[leading_fj.isclean]
        selection.add('iszeroL', (e_nloose==0)&(mu_nloose==0)&(tau_nloose==0)&(pho_nloose==0))
        selection.add('isoneM', (e_nloose==0)&(mu_ntight==1)&(mu_nloose==1)&(tau_nloose==0)&(pho_nloose==0))
        selection.add('isoneE', (e_ntight==1)&(e_nloose==1)&(mu_nloose==0)&(tau_nloose==0)&(pho_nloose==0))
        selection.add('leading_e_pt',(e_loose.pt.max()>40))
        selection.add('noextrab', (j_ndflvL==0))
        selection.add('extrab', (j_ndflvL>0))
        selection.add('fatjet', (fj_nclean>0))
        selection.add('noHEMj', noHEMj)
        selection.add('noHEMmet', noHEMmet)
        selection.add('met120',(met.pt<120))
        selection.add('met100',(met.pt>100))
        selection.add('msd40',(leading_fj.msd_corr.sum()>40))
        selection.add('recoil_qcdcr', (u['qcdcr'].mag>250))
        selection.add('mindphi_qcdcr', (abs(u['qcdcr'].delta_phi(j_clean.T)).min()<0.1))
        selection.add('minDphi_qcdcr', (abs(u['qcdcr'].delta_phi(fj_clean.T)).min()>1.5))
        selection.add('calo_qcdcr', ( (abs(calomet.pt - met.pt) / u['qcdcr'].mag)<0.5))
            
        #selection.add('mindphimet',(abs(met.T.delta_phi(j_clean.T)).min())>0.7)

        regions = {
            #'sr': ['iszeroL','fatjet','noextrab','noHEMmet','met_filters','met_triggers','noHEMj'],
            'sr': ['msd40','fatjet', 'noHEMj','iszeroL','noextrab','met_filters','met_triggers','noHEMmet'],
            'wmcr': ['msd40','isoneM','fatjet','noextrab','noHEMj','met_filters','met_triggers'],
            'tmcr': ['msd40','isoneM','fatjet','extrab','noHEMj','met_filters','met_triggers'],
            'wecr': ['msd40','isoneE','fatjet','noextrab','noHEMj','met_filters','singleelectron_triggers','met100'],
            'tecr': ['msd40','isoneE','fatjet','extrab','noHEMj','met_filters','singleelectron_triggers','met100'],
            'qcdcr': ['recoil_qcdcr','mindphi_qcdcr','minDphi_qcdcr','calo_qcdcr','msd40','fatjet', 'noHEMj','iszeroL','noextrab','met_filters','met_triggers','noHEMmet'],
        }

        def normalize(val, cut):
            if cut is None:
                ar = ak.to_numpy(ak.fill_none(val, np.nan))
                return ar
            else:
                ar = ak.to_numpy(ak.fill_none(val[cut], np.nan))
                return ar
                
        def fill(region, systematic):
            cut = selection.all(*regions[region])
            sname = 'nominal' if systematic is None else systematic
            if systematic in weights.variations:
                weight = weights.weight(modifier=systematic)[cut]
            else:
                weight = weights.weight()[cut]
            output['template'].fill(
                  region=region,
                  systematic=sname,
                  recoil=normalize(u[region].mag, cut),
                  fjmass=normalize(leading_fj.msd_corr, cut),
                  ZHbbvsQCD=normalize(leading_fj.ZHbbvsQCD, cut),
                  weight=weight
            )
            if systematic is None:
                variables = {
                    'mindphirecoil':          abs(u[region].delta_phi(j_clean.T)).min(),
                    'minDphirecoil':          abs(u[region].delta_phi(fj_clean.T)).min(),
                    'CaloMinusPfOverRecoil':  abs(calomet.pt - met.pt) / u[region].mag,
                    'met':                    met.pt.flatten(),
                    'metphi':                 met.phi.flatten(),
                    'mindphimet':             abs(met.T.delta_phi(j_clean.T)).min(),
                    'minDphimet':             abs(met.T.delta_phi(fj_clean.T)).min(),
                    'j1pt':                   leading_j.pt.sum(),
                    'j1eta':                  leading_j.eta.sum(),
                    'j1phi':                  leading_j.phi.sum(),
                    'fj1pt':                  leading_fj.sd.pt.sum(),
                    'fj1eta':                 leading_fj.sd.eta.sum(),
                    'fj1phi':                 leading_fj.sd.phi.sum(),
                    'njets':                  j_nclean,
                    'ndflvL':                 j_ndflvL,
                    'nfjclean':               fj_nclean,
                }
                if region in mT:
                    variables['mT']           = mT[region]
                if 'e' in region:
                    variables['l1pt']      = leading_e.pt.sum()
                    variables['l1phi']     = leading_e.phi.sum()
                    variables['l1eta']     = leading_e.eta.sum()
                if 'm' in region:
                    variables['l1pt']      = leading_mu.pt.sum()
                    variables['l1phi']     = leading_mu.phi.sum()
                    variables['l1eta']     = leading_mu.eta.sum()
                for variable in output:
                    if variable not in variables:
                        continue
                    normalized_variable = {variable: normalize(variables[variable],cut)}
                    output[variable].fill(
                        region=region,
                        ZHbbvsQCD=normalize(leading_fj.ZHbbvsQCD,cut),
                        **normalized_variable,
                        weight=weight,
                    )
                output['ZHbbvsQCD'].fill(
                      region=region,
                      ZHbbvsQCD=normalize(leading_fj.ZHbbvsQCD, cut),
                      weight=weight
                )

        if shift_name is None:
            systematics = [None] + list(weights.variations)
        else:
            systematics = [shift_name]
            
        for region in regions:
            if region not in selected_regions: continue

            ###
            # Adding recoil and minDPhi requirements
            ###

            selection.add('recoil_'+region, (u[region].mag>250))
            selection.add('mindphi_'+region, (abs(u[region].delta_phi(j_clean.T)).min()>0.5))
            selection.add('minDphi_'+region, (abs(u[region].delta_phi(fj_clean.T)).min()>1.5))
            selection.add('calo_'+region, ( (abs(calomet.pt - met.pt) / u[region].mag) < 0.5))
            if 'qcd' not in region:
                regions[region].insert(0, 'recoil_'+region)
                regions[region].insert(3, 'mindphi_'+region)
                regions[region].insert(4, 'minDphi_'+region)
                regions[region].insert(5, 'calo_'+region)

            for systematic in systematics:
                if isData and systematic is not None:
                    continue
                fill(region, systematic)
                
        return dataset, output

    def postprocess(self, accumulator):
        dataset, output = accumulator
        print('Scaling:', dataset)
        print('Cross section:',self._xsec[dataset])

        scale = 1
        if self._xsec[dataset]!= -1: 
            scale = self._lumi*self._xsec[dataset]

        for key in output:
            if key=='sumw': 
                continue
            output[key] *= scale
            
        return output

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('-y', '--year', help='year', dest='year')
    parser.add_option('-m', '--metadata', help='metadata', dest='metadata')
    parser.add_option('-n', '--name', help='name', dest='name')
    (options, args) = parser.parse_args()


    with open('metadata/'+options.metadata+'.json') as fin:
        samplefiles = json.load(fin)
        xsec = {k: v['xs'] for k,v in samplefiles.items()}

    corrections = load('data/corrections.coffea')
    ids         = load('data/ids.coffea')
    common      = load('data/common.coffea')

    processor_instance=AnalysisProcessor(year=options.year,
                                         xsec=xsec,
                                         corrections=corrections,
                                         ids=ids,
                                         common=common)

    save(processor_instance, 'data/darkhiggs'+options.name+'.processor')

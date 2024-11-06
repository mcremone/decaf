#!/usr/bin/env python
import logging
import numpy as np
import awkward as ak
import json
import copy
from collections import defaultdict
from coffea import processor
import cachetools
import hist
from coffea.analysis_tools import Weights, PackedSelection
from coffea.lumi_tools import LumiMask
from coffea.util import load, save
from optparse import OptionParser
from coffea.nanoevents.methods import vector
import gzip

def update(events, collections):
    """Return a shallow copy of events array with some collections swapped out"""
    out = events
    for name, value in collections.items():
        out = ak.with_field(out, value, name)
    return out

class AnalysisProcessor(processor.ProcessorABC):

    lumis = { 
        #Values from https://twiki.cern.ch/twiki/bin/view/CMS/LumiRecommendationsRun2                                                      
        '2016postVFP': 16.81,
        '2016preVFP': 19.52,
        '2017': 41.48,
        '2018': 59.83
    }

    lumiMasks = {
        '2016postVFP': LumiMask("data/jsons/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt"),
        '2016preVFP': LumiMask("data/jsons/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt"),
        '2017': LumiMask("data/jsons/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt"),
        '2018': LumiMask("data/jsons/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt"),
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
        self._systematics = False
        self._skipJER = False

        self._samples = {
            'sr':('Z1Jets','Z2Jets','WJets','DY','TT','ST','WW','WZ','ZZ','QCD','MET','TPhiTo2Chi'),
            'wmcr':('WJets','DY','TT','ST','WW','WZ','ZZ','QCD','MET'),
            'tmcr':('WJets','DY','TT','ST','WW','WZ','ZZ','QCD','MET'),
            'wecr':('WJets','DY','TT','ST','WW','WZ','ZZ','QCD','SingleElectron','EGamma'),
            'tecr':('WJets','DY','TT','ST','WW','WZ','ZZ','QCD','SingleElectron','EGamma'),
            'zmcr':('DY','TT','ST','WW','WZ','ZZ','QCD','MET'),
            'zecr':('DY','TT','ST','WW','WZ','ZZ','QCD','SingleElectron','EGamma'),
            'gcr' :('G1Jet','QCD','SinglePhoton','EGamma')
        }
        
        self._TvsQCDwp = { ## for 2018 top-tagger 0.26
            '2016preVFP': 0.26,
            '2016postVFP': 0.26,
            '2017': 0.26,
            '2018': 0.26
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
        self._singlephoton_triggers = {
            '2016': [
                'Photon175',
                'Photon165_HE10'
            ],
            '2017': [
                'Photon200'
            ],
            '2018': [
                'Photon200'
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
                'Photon200'
            ],
            '2018': [
                'Ele32_WPTight_Gsf',
                'Photon200'
            ]
        }
        self._singlemuon_triggers = {
            '2016': [
                'IsoMu24',
                'IsoTkMu24',
            ],
            '2017':
                [
                'IsoMu27',
            ],
            '2018':
                [
                'IsoMu24',
            ]
        }
        self._corrections = corrections
        self._ids = ids
        self._common = common

        #ptbins=np.arange(250,1210,38)

        self.make_output = lambda: {
            'sumw': 0.,
            'cutflow': hist.Hist(
                hist.axis.StrCategory([], name='region', growth=True),
                hist.axis.StrCategory([], name='cutname', growth=True),
                hist.axis.Variable([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17], name='cutflow', label='cut'),
                storage=hist.storage.Weight(),
            ),
            'TvsQCD': hist.Hist(
                hist.axis.StrCategory([], name='region', growth=True),
                hist.axis.Regular(15,0,1, name='TvsQCD', label='TvsQCD'),
                storage=hist.storage.Weight(),
            ),
            'mindphirecoil': hist.Hist(
                hist.axis.StrCategory([], name='region', growth=True),
                hist.axis.Regular(30,0,3.5, name='mindphirecoil', label='Min dPhi(Recoil,AK4s)'),
                hist.axis.Variable([0, self._TvsQCDwp[self._year], 1], name='TvsQCD', label='TvsQCD', flow=False),
                storage=hist.storage.Weight(),
            ),
            'minDphirecoil': hist.Hist(
                hist.axis.StrCategory([], name='region', growth=True),
                hist.axis.Regular(30,0,3.5, name='minDphirecoil', label='Min dPhi(Recoil, leading AK15s)'),
                hist.axis.Variable([0, self._TvsQCDwp[self._year], 1], name='TvsQCD', label='TvsQCD', flow=False),
                storage=hist.storage.Weight(),
            ),
            'ut': hist.Hist(
                hist.axis.StrCategory([], name='region', growth=True),
                hist.axis.Regular(65,350,1000, name='ut', label='U_{T}'),
                hist.axis.Variable([0, self._TvsQCDwp[self._year], 1], name='TvsQCD', label='TvsQCD', flow=False),
                storage=hist.storage.Weight(),
            ),
            'uphi': hist.Hist(
                hist.axis.StrCategory([], name='region', growth=True),
                hist.axis.Regular(35,-3.5,3.5, name='uphi', label='U_{phi}'),
                hist.axis.Variable([0, self._TvsQCDwp[self._year], 1], name='TvsQCD', label='TvsQCD', flow=False),
                storage=hist.storage.Weight(),
            ),
            'met': hist.Hist(
                hist.axis.StrCategory([], name='region', growth=True),
                hist.axis.Regular(30,0,1200, name='met', label='MET $p_{T}$ [GeV]'),
                hist.axis.Variable([0, self._TvsQCDwp[self._year], 1], name='TvsQCD', label='TvsQCD', flow=False),
                storage=hist.storage.Weight(),
            ),
            'metphi': hist.Hist(
                hist.axis.StrCategory([], name='region', growth=True),
                hist.axis.Regular(35,-3.5,3.5, name='metphi', label='MET $\phi$'),
                hist.axis.Variable([0, self._TvsQCDwp[self._year], 1], name='TvsQCD', label='TvsQCD', flow=False),
                storage=hist.storage.Weight(),
            ),
            'j1pt': hist.Hist(
                hist.axis.StrCategory([], name='region', growth=True),
                hist.axis.Regular(38,250,1200, name='j1pt', label='AK4 Leading Jet $p_{T}$'),
                hist.axis.Variable([0, self._TvsQCDwp[self._year], 1], name='TvsQCD', label='TvsQCD', flow=False),
                storage=hist.storage.Weight(),
            ),
            'j1eta': hist.Hist(
                hist.axis.StrCategory([], name='region', growth=True),
                hist.axis.Regular(35,-3.5,3.5, name='j1eta', label='AK4 Leading Jet Eta'),
                hist.axis.Variable([0, self._TvsQCDwp[self._year], 1], name='TvsQCD', label='TvsQCD', flow=False),
                storage=hist.storage.Weight(),
            ),
            'j1phi': hist.Hist(
                hist.axis.StrCategory([], name='region', growth=True),
                hist.axis.Regular(35,-3.5,3.5, name='j1phi', label='AK4 Leading Jet Phi'),
                hist.axis.Variable([0, self._TvsQCDwp[self._year], 1], name='TvsQCD', label='TvsQCD', flow=False),
                storage=hist.storage.Weight(),
            ),
            'fj1pt': hist.Hist(
                hist.axis.StrCategory([], name='region', growth=True),
                hist.axis.Regular(38,250,1200, name='fj1pt', label='AK15 Leading SoftDrop Jet Pt'),
                hist.axis.Variable([0, self._TvsQCDwp[self._year], 1], name='TvsQCD', label='TvsQCD', flow=False),
                storage=hist.storage.Weight(),
            ),
            'fj1eta': hist.Hist(
                hist.axis.StrCategory([], name='region', growth=True),
                hist.axis.Regular(35,-3.5,3.5, name='fj1eta', label='AK15 Leading SoftDrop Jet Eta'),
                hist.axis.Variable([0, self._TvsQCDwp[self._year], 1], name='TvsQCD', label='TvsQCD', flow=False),
                storage=hist.storage.Weight(),
            ),
            'fj1phi': hist.Hist(
                hist.axis.StrCategory([], name='region', growth=True),
                hist.axis.Regular(35,-3.5,3.5, name='fj1phi', label='AK15 Leading SoftDrop Jet Phi'),
                hist.axis.Variable([0, self._TvsQCDwp[self._year], 1], name='TvsQCD', label='TvsQCD', flow=False),
                storage=hist.storage.Weight(),
            ),
            'njets': hist.Hist(
                hist.axis.StrCategory([], name='region', growth=True),
                hist.axis.IntCategory([0, 1, 2, 3, 4, 5, 6], name='njets', label='AK4 Number of Jets'),
                hist.axis.Variable([0, self._TvsQCDwp[self._year], 1], name='TvsQCD', label='TvsQCD', flow=False),
                storage=hist.storage.Weight(),
            ),
            'ndflvL': hist.Hist(
                hist.axis.StrCategory([], name='region', growth=True),
                hist.axis.IntCategory([0, 1, 2, 3, 4, 5, 6], name='ndflvL', label='AK4 Number of deepFlavor Loose Jets'),
                hist.axis.Variable([0, self._TvsQCDwp[self._year], 1], name='TvsQCD', label='TvsQCD', flow=False),
                storage=hist.storage.Weight(),
            ),
            'nfjclean': hist.Hist(
                hist.axis.StrCategory([], name='region', growth=True),
                hist.axis.IntCategory([0, 1, 2, 3, 4], name='nfjclean', label='AK15 Number of Cleaned Jets'),
                hist.axis.Variable([0, self._TvsQCDwp[self._year], 1], name='TvsQCD', label='TvsQCD', flow=False),
                storage=hist.storage.Weight(),
            ),
            'mT': hist.Hist(
                hist.axis.StrCategory([], name='region', growth=True),
                hist.axis.Regular(20,0,600, name='mT', label='Transverse Mass'),
                hist.axis.Variable([0, self._TvsQCDwp[self._year], 1], name='TvsQCD', label='TvsQCD', flow=False),
                storage=hist.storage.Weight(),
            ),
            'l1pt': hist.Hist(
                hist.axis.StrCategory([], name='region', growth=True),
                hist.axis.Regular(38,250,1200, name='l1pt', label='Leading Lepton Pt'),
                hist.axis.Variable([0, self._TvsQCDwp[self._year], 1], name='TvsQCD', label='TvsQCD', flow=False),
                storage=hist.storage.Weight(),
            ),
            'l1eta': hist.Hist(
                hist.axis.StrCategory([], name='region', growth=True),
                hist.axis.Regular(48,-2.4,2.4, name='l1eta', label='Leading Lepton Eta'),
                hist.axis.Variable([0, self._TvsQCDwp[self._year], 1], name='TvsQCD', label='TvsQCD', flow=False),
                storage=hist.storage.Weight(),
            ),
            'l1phi': hist.Hist(
                hist.axis.StrCategory([], name='region', growth=True),
                hist.axis.Regular(64,-3.2,3.2, name='l1phi', label='Leading Lepton Phi'),
                hist.axis.Variable([0, self._TvsQCDwp[self._year], 1], name='TvsQCD', label='TvsQCD', flow=False),
                storage=hist.storage.Weight(),
            ),
    }

    def process(self, events):
        isData = not hasattr(events, "genWeight")
        if isData:
            # Nominal JEC are already applied in data
            return self.process_shift(events, None)

        jet_factory              = self._corrections['jet_factory']
        fatjet_factory           = self._corrections['fatjet_factory']
        subjet_factory           = self._corrections['subjet_factory']
        met_factory              = self._corrections['met_factory']

        
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
        fatjets = fatjet_factory[thekey].build(add_jec_variables(events.AK15PFPuppiJet, events.fixedGridRhoFastjetAll), jec_cache)
        #subjets = subjet_factory[thekey].build(add_jec_variables(events.AK15PFPuppiSubJet, events.fixedGridRhoFastjetAll), jec_cache)
        met = met_factory.build(events.MET, jets, {})

        shifts = [({"Jet": jets, "AK15PFPuppiJet": fatjets, "MET": met}, None)]
        if self._systematics:
            shifts.extend([
                ({"Jet": jets.JES_jes.up, "AK15PFPuppiJet": fatjets.JES_jes.up, "MET": met.JES_jes.up}, "JESUp"),
                ({"Jet": jets.JES_jes.down, "AK15PFPuppiJet": fatjets.JES_jes.down, "MET": met.JES_jes.down}, "JESDown"),
                ({"Jet": jets, "AK15PFPuppiJet": fatjets, "MET": met.MET_UnclusteredEnergy.up}, "UESUp"),
                ({"Jet": jets, "AK15PFPuppiJet": fatjets, "MET": met.MET_UnclusteredEnergy.down}, "UESDown"),
            ])
            if not self._skipJER:
                shifts.extend([
                    ({"Jet": jets.JER.up, "AK15PFPuppiJet": fatjets.JER.up, "MET": met.JER.up}, "JERUp"),
                    ({"Jet": jets.JER.down, "AK15PFPuppiJet": fatjets.JER.down, "MET": met.JER.down}, "JERDown"),
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
        selection = PackedSelection(dtype="uint64")
        weights = Weights(len(events), storeIndividual=True)
        output = self.make_output()
        if shift_name is None and not isData:
            output['sumw'] = ak.sum(events.genWeight)

        ###
        #Getting corrections, ids from .coffea files
        ###

        get_met_trig_weight      = self._corrections['get_met_trig_weight']
        get_ele_loose_id_sf      = self._corrections['get_ele_loose_id_sf']
        get_ele_tight_id_sf      = self._corrections['get_ele_tight_id_sf']
        get_ele_trig_weight      = self._corrections['get_ele_trig_weight']
        get_ele_reco_sf_below20  = self._corrections['get_ele_reco_sf_below20']
        get_ele_reco_sf_above20  = self._corrections['get_ele_reco_sf_above20']
        get_mu_loose_id_sf       = self._corrections['get_mu_loose_id_sf']
        get_mu_tight_id_sf       = self._corrections['get_mu_tight_id_sf']
        get_mu_loose_iso_sf      = self._corrections['get_mu_loose_iso_sf']
        get_mu_tight_iso_sf      = self._corrections['get_mu_tight_iso_sf']
        #get_mu_rochester_sf      = self._corrections['get_mu_rochester_sf'][self._year]
        get_met_xy_correction    = self._corrections['get_met_xy_correction']
        get_pu_weight            = self._corrections['get_pu_weight']    
        get_nlo_ewk_weight       = self._corrections['get_nlo_ewk_weight']    
        get_nnlo_nlo_weight      = self._corrections['get_nnlo_nlo_weight']
        get_msd_corr             = self._corrections['get_msd_corr']
        get_btag_weight      = self._corrections['get_btag_weight']
        get_ttbar_weight     = self._corrections['get_ttbar_weight']
        
        isLooseElectron = self._ids['isLooseElectron'] 
        isTightElectron = self._ids['isTightElectron'] 
        isLooseMuon     = self._ids['isLooseMuon']     
        isTightMuon     = self._ids['isTightMuon']     
        #isLooseTau      = self._ids['isLooseTau']      
        isLoosePhoton   = self._ids['isLoosePhoton']   
        isTightPhoton   = self._ids['isTightPhoton']   
        isGoodAK4       = self._ids['isGoodAK4']       
        isGoodAK15    = self._ids['isGoodAK15']    
        isHEMJet        = self._ids['isHEMJet']        
        
        deepflavWPs = self._common['btagWPs']['deepflav'][self._year]
        deepcsvWPs = self._common['btagWPs']['deepcsv'][self._year]
        sigmoid = self._common['sigmoid'] #to calculate photon trigger efficiency

        ###
        #Initialize global quantities (MET ecc.)
        ###

        npv = events.PV.npvsGood
        run = events.run
        calomet = events.CaloMET
        met = events.MET
        met['pt'] , met['phi'] = get_met_xy_correction(self._year, npv, run, met.pt, met.phi, isData)

        ###
        #Initialize physics objects
        ###

        mu = events.Muon
        mu['isloose'] = isLooseMuon(mu)
        mu['id_sf'] = ak.where(
            mu.isloose, 
            get_mu_loose_id_sf(self._year, abs(mu.eta), mu.pt), 
            ak.ones_like(mu.pt)
        )
        mu['iso_sf'] = ak.where(
            mu.isloose, 
            get_mu_loose_iso_sf(self._year, abs(mu.eta), mu.pt), 
            ak.ones_like(mu.pt)
        )
        mu['istight'] = isTightMuon(mu)
        mu['id_sf'] = ak.where(
            mu.istight, 
            get_mu_tight_id_sf(self._year, abs(mu.eta), mu.pt), 
            mu.id_sf
        )
        mu['iso_sf'] = ak.where(
            mu.istight, 
            get_mu_tight_iso_sf(self._year, abs(mu.eta), mu.pt), 
            mu.iso_sf
        )
        mu['T'] = ak.zip(
            {
                "r": mu.pt,
                "phi": mu.phi,
            },
            with_name="PolarTwoVector",
            behavior=vector.behavior,
        )
        mu_loose=mu[mu.isloose]
        mu_tight=mu[mu.istight]
        mu_ntot = ak.num(mu, axis=1)
        mu_nloose = ak.num(mu_loose, axis=1)
        mu_ntight = ak.num(mu_tight, axis=1)
        leading_mu = ak.firsts(mu_tight)
        mu_pairs = ak.combinations(mu_loose, 2)
        dimu = mu_pairs["0"]+mu_pairs["1"]
        dimu['T'] = ak.zip(
            {
                "r": dimu.pt,
                "phi": dimu.phi,
            },
            with_name="PolarTwoVector",
            behavior=vector.behavior,
        )
        leading_mu_pair = ak.firsts(mu_pairs[ak.argmax(dimu.pt, axis=1, keepdims=True)])
        leading_dimu = ak.firsts(dimu[ak.argmax(dimu.pt, axis=1, keepdims=True)])

        e = events.Electron
        ## ReMOVE e cleanning from mu (eta-phi 0.3) 2024.4.1 (Not joke!)
        e['isclean'] = ak.all(e.metric_table(mu_loose) > -99999, axis=2)
        e['reco_sf'] = ak.where(
            (e.pt<20),
            get_ele_reco_sf_below20(self._year, e.eta+e.deltaEtaSC, e.pt), 
            get_ele_reco_sf_above20(self._year, e.eta+e.deltaEtaSC, e.pt)
        )
        e['isloose'] = isLooseElectron(e)
        e['id_sf'] = ak.where(
            e.isloose,
            get_ele_loose_id_sf(self._year, e.eta+e.deltaEtaSC, e.pt),
            ak.ones_like(e.pt)
        )
        e['istight'] = isTightElectron(e)
        e['id_sf'] = ak.where(
            e.istight,
            get_ele_tight_id_sf(self._year, e.eta+e.deltaEtaSC, e.pt),
            e.id_sf
        )
        e['T'] = ak.zip(
            {
                "r": e.pt,
                "phi": e.phi,
            },
            with_name="PolarTwoVector",
            behavior=vector.behavior,
        )
        e_clean = e[e.isclean]
        e_loose = e_clean[e_clean.isloose]
        e_tight = e_clean[e_clean.istight]
        e_ntot = ak.num(e, axis=1)
        e_nloose = ak.num(e_loose, axis=1)
        e_ntight = ak.num(e_tight, axis=1)
        leading_e = ak.firsts(e_tight)
        ele_pairs = ak.combinations(ele_loose, 2)
        diele = ele_pairs["0"]+ele_pairs["1"]
        diele['T'] = ak.zip(
            {
                "r": diele.pt,
                "phi": diele.phi,
            },
            with_name="PolarTwoVector",
            behavior=vector.behavior,
        )
        leading_ele_pair = ak.firsts(ele_pairs[ak.argmax(diele.pt, axis=1, keepdims=True)])
        leading_diele = ak.firsts(diele[ak.argmax(diele.pt, axis=1, keepdims=True)])
        
        pho = events.Photon
        pho['isclean']=(
            ak.all(pho.metric_table(mu_loose) > 0.5, axis=2)
            & ak.all(pho.metric_table(e_loose) > 0.5, axis=2)
        )
        pho['T'] = ak.zip(
            {
                "r": pho.pt,
                "phi": pho.phi,
            },
            with_name="PolarTwoVector",
            behavior=vector.behavior,
        )
        pho['isloose']=isLoosePhoton(pho)
        pho['id_sf'] = ak.where(
            pho.isloose,
            get_pho_loose_id_sf(self._year, pho.eta, pho.pt),
            ak.ones_like(pho.pt)
        )
        pho['istight']=isTightPhoton(pho)
        pho['id_sf'] = ak.where(
            pho.istight,
            get_pho_tight_id_sf(self._year, pho.eta, pho.pt),
            pho.id_sf
        )
        pho_clean=pho[pho.isclean]
        pho_loose=pho_clean[pho_clean.isloose]
        pho_tight=pho_clean[pho_clean.istight]
        pho_ntot=ak.num(pho, axis=1)
        pho_nloose=ak.num(pho_loose, axis=1)
        pho_ntight=ak.num(pho_tight, axis=1)
        leading_pho = ak.firsts(pho_tight)

        fj = events.AK15PFPuppiJet
        #fj['pt'] = fj.subjets.sum().pt
        fj['msd_corr'] = get_msd_corr(fj)
        fj['isclean'] = (
            ak.all(fj.metric_table(mu_loose) > 1.5, axis=2)
            & ak.all(fj.metric_table(e_loose) > 1.5, axis=2)
            & ak.all(fj.metric_table(pho_loose) > 1.5, axis=2)
        )
        fj['isgood'] = isGoodAK15(fj)
        fj['T'] = ak.zip(
            {
                "r": fj.pt,
                "phi": fj.phi,
            },
            with_name="PolarTwoVector",
            behavior=vector.behavior,
        )
        probQCD=fj.particleNetAK15_QCDbb+fj.particleNetAK15_QCDcc+fj.particleNetAK15_QCDb+fj.particleNetAK15_QCDc+fj.particleNetAK15_QCDothers
        probT=fj.particleNetAK15_Tbqq+fj.particleNetAK15_Tbcq
        fj['TvsQCD'] = probT/(probT+probQCD)
        fj_good = fj[fj.isgood]
        fj_clean = fj_good[fj_good.isclean]
        fj_ntot = ak.num(fj, axis=1)
        fj_ngood = ak.num(fj_good, axis=1)
        fj_nclean = ak.num(fj_clean, axis=1)
        leading_fj = ak.firsts(fj_clean)

        j = events.Jet
        j['isgood'] = isGoodAK4(j, self._year)
        j['isHEM'] = isHEMJet(j)
        j['isclean'] = (
            ak.all(j.metric_table(mu_loose) > 0.4, axis=2)
            & ak.all(j.metric_table(e_loose) > 0.4, axis=2)
            & ak.all(j.metric_table(pho_loose) > 0.4, axis=2)
        )
        j['isiso'] = ak.all(j.metric_table(leading_fj) > 1.5, axis=2)
        j['isdcsvL'] = (j.btagDeepB>deepcsvWPs['loose'])
        j['isdflvL'] = (j.btagDeepFlavB>deepflavWPs['loose'])
        j['T'] = ak.zip(
            {
                "r": j.pt,
                "phi": j.phi,
            },
            with_name="PolarTwoVector",
            behavior=vector.behavior,
        )
        j_good = j[j.isgood]
        j_clean = j_good[j_good.isclean]
        j_iso = j_clean[j_clean.isiso]
        j_dcsvL = j_iso[j_iso.isdcsvL]
        j_dflvL = j_iso[j_iso.isdflvL]
        j_HEM = j[j.isHEM]
        j_ntot=ak.num(j, axis=1)
        j_ngood=ak.num(j_good, axis=1)
        j_nclean=ak.num(j_clean, axis=1)
        j_niso=ak.num(j_iso, axis=1)
        j_ndcsvL=ak.num(j_dcsvL, axis=1)
        j_ndflvL=ak.num(j_dflvL, axis=1)
        j_nHEM = ak.num(j_HEM, axis=1)
        leading_j = ak.firsts(j_clean)

        ###
        # Calculate recoil and transverse mass
        ###

        u = {
            'sr'    : met,
            'wecr'  : met+leading_e.T,
            'tecr'  : met+leading_e.T,
            'wmcr'  : met+leading_mu.T,
            'tmcr'  : met+leading_mu.T,
            'zecr'  : met+leading_diele.T,
            'zmcr'  : met+leading_dimu.T,
            'gcr'   : met+leading_pho.T,
        }

        mT = {
            'wecr'  : np.sqrt(2*leading_e.pt*met.pt*(1-np.cos(met.delta_phi(leading_e.T)))),
            'tecr'  : np.sqrt(2*leading_e.pt*met.pt*(1-np.cos(met.delta_phi(leading_e.T)))),
            'wmcr'  : np.sqrt(2*leading_mu.pt*met.pt*(1-np.cos(met.delta_phi(leading_mu.T)))),
            'tmcr'  : np.sqrt(2*leading_mu.pt*met.pt*(1-np.cos(met.delta_phi(leading_mu.T)))),
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
            nlo = np.ones(len(events), dtype='float')
            if('TTTo' in dataset): 
                nlo = np.sqrt(get_ttbar_weight(genTops[:,0].pt) * get_ttbar_weight(genTops[:,1].pt))
                
            gen['isW'] = (abs(gen.pdgId)==24)&gen.hasFlags(['fromHardProcess', 'isLastCopy'])
            gen['isZ'] = (abs(gen.pdgId)==23)&gen.hasFlags(['fromHardProcess', 'isLastCopy'])
            gen['isA'] = (abs(gen.pdgId)==22) & gen.hasFlags(['isPrompt', 'fromHardProcess', 'isLastCopy']) & (gen.status == 1)
            
            ###
            # Calculating gen photon dynamic isolation as in https://arxiv.org/pdf/1705.04664.pdf
            ###

            epsilon_0_dyn = 0.1
            n_dyn = 1
            gen['R_dyn'] = ak.where(
                gen.isA,
                91.1876/(gen.pt * np.sqrt(epsilon_0_dyn)),
                -999
            )
            gen['R_0_dyn'] = ak.where(
                (gen.R_dyn<1.0),
                gen.R_dyn,
                1.0
            )

            def isolation(R):
                hadrons = gen[ #Stable hadrons not in NanoAOD, using quarks/glouns instead
                    ((abs(gen.pdgId)<=5)|(abs(gen.pdgId)==21)) &
                    gen.hasFlags(['fromHardProcess', 'isFirstCopy'])
                ]
                pair_gen, pair_hadrons = ak.unzip(ak.cartesian([gen, hadrons], nested=True))
                distance = pair_hadrons.metric_table(pair_gen)
                hadronic_et = pair_hadrons[(distance <= R)].pt.sum()
                condition = ak.where(
                    ak.num(hadrons, axis=1)>0,
                    (hadronic_et<=(epsilon_0_dyn * gen.pt * np.power((1 - np.cos(R)) / (1 - np.cos(gen.R_0_dyn)), n_dyn))),
                    True
                )
                return condition

            isIsoA=gen.isA
            iterations = 5.
            for i in range(1, int(iterations) + 1):
                isIsoA=isIsoA&isolation(gen.R_0_dyn*i/iterations)
            gen['isIsoA']=isIsoA
            
            genWs = gen[gen.isW] 
            genZs = gen[gen.isZ]
            genDYs = gen[gen.isZ&(gen.mass>30)]
            genIsoAs = gen[gen.isIsoA] 
            
            nnlo_nlo = {}
            nlo_ewk = np.ones(len(events), dtype='float')
            if('G1Jet' in dataset): 
                nlo_ewk = get_nlo_ewk_weight('a', ak.firsts(genIsoAs).pt)
                for systematic in get_nnlo_nlo_weight(self._year, 'a', ak.firsts(genIsoAs).pt):
                    nnlo_nlo[systematic] = ak.where(
                        ((ak.num(genIsoAs, axis=1)>0)&(ak.firsts(genIsoAs).pt>=100)),
                        get_nnlo_nlo_weight(self._year, 'a', ak.firsts(genIsoAs).pt)[systematic],
                        np.ones(len(events), dtype='float')
                    )
            if('WJets' in dataset): 
                nlo_ewk = get_nlo_ewk_weight('w', ak.firsts(genWs).pt)
                for systematic in get_nnlo_nlo_weight(self._year, 'w', ak.firsts(genWs).pt):
                    nnlo_nlo[systematic] = ak.where(
                        ((ak.num(genWs, axis=1)>0)&(ak.firsts(genWs).pt>=100)),
                        get_nnlo_nlo_weight(self._year, 'w', ak.firsts(genWs).pt)[systematic],
                        np.ones(len(events), dtype='float')
                    )
            elif('DY' in dataset): 
                nlo_ewk = get_nlo_ewk_weight('dy', ak.firsts(genDYs).pt)
                for systematic in get_nnlo_nlo_weight(self._year, 'dy', ak.firsts(genDYs).pt):
                    nnlo_nlo[systematic] = ak.where(
                        ((ak.num(genDYs, axis=1)>0)&(ak.firsts(genDYs).pt>=100)),
                        get_nnlo_nlo_weight(self._year, 'dy', ak.firsts(genDYs).pt)[systematic],
                        np.ones(len(events), dtype='float')
                    )
            elif('Z1Jets' in dataset or 'Z2Jets' in dataset): 
                nlo_ewk = get_nlo_ewk_weight('z', ak.firsts(genZs).pt)
                for systematic in get_nnlo_nlo_weight(self._year, 'z', ak.firsts(genZs).pt):
                    nnlo_nlo[systematic] = ak.where(
                        ((ak.num(genZs, axis=1)>0)&(ak.firsts(genZs).pt>=100)),
                        get_nnlo_nlo_weight(self._year, 'z', ak.firsts(genZs).pt)[systematic],
                        np.ones(len(events), dtype='float')
                    )

            ###
            # Calculate PU weight and systematic variations
            ###
            pu = get_pu_weight(self._year, events.Pileup.nTrueInt)

            ###
            # Trigger efficiency weight
            ###

            e1sf = ak.where(
                (leading_ele_pair["0"].pt > 40),
                get_ele_trig_weight(self._year, leading_ele_pair["0"].eta+leading_ele_pair["0"].deltaEtaSC, leading_ele_pair["0"].pt),
                0.
            )
            e1sf = ak.where(
                (leading_ele_pair["1"].pt > 40),
                get_ele_trig_weight(self._year, leading_ele_pair["1"].eta+leading_ele_pair["1"].deltaEtaSC, leading_ele_pair["1"].pt),
                0.
            )
            
            if self._year == '2016':
                sf =  get_pho_trig_weight(leading_pho.pt)
            elif self._year == '2017': #Sigmoid used for 2017 and 2018, values from monojet
                sf = sigmoid(leading_pho.pt.sum(),0.335,217.91,0.065,0.996) / sigmoid(leading_pho.pt.sum(),0.244,212.34,0.050,1.000)
                sf[np.isnan(sf) | np.isinf(sf)] == 1
            elif self._year == '2018':
                sf = sigmoid(leading_pho.pt.sum(),1.022, 218.39, 0.086, 0.999) / sigmoid(leading_pho.pt.sum(), 0.301,212.83,0.062,1.000)
                sf[np.isnan(sf) | np.isinf(sf)] == 1
            
            trig = {
                'sr':   get_met_trig_weight(self._year, met.pt),
                'wmcr': get_met_trig_weight(self._year, u['wmcr'].r),
                'tmcr': get_met_trig_weight(self._year, u['tmcr'].r),
                'wecr': get_ele_trig_weight(self._year, leading_e.eta+leading_e.deltaEtaSC, leading_e.pt),
                'tecr': get_ele_trig_weight(self._year, leading_e.eta+leading_e.deltaEtaSC, leading_e.pt),
                'zmcr': get_met_trig_weight(self._year, u['zmcr'].r),
                'zecr': 1 - (1 - e1sf)*(1 - e2sf),
                'gcr': sf
            }

            ### 
            # Calculating electron and muon ID weights
            ###
            ids ={
                'sr':  np.ones(len(events), dtype='float'),
                'wmcr': leading_mu.id_sf,
                'tmcr': leading_mu.id_sf,
                'wecr': leading_e.id_sf,
                'tecr': leading_e.id_sf,
                'zmcr': leading_mu_pair["0"].id_sf*leading_mu_pair["1"].id_sf,
                'zecr': leading_ele_pair["0"].id_sf*leading_ele_pair["1"].id_sf,
                'gcr':  leading_pho.id_sf,
            }
           
            ###
            # Reconstruction weights for electrons
            ###                                       
            reco = {
                'sr': np.ones(len(events), dtype='float'),
                'wmcr': np.ones(len(events), dtype='float'),
                'tmcr': np.ones(len(events), dtype='float'),
                'wecr': leading_e.reco_sf,
                'tecr': leading_e.reco_sf,
                'zmcr': np.ones(len(events), dtype='float'),
                'zecr': leading_e.reco_sf,
                'gcr':  np.ones(len(events), dtype='float'),
            }

            ###
            # Isolation weights for muons
            ###
            isolation = {
                'sr': np.ones(len(events), dtype='float'),
                'wmcr': leading_mu.iso_sf,
                'tmcr': leading_mu.iso_sf,
                'wecr': np.ones(len(events), dtype='float'),
                'tecr': np.ones(len(events), dtype='float'),
                'zmcr': leading_mu_pair["0"].iso_sf*leading_mu_pair["1"].iso_sf,
                'zecr': np.ones(len(events), dtype='float'),
                'gcr':  np.ones(len(events), dtype='float'),
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
            btagSFlight_uncorrelatedDown  = get_btag_weight('deepflav',self._year,'loose').btag_weight(
                j_iso.pt,
                j_iso.eta,
                j_iso.hadronFlavour,
                j_iso.isdflvL
            )

            if hasattr(events, "L1PreFiringWeight"): 
                weights.add('prefiring', events.L1PreFiringWeight.Nom, events.L1PreFiringWeight.Up, events.L1PreFiringWeight.Dn)
            weights.add('genw',events.genWeight)
            weights.add('nlo_ewk',nlo_ewk)
            if 'cen' in nnlo_nlo:
                #weights.add('nnlo_nlo',nnlo_nlo['cen'])
                weights.add('qcd1',np.ones(len(events), dtype='float'), nnlo_nlo['qcd1up']/nnlo_nlo['cen'], nnlo_nlo['qcd1do']/nnlo_nlo['cen'])
                weights.add('qcd2',np.ones(len(events), dtype='float'), nnlo_nlo['qcd2up']/nnlo_nlo['cen'], nnlo_nlo['qcd2do']/nnlo_nlo['cen'])
                weights.add('qcd3',np.ones(len(events), dtype='float'), nnlo_nlo['qcd3up']/nnlo_nlo['cen'], nnlo_nlo['qcd3do']/nnlo_nlo['cen'])
                weights.add('ew1',np.ones(len(events), dtype='float'), nnlo_nlo['ew1up']/nnlo_nlo['cen'], nnlo_nlo['ew1do']/nnlo_nlo['cen'])
                weights.add('ew2G',np.ones(len(events), dtype='float'), nnlo_nlo['ew2Gup']/nnlo_nlo['cen'], nnlo_nlo['ew2Gdo']/nnlo_nlo['cen'])
                weights.add('ew3G',np.ones(len(events), dtype='float'), nnlo_nlo['ew3Gup']/nnlo_nlo['cen'], nnlo_nlo['ew3Gdo']/nnlo_nlo['cen'])
                weights.add('ew2W',np.ones(len(events), dtype='float'), nnlo_nlo['ew2Wup']/nnlo_nlo['cen'], nnlo_nlo['ew2Wdo']/nnlo_nlo['cen'])
                weights.add('ew3W',np.ones(len(events), dtype='float'), nnlo_nlo['ew3Wup']/nnlo_nlo['cen'], nnlo_nlo['ew3Wdo']/nnlo_nlo['cen'])
                weights.add('ew2Z',np.ones(len(events), dtype='float'), nnlo_nlo['ew2Zup']/nnlo_nlo['cen'], nnlo_nlo['ew2Zdo']/nnlo_nlo['cen'])
                weights.add('ew3Z',np.ones(len(events), dtype='float'), nnlo_nlo['ew3Zup']/nnlo_nlo['cen'], nnlo_nlo['ew3Zdo']/nnlo_nlo['cen'])
                weights.add('mix',np.ones(len(events), dtype='float'), nnlo_nlo['mixup']/nnlo_nlo['cen'], nnlo_nlo['mixdo']/nnlo_nlo['cen'])
                weights.add('muF',np.ones(len(events), dtype='float'), nnlo_nlo['muFup']/nnlo_nlo['cen'], nnlo_nlo['muFdo']/nnlo_nlo['cen'])
                weights.add('muR',np.ones(len(events), dtype='float'), nnlo_nlo['muRup']/nnlo_nlo['cen'], nnlo_nlo['muRdo']/nnlo_nlo['cen'])
            weights.add('pileup',pu)
            weights.add('trig', trig[region])
            weights.add('ids', ids[region])
            weights.add('reco', reco[region])
            weights.add('isolation', isolation[region])
            weights.add('btagSF',btagSF)
            weights.add('btagSFbc_correlated',np.ones(len(events), dtype='float'), btagSFbc_correlatedUp/btagSF, btagSFbc_correlatedDown/btagSF)
            weights.add('btagSFbc_uncorrelated',np.ones(len(events), dtype='float'), btagSFbc_uncorrelatedUp/btagSF, btagSFbc_uncorrelatedDown/btagSF)
            weights.add('btagSFlight_correlated',np.ones(len(events), dtype='float'), btagSFlight_correlatedUp/btagSF, btagSFlight_correlatedDown/btagSF)
            weights.add('btagSFlight_uncorrelated',np.ones(len(events), dtype='float'), btagSFlight_uncorrelatedUp/btagSF, btagSFlight_uncorrelatedDown/btagSF)
            
        ###
        # Selections
        ###

        lumimask = np.ones(len(events), dtype='bool')
        if isData:
            lumimask = AnalysisProcessor.lumiMasks[self._year](events.run, events.luminosityBlock)
        selection.add('lumimask', lumimask)

        met_filters =  np.ones(len(events), dtype='bool')
        #if isData: met_filters = met_filters & events.Flag['eeBadScFilter']#this filter is recommended for data only
        for flag in AnalysisProcessor.met_filters[self._year]:
            met_filters = met_filters & events.Flag[flag]
        selection.add('met_filters',met_filters)

        triggers = np.zeros(len(events), dtype='bool')
        for path in self._met_triggers[self._year]:
            if not hasattr(events.HLT, path): continue
            triggers = triggers | events.HLT[path]
        selection.add('met_triggers', triggers)

        triggers = np.zeros(len(events), dtype='bool')
        for path in self._singleelectron_triggers[self._year]:
            if not hasattr(events.HLT, path): continue
            triggers = triggers | events.HLT[path]
        selection.add('singleelectron_triggers', triggers)

        triggers = np.zeros(len(events), dtype='bool')
        for path in self._singlephoton_triggers[self._year]:
            if not hasattr(events.HLT, path): continue
            triggers = triggers | events.HLT[path]
        selection.add('single_photon_triggers', ak.to_numpy(triggers))

        triggers = np.zeros(len(events), dtype='bool')
        for path in self._singlemuon_triggers[self._year]:
            if path not in events.HLT.fields:
                continue
            triggers = triggers | events.HLT[path]
        selection.add('single_muon_triggers', ak.to_numpy(triggers))

        noHEMj = np.ones(len(events), dtype='bool')
        if self._year=='2018':
            noHEMj = (j_nHEM==0)

        noHEMmet = np.ones(len(events), dtype='bool')
        if self._year=='2018':
            noHEMmet = (met.pt>470)|(met.phi>-0.62)|(met.phi<-1.62)

        if ('WJetsToLNu' in dataset) & ('Pt' in dataset):
            remove_overlap = (gen[gen.hasFlags(['fromHardProcess', 'isFirstCopy', 'isPrompt']) & ((abs(gen.pdgId) == 24))].pt >120) ## W
            selection.add("exclude_wjets_greater_400", ak.to_numpy(ak.all(remove_overlap, axis=1)))
        else:
            selection.add("exclude_wjets_greater_400", np.full(len(events), True))

        if ('WJetsToLNu' in dataset) & (not ('Pt' in dataset)):
            remove_overlap = (gen[gen.hasFlags(['fromHardProcess', 'isFirstCopy', 'isPrompt']) & ((abs(gen.pdgId) == 24))].pt <= 120) ## w
            selection.add("exclude_wjets_less_400", ak.to_numpy(ak.all(remove_overlap, axis=1)))
        else:
            selection.add("exclude_wjets_less_400", np.full(len(events), True))

        selection.add('iszeroL', (e_nloose==0)&(mu_nloose==0)&(pho_nloose==0))
        selection.add('isoneM', (e_nloose==0)&(mu_ntight==1)&(mu_nloose==1)&(pho_nloose==0))
        selection.add('isoneE', (e_ntight==1)&(e_nloose==1)&(mu_nloose==0)&(pho_nloose==0))
        selection.add('isoneG', (e_nloose==0)&(mu_nloose==0)&(pho_nloose==1)&(pho_ntight==1))
        selection.add('istwoM', (e_nloose==0)&(mu_nloose==2)&(pho_nloose==0))
        selection.add('istwoE', (e_nloose==2)&(mu_nloose==0)&(pho_nloose==0))
        selection.add('one_ak4', (j_nclean>0))
        selection.add('one_ak15', (fj_nclean>0))
        selection.add('leading_fj250', (leading_fj.pt>250))

        selection.add('noextrab', (j_ndflvL==0))
        selection.add('extrab', (j_ndflvL>0))
        selection.add('oneb', (j_ndflvL==1))
        selection.add('noHEMj', noHEMj)
        selection.add('noHEMmet', noHEMmet)

        #selection.add('msd40',(leading_fj.msd_corr>40)) ## 2024.5.9

        selection.add('met120',(met.pt<120))
        selection.add('met150',(met.pt>150))
        selection.add('diele60',(leading_diele.mass>60))
        selection.add('diele120',(leading_diele.mass<120))
        selection.add('dimu60',(leading_dimu.mass>60))
        selection.add('dimu120',(leading_dimu.mass<120))
        selection.add('leading_ele40',(leading_diele_pair["0"].pt>40)|(leading_diele_pair["1"].pt>40))

        regions = {
            'sr': [
                    'met_filters', 'met_triggers',
                    'exclude_wjets_greater_400', 'exclude_wjets_less_400',
                    'one_ak15',
                    'leading_fj250',
                    'iszeroL',
                    'noextrab',
                    'noHEMj', 'noHEMmet',
            ],
            'wmcr': [
                    'met_filters', 'met_triggers',
                    'exclude_wjets_greater_400', 'exclude_wjets_less_400',
                    'one_ak15',
                    'leading_fj250',
                    'isoneM',
                    'met150',
                    'noextrab',
                    'noHEMj',
            ],
            'wecr': [
                    'met_filters', 'singleelectron_triggers',
                    'exclude_wjets_greater_400', 'exclude_wjets_less_400',
                    'one_ak15',
                    'leading_fj250',
                    'isoneE',
                    'met150',
                    'noextrab',
                    'noHEMj',
            ],
            'tmcr': [
                    'met_filters', 'met_triggers',
                    'exclude_wjets_greater_400', 'exclude_wjets_less_400',
                    'one_ak15',
                    'leading_fj250',
                    'isoneM',
                    'met150',
                    'extrab',
                    'noHEMj',
            ],
            'tecr': [
                    'met_filters', 'singleelectron_triggers',
                    'exclude_wjets_greater_400', 'exclude_wjets_less_400',
                    'one_ak15',
                    'leading_fj250',
                    'isoneE',
                    'met150',
                    'extrab',
                    'noHEMj',
            ],
            'zmcr': [
                    'met_filters', 'met_triggers',
                    'exclude_wjets_greater_400', 'exclude_wjets_less_400',
                    'one_ak15',
                    'leading_fj250',
                    'istwoM',
                    'met120',
                    'dimu60', 'dimu120',
                    'noHEMj',
            ],
            'zecr': [
                    'met_filters', 'singleelectron_triggers',
                    'exclude_wjets_greater_400', 'exclude_wjets_less_400',
                    'one_ak15',
                    'leading_fj250',
                    'istwoE',
                    'met120',
                    'diele60', 'diele120',
                    'leading_ele40',
                    'noHEMj',
            ],
            'gcr': [
                    'met_filters', 'single_photon_triggers',
                    'exclude_wjets_greater_400', 'exclude_wjets_less_400',
                    'one_ak15',
                    'leading_fj250',
                    'isoneG',
                    'noextrab',
                    'noHEMj',
            ],
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
            
            if systematic is None:
                variables = {
                    'mindphirecoil':          ak.min(abs(u[region].delta_phi(j_clean.T)), axis=1,mask_identity=False),
                    #'minDphirecoil':          ak.min(abs(u[region].delta_phi(leading_fj.T)), axis=1,mask_identity=False),
                    'minDphirecoil':           abs(u[region].delta_phi(leading_fj.T)),
                    'ut' :                   u[region].r,
                    'uphi' :                 u[region].phi,
                    'met':                    met.pt,
                    'metphi':                 met.phi,
                    'j1pt':                   leading_j.pt,
                    'j1eta':                  leading_j.eta,
                    'j1phi':                  leading_j.phi,
                    'fj1pt':                  leading_fj.pt,
                    'fj1eta':                 leading_fj.eta,
                    'fj1phi':                 leading_fj.phi,
                    'njets':                  j_nclean,
                    'ndflvL':                 j_ndflvL,
                    'nfjclean':               fj_nclean,
                }
                if region in mT:
                    variables['mT']           = mT[region]
                if 'e' in region:
                    if 'z' in region:
                        variables['l1pt']      = leading_diele.pt
                        variables['l1phi']     = leading_diele.phi
                        variables['l1eta']     = leading_diele.eta
                    else:
                        variables['l1pt']      = leading_e.pt
                        variables['l1phi']     = leading_e.phi
                        variables['l1eta']     = leading_e.eta
                if 'm' in region:
                    if 'z' in region:
                        variables['l1pt']      = leading_dimu.pt
                        variables['l1phi']     = leading_dimu.phi
                        variables['l1eta']     = leading_dimu.eta
                    else:
                        variables['l1pt']      = leading_mu.pt
                        variables['l1phi']     = leading_mu.phi
                        variables['l1eta']     = leading_mu.eta
                if 'g' in region:
                    variables['l1pt']      = leading_pho.pt
                    variables['l1phi']     = leading_pho.phi
                    variables['l1eta']     = leading_pho.eta
                for variable in output:
                    if variable not in variables:
                        continue
                    normalized_variable = {variable: normalize(variables[variable],cut)}
                    output[variable].fill(
                        region=region,
                        TvsQCD=normalize(leading_fj.TvsQCD,cut),
                        **normalized_variable,
                        weight=weight,
                    )
                output['TvsQCD'].fill(
                      region=region,
                      TvsQCD=normalize(leading_fj.TvsQCD, cut),
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

            if 'qcd' not in region:
                selection.add('recoil_'+region, (u[region].r>350))
                selection.add('mindphi_'+region, (ak.min(abs(u['sr'].delta_phi(j_clean.T)), axis=1, mask_identity=False) > 0.5))
                selection.add('minDphi_'+region, (abs(u[region].delta_phi(leading_fj.T)) > 1.5))
                regions[region].insert(7, 'recoil_'+region)
                regions[region].insert(9, 'mindphi_'+region)
                regions[region].insert(10, 'minDphi_'+region)
                if region in mT:
                    selection.add('mT_'+region, (mT[region]<150))

            for systematic in systematics:
                if isData and systematic is not None:
                    continue
                fill(region, systematic)
            vcut=np.zeros(len(events), dtype=np.int)
            output['cutflow'].fill(region=region,cutname='Initial', cutflow=vcut, weight=weights.weight())
            cuts = regions[region]
            allcuts = set()
            for i, icut in enumerate(cuts):
                allcuts.add(icut)
                jcut = selection.all(*allcuts)
                vcut = (i+1)*jcut
                output['cutflow'].fill(region=region,cutname=icut, cutflow=vcut, weight=weights.weight()*jcut)


        scale = 1
        if self._xsec[dataset]!= -1: 
            scale = self._lumi*self._xsec[dataset]

        for key in output:
            if key=='sumw': 
                continue
            output[key] *= scale
                
        return output

    def postprocess(self, accumulator):

        return accumulator

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('-y', '--year', help='year', dest='year')
    parser.add_option('-m', '--metadata', help='metadata', dest='metadata')
    parser.add_option('-n', '--name', help='name', dest='name')
    (options, args) = parser.parse_args()


    with gzip.open('metadata/'+options.metadata+'.json.gz') as fin:
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

    save(processor_instance, 'data/hadmonotop'+options.name+'.processor')

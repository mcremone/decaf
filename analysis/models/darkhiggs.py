from __future__ import print_function, division
from optparse import OptionParser
from collections import defaultdict, OrderedDict
import concurrent.futures
import sys
import os
import rhalphalib as rl
from libs.myrhalphalib import TransferFactorSample
import numpy as np
import pickle
import gzip
import json
import uncertainties.unumpy as unumpy 
from coffea import hist, processor
from coffea.util import load, save
from scipy import stats
import ROOT

rl.util.install_roofit_helpers()
rl.ParametericSample.PreferRooParametricHist = False

mass_binning = [40., 50., 60., 70., 80., 90., 100., 120., 150., 180., 240., 300.,]
#mass_binning = [40., 50., 60., 70., 80., 90., 100., 120., 160., 200., 300.,] 
recoil_binning = [250., 310., 370., 470., 590., 3000.]
category_map = {"pass": 1, "fail": 0}

def template(dictionary, process, systematic, recoil, region, category, mass, min_value=1e-5, read_sumw2=False):
    histogram = dictionary[region].integrate("process", process)
    nominal, sumw2 = histogram.integrate("systematic", "nominal").values(sumw2=True)[()]
    nominal=nominal[recoil, :, category_map[category]]
    sumw2=sumw2[recoil, :, category_map[category]]
    zerobins = nominal <= 0.
    output = nominal
    if "data" not in systematic:
        output[zerobins] = min_value
        sumw2[zerobins] = 0.
    if "nominal" not in systematic and "data" not in systematic:
        output = histogram.integrate("systematic", systematic).values()[()][recoil, :, category_map[category]]
        output[zerobins] = 1.
        output[~zerobins] /= nominal[~zerobins]
        output[~zerobins] = np.maximum(output[~zerobins], 1e-5)
        output[np.isnan(output)] = 1.
    binning = (
        dictionary[region]
        .integrate("process", process)
        .integrate("systematic", systematic)
        .axis("fjmass")
        .edges()
    )
    if read_sumw2:
        return (output, binning, 'fjmass'+mass, sumw2)
    return (output, binning, 'fjmass'+mass)

def remap_histograms(hists):
    data_hists = {}
    bkg_hists = {}
    signal_hists = {}
    fakedata_map = OrderedDict()
  
    process = hist.Cat("process", "Process", sorting="placement")
    cats = ("process",)
    sig_map = OrderedDict()
    bkg_map = OrderedDict()
    data_map = OrderedDict()
    bkg_map["Hbb"] = ("Hbb*",)
    bkg_map["DY+jets"] = ("DY+*",)
    bkg_map["VV"] = (["WW", "WZ", "ZZ"],)
    bkg_map["ST"] = ("ST*",)
    bkg_map["TT"] = ("TT*",)
    bkg_map["W+jets"] = ("W+*",)
    bkg_map["W+HF"] = ("W+HF",)
    bkg_map["W+LF"] = ("W+LF",)
    bkg_map["Z+jets"] = ("Z+*",)
    bkg_map["Z+HF"] = ("Z+HF",)
    bkg_map["Z+LF"] = ("Z+LF",)
    bkg_map["G+jets"] = ("G+*",)
    bkg_map["QCD"] = ("QCD*",)
    data_map["MET"] = ("MET",)
    data_map["SingleElectron"] = ("SingleElectron",)
    data_map["SinglePhoton"] = ("SinglePhoton",)
    data_map["EGamma"] = ("EGamma",)
    
    for signal in hists['sig']['template'].identifiers('process'):
        if 'mhs' not in str(signal): continue
        sig_map[str(signal)] = (str(signal),)  ## signals
        
    fakedata_list = []
    for bkg in hists['bkg']['template'].identifiers('process'):
        fakedata_list.append(str(bkg))
    fakedata_map['FakeData'] = (fakedata_list,)
    
    for key in hists["data"].keys():
        bkg_hists[key] = hists["bkg"][key].group(cats, process, bkg_map)
        signal_hists[key] = hists["sig"][key].group(cats, process, sig_map)
        data_hists[key] = hists["data"][key].group(cats, process, data_map)
        data_hists[key] += hists["bkg"][key].group(cats, process, fakedata_map)
    
    bkg_hists["template"] = bkg_hists["template"].rebin(
        "fjmass", hist.Bin("fjmass", "Mass", mass_binning)
    )
    signal_hists["template"] = signal_hists["template"].rebin(
        "fjmass", hist.Bin("fjmass", "Mass", mass_binning)
    )
    data_hists["template"] = data_hists["template"].rebin(
        "fjmass", hist.Bin("fjmass", "Mass", mass_binning)
    )

    bkg_hists["template"] = bkg_hists["template"].rebin(
        "recoil", hist.Bin("recoil", "Recoil", recoil_binning)
    )
    signal_hists["template"] = signal_hists["template"].rebin(
        "recoil", hist.Bin("recoil", "Recoil", recoil_binning)
    )
    data_hists["template"] = data_hists["template"].rebin(
        "recoil", hist.Bin("recoil", "Recoil", recoil_binning)
    )

    hists = {"bkg": bkg_hists, "sig": signal_hists, "data": data_hists}

    return hists

def makeTF(num, den):

    tf = num.getExpectation()/den.getExpectation()
    num=unumpy.uarray(( num._nominal, np.minimum(np.sqrt(num._sumw2),num._nominal) ))  
    den=unumpy.uarray(( den._nominal, np.minimum(np.sqrt(den._sumw2),den._nominal) ))  
    ratio=num/den
    unc = unumpy.std_devs(ratio)/unumpy.nominal_values(ratio)
    
    return tf, unc


def addBBLiteSyst(channel, epsilon=1e-5, effect_threshold=0.01, threshold=0, include_signal=0, channel_name=None):
    """
    Barlow-Beeston-lite method i.e. single stats parameter for all processes per bin.
    Same general algorithm as described in
    https://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/part2/bin-wise-stats/
    but *without the analytic minimisation*.
    `include_signal` only refers to whether signal stats are included in the *decision* to use bb-lite or not.
    """
    if not len(channel._samples):
        raise RuntimeError("Channel %r has no samples for which to run autoMCStats" % (channel))
    
    name = channel._name if channel_name is None else channel_name
    
    first_sample = channel._samples[list(channel._samples.keys())[0]]

    ntot_bb, etot2_bb = np.zeros_like(first_sample._nominal), np.zeros_like(first_sample._sumw2)
    ntot, etot2 = np.zeros_like(first_sample._nominal), np.zeros_like(first_sample._sumw2)
    
    for sample in channel._samples.values():
        if isinstance(sample, TransferFactorSample):
            ntot += sample._nominal_values
            etot2 += (sample._stat_unc*sample._nominal_values)**2
        elif isinstance(sample, rl.TemplateSample):
            ntot += sample._nominal
            etot2 += sample._sumw2
            if not include_signal and sample._sampletype == rl.Sample.SIGNAL:
                continue
            ntot_bb += sample._nominal
            etot2_bb += sample._sumw2
        else:
            continue

    for sample in channel._samples.values():
        if not isinstance(sample, TransferFactorSample):
            continue
        channel._samples[sample.name] = TransferFactorSample(sample.name, 
                                                             rl.Sample.BACKGROUND, 
                                                             sample.transferfactor, 
                                                             sample.dependentsample, 
                                                             nominal_values=sample.nominal_values, 
                                                             stat_unc=np.sqrt(etot2)/ntot, 
                                                             channel_name=name)
        
    for i in range(first_sample.observable.nbins):
        if etot2[i] <= 0.0:
            continue
        elif etot2_bb[i] <= 0:
            # this means there is signal but no background, so create stats unc. for signal only
            for sample in channel._samples.values():
                if sample._sampletype == rl.Sample.SIGNAL:
                    sample_name = None if channel_name is None else channel_name + "_" + sample._name[sample._name.find("_") + 1 :]
                    sample.autoMCStats(epsilon=epsilon, sample_name=sample_name, bini=i)
    
            continue

        neff_bb = ntot_bb[i]**2 / etot2_bb[i]
        if neff_bb <= threshold:
            for sample in channel._samples.values():
                sample_name = None if channel_name is None else channel_name + "_" + sample._name[sample._name.find("_") + 1 :]
                sample.autoMCStats(epsilon=epsilon, sample_name=sample_name, bini=i)
        else:
            effect_up = np.ones_like(first_sample._nominal)
            effect_down = np.ones_like(first_sample._nominal)

            effect = np.sqrt(etot2[i])/ntot[i]
            if effect < effect_threshold:
                continue
            effect_up[i] = 1.0 + min(1.0, effect)
            effect_down[i] = max(epsilon, 1.0 - min(1.0, effect))
    
            param = rl.NuisanceParameter(name + "_mcstat_bin%i" % i, combinePrior="shape")
    
            for sample in channel._samples.values():
                if not isinstance(sample, rl.TemplateSample):
                    continue
                if sample._nominal[i] <= 1e-5:
                    continue
                #print(sample._name, i, effect, name + "_mcstat_bin%i" % i)
                sample.setParamEffect(param, effect_up, effect_down)
                
def addBtagSyst(dictionary, recoil, process, region, templ, category, mass):
    btagUp = template(dictionary, process, "btagSFbc_correlatedUp", recoil, region, category, mass)[0]
    btagDown = template(dictionary, process, "btagSFbc_correlatedDown", recoil, region, category, mass)[0]
    templ.setParamEffect(btagSFbc_correlated, btagUp, btagDown)
    btagUp = template(dictionary, process, "btagSFbc_uncorrelatedUp", recoil, region, category, mass)[0]
    btagDown = template(dictionary, process, "btagSFbc_uncorrelatedDown", recoil, region, category, mass)[0]
    templ.setParamEffect(btagSFbc, btagUp, btagDown)
    btagUp = template(dictionary, process, "btagSFlight_correlatedUp", recoil, region, category, mass)[0]
    btagDown = template(dictionary, process, "btagSFlight_correlatedDown", recoil, region, category, mass)[0]
    templ.setParamEffect(btagSFlight_correlated, btagUp, btagDown)
    btagUp = template(dictionary, process, "btagSFlight_uncorrelatedUp", recoil, region, category, mass)[0]
    btagDown = template(dictionary, process, "btagSFlight_uncorrelatedDown", recoil, region, category, mass)[0]
    templ.setParamEffect(btagSFlight, btagUp, btagDown)

#def addDoubleBtagSyst(dictionary, recoil, process, region, templ, category, mass):
#    doublebtagUp = template(dictionary, process, 'doublebtagUp', recoil, region, category, mass)[0]
#    doublebtagDown = template(dictionary, process, 'doublebtagDown', recoil, region, category, mass)[0]
#    templ.setParamEffect(doublebtag, doublebtagUp, doublebtagDown)

def addDoubleBtagSyst(dictionary, process, region, templ, category): 
    def addSyst(dictionary, process, region, templ, category, syst, string):
        histogram = dictionary[region].integrate("process", process)
        nominal=histogram.integrate("systematic", "nominal").sum('recoil','fjmass').values()[()][category_map[category]]
        up=histogram.integrate("systematic", string+"Up").sum('recoil','fjmass').values()[()][category_map[category]]
        systUp = up / nominal
        systUp = np.nan_to_num(systUp, nan=1.)
        print(templ._name,systUp)
        templ.setParamEffect(syst, systUp)
    addSyst(dictionary, process, region, templ, category, doublebtag, "doublebtag")

#def addDoubleBtagSyst(dictionary, recoil, process, region, templ, category): 
#    def addSyst(dictionary, recoil, process, region, templ, category, syst, string):
#        histogram = dictionary[region].integrate("process", process)
#        nominal=histogram.integrate("systematic", "nominal").values()[()][recoil, :, category_map[category]]
#        up=histogram.integrate("systematic", string+"Up").values()[()][recoil, :, category_map[category]]
#        down=histogram.integrate("systematic",string+"Down").values()[()][recoil, :, category_map[category]]
#        systUp = np.array( up.sum() / nominal.sum() )
#        systUp[np.isnan(systUp)] = 1.
#        systUp = systUp.sum()
#        templ.setParamEffect(syst, systUp)
#    addSyst(dictionary, recoil, process, region, templ, category, doublebtag, "doublebtag")

def addVJetsSyst(dictionary, recoil, process, region, templ, category):
    def addSyst(dictionary, recoil, process, region, templ, category, syst, string):
        histogram = dictionary[region].integrate("process", process)
        nominal=histogram.integrate("systematic", "nominal").values()[()][recoil, :, category_map[category]]
        up=histogram.integrate("systematic", string+"Up").values()[()][recoil, :, category_map[category]]
        down=histogram.integrate("systematic",string+"Down").values()[()][recoil, :, category_map[category]]
        systUp = np.array( up.sum() / nominal.sum() )
        systUp[np.isnan(systUp)] = 1.
        systUp = systUp.sum()
        templ.setParamEffect(syst, systUp)
    addSyst(dictionary, recoil, process, region, templ, category, ew1, "ew1")
    addSyst(dictionary, recoil, process, region, templ, category, ew2W, "ew2W")
    addSyst(dictionary, recoil, process, region, templ, category, ew2Z, "ew2Z")
    addSyst(dictionary, recoil, process, region, templ, category, ew3W, "ew3W")
    addSyst(dictionary, recoil, process, region, templ, category, ew3Z, "ew3Z")
    addSyst(dictionary, recoil, process, region, templ, category, mix, "mix")
    addSyst(dictionary, recoil, process, region, templ, category, qcd1, "qcd1")
    addSyst(dictionary, recoil, process, region, templ, category, qcd2, "qcd2")
    addSyst(dictionary, recoil, process, region, templ, category, qcd3, "qcd3")
    #addSyst(dictionary, recoil, process, region, templ, category, muF, "muF")
    #addSyst(dictionary, recoil, process, region, templ, category, muR, "muR")

def addLumiSyst(templ, year):
    vlumi={
        '2016': 1.01,
        '2017': 1.02,
        '2018': 1.015,
    }
    vlumi_corr={
        '2016': 1.006,
        '2017': 1.009,
        '2018': 1.02,
    }
    vlumi_1718={
        '2017': 1.006,
        '2018': 1.002,
    }
    templ.setParamEffect(lumi, vlumi[year])
    templ.setParamEffect(lumi_corr, vlumi_corr[year])
    if '2016' not in year: templ.setParamEffect(lumi_1718, vlumi_1718[year])

def addPileupSyst(templ):
    templ.setParamEffect(pu, 1.01)

def addMETSyst(templ):
    templ.setParamEffect(met, 1.05)

def addJESSyst(templ):
    templ.setParamEffect(jes, 1.04)

def addPrefiringSyst(templ, year):
    if '2018' not in year: templ.setParamEffect(prefiring, 1.01)

def addSingleTopNormSyst(templ):
    templ.setParamEffect(st_norm, 1.1)

def addTTbarNormSyst(templ):
    templ.setParamEffect(ttMC_norm, 1.1)

def addHbbNormSyst(templ):
    templ.setParamEffect(hbb_norm, 1.2)

def addDibosonNormSyst(templ):
    templ.setParamEffect(vv_norm, 1.2)

def addDrellYanNormSyst(templ):
    templ.setParamEffect(zjetsMC_norm, 1.2)

def addWjetsNormSyst(templ):
    templ.setParamEffect(wjetsMC_norm, 1.4)

def addMETTrigSyst(templ, year):
    if '2016' in year:
        templ.setParamEffect(trig_met, 1.02)
    else:
        templ.setParamEffect(trig_met, 1.01)

def addEleIDSyst(templ, year):
    if '2016' in year:
        templ.setParamEffect(id_e, 1.02)
    else:
        templ.setParamEffect(id_e, 1.03)


    
def add_crs(year, mass, recoil, model, category):

    model_id = model.name

    ###
    ###
    # Single muon W control region
    ###
    ###

    ch_name = "wmcr" + model_id
    wmcr = rl.Channel(ch_name)
    model.addChannel(wmcr)

    ###
    # Add data distribution to the channel
    ###

    dataTemplate = template(data, "MET", "data", recoil, "wmcr", category, mass)
    wmcr.setObservation(dataTemplate)

    
    ###
    # Other MC-driven processes
    ###

    if iswjetsMC:
        wmcr_wjetsTemplate = template(background, "W+jets", "nominal", recoil, "wmcr", category, mass, read_sumw2=True)
        wmcr_wjets = rl.TemplateSample("wmcr" + model_id + "_wjetsMC", rl.Sample.BACKGROUND, wmcr_wjetsTemplate)
        addLumiSyst(wmcr_wjets, year)
        addPileupSyst(wmcr_wjets)
        addPrefiringSyst(wmcr_wjets, year)
        addJESSyst(wmcr_wjets)
        addMETTrigSyst(wmcr_wjets, year)
        wmcr_wjets.setParamEffect(veto_tau, nveto_tau)
        addWjetsNormSyst(wmcr_wjets)
        wmcr_wjets.setParamEffect(id_mu, nlepton)
        wmcr_wjets.setParamEffect(iso_mu, nlepton)
        addBtagSyst(background, recoil, "W+jets", "wmcr", wmcr_wjets, category, mass)
        addVJetsSyst(background, recoil, "W+jets", "wmcr", wmcr_wjets, category)
        wmcr.addSample(wmcr_wjets)

    if isttMC: 
        wmcr_ttTemplate = template(background, "TT", "nominal", recoil, "wmcr", category, mass, read_sumw2=True)
        wmcr_tt = rl.TemplateSample( "wmcr" + model_id + "_ttMC", rl.Sample.BACKGROUND, wmcr_ttTemplate)
        addLumiSyst(wmcr_tt, year)
        addPileupSyst(wmcr_tt)
        addPrefiringSyst(wmcr_tt, year)
        addJESSyst(wmcr_tt)
        addMETTrigSyst(wmcr_tt, year)
        wmcr_tt.setParamEffect(veto_tau, nveto_tau)
        wmcr_tt.setParamEffect(id_mu, nlepton)
        wmcr_tt.setParamEffect(iso_mu, nlepton)
        addTTbarNormSyst(wmcr_tt)
        addBtagSyst(background, recoil, "TT", "wmcr", wmcr_tt, category, mass)
        wmcr.addSample(wmcr_tt)
                    
    wmcr_stTemplate = template(background, "ST", "nominal", recoil, "wmcr", category, mass, read_sumw2=True)
    wmcr_st = rl.TemplateSample(ch_name + "_stMC", rl.Sample.BACKGROUND, wmcr_stTemplate)
    addLumiSyst(wmcr_st, year)
    addPileupSyst(wmcr_st)
    addPrefiringSyst(wmcr_st, year)
    addJESSyst(wmcr_st)
    addMETTrigSyst(wmcr_st, year)
    wmcr_st.setParamEffect(veto_tau, nveto_tau)
    addSingleTopNormSyst(wmcr_st)
    wmcr_st.setParamEffect(id_mu, nlepton)
    wmcr_st.setParamEffect(iso_mu, nlepton)
    addBtagSyst(background, recoilbin, "ST", "wmcr", wmcr_st, category, mass)
    wmcr.addSample(wmcr_st)

    wmcr_dyjetsTemplate = template(background, "DY+jets", "nominal", recoil, "wmcr", category, mass, read_sumw2=True)
    wmcr_dyjets = rl.TemplateSample(ch_name + "_dyjetsMC", rl.Sample.BACKGROUND, wmcr_dyjetsTemplate)
    addLumiSyst(wmcr_dyjets, year)
    addPileupSyst(wmcr_dyjets)
    addPrefiringSyst(wmcr_dyjets, year)
    addJESSyst(wmcr_dyjets)
    addMETTrigSyst(wmcr_dyjets, year)
    wmcr_dyjets.setParamEffect(veto_tau, nveto_tau)
    addDrellYanNormSyst(wmcr_dyjets)
    wmcr_dyjets.setParamEffect(id_mu, nlepton)
    wmcr_dyjets.setParamEffect(iso_mu, nlepton)
    addBtagSyst(background, recoilbin, "DY+jets", "wmcr", wmcr_dyjets, category, mass)
    addVJetsSyst(background, recoil, "DY+jets", "wmcr", wmcr_dyjets, category)
    wmcr.addSample(wmcr_dyjets)

    wmcr_vvTemplate = template(background, "VV", "nominal", recoil, "wmcr", category, mass, read_sumw2=True)
    wmcr_vv = rl.TemplateSample(ch_name + "_vvMC", rl.Sample.BACKGROUND, wmcr_vvTemplate)
    addLumiSyst(wmcr_vv, year)
    addPileupSyst(wmcr_vv)
    addPrefiringSyst(wmcr_vv, year)
    addJESSyst(wmcr_vv)
    addMETTrigSyst(wmcr_vv, year)
    wmcr_vv.setParamEffect(veto_tau, nveto_tau)
    addDibosonNormSyst(wmcr_vv)
    wmcr_vv.setParamEffect(id_mu, nlepton)
    wmcr_vv.setParamEffect(iso_mu, nlepton)
    addBtagSyst(background, recoilbin, "VV", "wmcr", wmcr_vv, category, mass)
    wmcr.addSample(wmcr_vv)

    wmcr_hbbTemplate = template(background, "Hbb", "nominal", recoil, "wmcr", category, mass, read_sumw2=True)
    wmcr_hbb = rl.TemplateSample(ch_name + "_hbbMC", rl.Sample.BACKGROUND, wmcr_hbbTemplate)
    addLumiSyst(wmcr_hbb, year)
    addPileupSyst(wmcr_hbb)
    addPrefiringSyst(wmcr_hbb, year)
    addJESSyst(wmcr_hbb)
    addMETTrigSyst(wmcr_hbb, year)
    wmcr_hbb.setParamEffect(veto_tau, nveto_tau)
    addHbbNormSyst(wmcr_hbb)
    wmcr_hbb.setParamEffect(id_mu, nlepton)
    wmcr_hbb.setParamEffect(iso_mu, nlepton)
    addBtagSyst(background, recoilbin, "Hbb", "wmcr", wmcr_hbb, category, mass)
    wmcr.addSample(wmcr_hbb)

    wmcr_qcdTemplate = template(background, "QCD", "nominal", recoil, "wmcr", category, mass, read_sumw2=True)
    wmcr_qcd = rl.TemplateSample(ch_name + "_qcdMC", rl.Sample.BACKGROUND, wmcr_qcdTemplate)
    addLumiSyst(wmcr_qcd, year)
    addPileupSyst(wmcr_qcd)
    addPrefiringSyst(wmcr_qcd, year)
    addJESSyst(wmcr_qcd)
    addMETTrigSyst(wmcr_qcd, year)
    wmcr_qcd.setParamEffect(veto_tau, nveto_tau)
    wmcr_qcd.setParamEffect(qcdmu_norm, nqcd_norm)
    wmcr_qcd.setParamEffect(id_mu, nlepton)
    wmcr_qcd.setParamEffect(iso_mu, nlepton)
    addBtagSyst(background, recoilbin, "QCD", "wmcr", wmcr_qcd, category, mass)
    wmcr.addSample(wmcr_qcd)

    ###
    # W(->lnu)+jets data-driven model
    ###

    if not iswjetsMC:
        wmcr_wjetsTemplate = template(background, "W+jets", "nominal", recoil, "wmcr", category, mass, min_value=1., read_sumw2=True)
        wmcr_wjetsMC = rl.TemplateSample("wmcr" + model_id + "_wjetsMC", rl.Sample.BACKGROUND, wmcr_wjetsTemplate)
        addMETTrigSyst(wmcr_wjetsMC, year)
        wmcr_wjetsMC.setParamEffect(id_mu, nlepton)
        wmcr_wjetsMC.setParamEffect(iso_mu, nlepton)
        addVJetsSyst(background, recoil, "W+jets", "wmcr", wmcr_wjetsMC, category)
        
        tf, unc = makeTF(wmcr_wjetsMC, sr_wjetsMC)
        wmcr_wjets = TransferFactorSample(ch_name + "_wjets", rl.Sample.BACKGROUND, tf, sr._samples[sr_wjets.name], nominal_values=wmcr_wjetsMC._nominal, stat_unc=unc)#stat_unc=None)
        wmcr.addSample(wmcr_wjets)

    ###
    # top-antitop data-driven model
    ###

    if not isttMC:
        wmcr_ttTemplate = template(background, "TT", "nominal", recoil, "wmcr", category, mass, min_value=1., read_sumw2=True)
        wmcr_ttMC = rl.TemplateSample( "wmcr" + model_id + "_ttMC", rl.Sample.BACKGROUND, wmcr_ttTemplate)
        addMETTrigSyst(wmcr_ttMC, year)
        wmcr_ttMC.setParamEffect(id_mu, nlepton)
        wmcr_ttMC.setParamEffect(iso_mu, nlepton)
        addBtagSyst(background, recoil, "TT", "wmcr", wmcr_ttMC, category, mass)

        tf, unc = makeTF(wmcr_ttMC, sr_ttMC)
        wmcr_tt = TransferFactorSample(ch_name + "_tt", rl.Sample.BACKGROUND, tf, sr._samples[sr_tt.name], nominal_values=wmcr_ttMC._nominal, stat_unc=unc)#stat_unc=None)
        wmcr.addSample(wmcr_tt)

    ###
    # Add BB-lite
    ###

    addBBLiteSyst(wmcr)

    
    
    ###
    # End of single muon W control region
    ###

    ###
    ###
    # Single electron W control region
    ###
    ###

    ch_name = "wecr" + model_id
    wecr = rl.Channel(ch_name)
    model.addChannel(wecr)

    ###
    # Add data distribution to the channel
    ###

    if year == "2018":
        dataTemplate = template(data, "EGamma", "data", recoil, "wecr", category, mass)
    else:
        dataTemplate = template(data, "SingleElectron", "data", recoil, "wecr", category, mass)
    wecr.setObservation(dataTemplate)
    
    ###
    # Other MC-driven processes
    ###

    if iswjetsMC:
        wecr_wjetsTemplate = template(background, "W+jets", "nominal", recoil, "wecr", category, mass, read_sumw2=True)
        wecr_wjets = rl.TemplateSample("wecr" + model_id + "_wjetsMC", rl.Sample.BACKGROUND, wecr_wjetsTemplate)
        addLumiSyst(wecr_wjets, year)
        addPileupSyst(wecr_wjets)
        addPrefiringSyst(wecr_wjets, year)
        addJESSyst(wecr_wjets)
        addEleIDSyst(wecr_wjets, year)
        wecr_wjets.setParamEffect(veto_tau, nveto_tau)
        addWjetsNormSyst(wecr_wjets)
        wecr_wjets.setParamEffect(id_e, nlepton)
        wecr_wjets.setParamEffect(reco_e, nlepton)
        addBtagSyst(background, recoil, "W+jets", "wecr", wecr_wjets, category, mass)
        addVJetsSyst(background, recoil, "W+jets", "wecr", wecr_wjets, category)
        wecr.addSample(wecr_wjets)

    if isttMC: 
        wecr_ttTemplate = template(background, "TT", "nominal", recoil, "wecr", category, mass, read_sumw2=True)
        wecr_tt = rl.TemplateSample("wecr" + model_id + "_ttMC", rl.Sample.BACKGROUND, wecr_ttTemplate)
        addLumiSyst(wecr_tt, year)
        addPileupSyst(wecr_tt)
        addPrefiringSyst(wecr_tt, year)
        addJESSyst(wecr_tt)
        wecr_tt.setParamEffect(trig_e, ntrig_e)
        wecr_tt.setParamEffect(veto_tau, nveto_tau)
        addEleIDSyst(wecr_tt, year)
        wecr_tt.setParamEffect(reco_e, nlepton)
        addTTbarNormSyst(wecr_tt)
        addBtagSyst(background, recoil, "TT", "wecr", wecr_tt, category, mass)
        wecr.addSample(wecr_tt)

    wecr_stTemplate = template(background, "ST", "nominal", recoil, "wecr", category, mass, read_sumw2=True)
    wecr_st = rl.TemplateSample(ch_name + "_stMC", rl.Sample.BACKGROUND, wecr_stTemplate)
    addLumiSyst(wecr_st, year)
    addPileupSyst(wecr_st)
    addPrefiringSyst(wecr_st, year)
    addJESSyst(wecr_st)
    wecr_st.setParamEffect(trig_e, ntrig_e)
    wecr_st.setParamEffect(veto_tau, nveto_tau)
    addSingleTopNormSyst(wecr_st)
    addEleIDSyst(wecr_st, year)
    wecr_st.setParamEffect(reco_e, nlepton)
    addBtagSyst(background, recoilbin, "ST", "wecr", wecr_st, category, mass)
    wecr.addSample(wecr_st)

    wecr_dyjetsTemplate = template(background, "DY+jets", "nominal", recoil, "wecr", category, mass, read_sumw2=True)
    wecr_dyjets = rl.TemplateSample(ch_name + "_dyjetsMC", rl.Sample.BACKGROUND, wecr_dyjetsTemplate)
    addLumiSyst(wecr_dyjets, year)
    addPileupSyst(wecr_dyjets)
    addPrefiringSyst(wecr_dyjets, year)
    addJESSyst(wecr_dyjets)
    wecr_dyjets.setParamEffect(trig_e, ntrig_e)
    wecr_dyjets.setParamEffect(veto_tau, nveto_tau)
    addDrellYanNormSyst(wecr_dyjets)
    addEleIDSyst(wecr_dyjets, year)
    wecr_dyjets.setParamEffect(reco_e, nlepton)
    addBtagSyst(background, recoilbin, "DY+jets", "wecr", wecr_dyjets, category, mass)
    addVJetsSyst(background, recoil, "DY+jets", "wecr", wecr_dyjets, category)
    wecr.addSample(wecr_dyjets)

    wecr_vvTemplate = template(background, "VV", "nominal", recoil, "wecr", category, mass, read_sumw2=True)
    wecr_vv = rl.TemplateSample(ch_name + "_vvMC", rl.Sample.BACKGROUND, wecr_vvTemplate)
    addLumiSyst(wecr_vv, year)
    addPileupSyst(wecr_vv)
    addPrefiringSyst(wecr_vv, year)
    addJESSyst(wecr_vv)
    wecr_vv.setParamEffect(trig_e, ntrig_e)
    wecr_vv.setParamEffect(veto_tau, nveto_tau)
    addDibosonNormSyst(wecr_vv)
    addEleIDSyst(wecr_vv, year)
    wecr_vv.setParamEffect(reco_e, nlepton)
    addBtagSyst(background, recoilbin, "VV", "wecr", wecr_vv, category, mass)
    wecr.addSample(wecr_vv)

    wecr_hbbTemplate = template(background, "Hbb", "nominal", recoil, "wecr", category, mass, read_sumw2=True)
    wecr_hbb = rl.TemplateSample(ch_name + "_hbbMC", rl.Sample.BACKGROUND, wecr_hbbTemplate)
    addLumiSyst(wecr_hbb, year)
    addPileupSyst(wecr_hbb)
    addPrefiringSyst(wecr_hbb, year)
    addJESSyst(wecr_hbb)
    wecr_hbb.setParamEffect(trig_e, ntrig_e)
    wecr_hbb.setParamEffect(veto_tau, nveto_tau)
    addHbbNormSyst(wecr_hbb)
    addEleIDSyst(wecr_hbb, year)
    wecr_hbb.setParamEffect(reco_e, nlepton)
    addBtagSyst(background, recoilbin, "Hbb", "wecr", wecr_hbb, category, mass)
    wecr.addSample(wecr_hbb)

    wecr_qcdTemplate = template(background, "QCD", "nominal", recoil, "wecr", category, mass, read_sumw2=True)
    wecr_qcd = rl.TemplateSample(ch_name + "_qcdMC", rl.Sample.BACKGROUND, wecr_qcdTemplate)
    addLumiSyst(wecr_qcd, year)
    addPileupSyst(wecr_qcd)
    addPrefiringSyst(wecr_qcd, year)
    addJESSyst(wecr_qcd)
    wecr_qcd.setParamEffect(trig_e, ntrig_e)
    wecr_qcd.setParamEffect(veto_tau, nveto_tau)
    wecr_qcd.setParamEffect(qcde_norm, nqcd_norm)
    addEleIDSyst(wecr_qcd, year)
    wecr_qcd.setParamEffect(reco_e, nlepton)
    addBtagSyst(background, recoilbin, "QCD", "wecr", wecr_qcd, category, mass)
    wecr.addSample(wecr_qcd)


    ###
    # W(->lnu)+jets data-driven model
    ###

    if not iswjetsMC:
        wecr_wjetsTemplate = template(background, "W+jets", "nominal", recoil, "wecr", category, mass, min_value=1., read_sumw2=True)
        wecr_wjetsMC = rl.TemplateSample("wecr" + model_id + "_wjetsMC", rl.Sample.BACKGROUND, wecr_wjetsTemplate)
        wecr_wjetsMC.setParamEffect(trig_e, ntrig_e)
        addEleIDSyst(wecr_wjetsMC, year)
        wecr_wjetsMC.setParamEffect(reco_e, nlepton)
        addVJetsSyst(background, recoil, "W+jets", "wecr", wecr_wjetsMC, category)

        tf, unc = makeTF(wecr_wjetsMC, sr_wjetsMC)
        wecr_wjets = TransferFactorSample( ch_name + "_wjets", rl.Sample.BACKGROUND, tf, sr._samples[sr_wjets.name], nominal_values=wecr_wjetsMC._nominal, stat_unc=unc)#stat_unc=None)
        wecr.addSample(wecr_wjets)
        

    ###
    # top-antitop data-driven model
    ###

    if not isttMC: 
        wecr_ttTemplate = template(background, "TT", "nominal", recoil, "wecr", category, mass, min_value=1., read_sumw2=True)
        wecr_ttMC = rl.TemplateSample("wecr" + model_id + "_ttMC", rl.Sample.BACKGROUND, wecr_ttTemplate)
        wecr_ttMC.setParamEffect(trig_e, ntrig_e)
        addEleIDSyst(wecr_ttMC, year)
        wecr_ttMC.setParamEffect(reco_e, nlepton)
        addBtagSyst(background, recoil, "TT", "wecr", wecr_ttMC, category, mass)

        tf, unc = makeTF(wecr_ttMC, sr_ttMC)
        wecr_tt = TransferFactorSample( ch_name + "_tt", rl.Sample.BACKGROUND, tf, sr._samples[sr_tt.name], nominal_values=wecr_ttMC._nominal, stat_unc=unc)#stat_unc=None)
        wecr.addSample(wecr_tt)

    ###
    # Add BB-lite
    ###

    addBBLiteSyst(wecr) 

    
    ###
    # End of single electron W control region
    ###

    if category=="fail": return model

    ###
    ###
    # Single muon top control region
    ###
    ###

    ch_name = "tmcr" + model_id
    tmcr = rl.Channel(ch_name)
    model.addChannel(tmcr)

    ###
    # Add data distribution to the channel
    ###

    dataTemplate = template(data, "MET", "data", recoil, "tmcr", category, mass)
    tmcr.setObservation(dataTemplate)
    
    ###
    # Other MC-driven processes
    ###

    if isttMC:
        tmcr_ttTemplate = template(background, "TT", "nominal", recoil, "tmcr", category, mass, read_sumw2=True)
        tmcr_tt = rl.TemplateSample("tmcr" + model_id + "_ttMC", rl.Sample.BACKGROUND, tmcr_ttTemplate)
        addLumiSyst(tmcr_tt, year)
        addPileupSyst(tmcr_tt)
        addPrefiringSyst(tmcr_tt, year)
        addJESSyst(tmcr_tt)
        addMETTrigSyst(tmcr_tt, year)
        tmcr_tt.setParamEffect(veto_tau, nveto_tau)
        addTTbarNormSyst(tmcr_tt)
        tmcr_tt.setParamEffect(id_mu, nlepton)
        tmcr_tt.setParamEffect(iso_mu, nlepton)
        addBtagSyst(background, recoil, "TT", "tmcr", tmcr_tt, category, mass)
        tmcr.addSample(tmcr_tt)

    tmcr_wjetsTemplate = template(background, "W+jets", "nominal", recoil, "tmcr", category, mass, read_sumw2=True)
    tmcr_wjets = rl.TemplateSample(ch_name + "_wjetsMC", rl.Sample.BACKGROUND, tmcr_wjetsTemplate)
    addLumiSyst(tmcr_wjets, year)
    addPileupSyst(tmcr_wjets)
    addPrefiringSyst(tmcr_wjets, year)
    addJESSyst(tmcr_wjets)
    addMETTrigSyst(tmcr_wjets, year)
    tmcr_wjets.setParamEffect(veto_tau, nveto_tau)
    addWjetsNormSyst(tmcr_wjets)
    tmcr_wjets.setParamEffect(id_mu, nlepton)
    tmcr_wjets.setParamEffect(iso_mu, nlepton)
    addBtagSyst(background, recoilbin, "W+jets", "tmcr", tmcr_wjets, category, mass)
    addVJetsSyst(background, recoil, "W+jets", "tmcr", tmcr_wjets, category)
    tmcr.addSample(tmcr_wjets)

    tmcr_stTemplate = template(background, "ST", "nominal", recoil, "tmcr", category, mass, read_sumw2=True)
    tmcr_st = rl.TemplateSample(ch_name + "_stMC", rl.Sample.BACKGROUND, tmcr_stTemplate)
    addLumiSyst(tmcr_st, year)
    addPileupSyst(tmcr_st)
    addPrefiringSyst(tmcr_st, year)
    addJESSyst(tmcr_st)
    addMETTrigSyst(tmcr_st, year)
    tmcr_st.setParamEffect(veto_tau, nveto_tau)
    addSingleTopNormSyst(tmcr_st)
    tmcr_st.setParamEffect(id_mu, nlepton)
    tmcr_st.setParamEffect(iso_mu, nlepton)
    addBtagSyst(background, recoilbin, "ST", "tmcr", tmcr_st, category, mass)
    tmcr.addSample(tmcr_st)

    tmcr_dyjetsTemplate = template(background, "DY+jets", "nominal", recoil, "tmcr", category, mass, read_sumw2=True)
    tmcr_dyjets = rl.TemplateSample(ch_name + "_dyjetsMC", rl.Sample.BACKGROUND, tmcr_dyjetsTemplate)
    addLumiSyst(tmcr_dyjets, year)
    addPileupSyst(tmcr_dyjets)
    addPrefiringSyst(tmcr_dyjets, year)
    addJESSyst(tmcr_dyjets)
    addMETTrigSyst(tmcr_dyjets, year)
    tmcr_dyjets.setParamEffect(veto_tau, nveto_tau)
    addDrellYanNormSyst(tmcr_dyjets)
    tmcr_dyjets.setParamEffect(id_mu, nlepton)
    tmcr_dyjets.setParamEffect(iso_mu, nlepton)
    addBtagSyst(background, recoilbin, "DY+jets", "tmcr", tmcr_dyjets, category, mass)
    addVJetsSyst(background, recoil, "DY+jets", "tmcr", tmcr_dyjets, category)
    tmcr.addSample(tmcr_dyjets)

    tmcr_vvTemplate = template(background, "VV", "nominal", recoil, "tmcr", category, mass, read_sumw2=True)
    tmcr_vv = rl.TemplateSample(ch_name + "_vvMC", rl.Sample.BACKGROUND, tmcr_vvTemplate)
    addLumiSyst(tmcr_vv, year)
    addPileupSyst(tmcr_vv)
    addPrefiringSyst(tmcr_vv, year)
    addJESSyst(tmcr_vv)
    addMETTrigSyst(tmcr_vv, year)
    tmcr_vv.setParamEffect(veto_tau, nveto_tau)
    addDibosonNormSyst(tmcr_vv)
    tmcr_vv.setParamEffect(id_mu, nlepton)
    tmcr_vv.setParamEffect(iso_mu, nlepton)
    addBtagSyst(background, recoilbin, "VV", "tmcr", tmcr_vv, category, mass)
    tmcr.addSample(tmcr_vv)

    tmcr_hbbTemplate = template(background, "Hbb", "nominal", recoil, "tmcr", category, mass, read_sumw2=True)
    tmcr_hbb = rl.TemplateSample(ch_name + "_hbbMC", rl.Sample.BACKGROUND, tmcr_hbbTemplate)
    addLumiSyst(tmcr_hbb, year)
    addPileupSyst(tmcr_hbb)
    addPrefiringSyst(tmcr_hbb, year)
    addJESSyst(tmcr_hbb)
    addMETTrigSyst(tmcr_hbb, year)
    tmcr_hbb.setParamEffect(veto_tau, nveto_tau)
    addHbbNormSyst(tmcr_hbb)
    tmcr_hbb.setParamEffect(id_mu, nlepton)
    tmcr_hbb.setParamEffect(iso_mu, nlepton)
    addBtagSyst(background, recoilbin, "Hbb", "tmcr", tmcr_hbb, category, mass)
    tmcr.addSample(tmcr_hbb)

    tmcr_qcdTemplate = template(background, "QCD", "nominal", recoil, "tmcr", category, mass, read_sumw2=True)
    tmcr_qcd = rl.TemplateSample(ch_name + "_qcdMC", rl.Sample.BACKGROUND, tmcr_qcdTemplate)
    addLumiSyst(tmcr_qcd, year)
    addPileupSyst(tmcr_qcd)
    addPrefiringSyst(tmcr_qcd, year)
    addJESSyst(tmcr_qcd)
    addMETTrigSyst(tmcr_qcd, year)
    tmcr_qcd.setParamEffect(veto_tau, nveto_tau)
    tmcr_qcd.setParamEffect(qcdmu_norm, nqcd_norm)
    tmcr_qcd.setParamEffect(id_mu, nlepton)
    tmcr_qcd.setParamEffect(iso_mu, nlepton)
    addBtagSyst(background, recoilbin, "QCD", "tmcr", tmcr_qcd, category, mass)
    tmcr.addSample(tmcr_qcd)

    ###
    # top-antitop data-driven model
    ###

    if not isttMC:
        tmcr_ttTemplate = template(background, "TT", "nominal", recoil, "tmcr", category, mass, min_value=1., read_sumw2=True)
        tmcr_ttMC = rl.TemplateSample("tmcr" + model_id + "_ttMC", rl.Sample.BACKGROUND, tmcr_ttTemplate)
        addMETTrigSyst(tmcr_ttMC, year)
        tmcr_ttMC.setParamEffect(id_mu, nlepton)
        tmcr_ttMC.setParamEffect(iso_mu, nlepton)
        addBtagSyst(background, recoil, "TT", "tmcr", tmcr_ttMC, category, mass)

        tf, unc = makeTF(tmcr_ttMC, sr_ttMC)
        tmcr_tt = TransferFactorSample(ch_name + "_tt", rl.Sample.BACKGROUND, tf, sr._samples[sr_tt.name], nominal_values=tmcr_ttMC._nominal, stat_unc=unc)#stat_unc=None)
        tmcr.addSample(tmcr_tt)

    ###
    # Add BB-lite
    ###

    addBBLiteSyst(tmcr)

    

    ###
    # End of single muon top control region
    ###

    ###
    ###
    # Single electron top control region
    ###
    ###

    ch_name = "tecr" + model_id
    tecr = rl.Channel(ch_name)
    model.addChannel(tecr)

    ###
    # Add data distribution to the channel
    ###

    if year == "2018":
        dataTemplate = template(data, "EGamma", "data", recoil, "tecr", category, mass)
    else:
        dataTemplate = template(data, "SingleElectron", "data", recoil, "tecr", category, mass)
    tecr.setObservation(dataTemplate)
    
    ###
    # Other MC-driven processes
    ###

    if isttMC:
        tecr_ttTemplate = template(background, "TT", "nominal", recoil, "tecr", category, mass, read_sumw2=True)
        tecr_tt = rl.TemplateSample("tecr" + model_id + "_ttMC", rl.Sample.BACKGROUND, tecr_ttTemplate)
        addLumiSyst(tecr_tt, year)
        addPileupSyst(tecr_tt)
        addPrefiringSyst(tecr_tt, year)
        addJESSyst(tecr_tt)
        tecr_tt.setParamEffect(trig_e, ntrig_e)
        tecr_tt.setParamEffect(veto_tau, nveto_tau)
        addTTbarNormSyst(tecr_tt)
        addEleIDSyst(tecr_tt, year)
        tecr_tt.setParamEffect(reco_e, nlepton)
        addBtagSyst(background, recoil, "TT", "tecr", tecr_tt, category, mass)
        tecr.addSample(tecr_tt)

    tecr_wjetsTemplate = template(background, "W+jets", "nominal", recoil, "tecr", category, mass, read_sumw2=True)
    tecr_wjets = rl.TemplateSample(ch_name + "_wjetsMC", rl.Sample.BACKGROUND, tecr_wjetsTemplate)
    addLumiSyst(tecr_wjets, year)
    addPileupSyst(tecr_wjets)
    addPrefiringSyst(tecr_wjets, year)
    addJESSyst(tecr_wjets)
    tecr_wjets.setParamEffect(trig_e, ntrig_e)
    tecr_wjets.setParamEffect(veto_tau, nveto_tau)
    addWjetsNormSyst(tecr_wjets)
    addEleIDSyst(tecr_wjets, year)
    tecr_wjets.setParamEffect(reco_e, nlepton)
    addBtagSyst(background, recoilbin, "W+jets", "tecr", tecr_wjets, category, mass)
    addVJetsSyst(background, recoil, "W+jets", "tecr", tecr_wjets, category)
    tecr.addSample(tecr_wjets)

    tecr_stTemplate = template(background, "ST", "nominal", recoil, "tecr", category, mass, read_sumw2=True)
    tecr_st = rl.TemplateSample(ch_name + "_stMC", rl.Sample.BACKGROUND, tecr_stTemplate)
    addLumiSyst(tecr_st, year)
    addPileupSyst(tecr_st)
    addPrefiringSyst(tecr_st, year)
    addJESSyst(tecr_st)
    tecr_st.setParamEffect(trig_e, ntrig_e)
    tecr_st.setParamEffect(veto_tau, nveto_tau)
    addSingleTopNormSyst(tecr_st)
    addEleIDSyst(tecr_st, year)
    tecr_st.setParamEffect(reco_e, nlepton)
    addBtagSyst(background, recoilbin, "ST", "tecr", tecr_st, category, mass)
    tecr.addSample(tecr_st)

    tecr_dyjetsTemplate = template(background, "DY+jets", "nominal", recoil, "tecr", category, mass, read_sumw2=True)
    tecr_dyjets = rl.TemplateSample(ch_name + "_dyjetsMC", rl.Sample.BACKGROUND, tecr_dyjetsTemplate)
    addLumiSyst(tecr_dyjets, year)
    addPileupSyst(tecr_dyjets)
    addPrefiringSyst(tecr_dyjets, year)
    addJESSyst(tecr_dyjets)
    tecr_dyjets.setParamEffect(trig_e, ntrig_e)
    tecr_dyjets.setParamEffect(veto_tau, nveto_tau)
    addDrellYanNormSyst(tecr_dyjets)
    addEleIDSyst(tecr_dyjets, year)
    tecr_dyjets.setParamEffect(reco_e, nlepton)
    addBtagSyst(background, recoilbin, "DY+jets", "tecr", tecr_dyjets, category, mass)
    addVJetsSyst(background, recoil, "DY+jets", "tecr", tecr_dyjets, category)
    tecr.addSample(tecr_dyjets)

    tecr_vvTemplate = template(background, "VV", "nominal", recoil, "tecr", category, mass, read_sumw2=True)
    tecr_vv = rl.TemplateSample(ch_name + "_vvMC", rl.Sample.BACKGROUND, tecr_vvTemplate)
    addLumiSyst(tecr_vv, year)
    addPileupSyst(tecr_vv)
    addPrefiringSyst(tecr_vv, year)
    addJESSyst(tecr_vv)
    tecr_vv.setParamEffect(trig_e, ntrig_e)
    tecr_vv.setParamEffect(veto_tau, nveto_tau)
    addDibosonNormSyst(tecr_vv)
    addEleIDSyst(tecr_vv, year)
    tecr_vv.setParamEffect(reco_e, nlepton)
    addBtagSyst(background, recoilbin, "VV", "tecr", tecr_vv, category, mass)
    tecr.addSample(tecr_vv)

    tecr_hbbTemplate = template(background, "Hbb", "nominal", recoil, "tecr", category, mass, read_sumw2=True)
    tecr_hbb = rl.TemplateSample(ch_name + "_hbbMC", rl.Sample.BACKGROUND, tecr_hbbTemplate)
    addLumiSyst(tecr_hbb, year)
    addPileupSyst(tecr_hbb)
    addPrefiringSyst(tecr_hbb, year)
    addJESSyst(tecr_hbb)
    tecr_hbb.setParamEffect(trig_e, ntrig_e)
    tecr_hbb.setParamEffect(veto_tau, nveto_tau)
    addHbbNormSyst(tecr_hbb)
    addEleIDSyst(tecr_hbb, year)
    tecr_hbb.setParamEffect(reco_e, nlepton)
    addBtagSyst(background, recoilbin, "Hbb", "tecr", tecr_hbb, category, mass)
    tecr.addSample(tecr_hbb)

    tecr_qcdTemplate = template(background, "QCD", "nominal", recoil, "tecr", category, mass, read_sumw2=True)
    tecr_qcd = rl.TemplateSample(ch_name + "_qcdMC", rl.Sample.BACKGROUND, tecr_qcdTemplate)
    addLumiSyst(tecr_qcd, year)
    addPileupSyst(tecr_qcd)
    addPrefiringSyst(tecr_qcd, year)
    addJESSyst(tecr_qcd)
    tecr_qcd.setParamEffect(trig_e, ntrig_e)
    tecr_qcd.setParamEffect(veto_tau, nveto_tau)
    tecr_qcd.setParamEffect(qcde_norm, nqcd_norm)
    addEleIDSyst(tecr_qcd, year)
    tecr_qcd.setParamEffect(reco_e, nlepton)
    addBtagSyst(background, recoilbin, "QCD", "tecr", tecr_qcd, category, mass)
    tecr.addSample(tecr_qcd)

    ###
    # top-antitop data-driven model
    ###

    if not isttMC:
        tecr_ttTemplate = template(background, "TT", "nominal", recoil, "tecr", category, mass, min_value=1., read_sumw2=True)
        tecr_ttMC = rl.TemplateSample("tecr" + model_id + "_ttMC", rl.Sample.BACKGROUND, tecr_ttTemplate)
        tecr_ttMC.setParamEffect(trig_e, ntrig_e)
        addEleIDSyst(tecr_ttMC, year)
        tecr_ttMC.setParamEffect(reco_e, nlepton)
        addBtagSyst(background, recoil, "TT", "tecr", tecr_ttMC, category, mass)

        tf, unc = makeTF(tecr_ttMC, sr_ttMC)
        tecr_tt = TransferFactorSample(ch_name + "_tt", rl.Sample.BACKGROUND, tf, sr._samples[sr_tt.name], nominal_values=tecr_ttMC._nominal, stat_unc=unc)#stat_unc=None)
        tecr.addSample(tecr_tt)

    ###
    # Add BB-lite
    ###

    addBBLiteSyst(tecr)

    
    ###
    # End of single electron top control region
    ###

    return model


if __name__ == "__main__":
    #if not os.path.exists("datacards"):
    #    os.mkdir("datacards")
    parser = OptionParser()
    parser.add_option("-y", "--year", help="year", dest="year", default="2018")
    parser.add_option("-m", "--mass", help="mass", dest="mass", default="40to300")
    parser.add_option("-f", "--fakedata", help="replace data to sum of backgrounds", action="store_true", dest="fakedata")
    (options, args) = parser.parse_args()
    year = options.year
    mass = options.mass

    #####
    ###
    # Preparing Rhalphabeth
    ###
    #####

    ###
    # Extract histograms from input file and remap
    ###

    hists = load("hists/darkhiggs" + year + ".scaled")
    hists = remap_histograms(hists)

    ###
    # Preparing histograms for Rhalphabeth
    ###

    background = {}
    for r in hists["bkg"]["template"].identifiers("region"):
        background[str(r)] = hists["bkg"]["template"].integrate("region", r)
    
    ###
    # Establishing 2D binning
    ###
    
    recoilbins = np.array(recoil_binning)
    nrecoil = len(recoilbins) - 1
    msdbins = np.array(mass_binning)
    msd = rl.Observable('fjmass', msdbins)
    # here we derive these all at once with 2D array
    ptpts, msdpts = np.meshgrid(recoilbins[:-1] + 0.3 * np.diff(recoilbins), msdbins[:-1] + 0.5 * np.diff(msdbins), indexing='ij')
    recoilscaled = (ptpts - recoil_binning[0]) / (recoil_binning[-1] - recoil_binning[0])
    msdscaled = (msdpts - mass_binning[0]) / (mass_binning[-1] - mass_binning[0])

    ###
    # Calculating average pass-to-fail ratio
    ###
    
    def efficiency(pass_templ, fail_templ, qcdmodel):
        qcdpass, qcdfail = 0., 0.
        for recoilbin in range(nrecoil):
            failCh = rl.Channel("recoilbin%d%s" % (recoilbin, 'fail'))
            passCh = rl.Channel("recoilbin%d%s" % (recoilbin, 'pass'))
            qcdmodel.addChannel(failCh)
            qcdmodel.addChannel(passCh)
            failCh.setObservation(fail_templ[recoilbin])
            passCh.setObservation(pass_templ[recoilbin])
            qcdfail += failCh.getObservation().sum()
            qcdpass += passCh.getObservation().sum()
            '''
            try:
                qcdfail = np.append(qcdfail, [failCh.getObservation()], axis = 0)
                qcdpass = np.append(qcdpass, [passCh.getObservation()], axis = 0)
            except:
                qcdfail = np.array([failCh.getObservation()])
                qcdpass = np.array([passCh.getObservation()])
            ''' 

        return qcdpass / qcdfail

    ###
    # Creating first Bernstein polynomial that represents the MC pass-to-fail ratio
    # Incorporates the dependence on mass/recoil (residual tagger correlation, HF fraction)
    # Includes statistical uncertainty by fitting the by-by-bin MC pass-to-fail ratio
    ###

    zjetspass_templ = []
    zjetsfail_templ = []
    for recoilbin in range(nrecoil):
        zjetspass_templ.append(template(background, "Z+jets", "nominal", recoilbin, "sr", "pass", ""))
        zjetsfail_templ.append(template(background, "Z+jets", "nominal", recoilbin, "sr", "fail", ""))

    zjetsmodel = rl.Model("zjetsmodel")
    zjetseff = efficiency(zjetspass_templ, zjetsfail_templ, zjetsmodel)
    tf_MCtemplZ = rl.BernsteinPoly("tf_MCtemplZ"+year, (1, 1), ['recoil', 'fjmass'], limits=(1e-5, 10))
    tf_MCtemplZ_params = zjetseff * tf_MCtemplZ(recoilscaled, msdscaled)

    wjetspass_templ = []
    wjetsfail_templ = []
    for recoilbin in range(nrecoil):
        wjetspass_templ.append(template(background, "W+jets", "nominal", recoilbin, "sr", "pass", ""))
        wjetsfail_templ.append(template(background, "W+jets", "nominal", recoilbin, "sr", "fail", ""))

    wjetsmodel = rl.Model("wjetsmodel")
    wjetseff = efficiency(wjetspass_templ, wjetsfail_templ, wjetsmodel)
    tf_MCtemplW = rl.BernsteinPoly("tf_MCtemplW"+year, (1, 1), ['recoil', 'fjmass'], limits=(1e-5, 10))
    tf_MCtemplW_params = wjetseff * tf_MCtemplW(recoilscaled, msdscaled)
    
    ###
    # Prepare model for the MC ratio fit
    ##

    def rhalphabeth(pass_templ, fail_templ, qcdmodel, tf_MCtempl_params):

        for recoilbin in range(nrecoil):
            failCh = qcdmodel['recoilbin%dfail' % recoilbin]
            passCh = qcdmodel['recoilbin%dpass' % recoilbin]
            failObs = failCh.getObservation()
            qcdparams = np.array([rl.IndependentParameter('qcdparam_ptbin%d_msdbin%d' % (recoilbin, i), 0) for i in range(msd.nbins)])
            sigmascale = 10.
            scaledparams = failObs * (1 + sigmascale/np.maximum(1., np.sqrt(failObs)))**qcdparams
            fail_qcd = rl.ParametericSample('recoilbin'+str(recoilbin)+'fail_'+qcdmodel.name, rl.Sample.BACKGROUND, msd, scaledparams)
            failCh.addSample(fail_qcd)
            pass_qcd = rl.TransferFactorSample('recoilbin'+str(recoilbin)+'pass_'+qcdmodel.name, rl.Sample.BACKGROUND, tf_MCtempl_params[recoilbin, :], fail_qcd)
            passCh.addSample(pass_qcd)

        return qcdmodel

    zjetsmodel = rhalphabeth(zjetspass_templ, zjetsfail_templ, zjetsmodel, tf_MCtemplZ_params)
    wjetsmodel = rhalphabeth(wjetspass_templ, wjetsfail_templ, wjetsmodel, tf_MCtemplW_params)
    
    ###
    # Perform the fit to the bin-by-bin MC ratio
    ###

    def fit(model):
        qcdfit_ws = ROOT.RooWorkspace('qcdfit_ws')
        simpdf, obs = model.renderRoofit(qcdfit_ws)
        qcdfit = simpdf.fitTo(obs,
                            ROOT.RooFit.Extended(True),
                            ROOT.RooFit.SumW2Error(True),
                            ROOT.RooFit.Strategy(2),
                            ROOT.RooFit.Save(),
                            ROOT.RooFit.Minimizer('Minuit2', 'migrad'),
                            ROOT.RooFit.PrintLevel(-1),
                            )
        qcdfit_ws.add(qcdfit)
        #if "pytest" not in sys.modules:
        #    qcdfit_ws.writeToFile(os.path.join(str(tmpdir), 'testModel_qcdfit.root'))
        if qcdfit.status() != 0:
            raise RuntimeError('Could not fit qcd')

        return qcdfit

    zjetsfit = fit(zjetsmodel)
    wjetsfit = fit(wjetsmodel)

    ###
    # Use the post-fit values of the Bernstein polynomial coefficients
    ###
    
    def shape(fit, tf_MCtempl):
        param_names = [p.name for p in tf_MCtempl.parameters.reshape(-1)]
        decoVector = rl.DecorrelatedNuisanceVector.fromRooFitResult(tf_MCtempl.name + '_deco', fit, param_names)
        tf_MCtempl.parameters = decoVector.correlated_params.reshape(tf_MCtempl.parameters.shape)
        tf_MCtempl_params_final = tf_MCtempl(recoilscaled, msdscaled)

        return tf_MCtempl_params_final

    tf_MCtemplW_params_final = shape(wjetsfit, tf_MCtemplW)
    tf_MCtemplZ_params_final = shape(zjetsfit, tf_MCtemplZ)
    
    ###
    # Create Bernstein polynomials that represent the correction to the MC ratio
    ###
    
    tf_dataResidualW = rl.BernsteinPoly("tf_dataResidualW"+year, (0, 1), ['recoil', 'fjmass'], limits=(1e-5, 10))
    tf_dataResidualW_params = tf_dataResidualW(recoilscaled, msdscaled)
    tf_dataResidualZ = rl.BernsteinPoly("tf_dataResidualZ"+year, (0, 1), ['recoil', 'fjmass'], limits=(1e-5, 10))
    tf_dataResidualZ_params = tf_dataResidualZ(recoilscaled, msdscaled)

    #####
    ###
    # End of Rhalphabeth preparation
    ###
    #####
    
    ###
    ###
    # Prepare histograms for the fit
    ###
    ###
    
    ###
    # Split mass range
    ###
    
    
    if '40to' in mass:
        cut = mass.split('40to')[1]
        index = mass_binning.index(int(cut))
        mass_binning = mass_binning[:(index+1)]
        nmass = len(mass_binning) - 1
        tf_MCtemplZ_params_final = tf_MCtemplZ_params_final[:, :nmass]
        tf_dataResidualZ_params = tf_dataResidualZ_params[:, :nmass]
        tf_MCtemplW_params_final = tf_MCtemplW_params_final[:, :nmass]
        tf_dataResidualW_params = tf_dataResidualW_params[:, :nmass]
    if 'to300' in mass:
        nmass = len(mass_binning) - 1
        cut = mass.split('to300')[0]
        index = mass_binning.index(int(cut))
        mass_binning = mass_binning[index:]
        nmass = nmass - (len(mass_binning) - 1)
        tf_MCtemplZ_params_final = tf_MCtemplZ_params_final[:, nmass:]
        tf_dataResidualZ_params = tf_dataResidualZ_params[:, nmass:]
        tf_MCtemplW_params_final = tf_MCtemplW_params_final[:, nmass:]
        tf_dataResidualW_params = tf_dataResidualW_params[:, nmass:]
        
    ###
    # Reload and remap histograms 
    ###

    hists = load("hists/darkhiggs" + year + ".scaled")
    hists = remap_histograms(hists)

    ###
    # Manipulate histograms to be fed to the model
    ###

    signal_hists = hists["sig"]
    signal = {}
    for r in signal_hists["template"].identifiers("region"):
        signal[str(r)] = signal_hists["template"].integrate("region", r)
        
    bkg_hists = hists["bkg"]
    background = {}
    for r in bkg_hists["template"].identifiers("region"):
        background[str(r)] = bkg_hists["template"].integrate("region", r)

    data_hists = hists["data"]
    data = {}
    for r in data_hists["template"].identifiers("region"):
        data[str(r)] = data_hists["template"].integrate("region", r)

    ###
    ###
    # Set up other systematics
    ###
    ###
    
    lumi = rl.NuisanceParameter("lumi_" + year, "lnN")
    lumi_corr = rl.NuisanceParameter("lumi_13TeV_correlated", "lnN")
    lumi_1718 = rl.NuisanceParameter("lumi_13TeV_1718", "lnN")
    pu = rl.NuisanceParameter("CMS_pileup_" + year, "lnN")
    prefiring = rl.NuisanceParameter("CMS_l1_ecal_prefiring_" + year, "lnN")
    id_e = rl.NuisanceParameter("CMS_eff_e_id_" + year, "lnN")
    id_mu = rl.NuisanceParameter("CMS_eff_m_id_" + year, "lnN")
    id_pho = rl.NuisanceParameter("CMS_eff_g_IDMVA_" + year, "lnN")
    reco_e = rl.NuisanceParameter("CMS_eff_e_reco_" + year, "lnN")
    iso_mu = rl.NuisanceParameter("CMS_eff_m_iso_" + year, "lnN")
    trig_e = rl.NuisanceParameter("CMS_eff_e_trigger_" + year, "lnN")
    trig_met = rl.NuisanceParameter("trig_met" + year, "lnN")
    trig_pho = rl.NuisanceParameter("trig_pho" + year, "lnN")
    veto_tau = rl.NuisanceParameter("veto_tau" + year, "lnN")
    jes = rl.NuisanceParameter("CMS_scale_j_" + year, "lnN")
    met = rl.NuisanceParameter("CMS_scale_met" + year, "lnN")
    btagSFbc_correlated = rl.NuisanceParameter("CMS_btag_heavy_correlated", "shape")  # AK4 btag
    btagSFbc = rl.NuisanceParameter("CMS_btag_heavy" + year, "shape")  # AK4 btag
    btagSFlight_correlated = rl.NuisanceParameter("CMS_btag_light_correlated", "shape")  # AK4 btag
    btagSFlight = rl.NuisanceParameter("CMS_btag_light" + year, "shape")  # AK4 btag 
    ew1 = rl.NuisanceParameter("ew1", "lnN") #Effects of unknown Sudakov logs
    #ew2G = rl.NuisanceParameter("ew2G", "lnN")
    ew2W = rl.NuisanceParameter("ew2W", "lnN") #Missing NNLO effects for W boson
    ew2Z = rl.NuisanceParameter("ew2Z", "lnN") #Missing NNLO effects for Z boson
    #ew3G = rl.NuisanceParameter("ew3G", "lnN")
    ew3W = rl.NuisanceParameter("ew3W", "lnN") #Effects of NLL Sudakov approx. w W boson
    ew3Z = rl.NuisanceParameter("ew3Z", "lnN") #Effects of NLL Sudakov approx. w Z boson
    mix = rl.NuisanceParameter("mix", "lnN")
    #muF = rl.NuisanceParameter("muF", "lnN")
    #muR = rl.NuisanceParameter("muR", "lnN")
    qcd1 = rl.NuisanceParameter("qcd1", "lnN")
    qcd2 = rl.NuisanceParameter("qcd2", "lnN") #pT shape dependence
    qcd3 = rl.NuisanceParameter("qcd3", "lnN") #Process dependence
    #doublebtag = rl.NuisanceParameter("doublebtag_" + year, "shape")
    doublebtag = rl.NuisanceParameter("doublebtag_" + year, "lnN")
        
    ###
    # Set lnN numbers
    ###

    ntrig_e = 1.01
    nveto_tau = 1.03
    nlepton = 1.01 ## id_mu, iso_mu, reco_e
    nqcd_norm = 2.0 ## qcdsig_norm, qcde_norm, qcdmu_norm

    ###
    ###
    # End of systematics setup
    ###
    ###

    for recoil in range(nrecoil):

        #######
        #####
        ###
        # Building the "fail" model
        ###  
        #####
        #######
        
        isttMC = True
        iswjetsMC = False
        
        qcdpho_norm = rl.NuisanceParameter("qcdpho_norm" + year + "fail", "lnN")
        qcde_norm = rl.NuisanceParameter("qcde_norm" + year + "fail", "lnN")
        qcdmu_norm = rl.NuisanceParameter("qcdmu_norm" + year + "fail", "lnN")
        qcdsig_norm = rl.NuisanceParameter("qcdsig_norm" + year + "fail", "lnN")
        st_norm = rl.NuisanceParameter("st_norm" + year + "fail", "lnN")
        ttMC_norm = rl.NuisanceParameter("tt_norm" + year + "fail", "lnN")
        vv_norm = rl.NuisanceParameter("vv_norm" + year + "fail", "lnN")
        hbb_norm = rl.NuisanceParameter("hbb_norm" + year + "fail", "lnN")
        wjetsMC_norm = rl.NuisanceParameter("wjets_norm" + year + "fail", "lnN")
        zjetsMC_norm = rl.NuisanceParameter("zjets_norm" + year + "fail", "lnN")

        model_id = year + "fail" + "mass" + mass+ "recoil" + str(recoil)
        model_fail = rl.Model(model_id)
    
        #####
        ###
        # Signal region "fail"
        ###  
        #####

        ch_name = "sr" + model_id
        sr_fail = rl.Channel(ch_name)
        model_fail.addChannel(sr_fail)
    
        ###
        # Add data distribution to the channel
        ###
    
        dataTemplate = template(data, "MET", "data", recoil, "sr", "fail", mass)
        sr_fail.setObservation(dataTemplate)
    
        ###
        # Other MC-driven processes
        ###
        
        sr_fail_ttTemplate = template(background, "TT", "nominal", recoil, "sr", "fail", mass, read_sumw2=True)
        sr_fail_tt = rl.TemplateSample("sr" + model_id + "_ttMC",rl.Sample.BACKGROUND, sr_fail_ttTemplate)
        addLumiSyst(sr_fail_tt, year)
        addPileupSyst(sr_fail_tt)
        addPrefiringSyst(sr_fail_tt, year)
        addMETSyst(sr_fail_tt)
        addJESSyst(sr_fail_tt)
        addMETTrigSyst(sr_fail_tt, year)
        sr_fail_tt.setParamEffect(veto_tau, nveto_tau)
        addTTbarNormSyst(sr_fail_tt)
        addBtagSyst(background, recoil, "TT", "sr", sr_fail_tt, "fail", mass)
        sr_fail.addSample(sr_fail_tt)
        
        sr_fail_stTemplate = template(background, "ST", "nominal", recoil, "sr", "fail", mass, read_sumw2=True)
        sr_fail_st = rl.TemplateSample(ch_name + "_stMC", rl.Sample.BACKGROUND, sr_fail_stTemplate)
        addLumiSyst(sr_fail_st, year)
        addPileupSyst(sr_fail_st)
        addPrefiringSyst(sr_fail_st, year)
        addMETSyst(sr_fail_st)
        addJESSyst(sr_fail_st)
        addMETTrigSyst(sr_fail_st, year)
        sr_fail_st.setParamEffect(veto_tau, nveto_tau)
        addSingleTopNormSyst(sr_fail_st)
        addBtagSyst(background, recoil, "ST", "sr", sr_fail_st, "fail", mass)
        sr_fail.addSample(sr_fail_st)
    
        sr_fail_dyjetsTemplate = template(background, "DY+jets", "nominal", recoil, "sr", "fail", mass, read_sumw2=True)
        sr_fail_dyjets = rl.TemplateSample(ch_name + "_dyjetsMC", rl.Sample.BACKGROUND, sr_fail_dyjetsTemplate)
        addLumiSyst(sr_fail_dyjets, year)
        addPileupSyst(sr_fail_dyjets)
        addPrefiringSyst(sr_fail_dyjets, year)
        addMETSyst(sr_fail_dyjets)
        addJESSyst(sr_fail_dyjets)
        addMETTrigSyst(sr_fail_dyjets, year)
        sr_fail_dyjets.setParamEffect(veto_tau, nveto_tau)
        addDrellYanNormSyst(sr_fail_dyjets)
        addBtagSyst(background, recoil, "DY+jets", "sr", sr_fail_dyjets, "fail", mass)
        addVJetsSyst(background, recoil, "DY+jets", "sr", sr_fail_dyjets, "fail")
        sr_fail.addSample(sr_fail_dyjets)
    
        sr__fail_vvTemplate = template(background, "VV", "nominal", recoil, "sr", "fail", mass, read_sumw2=True)
        sr_fail_vv = rl.TemplateSample(ch_name + "_vvMC", rl.Sample.BACKGROUND, sr__fail_vvTemplate)
        addLumiSyst(sr_fail_vv, year)
        addPileupSyst(sr_fail_vv)
        addPrefiringSyst(sr_fail_vv, year)
        addMETSyst(sr_fail_vv)
        addJESSyst(sr_fail_vv)
        addMETTrigSyst(sr_fail_vv, year)
        sr_fail_vv.setParamEffect(veto_tau, nveto_tau)
        addDibosonNormSyst(sr_fail_vv)
        addBtagSyst(background, recoil, "VV", "sr", sr_fail_vv, "fail", mass)
        sr_fail.addSample(sr_fail_vv)
    
        sr_fail_hbbTemplate = template(background, "Hbb", "nominal", recoil, "sr", "fail", mass, read_sumw2=True)
        sr_fail_hbb = rl.TemplateSample(ch_name + "_hbbMC", rl.Sample.BACKGROUND, sr_fail_hbbTemplate)
        addLumiSyst(sr_fail_hbb, year)
        addPileupSyst(sr_fail_hbb)
        addPrefiringSyst(sr_fail_hbb, year)
        addMETSyst(sr_fail_hbb)
        addJESSyst(sr_fail_hbb)
        addMETTrigSyst(sr_fail_hbb, year)
        sr_fail_hbb.setParamEffect(veto_tau, nveto_tau)
        addHbbNormSyst(sr_fail_hbb)
        addBtagSyst(background, recoil, "Hbb", "sr", sr_fail_hbb, "fail", mass)
        sr_fail.addSample(sr_fail_hbb)
    
        sr_fail_qcdTemplate = template(background, "QCD", "nominal", recoil, "sr", "fail", mass, read_sumw2=True)
        sr_fail_qcd = rl.TemplateSample(ch_name + "_qcdMC", rl.Sample.BACKGROUND, sr_fail_qcdTemplate)
        addLumiSyst(sr_fail_qcd, year)
        addPileupSyst(sr_fail_qcd)
        addPrefiringSyst(sr_fail_qcd, year)
        addMETSyst(sr_fail_qcd)
        addJESSyst(sr_fail_qcd)
        addMETTrigSyst(sr_fail_qcd, year)
        sr_fail_qcd.setParamEffect(veto_tau, nveto_tau)
        sr_fail_qcd.setParamEffect(qcdsig_norm, nqcd_norm)
        addBtagSyst(background, recoil, "QCD", "sr", sr_fail_qcd, "fail", mass)
        sr_fail.addSample(sr_fail_qcd)
    
        ###
        # Z+jets data-driven model
        ###
    
        sr_fail_zjetsMCTemplate = template(background, "Z+jets", "nominal", recoil, "sr", "fail", mass, min_value=1., read_sumw2=True)
        sr_fail_zjetsMC = rl.TemplateSample(
            "sr" + year + "fail" + "mass" + mass + "recoil" + str(recoil) + "_zjetsMC",
            rl.Sample.BACKGROUND,
            sr_fail_zjetsMCTemplate
        )
        addMETTrigSyst(sr_fail_zjetsMC, year)
        addVJetsSyst(background, recoil, "Z+jets", "sr", sr_fail_zjetsMC, "fail")

        sr_fail_zjetsObservable = rl.Observable("fjmass"+mass, sr_fail_zjetsMCTemplate[1])
        sr_fail_zjetsBinYields = np.array([rl.IndependentParameter("sr" + year + "fail" + "mass" + mass + "recoil" + str(recoil) + "_zjets_mu"+str(b), 
                                                                   sr_fail_zjetsMCTemplate[0][b], 
                                                                   1e-5, 
                                                                   sr_fail_zjetsMCTemplate[0].max()*2) for b in range(len(sr_fail_zjetsMCTemplate[0]))])

        sr_fail_zjets = rl.ParametericSample(
            "sr" + year + "fail" + "mass" + mass + "recoil" + str(recoil) + "_zjets",
            rl.Sample.BACKGROUND,
            sr_fail_zjetsObservable,
            sr_fail_zjetsBinYields
        )
        sr_fail.addSample(sr_fail_zjets)

        ###
        # W+jets data-driven model
        ###
      
        sr_fail_wjetsMCTemplate = template(background, "W+jets", "nominal", recoil, "sr", "fail", mass, min_value=1., read_sumw2=True)
        sr_fail_wjetsMC = rl.TemplateSample(
            "sr" + year + "fail" + "mass" + mass + "recoil" + str(recoil) + "_wjetsMC",
            rl.Sample.BACKGROUND,
            sr_fail_wjetsMCTemplate
        )
        addMETTrigSyst(sr_fail_wjetsMC, year)
        addVJetsSyst(background, recoil, "W+jets", "sr", sr_fail_wjetsMC, "fail")

        tf, unc = makeTF(sr_fail_wjetsMC, sr_fail_zjetsMC)
        sr_fail_wjets = TransferFactorSample(
            "sr" + year + "fail" + "mass" + mass + "recoil" + str(recoil) + "_wjets",
            rl.Sample.BACKGROUND,
            tf,
            sr_fail._samples[sr_fail_zjets.name],
            nominal_values=sr_fail_wjetsMC._nominal,
            stat_unc=unc#stat_unc=None
        )
        sr_fail.addSample(sr_fail_wjets)

        ###
        # Add BB-lite
        ###

        addBBLiteSyst(sr_fail)

        sr_zjets = sr_fail_zjets
        sr_wjets = sr_fail_wjets
        sr_wjetsMC = sr_fail_wjetsMC

        sr = sr_fail
            
        with open(
                "data/models/"
                + "darkhiggs"
                + "-"
                + year
                + "-"
                + "fail"
                + "-mass"
                + mass
                + "-recoil"
                + str(recoil)
                + ".model",
                "wb",
        ) as fout:
            pickle.dump(add_crs(year, mass, recoil, model_fail, "fail"), fout, protocol=2)

        #######
        #####
        ###
        # Building the "pass" model
        ###  
        #####
        #######

        isttMC = ('40to' in mass and not 'to300' in mass) | (recoil==4)
        iswjetsMC = (recoil==4)

        qcdpho_norm = rl.NuisanceParameter("qcdpho_norm" + year + "pass", "lnN")
        qcde_norm = rl.NuisanceParameter("qcde_norm" + year + "pass", "lnN")
        qcdmu_norm = rl.NuisanceParameter("qcdmu_norm" + year + "pass", "lnN")
        qcdsig_norm = rl.NuisanceParameter("qcdsig_norm" + year + "pass", "lnN")
        st_norm = rl.NuisanceParameter("st_norm" + year + "pass", "lnN")
        ttMC_norm = rl.NuisanceParameter("tt_norm" + year + "pass", "lnN")
        vv_norm = rl.NuisanceParameter("vv_norm" + year + "pass", "lnN")
        hbb_norm = rl.NuisanceParameter("hbb_norm" + year + "pass", "lnN")
        wjetsMC_norm = rl.NuisanceParameter("wjets_norm" + year + "pass", "lnN")
        zjetsMC_norm = rl.NuisanceParameter("zjets_norm" + year + "pass", "lnN")

        model_id = year + "pass" + "mass" + mass+ "recoil" + str(recoil)
        model_pass = rl.Model(model_id)
    

        #####
        ###
        # Signal region "pass"
        ###  
        #####

        ch_name = "sr" + model_id
        sr_pass = rl.Channel(ch_name)
        model_pass.addChannel(sr_pass)
    
        ###
        # Add data distribution to the channel
        ###

        if options.fakedata:
            dataTemplate = template(data, "FakeData", "data", recoil, "sr", "pass", mass)
        else:
            dataTemplate = template(data, "MET", "data", recoil, "sr", "pass", mass)
        sr_pass.setObservation(dataTemplate)
    
        ###
        # Other MC-driven processes
        ###
        
        if iswjetsMC: 
            sr_pass_wjetsTemplate = template(background, "W+jets", "nominal", recoil, "sr", "pass", mass, read_sumw2=True)
            sr_pass_wjets = rl.TemplateSample( "sr" + model_id + "_wjetsMC", rl.Sample.BACKGROUND, sr_pass_wjetsTemplate)
            addLumiSyst(sr_pass_wjets, year)
            addPileupSyst(sr_pass_wjets)
            addPrefiringSyst(sr_pass_wjets, year)
            addMETSyst(sr_pass_wjets)
            addJESSyst(sr_pass_wjets)
            addMETTrigSyst(sr_pass_wjets, year)
            sr_pass_wjets.setParamEffect(veto_tau, nveto_tau)
            addWjetsNormSyst(sr_pass_wjets)
            addBtagSyst(background, recoil, "W+jets", "sr", sr_pass_wjets, "pass", mass)
            addVJetsSyst(background, recoil, "W+jets", "sr", sr_pass_wjets, "pass")
            sr_pass.addSample(sr_pass_wjets)
    
        if isttMC: 
            sr_pass_ttTemplate = template(background, "TT", "nominal", recoil, "sr", "pass", mass, read_sumw2=True)
            sr_pass_tt = rl.TemplateSample("sr" + model_id + "_ttMC",rl.Sample.BACKGROUND, sr_pass_ttTemplate)
            addLumiSyst(sr_pass_tt, year)
            addPileupSyst(sr_pass_tt)
            addPrefiringSyst(sr_pass_tt, year)
            addMETSyst(sr_pass_tt)
            addJESSyst(sr_pass_tt)
            addMETTrigSyst(sr_pass_tt, year)
            sr_pass_tt.setParamEffect(veto_tau, nveto_tau)
            addTTbarNormSyst(sr_pass_tt)
            addBtagSyst(background, recoil, "TT", "sr", sr_pass_tt, "pass", mass)
            sr_pass.addSample(sr_pass_tt)
        
        sr_pass_stTemplate = template(background, "ST", "nominal", recoil, "sr", "pass", mass, read_sumw2=True)
        sr_pass_st = rl.TemplateSample(ch_name + "_stMC", rl.Sample.BACKGROUND, sr_pass_stTemplate)
        addLumiSyst(sr_pass_st, year)
        addPileupSyst(sr_pass_st)
        addPrefiringSyst(sr_pass_st, year)
        addMETSyst(sr_pass_st)
        addJESSyst(sr_pass_st)
        addMETTrigSyst(sr_pass_st, year)
        sr_pass_st.setParamEffect(veto_tau, nveto_tau)
        addSingleTopNormSyst(sr_pass_st)
        addBtagSyst(background, recoil, "ST", "sr", sr_pass_st, "pass", mass)
        sr_pass.addSample(sr_pass_st)
    
        sr_pass_dyjetsTemplate = template(background, "DY+jets", "nominal", recoil, "sr", "pass", mass, read_sumw2=True)
        sr_pass_dyjets = rl.TemplateSample(ch_name + "_dyjetsMC", rl.Sample.BACKGROUND, sr_pass_dyjetsTemplate)
        addLumiSyst(sr_pass_dyjets, year)
        addPileupSyst(sr_pass_dyjets)
        addPrefiringSyst(sr_pass_dyjets, year)
        addMETSyst(sr_pass_dyjets)
        addJESSyst(sr_pass_dyjets)
        addMETTrigSyst(sr_pass_dyjets, year)
        sr_pass_dyjets.setParamEffect(veto_tau, nveto_tau)
        addDrellYanNormSyst(sr_pass_dyjets)
        addBtagSyst(background, recoil, "DY+jets", "sr", sr_pass_dyjets, "pass", mass)
        addVJetsSyst(background, recoil, "DY+jets", "sr", sr_pass_dyjets, "pass")
        sr_pass.addSample(sr_pass_dyjets)
    
        sr_pass_vvTemplate = template(background, "VV", "nominal", recoil, "sr", "pass", mass, read_sumw2=True)
        sr_pass_vv = rl.TemplateSample(ch_name + "_vvMC", rl.Sample.BACKGROUND, sr_pass_vvTemplate)
        addLumiSyst(sr_pass_vv, year)
        addPileupSyst(sr_pass_vv)
        addPrefiringSyst(sr_pass_vv, year)
        addMETSyst(sr_pass_vv)
        addJESSyst(sr_pass_vv)
        addMETTrigSyst(sr_pass_vv, year)
        sr_pass_vv.setParamEffect(veto_tau, nveto_tau)
        addDibosonNormSyst(sr_pass_vv)
        addBtagSyst(background, recoil, "VV", "sr", sr_pass_vv, "pass", mass)
        sr_pass.addSample(sr_pass_vv)
    
        sr_pass_hbbTemplate = template(background, "Hbb", "nominal", recoil, "sr", "pass", mass, read_sumw2=True)
        sr_pass_hbb = rl.TemplateSample(ch_name + "_hbbMC", rl.Sample.BACKGROUND, sr_pass_hbbTemplate)
        addLumiSyst(sr_pass_hbb, year)
        addPileupSyst(sr_pass_hbb)
        addPrefiringSyst(sr_pass_hbb, year)
        addMETSyst(sr_pass_hbb)
        addJESSyst(sr_pass_hbb)
        addMETTrigSyst(sr_pass_hbb, year)
        sr_pass_hbb.setParamEffect(veto_tau, nveto_tau)
        addHbbNormSyst(sr_pass_hbb)
        addBtagSyst(background, recoil, "Hbb", "sr", sr_pass_hbb, "pass", mass)
        sr_pass.addSample(sr_pass_hbb)
    
        sr_pass_qcdTemplate = template(background, "QCD", "nominal", recoil, "sr", "pass", mass, read_sumw2=True)
        sr_pass_qcd = rl.TemplateSample(ch_name + "_qcdMC", rl.Sample.BACKGROUND, sr_pass_qcdTemplate)
        addLumiSyst(sr_pass_qcd, year)
        addPileupSyst(sr_pass_qcd)
        addPrefiringSyst(sr_pass_qcd, year)
        addMETSyst(sr_pass_qcd)
        addJESSyst(sr_pass_qcd)
        addMETTrigSyst(sr_pass_qcd, year)
        sr_pass_qcd.setParamEffect(veto_tau, nveto_tau)
        sr_pass_qcd.setParamEffect(qcdsig_norm, nqcd_norm)
        addBtagSyst(background, recoil, "QCD", "sr", sr_pass_qcd, "pass", mass)
        sr_pass.addSample(sr_pass_qcd)
    
        ###
        # top-antitop data-driven model
        ###
    
        if not isttMC:
            sr_pass_ttTemplate = template(background, "TT", "nominal", recoil, "sr", "pass", mass, min_value=1., read_sumw2=True)
            sr_pass_ttMC = rl.TemplateSample("sr" + model_id + "_ttMC",rl.Sample.BACKGROUND, sr_pass_ttTemplate)
            addMETTrigSyst(sr_pass_ttMC, year)
            addBtagSyst(background, recoil, "TT", "sr", sr_pass_ttMC, "pass", mass)
            
            sr_pass_ttObservable = rl.Observable("fjmass"+mass, sr_pass_ttTemplate[1])
            sr_pass_ttBinYields = np.array([rl.IndependentParameter(ch_name + "_tt_mu"+str(b), sr_pass_ttTemplate[0][b], 
                                                                    1e-5, 
                                                                    sr_pass_ttTemplate[0].max()*2) for b in range(len(sr_pass_ttTemplate[0]))])
            sr_pass_tt = rl.ParametericSample(ch_name + "_tt", rl.Sample.BACKGROUND, sr_pass_ttObservable, sr_pass_ttBinYields)
            sr_pass.addSample(sr_pass_tt)

        if not iswjetsMC:
            sr_pass_wjetsMCTemplate = template(background, "W+jets", "nominal", recoil, "sr", "pass", mass, min_value=1., read_sumw2=True)
            sr_pass_wjetsMC = rl.TemplateSample(
                "sr" + year + "pass" + "mass" + mass + "recoil" + str(recoil) + "_wjetsMC",
                rl.Sample.BACKGROUND,
                sr_pass_wjetsMCTemplate
            )
            addMETTrigSyst(sr_pass_wjetsMC, year)
            addVJetsSyst(background, recoil, "W+jets", "sr", sr_pass_wjetsMC, "pass")

            tf, unc = makeTF(sr_pass_wjetsMC, sr_fail_wjetsMC)
            tf_paramsW = tf * tf_dataResidualW_params[recoil, :]
            #tf_paramsW = wjetseff * tf_MCtemplW_params_final[recoil, :] * tf_dataResidualW_params[recoil, :]
            sr_pass_wjets = TransferFactorSample(
                "sr" + year + "pass" + "mass" + mass + "recoil" + str(recoil) + "_wjets",
                rl.Sample.BACKGROUND,
                tf_paramsW,
                sr_fail._samples[sr_fail_wjets.name],
                nominal_values=sr_pass_wjetsMC._nominal,
                stat_unc=unc#stat_unc=None
            )
            sr_pass.addSample(sr_pass_wjets)

        #####
        ###
        # Z+jets "pass"
        ###  
        ##### 

        sr_pass_zjetsMCTemplate = template(background, "Z+jets", "nominal", recoil, "sr", "pass", mass, min_value=1., read_sumw2=True)
        sr_pass_zjetsMC = rl.TemplateSample(
            "sr" + year + "pass" + "mass" + mass + "recoil" + str(recoil) + "_zjetsMC",
            rl.Sample.BACKGROUND,
            sr_pass_zjetsMCTemplate
        )
        addMETTrigSyst(sr_pass_zjetsMC, year)
        addVJetsSyst(background, recoil, "Z+jets", "sr", sr_pass_zjetsMC, "pass")

        tf, unc = makeTF(sr_pass_zjetsMC, sr_fail_zjetsMC)
        tf_paramsZ = tf * tf_dataResidualZ_params[recoil, :]
        #tf_paramsZ = zjetseff *tf_MCtemplZ_params_final[recoil, :] * tf_dataResidualZ_params[recoil, :]
        sr_pass_zjets = TransferFactorSample(
            "sr" + year + "pass" + "mass" + mass + "recoil" + str(recoil) + "_zjets",
            rl.Sample.BACKGROUND,
            tf_paramsZ,
            sr_fail._samples[sr_fail_zjets.name],
            nominal_values=sr_pass_zjetsMC._nominal,
            stat_unc=unc#stat_unc=None
        )
        sr_pass.addSample(sr_pass_zjets)

        ###
        # Add BB-lite
        ###
    
        addBBLiteSyst(sr_pass)
    
    
        ###
        # Signal
        ###
    
        
        for s in signal["sr"].identifiers("process"):
            sr_pass_signalTemplate = template(signal, s, "nominal", recoil, "sr", "pass", mass, read_sumw2=True)
            sr_pass_signal = rl.TemplateSample(ch_name + "_" + str(s), rl.Sample.SIGNAL, sr_pass_signalTemplate)
            addLumiSyst(sr_pass_signal, year)
            addPileupSyst(sr_pass_signal)
            addPrefiringSyst(sr_pass_signal, year)
            addMETSyst(sr_pass_signal)
            addJESSyst(sr_pass_signal)
            addMETTrigSyst(sr_pass_signal, year)
            sr_pass_signal.setParamEffect(veto_tau, nveto_tau)
            addBtagSyst(signal, recoil, str(s), "sr", sr_pass_signal, "pass", mass)
            #addDoubleBtagSyst(signal, recoil, str(s), "sr", sr_pass_signal, "pass")
            addDoubleBtagSyst(signal, str(s), "sr", sr_pass_signal, "pass")
            sr_pass.addSample(sr_pass_signal)
        
        if not iswjetsMC:
            sr_wjets = sr_pass_wjets
            sr_wjetsMC = sr_pass_wjetsMC
        if not isttMC:
            sr_tt = sr_pass_tt
            sr_ttMC = sr_pass_ttMC
        
        sr = sr_pass
            
        with open(
                "data/models/"
                + "darkhiggs"
                + "-"
                + year
                + "-"
                + "pass"
                + "-mass"
                + mass
                + "-recoil"
                + str(recoil)
                + ".model",
                "wb",
        ) as fout:
            pickle.dump(add_crs(year, mass, recoil, model_pass, "pass"), fout, protocol=2)


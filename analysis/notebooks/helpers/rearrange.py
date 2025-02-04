from __future__ import print_function, division
from collections import defaultdict, OrderedDict
import gzip
import pickle
import json
import os
import uproot
import matplotlib.pyplot as plt
import numpy as np
from coffea import processor 
from coffea.util import load, save
import hist
import mplhep as hep

def make_process_axis(histo):
    new_histo = {}
    for variable in histo.keys():
        if variable=='sumw': continue
        for process in histo[variable]:
            hnew = hist.Hist(
                hist.axis.StrCategory([process], name = 'process'), 
                *(ax for ax in histo[variable][process].axes), 
                    storage=hist.storage.Weight())
            # fill the new histogram with the old one except for the process



            new_histo[variable] = hnew
    
    return new_histo


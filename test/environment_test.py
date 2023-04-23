#!/usr/bin/env python
# coding: utf-8

# # Test of environment setup in SLURM job (remote file access, coffea, pytorch /w GPU support)
# 
# Uses PFNano with all PF candidates, SVs.


import torch
print(torch.cuda.is_available())

import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
import awkward as ak
import uproot
from coffea.nanoevents import NanoEventsFactory, PFNanoAODSchema


events = NanoEventsFactory.from_root(
        'root://grid-cms-xrootd.physik.rwth-aachen.de:1094//store/user/anstein/PFNano/ZH_HToBB_ZToLL_M-125_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1_PFNanoFromMiniV2/221125_132018/0000/nano_mc2017_ULv2_allPF_1-1.root',
        schemaclass=PFNanoAODSchema,
        metadata={"dataset": 'ZH_HToBB_ZToLL'},
    ).events()



gen_part = events.GenPart
#!/usr/bin/env python
# coding: utf-8

import glob
import awkward as ak
from tqdm import tqdm
import os

dir_to_parquet = '/hpcwork/rwth1377/OpenData/OpenNano/ttsemi_tests/custom_arrays_34738307'


parquet_files = sorted(glob.glob(f"{dir_to_parquet}/*.parquet"))


total_jets = 0
for a in tqdm(parquet_files):
    try:
        Njets = len(ak.flatten(ak.from_parquet(a)['Jet_CustomTagger_Cpfcan_puppiw_0'], axis = -1))
        total_jets += Njets
    except:
        print('The bad file:', a)
        os.remove(a) # boo! we only want good files to continue working with, so if there is some problem accessing information, throw it out now
        
print(total_jets)
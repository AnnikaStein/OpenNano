#!/usr/bin/env python
# coding: utf-8

# # Read OpenNano file and store interesting data for machine learning

import argparse
import awkward as ak
import numpy as np
import uproot


print('>>>> Read in file created in earlier step')
file = uproot.open('custom_tt_nanoaod.root')
customtagger_keys = file['Events'].keys(['Jet_CustomTagger*'])
customtagger_keys.extend(['Jet_pt', 'Jet_eta', 'Jet_nBHadrons', 'Jet_nCHadrons', 'Jet_hadronFlavour', 'Jet_partonFlavour', 'Jet_FlavSplit'])
print('<<<< General steps done')



print('>>>> START Variant A')
customtagger_arrays = file['Events'].arrays(customtagger_keys)
ak.to_parquet(customtagger_arrays,'arrays.parquet')
print('<<<< END Variant A: save awkward arrays directly to parquet')





parser = argparse.ArgumentParser(description="Setup input and output pattern")
parser.add_argument('-i',"--fname_in", help="Path to input file to run over")
parser.add_argument('-o',"--fname_out_pattern", help="Output pattern (everything in front of .parquet)")
args = parser.parse_args()

fname_in = args.fname_in
fname_out_pattern = args.fname_out_pattern

print('>>>> READ in nanoaod file created in earlier step')
file = uproot.open(fname_in)
customtagger_keys = file['Events'].keys(['Jet*'])
# customtagger_keys.extend(file['Events'].keys(['Fat*']))
# customtagger_keys.extend(file['Events'].keys(['Muon*']))
# customtagger_keys.extend(file['Events'].keys(['Electron*']))
# customtagger_keys.extend(file['Events'].keys(['Tau*']))
# customtagger_keys.extend(file['Events'].keys(['HLT*']))
# customtagger_keys.extend(file['Events'].keys(['*MET*']))
# customtagger_keys.extend(file['Events'].keys(['PV*']))
# customtagger_keys.extend(file['Events'].keys(['SV*']))
# customtagger_keys.extend(file['Events'].keys(['Photon*']))
# customtagger_keys.extend(file['Events'].keys(['PFCands*']))
# customtagger_keys.extend(file['Events'].keys(['Pileup*']))
# customtagger_keys.extend(file['Events'].keys(['event']))
# customtagger_keys.extend(file['Events'].keys(['luminosityBlock']))
# customtagger_keys.extend(file['Events'].keys(['run']))
customtagger_keys.extend(['event', 'luminosityBlock', 'run'])
print('<<<< General steps done')

print(customtagger_keys)
print('>>>> START actual parquet saving')
customtagger_arrays = file['Events'].arrays(customtagger_keys)
ak.to_parquet(customtagger_arrays,f'{fname_out_pattern}.parquet')
print('<<<< END: save awkward arrays directly to parquet')


print('>>>> START Variant B')
for data_pd in file['Events'].iterate(customtagger_keys, file['Events'].num_entries, library='pd'):
    pass
print(type(data_pd))
for k in range(len(data_pd)):
    data_pd[k].to_pickle(f'{fname_out_pattern}_{k}.pkl')
print('<<<< END Variant B: save pandas dataframe to pkl')
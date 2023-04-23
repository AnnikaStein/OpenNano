#!/usr/bin/env python
# coding: utf-8

# # Read OpenNano file and store interesting data for machine learning


# import torch
# print(torch.cuda.is_available())

# from coffea.nanoevents import NanoEventsFactory, PFNanoAODSchema

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



print('>>>> START Variant B')
for data_pd in file['Events'].iterate(customtagger_keys, file['Events'].num_entries, library='pd'):
    pass
data_pd.to_pickle('df.pkl')
print('<<<< END Variant B: save pandas dataframe to pkl')



print('>>>> START Variant C')
for data in file['Events'].iterate(customtagger_keys, file['Events'].num_entries):
    pass
njets = len(ak.flatten(data['Jet_pt']))
nkeys = len(customtagger_keys)
rectangular = np.zeros((nkeys, njets))
for ind,k in enumerate(customtagger_keys):
    col = ak.to_numpy(ak.flatten(data[k], axis=1))
    rectangular[ind] = col
np.save('arrays.npy',rectangular)
print('<<<< Variant C: save numpy arrays after flattening')
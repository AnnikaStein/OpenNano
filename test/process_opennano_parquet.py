#!/usr/bin/env python
# coding: utf-8

# # Read OpenNano file and store interesting data for machine learning

import argparse
import awkward as ak
import uproot

parser = argparse.ArgumentParser(description="Setup input and output pattern")
parser.add_argument('-i',"--fname_in", help="Path to input file to run over")
parser.add_argument('-o',"--fname_out_pattern", help="Output pattern (everything in front of .parquet)")
args = parser.parse_args()

fname_in = args.fname_in
fname_out_pattern = args.fname_out_pattern

print('>>>> READ in nanoaod file created in earlier step')
file = uproot.open(fname_in)
customtagger_keys = file['Events'].keys(['Jet*'])
customtagger_keys.extend(file['Events'].keys(['Fat*']))
customtagger_keys.extend(file['Events'].keys(['Muon*']))
customtagger_keys.extend(file['Events'].keys(['Electron*']))
customtagger_keys.extend(file['Events'].keys(['Tau*']))
customtagger_keys.extend(file['Events'].keys(['HLT*']))
customtagger_keys.extend(file['Events'].keys(['*MET*']))
customtagger_keys.extend(file['Events'].keys(['PV*']))
customtagger_keys.extend(file['Events'].keys(['SV*']))
customtagger_keys.extend(file['Events'].keys(['Photon*']))
customtagger_keys.extend(file['Events'].keys(['PFCands*']))
customtagger_keys.extend(file['Events'].keys(['Pileup*']))
customtagger_keys.extend(file['Events'].keys(['event']))
customtagger_keys.extend(file['Events'].keys(['luminosityBlock']))
customtagger_keys.extend(file['Events'].keys(['run']))
customtagger_keys.extend(['event', 'luminosityBlock', 'run'])
print('<<<< General steps done')


print('>>>> START actual parquet saving')
customtagger_arrays = file['Events'].arrays(customtagger_keys)
ak.to_parquet(customtagger_arrays,f'{fname_out_pattern}.parquet')
print('<<<< END: save awkward arrays directly to parquet')
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 21 12:06:01 2021

@author: omid
"""

import numpy as np
import pandas as pd

default_rnap = 'rnap'
default_rib = 'rib'

kdeg_enz = np.log(2)/20 # [h-1] # Normal value in ecoli
kdeg_mrna = 60*np.log(2)/5 # Normal value in ecoli

dna_nucleotides = {
                't': 'dttp_c',
                'g': 'dgtp_c',
                'a': 'datp_c',
                'c': 'dctp_c'}

ecoli_biomass_mets = {'lipid_c':'lipid', 
                'protein_c':'protein',
                'carbohydrate_c':'carbohydrate',
                'rna_c':'RNA',  
                'dna_c':'DNA',
                'cofactor_c':'cofactor',
                'ion_c':'ion',}

biomass_composition = {'protein': 0.65640138936118,
                     'RNA': 0.3241955179763448,
                     'DNA': 0.048793783529302,
                     'lipid': 0.09908925307962699,
                     'carbohydrate': 0.1277777399792,
                     'ion': 0.009438329829141999,
                     'cofactor': 0.008786582548793,
                     'total mass': 1.2744825963035888}

### Different valus for basal activity
tcpt_act = 0.01 # Basal value in ecoli
tnsl_act = 1.4e-6 # Basal value in ecoli

# finding average transcription activity in Ecoli (after solving the raw model):
# rnap_footprint_size = 40
# rnap_frac = []
# for gene in model.genes:
#     polysome_size = len(gene.sequence) / rnap_footprint_size
#     n_loci = gene.copy_number
#     scaled_conc = model.variables.DN_DNA.primal
#     try:
#         rnap_usage = model.get_variables_of_type(RNAPUsage).get_by_id(gene.id)
#         scaling_factor = model.dna.scaling_factor / rnap_usage.scaling_factor
#     except KeyError:
#         continue
#     try:
#         the_frac = rnap_usage.variable.primal/(scaling_factor*scaled_conc*polysome_size*n_loci)
#     except ZeroDivisionError:
#         the_frac = 0
#     rnap_frac.append(the_frac)
# average_act = sum(rnap_frac) / len(rnap_frac)
tcpt_ave_act = 0.011498327253910541

# finding average translation activity in Ecoli (after solving the raw model):
# ribo_footprint_size = 60
# rib_frac = []
# for mrna in model.mrnas:
#     polysome_size = len(mrna.gene.rna) / ribo_footprint_size
#     scaled_conc = model.get_variables_of_type(mRNAVariable).get_by_id(mrna.id)
#     try:
#         rib_usage = model.get_variables_of_type(RibosomeUsage).get_by_id(mrna.id)
#         scaling_factor = mrna.scaling_factor / rib_usage.scaling_factor
#     except KeyError:
#         continue
#     try:
#         the_frac = rib_usage.variable.primal/(scaling_factor*scaled_conc.variable.primal*polysome_size)
#     except ZeroDivisionError:
#         the_frac = 0
#     rib_frac.append(the_frac)
# average_act = sum(rib_frac) / len(rib_frac)
tnsl_ave_act = 0.028904996536918463

# determining the minimal activity for lac promoter
# Using Table VI in Peretti & Bailey --> find the ratio of RNAP~p to RNAP~c for copy_number=1
# The minimal activity then can be found by reverse engineering:
# The problem was solved for copy_number=100, adding the constraint 0.05*Ez_rnap - RM_p_beta_lac = 0
# Then, from the relation above, the fraction was found
# lac_tcpt_act = 0.030134586648487262

# determining the minimal activity for lac rib binding site
# Using Table VI in Peretti & Bailey --> find the ratio of Rib~p to Rib~c for copy_number=1
# The minimal activity then can be found by reverse engineering:
# The problem was solved for copy_number=100, adding the constraint 0.17*Ez_rib - RP_p_beta_lac = 0
# lac_tnsl_act = 0.005536975085192169

### Variable basal activity
lac_activity_data = pd.read_excel('../plasmid_data/lac_activity.xlsx', 
                                   header=0, usecols=[0,4,5])
ind = list(np.arange(0,1001,1))
lac_tcpt_act = np.interp(x=ind, xp=lac_activity_data['copy number'],
                        fp=lac_activity_data['rnap_act'])
lac_tnsl_act = np.interp(x=ind, xp=lac_activity_data['copy number'],
                        fp=lac_activity_data['rib_act'])
lac_tcpt_dict = {ind:x for ind,x in enumerate(lac_tcpt_act.tolist())}
lac_tnsl_dict = {ind:x for ind,x in enumerate(lac_tnsl_act.tolist())}


def read_seq(filename):
    with open(filename,'r') as fid:
        all_lines = fid.readlines()

    seq = ''.join([x for line in all_lines for x in line.split()
                   if not x.isdigit()])

    return seq
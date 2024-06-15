#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 21 18:23:59 2021

@author: omid
"""
import numpy as np
import pandas as pd
from etfl.io.json import load_json_model
from etfl.optim.variables import RNAPUsage, RibosomeUsage

COPY_NUMBER = 12


model_ref = 'models/iJO1366_iJO1366_cEFL_2007_enz_128_bins__20210621_133548_{}_copies_pMB1_plasmid.json'
model = load_json_model(model_ref.format(COPY_NUMBER))

gene_id = 'p_beta_Lac'
vector_id = 'pMB1'

phi_p = [0.2, 0.3, 0.4, 0.5]


rnap_usage = model.get_variables_of_type(RNAPUsage)
rib_usage = model.get_variables_of_type(RibosomeUsage)

ind = 0
ret = dict()

for this_phi in phi_p:
    model.change_vector_protein_alloc(this_phi, vector_id)
          
    ##
    sol = model.slim_optimize()
    if np.isnan(sol): # it is infeasible
        continue
    
    growth = model.growth_reaction.flux
    f_rnap = rnap_usage.get_by_id(gene_id).variable.primal/\
        (model.enzymes.rnap.variable.primal + model.enzymes.rnap_rrna.variable.primal)
    f_rib = rib_usage.get_by_id(gene_id).variable.primal/model.enzymes.rib.variable.primal
    ret[ind] = {'phi_p'      : this_phi,
                'growth'     : growth,
                'f_rnap'     : f_rnap,
                'f_rib'      : f_rib
                }
    
    ind = ind+1

data = pd.DataFrame.from_dict(ret,orient='index')
data.to_csv('outputs/growth_fit_{}_copies.csv'.format(COPY_NUMBER))
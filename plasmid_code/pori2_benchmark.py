#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  8 18:59:54 2021

@author: omid
"""
from etfl.io.json import load_json_model

import pandas as pd

vector_id = 'pOri2'

if __name__ == '__main__':
    
    recombinant = load_json_model('models/iJO1366_iJO1366_cEFL_2037_enz_128_bins__20210304_104349_410_copies_pOri2_plasmid.json')
    # Setting the phi_p and NGAM similar to the fitted for pMB1
    recombinant.change_vector_protein_alloc(0.2, vector_id)
    recombinant.reactions.ATPM.lower_bound = 15
        
    copy_number = recombinant.vector_copy_number[vector_id]
    flux_dict = dict()
    
    # In this study, the growth is set to differnt values h-1 and glucose uptake is minimized. Then, the boundary fluxes are investigated.
    upt_fix = {0 : -5.2, 50 : -7, 410 : -6.3} # Wang et al., Table 1
    upt = upt_fix[copy_number]
    recombinant.reactions.EX_glc__D_e.bounds = (upt, 0)
    
    recombinant.slim_optimize()
    flux_dict[copy_number] = {
        'Growth '           : recombinant.growth_reaction.flux,
        'Glucose uptake'    : recombinant.reactions.EX_glc__D_e.flux,
        'O2 uptake'         : recombinant.reactions.EX_o2_e.flux,
        'CO2 production'    : recombinant.reactions.EX_o2_e.flux,
        'Acetate secretion' : recombinant.reactions.EX_ac_e.flux,
        }
    
    data = (pd.DataFrame.from_dict(flux_dict, orient='index'))
    data.to_csv('outputs/pOri2_{}.csv'.format(copy_number))
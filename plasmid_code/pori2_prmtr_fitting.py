#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  8 18:59:54 2021

@author: omid
"""
from etfl.io.json import load_json_model

import pandas as pd
import numpy as np

vector_id = 'pOri2'
ngam = [15, 20, 25, 30, 35, 40]
phi_p = [0.2, 0.25, 0.3, 0.35]

ind = 0 

if __name__ == '__main__':
    recombinant = load_json_model('models/iJO1366_iJO1366_cEFL_2007_enz_128_bins__20210621_133548_410_copies_pOri2_plasmid.json')
    flux_dict = dict()
    
    for this_ngam in ngam:
        for this_phi in phi_p:
            # Setting the phi_p and NGAM similar to the fitted for pMB1
            recombinant.change_vector_protein_alloc(this_phi, vector_id)
            recombinant.reactions.ATPM.lower_bound = this_ngam
                
            copy_number = recombinant.vector_copy_number[vector_id]
            
            
            # In this study, the growth is set to differnt values h-1 and glucose uptake is minimized. Then, the boundary fluxes are investigated.
            upt_fix = {0 : -5.2, 50 : -7, 410 : -6.3} # Wang et al., Table 1
            upt = upt_fix[copy_number]
            recombinant.reactions.EX_glc__D_e.bounds = (upt, 0)
            
            sol = recombinant.slim_optimize()
            if np.isnan(sol):
                continue
            flux_dict[ind] = {
                'NGAM'              : this_ngam,
                'phi_p'             : this_phi,
                'Growth'            : recombinant.growth_reaction.flux,
                'Glucose uptake'    : recombinant.reactions.EX_glc__D_e.flux,
                'O2 uptake'         : recombinant.reactions.EX_o2_e.flux,
                'CO2 production'    : recombinant.reactions.EX_o2_e.flux,
                'Acetate secretion' : recombinant.reactions.EX_ac_e.flux,
                }
            ind += 1
    
    data = (pd.DataFrame.from_dict(flux_dict, orient='index'))
    data.to_csv('outputs/pOri2_{}.csv'.format(copy_number))
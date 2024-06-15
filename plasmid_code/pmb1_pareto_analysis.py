#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  8 18:46:40 2021

@author: omid
"""
from etfl.io.json import load_json_model

from pytfa.optim.utils import symbol_sum

import pandas as pd


def beneficial_vector(model, product_ids = [], biomass_weight = 1):
    '''
    Assume that the product of the vector is beneficial for the 

    Parameters
    ----------
    model : RecomModel
    product_ids : list, optional
       The enzyme IDs to porduce the proteins on the vector. The default is [].
   biomass_weight: float
       The coefficient for growth reaction in the objective function. (1-w) is set for the prodcts.

    Returns
    -------
    None.

    '''
    
    product_weight = 1 - biomass_weight
    
    if product_weight == 0 or len(product_ids) == 0:
        model.objective = model.growth_reaction
    else:
        products = [model.enzymes.get_by_id(id_) for id_ in product_ids]
        expr_terms = []
        for prot in products:
            the_complexation = model.reactions.get_by_id('{}_complex'.format(prot.id))
            # multiplying the MW in kDa and the flux in mmol/gDW/h results in g/gDW/h, similar unit to growth rate
            expr_terms.append(prot.molecular_weight * the_complexation.flux_expression)
        
        expr = product_weight * symbol_sum(expr_terms) + biomass_weight * model.growth_reaction.flux_expression
        model.objective = expr
    
    return

product_ids = ['beta_lac']
vector_id = 'pMB1'
copy_numbers = [12, 122, 408, 700] # From Seo & Bailey's paper

# Pareto front for max biomass and product level with different weights
if __name__ == '__main__':
    for copy_number in copy_numbers:
        pathfile = 'models/iJO1366_iJO1366_cEFL_2007_enz_128_bins__20210621_133548_{}_copies_pMB1_plasmid.json'.format(copy_number)
        recombinant = load_json_model(pathfile)
        # setting the fitted parameters: phi_p = 0.2 and NGAM = 15
        recombinant.change_vector_protein_alloc(0.2, vector_id)
        recombinant.reactions.ATPM.lower_bound = 15
        
        weights = [0.1, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
        results = dict()
        glu_uptake = 0.180 * recombinant.reactions.EX_glc__D_e.lower_bound # g/gDW/h
        for w in weights:
            beneficial_vector(recombinant, product_ids=product_ids, biomass_weight=w)
            recombinant.slim_optimize()
            growth = recombinant.growth_reaction.flux
            product = recombinant.reactions.beta_lac_complex.flux * recombinant.enzymes.beta_lac.molecular_weight
            results[w] ={'biomass yield': growth/glu_uptake,
                         'product yield': product/glu_uptake,
                         'growth'       : growth,
                         'product'      : recombinant.enzymes.beta_lac.variable.primal}
          
        data = abs(pd.DataFrame.from_dict(results, orient='index'))
        data.to_csv('outputs/pareto_front_{}.csv'.format(copy_number))
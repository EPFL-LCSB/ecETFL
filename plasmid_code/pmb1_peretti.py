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
copy_numbers = [12, 24, 60, 122, 408] # From Seo & Bailey's paper

# Pareto front for max biomass and product level with different weights
if __name__ == '__main__':
    results = dict()
    for copy_number in copy_numbers:
        pathfile = 'models/iJO1366_iJO1366_cEFL_2007_enz_128_bins__20210621_133548_{}_copies_pMB1_plasmid.json'.format(copy_number)
        recombinant = load_json_model(pathfile)
        # setting the fitted parameters: phi_p = 0.2 and NGAM = 15
        recombinant.change_vector_protein_alloc(0.2, vector_id)
        recombinant.reactions.ATPM.lower_bound = 15
        
        
        recombinant.objective = recombinant.growth_reaction
        recombinant.slim_optimize()
        growth = recombinant.growth_reaction.flux
        recombinant.growth_reaction.lower_bound = growth
        
        # min-max product
        recombinant.objective = recombinant.enzymes.beta_lac.variable
        recombinant.objective_direction = 'min'
        pr_lb = recombinant.slim_optimize()
        recombinant.objective_direction = 'max'
        pr_ub = recombinant.slim_optimize()
        results[copy_number] ={
                     'growth'       : growth,
                     'product_lb'   : pr_lb,
                     'product_ub'   : pr_ub,
                     }
          
    data = abs(pd.DataFrame.from_dict(results, orient='index'))
    data.to_csv('outputs/minmax_product.csv')
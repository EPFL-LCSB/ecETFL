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

COPY_NUMBER = 112


model_ref = 'models/iJO1366_iJO1366_cEFL_2007_enz_128_bins__20210621_133548_{}_copies_pMB1_plasmid.json'
model = load_json_model(model_ref.format(COPY_NUMBER))

growth_data = pd.read_csv('outputs/sensitivity_prmtrs_{}_copies.csv'.format(COPY_NUMBER))

gene_id = 'p_beta_Lac'
vector_id = 'pMB1'

phi_p = [0.3, 0.4, 0.5, 0.6, 0.7]
w_tcpt = [0.005, 0.01, 0.02, 0.03, 0.04]
w_tnsl = [0.0035, 0.0055, 0.0075, 0.0095, 0.0115, 0.0135]

rnap_usage = model.get_variables_of_type(RNAPUsage)
rib_usage = model.get_variables_of_type(RibosomeUsage)

ind = 0
ret = dict()

for this_phi in phi_p:
    model.change_vector_protein_alloc(this_phi, vector_id)
    for this_w_tcpt in w_tcpt:
        model.change_tcpt_basal_activity(this_w_tcpt, vector_id)
        for this_w_tnsl in w_tnsl:
            model.change_tnsl_basal_activity(this_w_tnsl, vector_id)
            
            # Set the lower bound for the growth to be maximum growth
            try:
                max_mu = float(growth_data[(growth_data.omega_tcpt == this_w_tcpt) & \
                                  (growth_data.omega_tnsl == this_w_tnsl) & \
                                      (growth_data.phi_p == this_phi)]['growth'])
                model.growth_reaction.lower_bound = max_mu
            except TypeError:
                continue
            
            # find the min f_rnap
            model.objective = rnap_usage.get_by_id(gene_id).variable
            model.objective_direction = 'min'
            sol = model.slim_optimize()
            if np.isnan(sol): # it is infeasible
                continue
            growth = model.growth_reaction.flux
            f_rnap = rnap_usage.get_by_id(gene_id).variable.primal/\
                (model.enzymes.rnap.variable.primal + model.enzymes.rnap_rrna.variable.primal)
            f_rib = rib_usage.get_by_id(gene_id).variable.primal/model.enzymes.rib.variable.primal
            ret['min_frnap_{}'.format(ind)] = {'phi_p'      : this_phi,
                        'omega_tcpt' : this_w_tcpt,
                        'omega_tnsl' : this_w_tnsl,
                        'growth'     : growth,
                        'f_rnap'     : f_rnap,
                        'f_rib'      : f_rib
                        }
            
            # find the max f_rib
            model.objective = rib_usage.get_by_id(gene_id).variable
            model.objective_direction = 'max'
            sol = model.slim_optimize()
            if np.isnan(sol): # it is infeasible
                continue
            growth = model.growth_reaction.flux
            f_rnap = rnap_usage.get_by_id(gene_id).variable.primal/\
                (model.enzymes.rnap.variable.primal + model.enzymes.rnap_rrna.variable.primal)
            f_rib = rib_usage.get_by_id(gene_id).variable.primal/model.enzymes.rib.variable.primal
            ret['max_frib_{}'.format(ind)] = {'phi_p'      : this_phi,
                        'omega_tcpt' : this_w_tcpt,
                        'omega_tnsl' : this_w_tnsl,
                        'growth'     : growth,
                        'f_rnap'     : f_rnap,
                        'f_rib'      : f_rib
                        }
            
            
            ind = ind+1

data = pd.DataFrame.from_dict(ret,orient='index')
data.to_csv('outputs/min_max_prmtrs_{}_copies.csv'.format(COPY_NUMBER))
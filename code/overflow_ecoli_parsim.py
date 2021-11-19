from collections import namedtuple
import pandas as pd
import numpy  as np

from etfl.io.json import load_json_model
from etfl.optim.config import standard_solver_config, growth_uptake_config

from etfl.optim.variables import GrowthActivation, BinaryActivator

from time import time
from copy import copy

from etfl.optim.utils import fix_growth, release_growth, \
                            get_active_growth_bounds, safe_optim
from etfl.optim.variables import EnzymeVariable, mRNAVariable
from etfl.optim.constraints import BackwardCatalyticConstraint, \
                                    ForwardCatalyticConstraint
                                    
from pytfa.analysis.chebyshev import chebyshev_center
from etfl.analysis.dynamic import compute_center_generic, compute_center

from etfl.optim.constraints import ModelConstraint
from pytfa.optim.utils import symbol_sum

from cobra.util.solver import set_objective

# flux_to_set = 'growth'
flux_to_set = 'glucose'

solver = 'optlang-gurobi'
#solver = 'optlang-cplex'

GLC_RXN_ID = 'EX_glc__D_e'
GROWTH_RXN_ID = 'BIOMASS_Ec_iJO1366_WT_53p95M'



def _va_sim(model):
    model.objective.direction = 'max'
    sol_max = safe_optim(model)

    model.objective.direction = 'min'
    sol_min = safe_optim(model)

    return sol_min, sol_max




def prep_sol(substrate_uptake, model):

    ret = {'obj':model.solution.objective_value,
           'mu':model.solution.fluxes.loc[model.growth_reaction.id],
           'available_substrate':-1*substrate_uptake,
           'uptake':-1*model.solution.fluxes[GLC_RXN_ID]
           }

#    for exch in model.exchanges:
#        ret[exch.id] = model.solution.fluxes.loc[exch.id]
    for rxn in model.reactions:
        ret[rxn.id] = model.solution.fluxes.loc[rxn.id]
    for enz in model.enzymes:
        ret['EZ_'+ enz.id] = model.solution.raw.loc['EZ_'+enz.id]

    return pd.Series(ret)

if __name__ == '__main__':
    # Do things

    
#    uptake_range = pd.Series(np.arange(-15,-16, -1))
    uptake_range = pd.Series(np.arange(-0,-9,-1))
#    growth_range = pd.Series(np.arange(0.025,0.3,0.025)).append(pd.Series(np.arange(0.3,0.4,0.01)))
    growth_range = pd.Series(np.arange(0.05,0.35,0.05)).append(pd.Series(np.arange(0.35,0.61,0.025)))


    model_files = {
        'cEFL':'iJO1366_cEFL_2037_enz_128_bins__20201103_191050.json',
        'vEFL':'iJO1366_vEFL_2037_enz_128_bins__20201103_190819.json',
        'cETFL':'SlackModel iJO1366_cETFL_2037_enz_128_bins__20201103_190927.json',
        'vETFL':'SlackModel iJO1366_vETFL_2037_enz_128_bins__20201103_190721.json', 
    }

    models = {k:load_json_model('models/'+v,solver=solver) for k,v in model_files.items()}
    data = {}
    sol = pd.Series()
    chebyshev_variables = None
    BIGM = 1000
    

    for name,model in models.items():
        # growth_uptake_config(model)
        model.warm_start = None
        model.logger.info('Simulating ...')
        start = time()
        if chebyshev_variables is None:
                chebyshev_variables = model.get_variables_of_type(EnzymeVariable) + \
                    model.get_variables_of_type(mRNAVariable)
        chebyshev_center(model, chebyshev_variables,
                     inplace=True,
                     big_m=BIGM,
                     include=[ForwardCatalyticConstraint,BackwardCatalyticConstraint],
                     exclude=[])
        
        if flux_to_set == 'growth':
            tol = 0.01    
            model.reactions.get_by_id(GLC_RXN_ID).upper_bound = 0
            model.reactions.get_by_id(GLC_RXN_ID).lower_bound = -1000
            for gr in growth_range:
                # minimize glucose uptake
                model.objective = GLC_RXN_ID
                model.objective_direction = 'max'
                
                model.reactions.get_by_id(GROWTH_RXN_ID).lower_bound = gr
                model.reactions.get_by_id(GROWTH_RXN_ID).upper_bound = gr
                temp_sol = safe_optim(model)
                if np.isnan(temp_sol.objective_value):
                    continue
                # fix substrate uptake
                upt = model.objective.value
                
                chebyshev_sol = compute_center_generic(model) # this also fixes growth with fix_growth function
                new_sol = prep_sol(upt, model)
                sol = pd.concat([sol,new_sol],axis =1)
                # revert the changes
                release_growth(model)
        elif flux_to_set == 'glucose':
            tol = 0.01    
            model.reactions.get_by_id(GROWTH_RXN_ID).upper_bound = 10
            model.reactions.get_by_id(GROWTH_RXN_ID).lower_bound = 0
            for upt in uptake_range:
                # minimize glucose uptake
                model.objective = GROWTH_RXN_ID
                model.objective_direction = 'max'
                
                model.reactions.get_by_id(GLC_RXN_ID).lower_bound = upt
                model.reactions.get_by_id(GLC_RXN_ID).upper_bound = upt
                temp_sol = safe_optim(model)
                # if np.isnan(temp_sol.objective_value):
                #     continue                
                
                # chebyshev_sol = compute_center(model) # this also fixes growth with fix_growth function
                new_sol = prep_sol(upt, model)
                sol = pd.concat([sol,new_sol],axis =1)
                # revert the changes
                release_growth(model)
                
                
        data[name] = sol
        stop = time()
        print('Elapsed time: {}'.format(stop - start))
        data[name].to_csv('outputs/overflow_parsi_cheby_{}.csv'.format(name))

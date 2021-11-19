from collections import namedtuple
import pandas as pd
import numpy  as np

from etfl.io.json import load_json_model
from etfl.optim.config import standard_solver_config, growth_uptake_config

from etfl.optim.variables import GrowthActivation, BinaryActivator, \
                                mRNAVariable, EnzymeVariable

from time import time
from copy import copy

from etfl.optim.utils import fix_growth, release_growth, \
                            get_active_growth_bounds, safe_optim

try:
    from gurobipy import GRB
except ModuleNotFoundError:
    pass

solver = 'optlang-gurobi'
# solver = 'optlang-cplex'


def _va_sim(model):
    model.objective.direction = 'max'
    sol_max = safe_optim(model)

    model.objective.direction = 'min'
    sol_min = safe_optim(model)

    return sol_min, sol_max


def simulate(available_uptake, model, variables, warm_start=None):

#    model.solver.problem.reset()
    model.logger.info('available_uptake = {}'.format(available_uptake))
    model.reactions.EX_glc__D_e.lower_bound = available_uptake
    model.growth_reaction.lower_bound = 0
    model.growth_reaction.upper_bound = 10

    model.objective = model.growth_reaction.id
    model.objective.direction = 'max'

    out = safe_optim(model)

    if model.solver.status == 'infeasible':
        ret = {'obj':np.nan,
               'mu': np.nan,
               'mu_lb':np.nan,
               'mu_ub':np.nan,
               'available_substrate':available_uptake,
               'uptake':np.nan,
               'prot_ratio':np.nan,
               'mrna_ratio':np.nan
               }
        for var in variables:
            ret[var + '_lb'] = np.nan
            ret[var + '_ub'] = np.nan
        print('INFEASIBLE SOLUTION AT q={}'.format(available_uptake))
        return pd.Series(ret)

    growth_solution = copy(model.solution)
    mu_i, mu_lb, mu_ub = get_active_growth_bounds(model)
    mu = model.growth_reaction.flux
    # release_warm_start(model)

    try:
        prot_ratio = model.interpolation_variable.prot_ggdw.variable.primal
        mrna_ratio = model.interpolation_variable.mrna_ggdw.variable.primal
    except AttributeError:
        # Model without neidhardt data
        prot = []
        mRNA = []
        RNAs = model.get_variables_of_type(mRNAVariable)
        Prots = model.get_variables_of_type(EnzymeVariable)
        for the_var in RNAs:
            mRNA += [model.solution.raw.loc[the_var.name]]
        for the_var in Prots:
            prot += [model.solution.raw.loc[the_var.name]]
        prot_ratio = sum(prot)
        mrna_ratio = sum(mRNA)


    ret = {'obj':model.solution.objective_value,
           'mu': mu,
           'mu_lb':mu_lb,
           'mu_ub':mu_ub,
           'available_substrate':-1*available_uptake,
           'uptake':-1*growth_solution.fluxes['EX_glc__D_e'],
           'prot_ratio':prot_ratio,
           'mrna_ratio':mrna_ratio
           }

    fix_growth(model, model.solution)

    for var in variables:
        # THIS WILL DO THE VA ON ETHANOL
        model.objective = model.variables.get(var)

        lb, ub = _va_sim(model)

        ret[var + '_lb'] = lb.objective_value#np.nan
        ret[var + '_ub'] = ub.objective_value#np.nan

    print(pd.Series(ret))

    release_growth(model)
    # apply_warm_start(model, growth_solution)
    
    # Add values of other secretions in the ret dictionnary
    for exch in model.reactions:
        ret[exch.id] = growth_solution.fluxes.loc[exch.id]

    return pd.Series(ret)

if __name__ == '__main__':
    # Do things

    variables = [
#                'EZ_rib',
#                'EZ_rnap',
                # 'EZ_dummy_enzyme',
                # 'MR_dummy_gene',
#                'DM_ac_e',    # acetate exchange
                 ]


    # uptake_range = pd.Series(np.arange(-1,-40, -1))
    # uptake_range = pd.Series(np.arange(-1,-30, -1))
    # uptake_range = pd.Series(np.arange(-1,-15, -0.5))
    # uptake_range = pd.Series(np.arange(-15,-20, -1))
    uptake_range = pd.Series(np.arange(-1/3,-5,-1/3)).append(pd.Series(np.arange(-6,-13,-1)))
    # uptake_range = pd.Series(np.arange(-0,-0.25,-0.25))

    model_files = {
        
#        'cEFL':'iJO1366_cEFL_2037_enz_128_bins__20201103_191050.json',
#        'vEFL':'iJO1366_vEFL_2037_enz_128_bins__20201103_190819.json',
        'cETFL':'SlackModel iJO1366_cETFL_2037_enz_128_bins__20201103_190927.json',
        'vETFL':'SlackModel iJO1366_vETFL_2037_enz_128_bins__20201103_190721.json',
    }

    models = {k:load_json_model('models/'+v,solver=solver) for k,v in model_files.items()}
    data = {}

    for name,model in models.items():
        # growth_uptake_config(model)
        model.warm_start = None
        model.logger.info('Simulating ...')
        start = time()
        data[name] = uptake_range.apply(simulate, args=[model,variables])
        stop = time()
        print('Elapsed time: {}'.format(stop - start))
        data[name].to_csv('outputs/benchmark_{}.csv'.format(name))

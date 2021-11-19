#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#       ecoli. Model
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

import logging

import pandas as pd

from pytfa.io import        load_thermoDB,                    \
                            read_lexicon, annotate_from_lexicon,            \
                            read_compartment_data, apply_compartment_data
from pytfa.optim.relaxation import relax_dgo

from etfl.core import Enzyme, Ribosome, RNAPolymerase, ThermoMEModel, MEModel

from etfl.io.json import save_json_model
from pytfa.utils.logger import get_timestr

from ecoli import   get_model, get_thermo_data, get_coupling_dict, \
                        get_mrna_dict, get_rnap, get_rnap_rrna, get_monomers_dict, \
                        get_nt_sequences, get_ratios, get_neidhardt_data, \
                        get_mrna_metrics, get_enz_metrics, get_trna_metrics,\
                        remove_from_biomass_equation, get_ecoli_gen_stats, \
                        get_essentials,  get_enzyme_fraction, get_rib, \
                        get_aa_sequences, remove_DNA, get_macromole_ratio, \
                        modify_GAM, update_copy_numbers, get_transcription_dict

from etfl.optim.config import standard_solver_config

from optlang.exceptions import SolverError

from multiprocessing import Pool
from etfl.debugging.debugging import relax_catalytic_constraints

from etfl.core.allocation import add_protein_mass_requirement, \
	    add_rna_mass_requirement, add_dna_mass_requirement, \
            add_lipid_mass_requirement, add_carbohydrate_mass_requirement, \
        add_ion_mass_requirement, fix_prot_ratio, fix_RNA_ratio, fix_DNA_ratio, \
            constrain_enzymes, constrain_trna
from etfl.analysis.summary import print_standard_sol

from ecoli_utils import find_fnct_pep

# Costrain further the enzymes just like gecko?
frac_proteome = True
# Constrain further the tRNA?
frac_trna = True

# Run model gen in parallel ?
PARALLEL = False # True

data_dir = '../organism_data/info_ecoli'

#solver='optlang-glpk'
#solver = 'optlang-cplex'
solver = 'optlang-gurobi'

# McCloskey2014 values
glc_uptake = 7.54
glc_uptake_std = 0.56
observed_growth_std = 0.01
observed_growth = 0.57

growth_reaction_id = 'BIOMASS_Ec_iJO1366_WT_53p95M'

def apply_bounds(model, bounds):
    for rid in bounds.index:
        rxn = model.reactions.get_by_id(rid)
        rxn.lower_bound = bounds.loc[rid].min()
        rxn.upper_bound = bounds.loc[rid].max()


def create_model(has_thermo, has_expression, var_allocation, 
				kcat_mode='kmax',
                infer_missing_enz=False,
                additional_enz = None,
				free_rib_ratio=0.2, # Medium-dependent control of the bacterial growth rate
				free_rnap_ratio=0.75, # Medium-dependent control of the bacterial growth rate
				add_displacement = False,
				n_mu_bins = 128):
    #------------------------------------------------------------
    # Initialisation
    #------------------------------------------------------------

    assert has_expression == True

    # this hack works because we are using the solver switch to update the var
    # names in the solver but really we should not do this
    # clean up model.sanitize_varnames
    vanilla_model = get_model('optlang-glpk')
    # D-glucose exchange
    vanilla_model.reactions.EX_glc__D_e.lower_bound = -1 * glc_uptake - glc_uptake_std
    vanilla_model.reactions.EX_glc__D_e.upper_bound = -1 * glc_uptake + glc_uptake_std
    
#    vanilla_model.reactions.r_4046.lower_bound = 0 # non-growth maintenance association
    
    vanilla_model.objective = growth_reaction_id
    fba_sol = vanilla_model.slim_optimize()
    
    
    # we need to load an unmodiifed FBa model and extract the biomass composition
    # this data is used to define the allocation constraints on the model for the basic case
    mass_ratios = get_macromole_ratio()
    modify_GAM(vanilla_model, growth_reaction_id, prot_rxn_id='Protrxn')

    mu_0 = fba_sol
    mu_range = [0, 1.5]
    n_mu_bins = n_mu_bins

    time_str = get_timestr()

    coupling_dict = get_coupling_dict(vanilla_model,
                                      mode = kcat_mode,
                                      atps_name = 'ATPS4rpp',
                                      infer_missing_enz = infer_missing_enz)

    pep_list = find_fnct_pep(coupling_dict)
    # if additional_enz is not None:
        



    aa_dict, rna_nucleotides, rna_nucleotides_mp, dna_nucleotides = get_monomers_dict()
    essentials = get_essentials()

    # Initialize the model
    model_name = 'ETFL' if has_thermo else 'EFL'
    model_name = ('v'+model_name) if var_allocation else ('c'+model_name)
    model_name = (model_name+'_infer') if bool(infer_missing_enz) else model_name

    name = 'iJO1366_{}_{}_enz_{}_bins_{}'.format(
        model_name,
        len(coupling_dict),
        n_mu_bins,
        time_str)

    if has_thermo:

        thermo_data, lexicon, compartment_data = get_thermo_data()

        ecoli = ThermoMEModel(thermo_data,model = vanilla_model,
                              growth_reaction = growth_reaction_id,
                              mu_range = mu_range,
                              n_mu_bins = n_mu_bins,
                              name = name,
                              )
    else:
        ecoli = MEModel(model = vanilla_model,
                        growth_reaction = growth_reaction_id,
                        mu_range = mu_range,
                        n_mu_bins = n_mu_bins,
                        name = name,
                        )
    ecoli.name = name
    ecoli.logger.setLevel(logging.WARNING)
    ecoli.sloppy = False
    # apply_bounds(ecoli,fva)

    ecoli.solver = solver
    standard_solver_config(ecoli, verbose=False)

    if has_thermo:
        # Annotate the cobra_model
        annotate_from_lexicon(ecoli, lexicon)
        apply_compartment_data(ecoli, compartment_data)

        # TFA conversion
        ecoli.prepare()
        
        # if we need to specifically relax a thermo constraint
#        relax_thermo(ecoli)
        
        ecoli.convert(add_displacement = add_displacement)


    mrna_dict = get_mrna_dict(ecoli)
    nt_sequences = get_nt_sequences()
    aa_sequences = get_aa_sequences()
    rnap = get_rnap()
    rnap_rrna = get_rnap_rrna()
    rib = get_rib()
    
    transcription_dict = get_transcription_dict()
    

    # Remove nucleotides and amino acids from biomass reaction as they will be
    # taken into account by the expression
    # in the case of ecoli8 this function does nothing!

    remove_from_biomass_equation(model = ecoli,
                                 nt_dict = rna_nucleotides,
                                 aa_dict = aa_dict,
                                 atp_id=essentials['atp'],
                                 adp_id=essentials['adp'],
                                 pi_id=essentials['pi'],
                                 h2o_id=essentials['h2o'],
                                 h_id=essentials['h'],
                                 )

    ##########################
    ##    MODEL CREATION    ##
    ##########################
            
    ecoli.add_nucleotide_sequences(nt_sequences)
    ecoli.add_peptide_sequences(aa_sequences)
    ecoli.add_essentials(  essentials=essentials,
                           aa_dict=aa_dict,
                           rna_nucleotides=rna_nucleotides,
                           rna_nucleotides_mp = rna_nucleotides_mp
                           )
    ecoli.add_mrnas(mrna_dict.values())
    ecoli.add_ribosome(rib,free_ratio=free_rib_ratio)
    ecoli.add_rnap(rnap, free_ratio=free_rnap_ratio)
    ecoli.add_rnap(rnap_rrna, free_ratio=free_rnap_ratio)
    # should be after adding ribosome, since the rrna genes are replaced
    ecoli.add_transcription_by(transcription_dict)
    
    ecoli.build_expression()
    ecoli.add_enzymatic_coupling(coupling_dict)
    
    
    # Dummy protein and mRNA should be added anyway
    nt_ratios, aa_ratios = get_ratios()
    chromosome_len, gc_ratio = get_ecoli_gen_stats()
    kdeg_mrna, mrna_length_avg  = get_mrna_metrics()
    kdeg_trna, trna_length_avg  = get_trna_metrics()
    kdeg_enz,  peptide_length_avg   = get_enz_metrics()
    ecoli.add_dummies(nt_ratios=nt_ratios,
                          mrna_kdeg=kdeg_mrna,
                          mrna_length=mrna_length_avg,
                          aa_ratios=aa_ratios,
                          enzyme_kdeg=kdeg_enz,
                          peptide_length=peptide_length_avg,
                          transcribed_by=[rnap.id],
                          translated_by=[rib.id])
    
    ecoli.add_dummy_trna(nt_ratios=nt_ratios,
                         trna_kdeg=kdeg_trna,
                         trna_length=trna_length_avg,
                         transcribed_by=[rnap.id])
        
            
    if var_allocation:

        remove_DNA(ecoli) 
        
        mu, rna, prot, dna, lipid, carbohydrate, ion = get_neidhardt_data()

        add_protein_mass_requirement(ecoli, mu, prot)
        add_rna_mass_requirement(ecoli, mu, rna)
        add_dna_mass_requirement(ecoli, mu_values=mu,
                                       dna_rel=dna,
                                       gc_ratio=gc_ratio,
                                       chromosome_len=chromosome_len,
                                       dna_dict=dna_nucleotides,
                                       # ppi='s_0633_c'
                                       )
        add_lipid_mass_requirement(ecoli, lipid_mets = ['lipid_c'],
                                   mass_ratio = mass_ratios['lipid'],
                                   mu_values=mu,
                                   lipid_rel=lipid,
                                   lipid_rxn='Liprxn')
        add_carbohydrate_mass_requirement(ecoli, carbohydrate_mets = ['carbohydrate_c'],
                                   mass_ratio = mass_ratios['carbohydrate'],
                                   mu_values=mu,
                                   carbohydrate_rel=carbohydrate,
                                   carbohydrate_rxn='Carbrxn')
        add_ion_mass_requirement(ecoli, ion_mets = ['ion_c'],
                                   mass_ratio = mass_ratios['ion'],
                                   mu_values=mu,
                                   ion_rel=ion,
                                   ion_rxn='Ionrxn')
    else:
        # constant allocation constraint to fix the sahere of RNA and protein
        fix_prot_ratio(ecoli, mass_ratios['protein'])
        fix_RNA_ratio(ecoli, mass_ratios['RNA'])
        # this is added to have a DNA variable anyway which is needed for the RNAP allocation
        # DNA should not be removed from the biomass in this case
        fix_DNA_ratio(ecoli, mass_ratios['DNA'], gc_ratio=gc_ratio, 
                      chromosome_len=chromosome_len)

    update_copy_numbers(ecoli) 
    
    
    # Need to put after, because dummy has to be taken into account if used, and DNA should be added
    ecoli.populate_expression()
    ecoli.add_trna_mass_balances()
    
    if frac_proteome:
        enz_ratio = get_enzyme_fraction(ecoli.peptides) 
        constrain_enzymes(ecoli, enz_ratio, mass_ratios['protein'])
    if frac_trna:
        constrain_trna(ecoli, 0.15, mass_ratios['RNA']) # in almost all organisms trna ratio is approxiamtely 15%

    ecoli.print_info()
    ecoli.growth_reaction.lower_bound = observed_growth # - 1*observed_growth_std
    need_relax = False

    ecoli.repair()

    try:
        ecoli.optimize()
    except (AttributeError, SolverError):
        need_relax = True

    
    if has_thermo and need_relax:
        # final_model, slack_model, relax_table = relax_dgo(ecoli)
        ecoli.solver.configuration.timeout = 10*3600 # increasing time limit as it's more difficult
        final_model, slack_model, relax_table = relax_dgo(ecoli, in_place = True)
    else:
        final_model = ecoli


    final_model.growth_reaction.lower_bound = 0
    final_model.solver.problem.Params.FeasibilityTol = 1e-9
#    activator_states, final_model, relaxation=relax_catalytic_constraints(final_model, min_growth= 0.5)
    # apply_bounds(ecoli, original_bounds)
    solution = final_model.optimize()
    flux_dict = {'  Glucose uptake    ' : 'EX_glc__D_e',
                 '  Growth            ' : final_model.growth_reaction.id}
    print_standard_sol(final_model, flux_dict = flux_dict)

    filepath = 'models/{}'.format(final_model.name)
    save_json_model(final_model, filepath)

    final_model.logger.info('Build complete for model {}'.format(final_model.name))

    return final_model

if __name__ == '__main__':

    

    models = dict()

    # Models defined by Thermo - Expression - Neidhardt
    model_calls = dict()
    
    
    model_calls[ 'cEFL'  ] = {'has_expression':True,
                              'has_thermo':False,
                              'var_allocation':False,
                              'kcat_mode':'kmax'}
#    model_calls[ 'cETFL' ] = {'has_expression':True,
#                               'has_thermo':True,
#                               'var_allocation':False,
#                               'kcat_mode':'kmax'}
#    model_calls['vEFL'  ] = {'has_expression':True,
#                               'has_thermo':False,
#                               'var_allocation':True,
#                               'kcat_mode':'kmax'}
#    model_calls['vETFL' ] = {'has_expression':True,
#                              'has_thermo':True,
#                              'var_allocation':True,
#                              'kcat_mode':'kmax'}


    if not PARALLEL:
        for mid,mc in model_calls.items():
            models[mid] = create_model(**mc)
    else:
        pool = Pool()

        for mid,mc in model_calls.items():
            def this_callback(result, mid=mid):
                models[mid] = result
            pool.apply_async(create_model, [], mc, callback=this_callback)

        pool.close()
        pool.join()

    print('Completed')


    # Make thermo model
    # make_thermo_model()

    # Save FBA model
    # make_fba_model()

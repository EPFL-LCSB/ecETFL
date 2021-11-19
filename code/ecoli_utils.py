#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 27 10:53:18 2020

@author: omid
"""
from cobra import Reaction, Metabolite
import numpy as np
from os.path import join as pjoin
import pandas as pd
import os

from etfl.core.enzyme import Enzyme
from etfl.utils.parsing import parse_gpr

import sympy


file_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = pjoin(file_dir,'../../organism_data/info_ecoli')

def infer_enzyme_from_gpr(reaction, default_kcat, default_kdeg):
    new_enzymes = list()
    compositions = compositions_from_gpr(reaction)
    for e,composition in enumerate(compositions):
        new_enzyme = Enzyme(id = reaction.id + '_inferred_{}'.format(e),
                            kcat=default_kcat,
                            kdeg=default_kdeg,
                            composition=composition)
        new_enzymes.append(new_enzyme)
    return new_enzymes

def compositions_from_gpr(reaction):
    """
    *Warning*: Use this function only if you have no information on the enzymes.
    Logically parses the GPR to automatically find isozymes ( logical OR )
    and subunits ( logical AND ), and creates the necessary complexation
    reactions: 1 per isozyme, requiring the peptides of each subunit

    :param reaction:
    :type reaction: cobra.Reaction
    :return:
    """

    model = reaction.model

    this_gpr = reaction.gene_reaction_rule

    sym_gpr = parse_gpr(this_gpr)

    if isinstance(sym_gpr, sympy.Symbol):
        # GPR of the type: '(gene0)'
        # Gene <=> Protein
        isozymes = [sym_gpr]
    elif isinstance(sym_gpr, sympy.And):
        # GPR of the type: '(gene0 & gene1)'
        # Subunits of one enzyme
        isozymes = [sym_gpr]
    elif isinstance(sym_gpr, sympy.Or):
        # GPR of the type: '(gene0 | gene1)', '((gene0 & gene1) | gene2)'
        # Two isozymes that are the arguments of the OR
        isozymes = sym_gpr.args

    compositions = []

    for e, this_isozyme in enumerate(isozymes):
        if isinstance(this_isozyme, sympy.And):
            # this is a GPR with several subunits
            peptides = {x.name: 1 \
                        for x in this_isozyme.args}
        elif isinstance(this_isozyme, sympy.Symbol):
            # there is only one subunit
            peptides = {this_isozyme.name: 1}
        else:
            # The GPR has been incorrectly parsed
            model.logger.error('Incorrect parsing of {}'.format(isozymes))
            raise TypeError

        compositions += [peptides]

    return compositions

def lump_biomass(model):
    # a function to lump Biomass Building Blocks based on different types of 
    # macromolecules, i.e. protein, carbohydrate, ... .
    biomass_rxn = model.reactions.BIOMASS_Ec_iJO1366_WT_53p95M
    # Categorizing the metabolites in biomass reaction
    ppi='ppi_c'
    
    aa_dict = {'A': 'ala__L_c',
               'R': 'arg__L_c',
               'N': 'asn__L_c',
               'D': 'asp__L_c',
               # 'B':'asx',
               'C': 'cys__L_c',
               'E': 'glu__L_c',
               'Q': 'gln__L_c',
               # 'Z':'glx',
               'G': 'gly_c',
               'H': 'his__L_c',
               'I': 'ile__L_c',
               'L': 'leu__L_c',
               'K': 'lys__L_c',
               'M': 'met__L_c',
               'F': 'phe__L_c',
               'P': 'pro__L_c',
               'S': 'ser__L_c',
               'T': 'thr__L_c',
               'W': 'trp__L_c',
               'Y': 'tyr__L_c',
               'V': 'val__L_c', }
    
    rna_nucleotides = {
                'u': 'utp_c',
                'g': 'gtp_c',
                'a': 'atp_c',
                'c': 'ctp_c'}
    
    dna_nucleotides = {
                't': 'dttp_c',
                'g': 'dgtp_c',
                'a': 'datp_c',
                'c': 'dctp_c'}
    
    cofactors = ['10fthf_c','2dmmql8_c', '5mthf_c', 'accoa_c', 'adocbl_c',
                 'amet_c', 'btn_c', 'bmocogdp_c', 'chor_c','coa_c','fad_c', 
                 'enter_c','gthrd_c','hemeO_c','lipopb_c', 'malcoa_c','mlthf_c',
                 'mococdp_c', 'mocogdp_c', 'mql8_c','nad_c', 'nadh_c','nadp_c',
                 'nadph_c', 'pheme_c','q8h2_c', 'ribflv_c', 'ptrc_c', 'pydx5p_c',
                 'sheme_c', 'succoa_c' , 'thf_c', 'udcpdp_c', 'thmpp_c', 'spmd_c',
                 ]
    
    carbohydrates = ['glycogen_c', 'murein3p3p_p', 'murein4p4p_p', 'murein3px4p_p',
                     'murein4px4p_p', 'murein4px4px4p_p',]
    
    lipids = ['clpn160_p','clpn161_p', 'clpn181_p','colipa_e','pe160_c', 'pe160_p',
              'pe161_c','pe161_p', 'pe181_c', 'pe181_p','pg160_c', 'pg160_p',
              'pg161_c','pg161_p', 'pg181_c', 'pg181_p']
    
    ions = ['2fe2s_c', '4fe4s_c','cl_c','ca2_c','cobalt2_c', 'cu2_c','fe2_c',
            'fe3_c','k_c', 'mg2_c', 'mn2_c', 'mobd_c', 'nh4_c','ni2_c','so4_c',
            'zn2_c']
    
    # Removing GAM from stoichiometric coefficient of ATP, later it will be added again
    adp = model.metabolites.get_by_id('adp_c')
    biomass_rxn.subtract_metabolites({
        rna_nucleotides['a']: - biomass_rxn.metabolites[adp]
            })
    
    # Defining the lumped metabolites and lumping reactions
    ##
    prot = Metabolite('protein_c',
                       formula='',
                       name='Protein',
                       compartment='c')
    prot_rxn = Reaction('Protrxn')
    model.add_reactions([prot_rxn])
    prot_rxn.add_metabolites({prot:1})
    biomass_rxn.add_metabolites({prot:-1})
    # adding metabolites to the lumping rxn
    prot_rxn.add_metabolites({k:v for k,v  \
          in biomass_rxn.metabolites.items() \
              if k.id in aa_dict.values()})
    # removing the metabolites from biomass rxn
    biomass_rxn.subtract_metabolites({k:v for k,v  \
          in biomass_rxn.metabolites.items() \
              if k.id in aa_dict.values()})
    ##
    rna = Metabolite('rna_c',
                       formula='',
                       name='RNA',
                       compartment='c')
    rna_rxn = Reaction('RNArxn')
    model.add_reactions([rna_rxn])
    rna_rxn.add_metabolites({rna:1})
    biomass_rxn.add_metabolites({rna:-1})
    # removing part of ppi from biomass and add it to the lumping reaction
    tot_mol_rna = abs(sum([v for k,v  \
          in biomass_rxn.metabolites.items() \
              if k.id in rna_nucleotides.values()]))
    rna_rxn.add_metabolites({ppi:tot_mol_rna})
    biomass_rxn.subtract_metabolites({ppi:tot_mol_rna})
    # adding metabolites to the lumping rxn
    rna_rxn.add_metabolites({k:v for k,v  \
          in biomass_rxn.metabolites.items() \
              if k.id in rna_nucleotides.values()})
    # removing the metabolites from biomass rxn
    biomass_rxn.subtract_metabolites({k:v for k,v  \
          in biomass_rxn.metabolites.items() \
              if k.id in rna_nucleotides.values()})
    ##
    dna = Metabolite('dna_c',
                       formula='',
                       name='DNA',
                       compartment='c')
    dna_rxn = Reaction('DNArxn')
    model.add_reactions([dna_rxn])
    dna_rxn.add_metabolites({dna:1})
    biomass_rxn.add_metabolites({dna:-1})
    # removing part of ppi from biomass and add it to the lumping reaction
    tot_mol_dna = abs(sum([v for k,v  \
          in biomass_rxn.metabolites.items() \
              if k.id in dna_nucleotides.values()]))
    dna_rxn.add_metabolites({ppi:tot_mol_dna})
    biomass_rxn.subtract_metabolites({ppi:tot_mol_dna})
    # adding metabolites to the lumping rxn
    dna_rxn.add_metabolites({k:v for k,v  \
          in biomass_rxn.metabolites.items() \
              if k.id in dna_nucleotides.values()})
    # removing the metabolites from biomass rxn
    biomass_rxn.subtract_metabolites({k:v for k,v  \
          in biomass_rxn.metabolites.items() \
              if k.id in dna_nucleotides.values()})
    ##
    lip = Metabolite('lipid_c',
                       formula='',
                       name='Lipid',
                       compartment='c')
    lip_rxn = Reaction('Liprxn')
    model.add_reactions([lip_rxn])
    lip_rxn.add_metabolites({lip:1})
    biomass_rxn.add_metabolites({lip:-1})
    # adding metabolites to the lumping rxn
    lip_rxn.add_metabolites({k:v for k,v  \
          in biomass_rxn.metabolites.items() \
              if k.id in lipids})
    # removing the metabolites from biomass rxn
    biomass_rxn.subtract_metabolites({k:v for k,v  \
          in biomass_rxn.metabolites.items() \
              if k.id in lipids})
    ##
    carb = Metabolite('carbohydrate_c',
                       formula='',
                       name='Carbohydrate',
                       compartment='c')
    carb_rxn = Reaction('Carbrxn')
    model.add_reactions([carb_rxn])
    carb_rxn.add_metabolites({carb:1})
    biomass_rxn.add_metabolites({carb:-1})
    # adding metabolites to the lumping rxn
    carb_rxn.add_metabolites({k:v for k,v  \
          in biomass_rxn.metabolites.items() \
              if k.id in carbohydrates})
    # removing the metabolites from biomass rxn
    biomass_rxn.subtract_metabolites({k:v for k,v  \
          in biomass_rxn.metabolites.items() \
              if k.id in carbohydrates})
    ##
    cofac = Metabolite('cofactor_c',
                       formula='',
                       name='Cofactor',
                       compartment='c')
    cofac_rxn = Reaction('Cofacrxn')
    model.add_reactions([cofac_rxn])
    cofac_rxn.add_metabolites({cofac:1})
    biomass_rxn.add_metabolites({cofac:-1})
    # adding metabolites to the lumping rxn
    cofac_rxn.add_metabolites({k:v for k,v  \
          in biomass_rxn.metabolites.items() \
              if k.id in cofactors})
    # removing the metabolites from biomass rxn
    biomass_rxn.subtract_metabolites({k:v for k,v  \
          in biomass_rxn.metabolites.items() \
              if k.id in cofactors})
    ##
    ion = Metabolite('ion_c',
                       formula='',
                       name='Ion',
                       compartment='c')
    ion_rxn = Reaction('Ionrxn')
    model.add_reactions([ion_rxn])
    ion_rxn.add_metabolites({ion:1})
    biomass_rxn.add_metabolites({ion:-1})
    # adding metabolites to the lumping rxn
    ion_rxn.add_metabolites({k:v for k,v  \
          in biomass_rxn.metabolites.items() \
              if k.id in ions})
    # removing the metabolites from biomass rxn
    biomass_rxn.subtract_metabolites({k:v for k,v  \
          in biomass_rxn.metabolites.items() \
              if k.id in ions})
    # Adding GAM to stoichiometric coefficient of ATP
    biomass_rxn.add_metabolites({
        rna_nucleotides['a']: - biomass_rxn.metabolites[adp]
            })
    return

def removing_excessive_mets_rxns(model):
    
    pseudo_rxns = ['Protrxn','RNArxn',]
    pseudo_mets = ['protein_c','rna_c',]
    
    rxn_list = pseudo_rxns
    met_list = pseudo_mets
    
    rxn_set = [model.reactions.get_by_id(x) for x in rxn_list]
    met_set = [model.metabolites.get_by_id(x) for x in met_list]
    
    for rxn in rxn_set:
        model.remove_reactions(rxn)
    for met in met_set:
        model.remove_metabolites(met)
    return

def find_bbb_ratio(model, BBB = ['all']):
    '''
    This function is designed for iJO1366!
    
    A function to find the mass ratio of each biomass building block in FBA biomass definition.
    Also can return total weight of biomass
    inputs:
       model: a cobra model before any modification
       BBB: A list of the building block(s) of interest [lipid, carbohydrate, protein, RNA, DNA, ion, cofactor, {all}]
    '''
    
    growth_rxn = 'BIOMASS_Ec_iJO1366_WT_53p95M'
    pseudo_rxns = {'protein':'Protrxn',
    'carbohydrate':'Carbrxn',
    'RNA':'RNArxn',
    'DNA':'DNArxn',
    'lipid':'Liprxn',
    'cofactor':'Cofacrxn',            
    'ion':'Ionrxn',}
    
    pseudo_mets = {'lipid_c':'lipid', 
                     'protein_c':'protein',
                     'carbohydrate_c':'carbohydrate',
                     'rna_c':'RNA',  
                     'dna_c':'DNA',
                     'cofactor_c':'cofactor',
                     'ion_c':'ion',}
    
    ratios = dict()
    
    if 'all' in BBB:
        BBB = ['lipid', 'carbohydrate', 'protein', 'RNA', 'DNA', 'ion', 'cofactor','all']
    if 'protein' in BBB:
        pseudo_id = pseudo_rxns['protein']
        pseudo_rxn = model.reactions.get_by_id(pseudo_id)
        # the protein psedoreaction is junk. I should find AAs first.
        AA_met = {'A': 'ala__L_c',
               'R': 'arg__L_c',
               'N': 'asn__L_c',
               'D': 'asp__L_c',
               # 'B':'asx',
               'C': 'cys__L_c',
               'E': 'glu__L_c',
               'Q': 'gln__L_c',
               # 'Z':'glx',
               'G': 'gly_c',
               'H': 'his__L_c',
               'I': 'ile__L_c',
               'L': 'leu__L_c',
               'K': 'lys__L_c',
               'M': 'met__L_c',
               'F': 'phe__L_c',
               'P': 'pro__L_c',
               'S': 'ser__L_c',
               'T': 'thr__L_c',
               'W': 'trp__L_c',
               'Y': 'tyr__L_c',
               'V': 'val__L_c', }
        AA_MWs = {k:model.metabolites.get_by_id(v).formula_weight \
                  for k,v in AA_met.items()}
        AA_stoic = {k:pseudo_rxn.get_coefficient(v) for k,v in AA_met.items()}
        prot_ratio = sum([-v * AA_MWs[k] for k,v in AA_stoic.items()])/1000
        ratios['protein'] = prot_ratio
    if 'RNA' in BBB:
        pseudo_id = pseudo_rxns['RNA']
        pseudo_rxn = model.reactions.get_by_id(pseudo_id)
        NMP_met = {'u': 'utp_c',
                'g': 'gtp_c',
                'a': 'atp_c',
                'c': 'ctp_c'}
        NMP_MWs = {k:model.metabolites.get_by_id(v).formula_weight \
                  for k,v in NMP_met.items()}
        NMP_stoic = {k:pseudo_rxn.get_coefficient(v) for k,v in NMP_met.items()}
        RNA_ratio = sum([-v * NMP_MWs[k] for k,v in NMP_stoic.items()])/1000
        ratios['RNA'] = RNA_ratio
    if 'DNA' in BBB:
        pseudo_id = pseudo_rxns['DNA']
        pseudo_rxn = model.reactions.get_by_id(pseudo_id)
        dNMP_met = {
                't': 'dttp_c',
                'g': 'dgtp_c',
                'a': 'datp_c',
                'c': 'dctp_c'}
        dNMP_MWs = {k:model.metabolites.get_by_id(v).formula_weight \
                  for k,v in dNMP_met.items()}
        dNMP_stoic = {k:pseudo_rxn.get_coefficient(v) for k,v in dNMP_met.items()}
        DNA_ratio = sum([-v * dNMP_MWs[k] for k,v in dNMP_stoic.items()])/1000
        ratios['DNA'] = DNA_ratio
    if 'lipid' in BBB:
        pseudo_id = pseudo_rxns['lipid']
        pseudo_rxn = model.reactions.get_by_id(pseudo_id)
        LC_MWs = {met.name:met.formula_weight for met in pseudo_rxn.metabolites.keys()}
        LC_stoic = {met.name:pseudo_rxn.get_coefficient(met.id) \
                       for met in pseudo_rxn.metabolites.keys()}
        lipid_ratio = sum([-v * LC_MWs[k] for k,v in LC_stoic.items() if v<0])/1000 
        ratios['lipid'] = lipid_ratio
    if 'carbohydrate' in BBB:
        pseudo_id = pseudo_rxns['carbohydrate']
        pseudo_rxn = model.reactions.get_by_id(pseudo_id)
        sugar_MWs = {met.name:met.formula_weight for met in pseudo_rxn.metabolites.keys()}
        sugar_stoic = {met.name:pseudo_rxn.get_coefficient(met.id) \
                       for met in pseudo_rxn.metabolites.keys()}
        carbohydrate_ratio = sum([-v * sugar_MWs[k] for k,v in sugar_stoic.items() if v<0])/1000
        ratios['carbohydrate'] = carbohydrate_ratio
    if 'ion' in BBB:
        pseudo_id = pseudo_rxns['ion']
        pseudo_rxn = model.reactions.get_by_id(pseudo_id)
        ions_MWs = {met.name:met.formula_weight for met in pseudo_rxn.metabolites.keys()}
        ions_stoic = {met.name:pseudo_rxn.get_coefficient(met.id) \
                       for met in pseudo_rxn.metabolites.keys()}
        ion_ratio = sum([-v * ions_MWs[k] for k,v in ions_stoic.items() if v<0])/1000
        ratios['ion'] = ion_ratio
    if 'cofactor' in BBB:
        pseudo_id = pseudo_rxns['cofactor']
        pseudo_rxn = model.reactions.get_by_id(pseudo_id)
        cofactor_MWs = {met.name:met.formula_weight for met in pseudo_rxn.metabolites.keys()}
        cofactor_stoic = {met.name:pseudo_rxn.get_coefficient(met.id) \
                       for met in pseudo_rxn.metabolites.keys()}
        cofactor_ratio = sum([-v * cofactor_MWs[k] for k,v in cofactor_stoic.items() if v<0])/1000
        ratios['cofactor'] = cofactor_ratio
    
    pseudo_rxn = model.reactions.get_by_id(growth_rxn)
    for k,v in pseudo_mets.items():
        ratios[v] *= abs(pseudo_rxn.get_coefficient(k))
    if 'all' in BBB:
        pseudo_rxn = model.reactions.get_by_id(growth_rxn)
        ratios['total mass'] = sum([x for x in ratios.values()])
    
    return ratios
    
def stoich_coeff_from_mass_ratio(model, alloc_data, macromole_ratio, \
                                 growth_id, biomass_comps, met_macromol_dict, \
                                 lumped_comp = True):
    '''
    This is a function to find new stoichiometric coefficients for biomass reaction
    based on experimental data for mass ratios.
    
    inputs:
        model: a cobra or pytfa model (not modified)
        alloc_data: a dict or DataFrame of different biomass components with
            their experimental mass ratios
        macromole_ratio: dict of mass ratios of different biomass components in
            GEM.
        growth_id: the rxn id for growth.
        biomass_comps: a list of biomass components (metabolite ids),
        lumped_comp: indicates if biomass metabolites are lumped into pseudometabolites
        met_macromol_dict: adictionary to relate metabolite ids with each type
            of macromolecules, e.g. protein, RNA, etc. Keys are met ids and values
            are macromolecule names compatible with tags in alloc_data & macromole_ratio
        
    ouputs:
        new model with modified stoichiometric coefficients in its growth reaction.
    '''
    
    # finding the original stoichiometric coefficients for each metabolite
    biomass_rxn = model.reactions.get_by_id(growth_id)
    
    if lumped_comp:
        # nothing to do :)
        org_coeffs = {x : biomass_rxn.get_coefficient(x) for x in biomass_comps}
    else:
        # first, we should lump different metabolites to relate them to macromolecules
        raise NotImplementedError()
        
    try:
        tot_mass = macromole_ratio['total mass']
    except KeyError:
        tot_mass = sum([x for x in macromole_ratio.values()])
        
    ratios = {k : v/tot_mass for k,v in macromole_ratio.items()}
    new_coeffs = {k : v * alloc_data[met_macromol_dict[k]] /\
                  ratios[met_macromol_dict[k]] for k,v in org_coeffs.items()}
    
    # the new stoichiometric coefficient is added to the previous one
    # to avoid redundancy I should first subtract the old coefficient
    change_coeffs = {k : v - org_coeffs[k] for k,v in new_coeffs.items()}
    biomass_rxn.add_metabolites(change_coeffs)
    
    return 


def get_enzyme_rna_dist():
    # the constant: tot_prot / tot_rna * K (see the presentation)
    # tot_prot in the cell is between 3e+6 and 4e+6 BNID 115702
    # tot_rna in the cell is 2400 BNID 112795
    constant = (3500000/2400) * (1/8275)
    # a function to find the fraction of proteome that is enzyme
    # peptides is the list of model peptides 
    # data for each related protein abundance from PaxDB
    prot_abundance =  pd.read_excel(pjoin(data_dir,'abundance_table.xlsx'),
                                   skiprows=range(0,11),
                                   header = 0,
                                   )
    # data for each related mrna abundance from PaxDB
    mrna_abundance =  pd.read_excel(pjoin(data_dir,'abundance_table_rna.xlsx'),
                                   usecols=[0,1],
                                   header = 0,
                                   )
    # data for each protein MWs (including enzymes) from KEGG
    prot_mws = pd.read_csv(pjoin(data_dir,'all_peptides_mw.csv'),header=0)
    # data for each protein MWs (including enzymes) from KEGG
    mrna_mws = pd.read_csv(pjoin(data_dir,'all_mrnas_mw.csv'),header=0)
    
    weight_dict_prot = dict()
    for _, pair in prot_abundance.iterrows():
        gene_id = pair[0]
        abundance = pair[1]
        weight_dict_prot[gene_id] = [y * abundance for y in \
                 prot_mws[prot_mws['gene name']==gene_id]['MW']]
    tot_prot = sum([x[0] for x in weight_dict_prot.values() if len(x)>0 and \
                   not np.isnan(x[0])])
    
    weight_dict_mrna = dict()
    for _, pair in mrna_abundance.iterrows():
        gene_id = pair[0]
        abundance = pair[1]
        weight_dict_mrna[gene_id] = [y * abundance for y in \
                 mrna_mws[mrna_mws['gene name']==gene_id]['MW']]
    tot_mrna = sum([x[0] for x in weight_dict_mrna.values() if len(x)>0 and \
                   not np.isnan(x[0])])
    
    frac_dict = dict()
    for k in weight_dict_mrna.keys():
        try:
            absol_prot = weight_dict_prot[k][0]/tot_prot
            absol_mrna = weight_dict_mrna[k][0]/tot_mrna
            if np.isclose(absol_mrna, 0):
                continue
            else:
                frac_dict[k] = constant * absol_prot/absol_mrna
        except KeyError:
            continue
    return frac_dict

def find_fnct_pep(coupling_dict):
    # to find which peptides are used to define the enzymes, i.e. associated to a function
    pep_list = []
    for _,v in coupling_dict.items():
        # v is a list of enzyme objects
        for enz in v:
            for pep in enz.composition.keys():
                if pep not in pep_list:
                    pep_list.append(pep)
    
    return pep_list
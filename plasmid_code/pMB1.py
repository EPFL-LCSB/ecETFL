#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 20 18:28:53 2021

@author: omid
"""
from plasmid_support import read_seq, kdeg_enz, kdeg_mrna, ecoli_biomass_mets, \
                lac_tcpt_dict, lac_tnsl_dict, default_rib, default_rnap, \
                    biomass_composition, dna_nucleotides

from etfl.core.genes import ExpressedGene, CodingGene
from etfl.core.rna import mRNA
from etfl.core.enzyme import Enzyme
from etfl.core.vector import Plasmid

from etfl.io.json import load_json_model, save_json_model
from etfl.core.recombmodel import RecombModel

COPY_NUMBER = 1


def get_pMB1():
    '''
    

    Returns
    -------
    vector 

    '''
    
    ### Defining the genes
    RNA1 = ExpressedGene(id = 'p_RNA_I',
                        name = 'RNA_I',
                        sequence = read_seq('../plasmid_data/RNA_1.txt')) # small RNA
    
    B_LAC = CodingGene(id = 'p_beta_Lac',
                        name = 'beta Lactamase',
                        sequence = read_seq('../plasmid_data/B_LAC.txt'))
    
    n_gene = COPY_NUMBER
    RNA1.min_tcpt_activity = lac_tcpt_dict[n_gene] # As a function of the gene copy number
    B_LAC.min_tcpt_activity = lac_tcpt_dict[n_gene] # As a function of the gene copy number
    B_LAC.min_tnsl_activity = lac_tnsl_dict[n_gene] # As a function of the gene copy number
    
    gene_list = [RNA1, B_LAC]
    
    ### Defining the proteins
    B_LAC_enz = Enzyme(
        id = 'beta_lac',
        kcat=0, # it is here a protein and not an enzyme
        kdeg=kdeg_enz,
        composition = {'p_beta_Lac':1}
    )
    protein_list = [B_LAC_enz]
    
    ### Defining the reactions and coupling
    rxn_list = [] # No reaction
    coupling_dict = {} 
    
    
    ### Defining the plasmid
    # the molecular weight of the vector is 5.5 MDa (Seo & Bailey, Table III)
    # Since the molecular weight is important at the end, the gc and length was approximated
    # based on the following equation to be compatible with the molecular weight:
    # Mw = ((1 - gc) * (ma + mt) + g * (mc + mg)) * length --> gc = 0.5, length = 7647 bp
    vector = Plasmid(id_ = 'pMB1',
                         length = 7647,
                         gc_ratio = 0.5,
                         genes = gene_list,
                         reactions = rxn_list,
                         coupling_dict = coupling_dict,
                         proteins = protein_list)
    
    vector.default_rib = default_rib
    vector.default_rnap = default_rnap
    vector.calibration_tcpt = lac_tcpt_dict
    vector.calibration_tnsl = lac_tnsl_dict
    
    ### Defining the mRNAs
    vector.build_default_mrna(kdeg_mrna)
    
    return vector

    
if __name__ == '__main__':
    
    model = load_json_model('../code/models/iJO1366_cEFL_2007_enz_128_bins__20210621_104409.json')
    recombinant = RecombModel(model, inplace=True)
    recombinant.biomass_metabolites  = ecoli_biomass_mets
    recombinant.biomass_composition = biomass_composition
    try:
        recombinant.dna_nucleotides
    except AttributeError:
        recombinant.dna_nucleotides = dna_nucleotides
        
    pMB1 = get_pMB1()
    pMB1.model = recombinant # This is needed for the allocations inside the add_vector
    recombinant.add_vector(pMB1, COPY_NUMBER)
    
    # saving
    name = 'iJO1366_{}_{}_copies_{}_plasmid'.format(
        recombinant.name,
        COPY_NUMBER,
        pMB1.id
        )
    filepath = 'models/{}'.format(name) 
    
    save_json_model(recombinant, filepath)
    # In this study the growth rate and the product level is benchmarked
    recombinant.slim_optimize()
    growth = recombinant.growth_reaction.flux
    product = recombinant.enzymes.beta_lac.variable.primal
    
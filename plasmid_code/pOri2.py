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


def get_pOri2():
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
                        sequence = read_seq('../plasmid_data/B_LAC.txt')) # as an ampicilin resistance gene
    
    ROM = CodingGene(id = 'p_rom',
                        name = 'RNA I inhibition modulator',
                        sequence = read_seq('../plasmid_data/ROM.txt')) # a regulatory protein
    
    n_gene = COPY_NUMBER
    RNA1.min_tcpt_activity = lac_tcpt_dict[n_gene] # As a function of the gene copy number
    B_LAC.min_tcpt_activity = lac_tcpt_dict[n_gene] # As a function of the gene copy number
    B_LAC.min_tnsl_activity = lac_tnsl_dict[n_gene] # As a function of the gene copy number
    ROM.min_tcpt_activity = lac_tcpt_dict[n_gene] # As a function of the gene copy number
    ROM.min_tnsl_activity = lac_tnsl_dict[n_gene] # As a function of the gene copy number
    
    gene_list = [RNA1, B_LAC, ROM]
    
    ### Defining the proteins
    B_LAC_enz = Enzyme(
        id = 'beta_lac',
        kcat=0, # it is here a protein and not an enzyme
        kdeg=kdeg_enz,
        composition = {'p_beta_Lac':1}
    )
    
    Rop = Enzyme(
        id = 'rop',
        kcat=0, # it is here a protein and not an enzyme
        kdeg=kdeg_enz,
        composition = {'p_rom':1}
    )
    
    protein_list = [B_LAC_enz, Rop]
    
    ### Defining the reactions and coupling
    rxn_list = [] # No reaction
    coupling_dict = {} 
    
    
    ### Defining the plasmid
    # The length of 2 versions of this plasmid is 4575 bp based on:
    # Effects of the presence of ColE1 plasmid DNA in Escherichia coli on the host cell metabolism, Wang et al.
    # GC ration is 0.506 based on Ow et al.
    vector = Plasmid(id_ = 'pOri2',
                         length = 4575,
                         gc_ratio = 0.506,
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
    
    model = load_json_model('../code/models/iJO1366_cEFL_2037_enz_128_bins__20210304_104349.json')
    recombinant = RecombModel(model, inplace=True)
    recombinant.biomass_metabolites  = ecoli_biomass_mets
    recombinant.biomass_composition = biomass_composition
    try:
        recombinant.dna_nucleotides
    except AttributeError:
        recombinant.dna_nucleotides = dna_nucleotides
        
    pOri2 = get_pOri2()
    pOri2.model = recombinant # This is needed for the allocations inside the add_vector
    recombinant.add_vector(pOri2, COPY_NUMBER)
    
    # recombinant.growth_reaction.bounds = (gr, gr)
    # recombinant.objective = 'DM_glc_e' # It is uptake, so the direction should be maximization
    
    # saving
    name = 'iJO1366_{}_{}_copies_{}_plasmid'.format(
        recombinant.name,
        COPY_NUMBER,
        pOri2.id
        )
    filepath = 'models/{}'.format(name) 
    
    save_json_model(recombinant, filepath)
    
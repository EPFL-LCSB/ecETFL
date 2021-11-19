# -*- coding: utf-8 -*-

from os.path import join as pjoin
import os

import cobra

import pandas as pd

import numpy as np

from pytfa.io import load_thermoDB,                    \
                            read_lexicon, annotate_from_lexicon,            \
                            read_compartment_data, apply_compartment_data


from etfl.core.expression import is_me_compatible
from etfl.core import Enzyme, Ribosome, RNAPolymerase, ThermoMEModel, MEModel
from etfl.core.rna import mRNA

from pytfa.optim.utils import symbol_sum

from collections import defaultdict
from numbers import Number

import re

from ecoli_utils import find_bbb_ratio, infer_enzyme_from_gpr

MEDIAN_KCAT = 65 * 3600 # the mean value of kcats in obrien 2013

def clean_string(s):

    s = s.replace('-','_')

    # Remove invalid characters
    s = re.sub('[^0-9a-zA-Z_]', '', s)

    # Remove leading characters until we find a letter or underscore
    s = re.sub('^[^a-zA-Z_]+', '', s)

    return s

data_dir = '../info_ecoli'

#########################
###     BASE MODEL    ###
#########################

def get_model(solver):
    vanilla_model = cobra.io.load_json_model(pjoin(data_dir,'iJO1366_mod.json'))
    # to remove the artificial gene s0001:
    cobra.manipulation.delete.remove_genes(vanilla_model, ['s0001'], remove_reactions=False)
    vanilla_model.slim_optimize()
    vanilla_model.solver = solver
    vanilla_model.slim_optimize()


    def sanitize_varnames(model):
        for met in model.metabolites:
            if met.id[0].isdigit():
                met.id = '_' + met.id
        for rxn in model.reactions:
            if rxn.id[0].isdigit():
                rxn.id = '_' + rxn.id
        model.repair()

        return model

    # Add cystein -> selenocystein transformation for convenience
    # selcys = cobra.Metabolite(id='selcys__L_c', compartment='c', formula='C3H7NO2Se')
    # selcys_rxn = cobra.Reaction(id='PSEUDO_selenocystein_synthase',
    #                             name='PSEUDO Selenocystein_Synthase')
    # selcys_rxn.add_metabolites(
    #     {vanilla_model.metabolites.cys__L_c: -1, selcys: +1})
    # vanilla_model.add_reactions([selcys_rxn])

    sanitize_varnames(vanilla_model)
    vanilla_model.slim_optimize()
    return vanilla_model

# ------------------------------------------------------------
# Thermo
# ------------------------------------------------------------

def get_thermo_data():
    # Load Thermo data
    thermo_data = load_thermoDB(pjoin(data_dir,'thermo_data','thermo_data.thermodb'))
    lexicon = read_lexicon(pjoin(data_dir,'thermo_data','iJO1366_lexicon.csv'))
    # lexicon = curate_lexicon(read_lexicon('thermo_data/iJO1366_lexicon.csv'))
    compartment_data = read_compartment_data(pjoin(data_dir,'thermo_data',
                                                   'iJO1366_compartment_data.json'))


    def curate_lexicon(lexicon):
        ix = pd.Series(lexicon.index)
        ix = ix.apply(lambda s: str.replace(s,'-','-'))
        ix = ix.apply(lambda s: '_'+s if s[0].isdigit() else s)
        lexicon.index = ix
        return lexicon

    lexicon = curate_lexicon(lexicon)

    return thermo_data, lexicon, compartment_data

#------------------------------------------------------------
# Data
#------------------------7.54--------------------------------

# Essentials
def get_essentials():
    return dict(atp='atp_c',
                adp='adp_c',
                amp='amp_c',
                gtp='gtp_c',
                gdp='gdp_c',
                pi='pi_c',
                ppi='ppi_c',
                h2o='h2o_c',
                h='h_c')

# Growth-related abundances
def get_neidhardt_data():
    # The data for RNA, DNA and protein is form Neidhardt book as ETFL.
    # The data for lipid, carbohydrate, ion and cofactor is derived from two sources:
    # iJO1366 itself for the lowest growth rate and for the other growth rates
    # from Bionumber: https://bionumbers.hms.harvard.edu/bionumber.aspx?s=n&v=16&id=104954
    neidhardt_data = pd.read_excel(pjoin(data_dir,'neidhardt_tab2.xlsx'),
                                   skiprows=range(0,6),
                                   skipfooter=39)
    mu_cols = ['mu=0.41','mu=0.7','mu=1.0','mu=1.38','mu=1.72']
    neidhardt_data.columns = ['parameter','symbol','units'] + mu_cols \
                             + ['observed_parameters','footnote']
    neidhardt_data.set_index('symbol', inplace=True)

    Pc = neidhardt_data.loc['Pc (μg)'][mu_cols] # μg/10^9 cells
    Rc = neidhardt_data.loc['Rc (μg)'][mu_cols] # μg/10^9 cells
    Dc = neidhardt_data.loc['Gc (μg)'][mu_cols] # μg/10^9 cells
    Mc = neidhardt_data.loc['Mc (μg)'][mu_cols] # μg dry weight/10^9 cells
    Cc = neidhardt_data.loc['Cc (%gDW)'][mu_cols] # μg/10^9 cells
    Lc = neidhardt_data.loc['Lc (%gDW)'][mu_cols] # μg/10^9 cells
    Ic = neidhardt_data.loc['Ic (%gDW)'][mu_cols] # μg/10^9 cells
    
    scaling_factor = get_macromole_ratio()['total mass']
    neidhardt_prel = (Pc/Mc).astype(float)*scaling_factor
    neidhardt_rrel = (Rc/Mc).astype(float)*scaling_factor
    neidhardt_drel = (Dc/Mc).astype(float)*scaling_factor
    crel = Cc.astype(float)
    lrel = Lc.astype(float)
    irel = Ic.astype(float)
    neidhardt_mu = pd.Series(Pc.index.str.replace('mu=','')).astype(float)

    return neidhardt_mu, neidhardt_rrel, neidhardt_prel, neidhardt_drel, \
        crel, lrel, irel



#------------------------------------------------------------
# Expression
#------------------------------------------------------------

# Data
# Sequences from KEGG
nt_sequences = pd.read_csv(pjoin(data_dir,'iJO1366_nt_seq_kegg.csv'),
                           index_col = 0,
                           header = None)[1]

def get_nt_sequences():
    return nt_sequences

# Sequences from KEGG
aa_sequences = pd.read_csv(pjoin(data_dir,'iJO1366_aa_seq_kegg.csv'),
                           index_col = 0,
                           header = None)[1]

def get_aa_sequences():
    return aa_sequences

# data for each related protein abundance from PaxDB
prot_abundance =  pd.read_excel(pjoin(data_dir,'abundance_table.xlsx'),
                                   skiprows=range(0,11),
                                   header = 0,
                                   )
# data for each protein MWs (including enzymes) from KEGG
mws = pd.read_csv(pjoin(data_dir,'all_peptides_mw.csv'),header=0)

# iJO kcat info from:
# Davidi, Dan, et al.
# "Global characterization of in vivo enzyme catalytic rates and their correspondence to in vitro kcat measurements."
# Proceedings of the National Academy of Sciences 113.12 (2016): 3401-3406.
kcat_info_milo = pd.read_excel(pjoin(data_dir,'pnas.1514240113.sd01.xlsx'),
                               sheet_name='kcat 1s',
                               header=2,
                               )
kmax_info_milo = pd.read_excel(pjoin(data_dir,'pnas.1514240113.sd01.xlsx'),
                               sheet_name='kmax 1s',
                               header=2,
                               )
kcat_info_aggregated    = pd.read_csv(pjoin(data_dir,'aggregated_kcats.csv'),
                                      index_col = 0)
ec_info_ecocyc          = pd.read_csv(pjoin(data_dir,'complex2ec.csv'),
                                      index_col = 0)
composition_info_ecocyc = pd.read_csv(pjoin(data_dir,'complex2genes.csv'),
                                      index_col = 0)
reaction2complexes_info_obrien = pd.read_excel(
    pjoin(data_dir, 'obrien2013_SI_tab10.xlsx'), index_col=0, usecols=[0, 1])
complexes2peptides_info_obrien = pd.read_excel(
    pjoin(data_dir, 'obrien2013_SI_tab1.xlsx'), index_col=0, usecols=[0, 1])
reaction2kcat = pd.read_excel(
    pjoin(data_dir, 'rxn2kcat.xlsx'), index_col=0,)

reaction2complexes_info_lloyd = pd.read_csv(
    pjoin(data_dir, 'lloyd_2018_enzyme_reaction_association.txt'),
    delimiter = '\t',
    index_col = 0,
    header=None)
reaction2complexes_info_lloyd.columns = ['Enzymes']
complexes2peptides_info_lloyd = pd.read_csv(
    pjoin(data_dir, 'lloyd_2018_protein_complexes.txt'),
    delimiter = '\t',
    index_col = 0,
    usecols=[0,2],
    header = None)
complexes2peptides_info_lloyd.columns = ['Gene composition']

gene_names = pd.read_csv(pjoin(data_dir,'gene2bname.txt'), delimiter='\t',
                         index_col=0)


# mRNA degardation rates from
# Bernstein et al. (2002) Proc. Natl. Acad. Sci. USA, 10.1073/pnas.112318199
# "Global analysis of mRNA decay and abundance in Escherichia coli at single-gene resolution using two-color fluorescent DNA microarrays"
bernstein_ecoli_deg_rates = pd.read_excel(
    pjoin(data_dir,'bernstein_2002_mrna_deg.xls'),
    skiprows=range(8),
    index_col=0)


# Bionumber
# http://bionumbers.hms.harvard.edu/bionumber.aspx?id=100528
# ID	100528
# Property	GC content of E. coli K12 chromosome
# Organism	Bacteria Escherichia coli
# Value	50.8
# Units	%
# Reference	Escherichia coli K12 (Escherichia coli K12 substr. MG1655) Genome Browser Gateway link scroll down to slightly above bottom of page
# Comments	GC content value given 50.79%. Length of chromosome 4,639,675bp. Retrieved March 29th 2017
# Assumed GC ratio of 0.5078
gc_ratio = 0.5078 # %
chromosome_len = 4639675 # bp

def get_ecoli_gen_stats():
    return chromosome_len, gc_ratio

def get_ratios():
    # Bionumbers
    # ID        104876
    # Property  Amino acid composition of the proteins from E. coli cell supernatant
    # http://bionumbers.hms.harvard.edu/bionumber.aspx?id=104876
    # Unit: per 100 moles aas
    aa_ratios = {
        'K':9.01/100,  # lys__L_c
        'H':1.91/100,  # his__L_c
        'R':7.30/100,  # arg__L_c
        'D':8.30/100,  # asp__L_c
        'E':10.08/100, # glu__L_c
        'G':8.18/100,  # gly_c
        'A':10.98/100, # ala__L_c
        'V':9.63/100,  # val__L_c
        'L':7.40/100,  # leu__L_c
        'I':5.51/100,  # ile__L_c
        'T':5.22/100,  # thr__L_c
        'S':4.38/100,  # ser__L_c
        'P':3.67/100,  # pro__L_c
        'Y':1.78/100,  # tyr__L_c
        'F':3.03/100,  # phe__L_c
        'W':0.69/100,  # trp__L_c
        'M':2.40/100,  # met__L_c
        'C':0.53/100,  # cys__L_c
    }


    nt_ratios = {'u':0.5*(1-gc_ratio),
                 'a':0.5*(1-gc_ratio),
                 'g':0.5*gc_ratio,
                 'c':0.5*gc_ratio,
                 }

    return nt_ratios,aa_ratios



def get_monomers_dict():

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
               # 'U': 'selcys__L_c',
               'W': 'trp__L_c',
               'Y': 'tyr__L_c',
               'V': 'val__L_c', }

    rna_nucleotides = {
                'u': 'utp_c',
                'g': 'gtp_c',
                'a': 'atp_c',
                'c': 'ctp_c'}

    rna_nucleotides_mp = {
                'u': 'ump_c',
                'g': 'gmp_c',
                'a': 'amp_c',
                'c': 'cmp_c'}

    dna_nucleotides = {
                't': 'dttp_c',
                'g': 'dgtp_c',
                'a': 'datp_c',
                'c': 'dctp_c'}


    return aa_dict, rna_nucleotides, rna_nucleotides_mp, dna_nucleotides


def remove_from_biomass_equation(model, nt_dict, aa_dict, atp_id, adp_id,
                                 h2o_id, h_id, pi_id):

    # According to discussions, should only remove GAM

    mets_to_rm = dict()

    old_total_stoich = abs(sum([x for x in
                                model.growth_reaction.metabolites.values() if
                                x<0]))

    expression_mets = list(nt_dict.values()) + list(aa_dict.values())

    for m,stoich in model.growth_reaction.metabolites.items():
        if m.id in expression_mets:
            mets_to_rm[m] = -1*stoich

    model.growth_reaction.add_metabolites(mets_to_rm)

    # We need to add back the ATP from the GAM that we just removed but should
    # still be taken into account. Indeed, nADP =/= nATP in the original
    # biomass reaction:
    # -54.119975 atp_c + .... --> 53.95 adp_c

    atp = model.metabolites.get_by_id(atp_id)
    adp = model.metabolites.get_by_id(adp_id)
    pi  = model.metabolites.get_by_id(pi_id)
    h2o = model.metabolites.get_by_id(h2o_id)
    h = model.metabolites.get_by_id(h_id)
    atp_recovery = model.growth_reaction.metabolites[adp]
    model.growth_reaction.add_metabolites({atp:-1*atp_recovery})



# Prot degradation
# BNID 111930
# Moran MA et al., Sizing up metatranscriptomics. ISME J. 2013 Feb7(2):237-43.
# A typical bacterial protein half-Life is ~20 h
# -------
# tau = 1/kdeg = t_0.5 /ln(2)
# kdeg = ln(2)/t_0.5
kdeg_enz = np.log(2)/20 # [h-1]

# The same as ribosome : Half life of a ribosome is 5 days
kdeg_rrna = np.log(2)/(5*24)

# From :
# http://book.bionumbers.org/how-fast-do-rnas-and-proteins-degrade/
# Figure 1: Measured half lives of mRNAs in E. coli, budding yeast and mouse NIH3T3 fibroblasts.
# (A, adapted from J. A. Bernstein et al., Proc. Natl Acad. Sci. USA 99:9697, 2002;
#  B, adapted from Y. Wang et al., Proc. Natl Acad. Sci. USA 99:5860, 2002;
#  C. adapted from B. Schwanhausser, Nature, 473:337, 2013).
# -------
# Mean half life of mrna is 5 minutes in ecoli
# tau = 1/kdeg = t_0.5 /ln(2)
# kdeg = ln(2)/t_0.5
kdeg_mrna = 60*np.log(2)/5

# Average mrna length from Bionumber 100023
# http://bionumbers.hms.harvard.edu/bionumber.aspx?&id=100023&ver=3
# mrna_length_avg = 370 # nm
mrna_length_avg = 1000

# Average trna length from Bionumber 109645
# https://bionumbers.hms.harvard.edu/bionumber.aspx?id=109645&ver=1&trm=average+length+of+trna&org=
# between 75-90 nucleotides
trna_length_avg = 80

# Average peptide length
peptide_length_avg = int(np.round(mrna_length_avg/3))

kdeg_trna = 0

def get_mrna_metrics():
    return kdeg_mrna, mrna_length_avg

def get_trna_metrics():
    return kdeg_trna, trna_length_avg

def get_enz_metrics():
    return kdeg_enz, peptide_length_avg

# Generate a coupling dict
def is_gpr(s):
    return bool(s) and s != '[]'


# Milo kcats
#############

def get_homomer_coupling_dict(model, mode = 'kcat'):
    if mode == 'kcat':
        k_info = kcat_info_milo
        k_column = 'kcat per active site [1/s]'
        n_column = 'catalytic sites per complex'
    elif mode == 'kmax':
        k_info = kmax_info_milo
        k_column = 'kmax per polypeptide chain [s-1]'
        n_column = 'polypeptides per complex'
    elif isinstance(mode, Number):
        k_info = kmax_info_milo
        k_column = 'kmax per polypeptide chain [s-1]'
        n_column = 'polypeptides per complex'

    else:
        raise Exception("Mode {} not understood. Should be 'kcat' or 'kmax' "
                        "or a number")

    coupling_dict = dict()

    for x in model.reactions:
        composition, kcat_bwd, kcat_fwd = get_rate_constant(x, k_info,
                                                            k_column,
                                                            n_column)
        if isinstance(mode, Number):
            # It is a number
            kcat_fwd = kcat_bwd = mode

        if kcat_fwd == 0 and kcat_bwd == 0:
            continue

        if kcat_bwd == 0:
            kcat_bwd = kcat_fwd

        if kcat_fwd == 0:
            kcat_fwd = kcat_bwd

        if not composition:
            # WE cannot make the enzyme. To infer a composition, use
            # infer_missing_enzymes in get_coupling_dict
            continue
        #FIXME several polypeptides per complex ??

        new_enzyme = Enzyme(x.id,
                            kcat_fwd=kcat_fwd,
                            kcat_bwd=kcat_bwd,
                            kdeg=kdeg_enz,
                            composition = composition)

        coupling_dict[x.id] = [new_enzyme]

    return coupling_dict


def get_rate_constant(reaction, k_info, k_column, n_column):
    k_data = k_info[k_info['reaction (model name)'] == reaction.id]
    k_data_reverse = k_info[k_info['reaction (model name)'] == reaction.id + '_reverse']
    data = kcat_info_milo[kcat_info_milo['reaction (model name)'] == reaction.id]
    data_reverse = kcat_info_milo[kcat_info_milo['reaction (model name)'] == reaction.id + '_reverse']
    kcat_fwd = 0
    kcat_bwd = 0
    composition = {}
    if len(k_data) > 0:
        candidate_k = k_data[k_column].iloc[0]
        n_peptides = data['polypeptides per complex'].iloc[0] \
            if len(data) > 0 else 1
        n_sites = data[n_column].iloc[0] \
            if len(data) > 0 else 1

        kcat_fwd = candidate_k \
                   * n_sites \
                   * 3600  # s/h
        composition = {k_data['bnumber'].iloc[0]: n_peptides}
    if len(k_data_reverse) > 0:
        candidate_k = k_data_reverse[k_column].iloc[0]
        n_peptides = data_reverse['polypeptides per complex'].iloc[0] \
            if len(data_reverse) > 0 else 1
        n_sites = data_reverse[n_column].iloc[0] \
            if len(data_reverse) > 0 else 1

        kcat_bwd = candidate_k \
                   * n_sites \
                   * 3600  # s/h
        composition = {
            k_data_reverse['bnumber'].iloc[0]: n_peptides} \
            if not composition else composition
    return composition, kcat_bwd, kcat_fwd


# Aggregated kcats
##################

def ec2ecocyc(ec_number):
    if not isinstance(ec_number, list):
        return ec_info_ecocyc[ec_info_ecocyc['ec'] == ec_number]
    else:
        return ec_info_ecocyc[ec_info_ecocyc['ec'].isin(ec_number)]

def score_against_genes(putative_genes, reaction_genes):
    score = 0

    putative_genes_list = putative_genes.split('" // "')
    reaction_gene_list  = [x.name for x in reaction_genes]
    s1 = len(set(putative_genes_list).intersection(reaction_gene_list))
    s2 = len(set(putative_genes_list).difference  (reaction_gene_list))
    s3 = len(set(reaction_gene_list) .difference  (putative_genes_list))

    score = 2*s1-s2-s3
    # print(putative_genes_list, reaction_gene_list, score)
    return score

def match_ec_genes_ecocyc(ecocyc, genes, threshold=0.5):
    this_data = composition_info_ecocyc[composition_info_ecocyc['complex'].isin(ecocyc)]
    scores = this_data['putative_genes'].apply(score_against_genes, args=[genes])
    selectable = this_data[scores>len(genes)*threshold]
    if len(selectable) == 0:
        return None, scores
    else:
        return selectable[scores == scores.max()].iloc[0], scores

def ecocyc2composition(ecocyc):
    ecocyc_comp = composition_info_ecocyc[composition_info_ecocyc['complex'] == ecocyc]
    # We do a left join to get the bnumbers that are used in the model
    composition_data = ecocyc_comp.merge(gene_names, right_index=True,
                                  how='left', left_on='gene')
    composition = defaultdict(int)
    for e,row in ecocyc_comp.iterrows():
        this_gene = row['gene']
        try:
            this_b_number = gene_names.loc[this_gene]['b#']
            composition[this_b_number] += composition_data['coeffs'].iloc[0]
        except:
            # ecoli.logger.warning('Could not find gene associated to {}'
            #                      .format(row['obj_ids']))
            # ecoli.logger.info(ecocyc_comp)
            pass

    return composition

comp_regex = re.compile(r'(b[0-9]{4})\((\d?)\)')

def complex2composition(complex_name):
    # Silence modifications
    if '_mod_' in complex_name:
        complex_name = complex_name[0:complex_name.index('_mod_')]

    try:
        composition_string = complexes2peptides_info_lloyd.loc[complex_name,'Gene composition']
    except KeyError:
        composition_string = complexes2peptides_info_obrien.loc[complex_name,'Gene composition']

    composition_dict = {}
    groups = comp_regex.findall(composition_string)
    for peptide, stoich in groups:
        if stoich == '':
            stoich = 1
        else:
            stoich = int(stoich)
        composition_dict[peptide] = stoich

    return composition_dict

def ec2kcat(ec_number):
    try:
        return kcat_info_aggregated['kcat'].loc[ec_number].max() * 3600  # s/h
    except KeyError:
        return None

def check_id_in_reaction_list(the_id, df):
    if the_id in df.index:
        return the_id
    elif the_id[0] == '_' and the_id[1:] in df.index:
        return the_id[1:]
    elif the_id + '1' in df.index:
        return the_id + '1'
    else:
        return ''

def rxn2kcat(rxn_id):
    kcat_fwd = reaction2kcat['kcat_fwd']\
        [reaction2kcat.index.get_loc(rxn_id)]
    kcat_bwd = reaction2kcat['kcat_bwd']\
        [reaction2kcat.index.get_loc(rxn_id)]
    # if kcat not found for the forward, set both to median value
    if np.isclose(kcat_fwd,0,atol=1e-6):
        kcat_fwd = MEDIAN_KCAT
        kcat_bwd = MEDIAN_KCAT
    # if kcat not found for the forward, set it equal to forward
    elif np.isclose(kcat_bwd,0,atol=1e-6):
        kcat_bwd = kcat_fwd
    return kcat_fwd, kcat_bwd

def get_generalized_coupling_dict(model):
    
    coupling_dict = dict()
    total_enz_dict = dict() # list of enzyme IDs
    for rxn in model.reactions:
        # the name of rxns has been sanisitized, reverse them for look up
        query = rxn.id[1:] if rxn.id[0] == '_' else rxn.id
        try:
            enz_set = reaction2complexes_info_obrien['Enzymes'] \
                [reaction2complexes_info_obrien.index.get_loc(query)]
        except KeyError:
            continue
        # there are some redundant rxn IDs, take the first one
        if not isinstance(enz_set,str):
            enz_set = enz_set.values[0]
        # at this stage we might have isozymes separated with OR
        isozymes = enz_set.split(' OR ')
        # the name of some complexes is contaminated with _mod
        isozymes = [x.split('_mod')[0] for x in isozymes]
        these_enz = []
        for enz in isozymes:
            if enz in total_enz_dict.keys(): # this enzyme is already defined?
                these_enz += [total_enz_dict[enz]]
            else:
                pep_set = complexes2peptides_info_obrien['Gene composition'] \
                    [complexes2peptides_info_obrien.index.get_loc(enz)]
                # pep_set includes the name of constituting peptides
                # seperated by AND and the stoichiometric coefficient in ()
                pep_list = pep_set.split(' AND ')
                pep_list = [x.replace(')','') for x in pep_list]
                composition = {x.split('(')[0] : int(x.split('(')[1]) \
                               for x in pep_list}
                kcat_fwd, kcat_bwd = rxn2kcat(query)
                # the names will be sanitized:
                id_ = '_'  + enz if enz[0].isdigit() else enz
                new_enz = Enzyme(id=id_,
                                 kcat_fwd=kcat_fwd,
                                 kcat_bwd=kcat_bwd, 
                                 kdeg=kdeg_enz,
                                 composition=composition)
                these_enz += [new_enz]
                total_enz_dict[enz] = new_enz

        coupling_dict[rxn.id] = these_enz
            
    return coupling_dict


def get_aggregated_coupling_dict(model, coupling_dict = dict()):
    aggregated_coupling_dict = defaultdict(list)

    for x in model.reactions:
        # reactions starting with a number have been sanitized to start with '_'
        if x.id in coupling_dict:
            # We already have info
            continue

        if 'ec_numbers' not in x.notes or x.notes['ec_numbers'] == ['nan']:
            # There is nothing we can do
            continue

        reaction_ecs = x.notes['ec_numbers']

        lloyd_id = check_id_in_reaction_list(x.id, reaction2complexes_info_lloyd)
        obrien_id = check_id_in_reaction_list(x.id, reaction2complexes_info_obrien)

        if lloyd_id:
            complex_names = reaction2complexes_info_lloyd.loc[lloyd_id,'Enzymes']\
                .split(' OR ')
        elif obrien_id:
            complex_names = reaction2complexes_info_obrien.loc[obrien_id,'Enzymes']\
                .split(' OR ')
        else:
            continue

        for e,this_complex_name in enumerate(complex_names):

            # Start with this:
            composition = complex2composition(this_complex_name)
            if not composition:
                # Skip this one
                continue

            this_ec = x.notes['ec_numbers'][0]
            kcat = ec2kcat(this_ec)

            if kcat is None:
                continue

            cleaned_cplx_name = clean_string(this_complex_name)

            new_enzyme = Enzyme('{}_{}_{}'.format(x.id,cleaned_cplx_name,e),
                                name='{}_{}: {}'.format(x.id, e, this_complex_name),
                                kcat=kcat,
                                kdeg=kdeg_enz,
                                composition = composition)

            new_enzyme.notes['EC'] = this_ec

            aggregated_coupling_dict[x.id].append(new_enzyme)

    return aggregated_coupling_dict


def get_lloyd_keffs():
    import json
    with open(pjoin(data_dir, 'lloyd_2018_keffs.json'), 'r') as fid:
        keffs = json.load(fid)

    new_keffs = dict()

    for key in keffs.keys():
        new_key = key.replace('keff_','')
        new_key = new_key.replace('_DASH_','-')
        new_keffs[new_key] = keffs[key] * 3600  # s/h
    return new_keffs


def get_coupling_dict(model, mode, atps_name = None, infer_missing_enz=False):
    coupling_dict = dict()
    
    homomer_coupling_dict = get_homomer_coupling_dict(model, mode=mode)
    aggregated_coupling_dict = get_aggregated_coupling_dict(model, homomer_coupling_dict)
    generalized_coupling_dict = get_generalized_coupling_dict(model)
    # Most important last
    # coupling_dict.update(lloyd_dict)
    coupling_dict.update(generalized_coupling_dict)
    coupling_dict.update(aggregated_coupling_dict)
    coupling_dict.update(homomer_coupling_dict)

    if infer_missing_enz:
        inferred_enz = dict()

        if isinstance(mode, Number):
            kcat = mode
        else:
            kcat = get_average_kcat()

        for r in model.reactions:
            if r.id in coupling_dict:
                continue

            if r.id not in coupling_dict \
                    and is_me_compatible(r):
                inferred_enz[r.id] = infer_enzyme_from_gpr(r,
                                                           default_kcat=kcat,
                                                           default_kdeg=kdeg_enz)

        coupling_dict.update(inferred_enz)

    # ATP Synthase bypasses numeric kcat
    if atps_name is not None:
        atps = get_atp_synthase_coupling(atps_name)
        coupling_dict.update(atps)

    # coupling_dict.update(get_transporters_coupling())

    return coupling_dict

def get_average_kcat():
    # return np.median(list(get_lloyd_keffs().values()))
    return np.mean(list(get_lloyd_keffs().values()))

def get_atp_synthase_coupling(atps_name):
    """
    ATP synthesis rate of F1F0 ATP synthase
    Range 	at room temperature ∼0.060-0.10 μmol/min/mg of membrane protein : at 37°C 0.20 μmol/min/mg of membrane protein
    Organism 	Bacteria Escherichia coli
    Reference 	Tomashek JJ, Glagoleva OB, Brusilow WS. The Escherichia coli F1F0 ATP synthase displays biphasic synthesis kinetics. J Biol Chem. 2004 Feb 6 279(6):4465-70 DOI: 10.1074/jbc.M310826200 p.4467 right column bottom paragraphPubMed ID14602713
    Primary Source 	[18] Etzold C, Deckers-Hebestreit G, Altendorf K. Turnover number of Escherichia coli F0F1 ATP synthase for ATP synthesis in membrane vesicles. Eur J Biochem. 1997 Jan 15 243(1-2):336-43.PubMed ID9030757
    Method 	Luciferase assay
    Comments 	P.4467 right column bottom paragraph: "Previously, Etzold et al. (primary source) used the luciferase assay to measure the turnover number of the ATP synthase during synthesis by membrane vesicles of E. coli. They measured ATP synthesis rates of ∼0.060-0.10 μmol/min/mg of membrane protein at room temperature and 0.20 μmol/min/mg of membrane protein at 37 °C."
    Entered by 	Uri M
    ID 	115175
    :return:
    """

    #      umol/(mg.min)  *  min/h  * mmol/umol * mg/g
    # kcat = 0.08           *  60     * 1
    kcat = 232           *  3600     * 1
    composition = {
        'b3733':1,
        'b3738':1,
        'b3732':3,
        'b3736':2,
        'b3731':1,
        'b3735':1,
        'b3737':10,
        'b3734':3,
        }

    atp_synthase = Enzyme(atps_name,
                        name='ATP Synthase',
                        kcat_fwd=kcat,
                        kcat_bwd=kcat,
                        kdeg=kdeg_enz,
                        composition=composition)

    return {atps_name:[atp_synthase]}



def get_mrna_dict(model):
    mrna_dict = dict()

    # Generate a mRNA dict

    for x in nt_sequences.index:
        try:
            the_gene = model.genes.get_by_id(x)
        except KeyError:
            model.genes += [cobra.Gene(id=x)]
            the_gene = model.genes.get_by_id(x)

        # Try to get half life from Bernstein et al.
        try:
            t_half = bernstein_ecoli_deg_rates.loc[x.upper()]['medium, min.1'] #M9 medium
            # Mean half life of mrna is 5 minutes in ecoli
            # tau = t_0.5 /ln(2)
            # kdeg = 1/tau [min^-1] * [min/h]
            this_kdeg_mrna = (60 * np.log(2) / t_half)
        except KeyError:
            if x in rrna_genes:
                this_kdeg_mrna = kdeg_rrna # Same as ribosome
            else:
                this_kdeg_mrna = kdeg_mrna # Average value of 5 mins

        if np.isnan(this_kdeg_mrna):
            this_kdeg_mrna = kdeg_mrna # Average value of 5 mins

        new_mrna = mRNA(x,
                        kdeg = this_kdeg_mrna,
                        gene_id = the_gene.id)
        mrna_dict[x] = new_mrna
    return mrna_dict


# Half life of a ribosome is 5 days
kdeg_rib = np.log(2)/(5*24)


#[ ribosomes and RNAP ]#
rpeptide_genes = pd.read_csv(pjoin(data_dir,'ribosomal_proteins_ecoli.tsv'),
                             delimiter='\t',
                             header=None)[0]
rpeptide_genes = rpeptide_genes.str.split(':').apply(lambda x:x[1])

rrna_genes = ['b3851', 'b3854', 'b3855']
def get_rib():
    """
    # Ribosome

    rRNA:
    b3851: K01977 16S ribosomal RNA | (RefSeq) rrsA; 16S ribosomal RNA of rrnA operon
    b3854: K01980 23S ribosomal RNA | (RefSeq) rrlA; 23S ribosomal RNA of rrnA operon
    b3855: K01985 5S ribosomal RNA | (RefSeq) rrfA; 5S ribosomal RNA of rrnA operon
    # rPeptides:
    See file ribosomal_proteins_ecoli.tsv

    :return:
    """
    # bionumber : 100059
    # between 12-20 aa/sec 
    rib = Ribosome(id='rib', name='Ribosome', kribo=12 * 3600, kdeg=kdeg_rib,
                   composition = rpeptide_genes, rrna=rrna_genes)
    return rib

# http://bionumbers.hms.harvard.edu/bionumber.aspx?&id=100060&ver=32
# Bionumber ID  100060
# Value 	    45 - 85 nt/sec
# Source        Bremer, H., Dennis, P. P. (1996) Modulation of chemical composition and other parameters of the cell by growth rate.
#               Neidhardt, et al. eds. Escherichia coli and Salmonella typhimurium: Cellular and Molecular Biology, 2nd ed. chapter 97 Table 3


rnap_genes = ['b3295','b3649','b3987','b3988']

def get_rnap():
    """
    # RNAP

    b3295: K03040 DNA-directed RNA polymerase subunit alpha [EC:2.7.7.6] | (RefSeq) rpoA; RNA polymerase, alpha subunit
    b3649: K03060 DNA-directed RNA polymerase subunit omega [EC:2.7.7.6] | (RefSeq) rpoZ; RNA polymerase, omega subunit
    b3987: K03043 DNA-directed RNA polymerase subunit beta [EC:2.7.7.6] | (RefSeq) rpoB; RNA polymerase, beta subunit
    b3988: K03046 DNA-directed RNA polymerase subunit beta' [EC:2.7.7.6] | (RefSeq) rpoC; RNA polymerase, beta prime subunit

    :return:
    """
    ktrans = 45
    rnap = RNAPolymerase(id='rnap',
                         name='RNA Polymerase',
                         ktrans = ktrans*3600,
                         # kdeg = 0.2,
                         kdeg = kdeg_enz,
                         composition = rnap_genes)


    return rnap

def get_rnap_rrna():
    """
    # RNAP

    b3295: K03040 DNA-directed RNA polymerase subunit alpha [EC:2.7.7.6] | (RefSeq) rpoA; RNA polymerase, alpha subunit
    b3649: K03060 DNA-directed RNA polymerase subunit omega [EC:2.7.7.6] | (RefSeq) rpoZ; RNA polymerase, omega subunit
    b3987: K03043 DNA-directed RNA polymerase subunit beta [EC:2.7.7.6] | (RefSeq) rpoB; RNA polymerase, beta subunit
    b3988: K03046 DNA-directed RNA polymerase subunit beta' [EC:2.7.7.6] | (RefSeq) rpoC; RNA polymerase, beta prime subunit

    :return:
    """
    ktrans = 85
    rnap = RNAPolymerase(id='rnap_rrna',
                         name='RNA Polymerase (rRNA)',
                         ktrans = ktrans*3600,
                         # kdeg = 0.2,
                         kdeg = kdeg_enz,
                         composition = rnap_genes)


    return rnap

gene_list =  pd.read_excel(pjoin(data_dir,'trans_list.xlsx'),
                                      sheet_name = 'Genes')
def get_transcription_dict():
    rnap = get_rnap()
    rnap_rrna = get_rnap_rrna()
    transcription_dict = dict()
    for the_gene_id in gene_list['Gene_id'].values:
        if the_gene_id in rrna_genes:
            transcription_dict[the_gene_id] = [rnap_rrna.id]
        else:
            transcription_dict[the_gene_id] = [rnap.id]
    
    return transcription_dict

def remove_DNA(model):
    # if allocation data is available we can remove DNA from biomass, at it 
    # is made growth dependent in the code.
    DNA = model.metabolites.get_by_id('dna_c')
    DNA_rxn = model.reactions.get_by_id('DNArxn')
    model.remove_metabolites(DNA)
    model.remove_reactions(DNA_rxn)
    return

# for some procedures we need to have unmodified model
cobra_model = cobra.io.load_json_model('../input_models/iJO1366.json')

def get_macromole_ratio():
    # model must be an unmodified model
    mass_ratios = find_bbb_ratio(cobra_model)
    return mass_ratios

def modify_GAM(model, growth_id, prot_rxn_id=None):
    '''
    This function modify GAM by removing the energetic cost of peptide synthesis.

    Parameters
    ----------
    model: cobra-model, with modified biomass
    growth_id: str
    prot_rxn_id : str, optional
        The id of protein pseudoreaction if aa are not taken into account in the
        biomass reaction. The default is None.

    Returns
    -------
    cobra-model with modified GAM

    '''
    # model must be an unmodified model
    growth_rxn = cobra_model.reactions.get_by_id(growth_id)
    
    if prot_rxn_id is None:
         aa_dict,_,_,_ = get_monomers_dict()
         aa_stoic = [abs(growth_rxn.get_coefficient(x)) \
                     for _,x in aa_dict.items()]  
    else:
         prot_rxn = cobra_model.reactions.get_by_id(prot_rxn_id)
         aa = [x.id for x in prot_rxn.reactants]
         aa_stoic = [abs(prot_rxn.get_coefficient(x)) \
                     for x in aa]
     
    tot_aa = sum(aa_stoic) # number of moles required to synthesize one unit of biomass
    ########### In the case of ecoli trna charging cost should be also considered
    gtp_expense = 3 * tot_aa # 2 moles of GTP needed per each mole of aa
    # now we should find GAM metabolites
    GAM_mets = {'atp':-1, 'h2o':-1,
                'adp':1,'h':1,'pi':1} # differentiating reactants & products
    essentials = get_essentials()
    # create a dict to subtract from metabolite ids
    mets = {essentials[k]:v*gtp_expense for k,v in GAM_mets.items()}
    # to apply changes on the modified model
    biomass = model.reactions.get_by_id(growth_id)
    biomass.subtract_metabolites(mets)
    
def get_enzyme_fraction(peptides):
    # a function to find the fraction of proteome that is enzyme
    # peptides is the list of model peptides 
    
    # it's possible it cannot find all the enzymes, so a flexibility is introduced
    alpha = 0 # 0% tolerance    
    weight_dict_enz = dict()
    for the_peptide in peptides:
        weight_dict_enz[the_peptide.id] = [the_peptide.molecular_weight * x for x in \
            prot_abundance[prot_abundance['genes']==the_peptide.id]\
                ['abundance']]
    tot_enz = sum([x[0] for x in weight_dict_enz.values() if len(x)>0 and \
                   not np.isnan(x[0])])
    tot_enz *=1000 #it is in KDa
    
    weight_dict_prot = dict()
    for _, pair in prot_abundance.iterrows():
        gene_id = pair[0]
        abundance = pair[1]
        weight_dict_prot[gene_id] = [y * abundance for y in \
                 mws[mws['gene name']==gene_id]['MW']]
    tot_prot = sum([x[0] for x in weight_dict_prot.values() if len(x)>0 and \
                   not np.isnan(x[0])])
    fraction = tot_enz/tot_prot
    return fraction + alpha*fraction
        
# Setting copy numbers for the genes:
def update_copy_numbers(model):
    # http://book.bionumbers.org/how-many-ribosomal-rna-gene-copies-are-in-the-genome/
    # For rRNAs we need to set the copy number equal to 7. 
    copy_dict = {'b3851':7, 'b3854':7, 'b3855':7 }
    for gene,copy_num in copy_dict.items():
        try:
            this_gene = model.genes.get_by_id(gene)
        except KeyError:
            continue
        this_gene.copy_number = copy_num
    return
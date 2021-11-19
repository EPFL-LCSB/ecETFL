#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 28 11:23:45 2020

@author: omid
"""


from Bio.SeqUtils import molecular_weight
import cobra
import pandas as pd
import numpy as np
import requests

from os.path import join as pjoin

data_dir = 'info_ecoli'

kegg_url = "http://rest.kegg.jp/get/{org}:{gene_name}/aaseq"

all_genes = pd.read_excel(pjoin(data_dir,'abundance_table.xlsx'),
                                   skiprows=range(0,11),
                                   header = 0,
                                   )
all_b_genes = pd.Series(['eco:{}'.format(x) for x in all_genes['genes']])




def get_from_kegg(gene_id):
    org,gene_name = gene_id.split(':')
    response = requests.post(kegg_url.format(org=org,
                                             gene_name=gene_name))
    if response.ok:
        eol_ix = response.text.find('\n')
        text = response.text[eol_ix+1:].replace('\n','')
        return text
    else:
        return np.nan

aa_sequences = all_b_genes.apply(get_from_kegg)
aa_sequences.index = all_b_genes.str.split(':').apply(lambda x:x[1])
aa_sequences.dropna(inplace = True)
mws = aa_sequences.apply(molecular_weight,seq_type='protein')/1000
mws.to_csv(pjoin(data_dir,'all_peptides_mw.csv'))
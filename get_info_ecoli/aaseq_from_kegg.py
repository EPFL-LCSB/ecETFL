#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 27 09:46:10 2020

@author: omid
"""


import cobra
import pandas as pd
import numpy as np
import requests

from os.path import join as pjoin

data_dir = 'info_ecoli'

kegg_url = "http://rest.kegg.jp/get/{org}:{gene_name}/aaseq"

ecoli = cobra.io.load_json_model('iJO1366.json')

all_b_genes = pd.Series(['eco:{}'.format(x.id) for x in ecoli.genes])

rnap_genes = pd.Series(['eco:b3295','eco:b3649','eco:b3987','eco:b3988'])
rrna_genes = pd.Series(['eco:b3851','eco:b3854','eco:b3855'])
rpeptide_genes = pd.read_csv(pjoin(data_dir,'ribosomal_proteins_ecoli.tsv'),
                             delimiter='\t',
                             header=None)[0]

# transporter_genes = pd.read_csv(pjoin(data_dir,'transporters_kcats.csv'),
#                                header=0, skiprows=[1,], # Units row
#                                index_col=0)['gene']

all_b_genes = pd.concat([all_b_genes, rnap_genes, rrna_genes, rpeptide_genes])
all_b_genes.drop_duplicates(inplace=True)


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
aa_sequences.to_csv(pjoin(data_dir,'iJO1366_aa_seq_kegg.csv'))
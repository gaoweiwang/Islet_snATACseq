#!/usr/bin/env python3

import os
import sys
import numpy as np
import pandas as pd
import subprocess
import scipy.stats


def run_imbalance(donor, cluster):
	
	snpinfo = pd.read_table('biobank.impute_TOPMed.lift_GRCh37.info.gz', sep='\t', header=0, dtype=str)
	snpinfo = snpinfo.loc[(snpinfo['REF(0)'].str.len()==1) & (snpinfo['ALT(1)'].str.len()==1)]
	snpinfo['CHR'] = snpinfo['SNP'].str.split(':', expand=True)[0].str.replace('chr','')
	snpinfo['CHR'] = 'chr' + snpinfo['CHR']
	snpinfo['POS'] = snpinfo['SNP'].str.split(':', expand=True)[1]
	snpinfo.index = snpinfo['CHR'] + ':' + snpinfo['POS'] + ':' + snpinfo['REF(0)']
	snpinfo = snpinfo.loc[~snpinfo.index.duplicated()]

	pileup = pd.read_table(os.path.join(donor, 'pileup', f'{donor}.{cluster}.pileup'), sep='\t', header=None)
	pileup.columns = ['CHR','POS','REF','TOTAL_COUNT','PILEUP','QUAL']
	pileup['REF'] = pileup['REF'].str.upper()
	pileup['ALT'] = (pileup['CHR'].astype(str) + ':' + pileup['POS'].astype(str) + ':' + pileup['REF']).map(snpinfo['ALT(1)'])
	pileup['IMP_R2'] = (pileup['CHR'].astype(str) + ':' + pileup['POS'].astype(str) + ':' + pileup['REF']).map(snpinfo['Rsq'])
	pileup['REF_COUNT'] = pileup['PILEUP'].str.count(r'\.') + pileup['PILEUP'].str.count(',')
	pileup = pileup.loc[~pd.isnull(pileup['ALT'])]
	pileup['ALT_COUNT'] = [row['PILEUP'].upper().count(row['ALT']) for i,row in pileup.iterrows()]
	pileup['TOTAL_COUNT'] = pileup['REF_COUNT'] + pileup['ALT_COUNT']
	pileup['SNPID'] = pileup['CHR'].astype(str) + ':' + pileup['POS'].astype(str) + ':' + pileup['REF'] + ':' + pileup['ALT']
	pileup = pileup[['SNPID','CHR','POS','REF','ALT','TOTAL_COUNT','REF_COUNT','ALT_COUNT','PILEUP','QUAL','IMP_R2']]
	pileup = pileup.loc[pileup['TOTAL_COUNT']>0]
	pileup['Z'] = (pileup['ALT_COUNT'] - 0.5*pileup['TOTAL_COUNT']) / np.sqrt(0.5*0.5*pileup['TOTAL_COUNT'])
	pileup['P'] = scipy.stats.norm.sf(abs(pileup['Z']))*2
	pileup.to_csv(os.path.join('/nfs/lab/joshchiou/biobank_imbalance/WASP_remap/', donor, 'imbalance', f'{donor}.{cluster}.imbalance.txt'), sep='\t', index=False)
	return

donor = sys.argv[1]

for cluster in ['b0','b1']:
	if not os.path.exists(os.path.join(donor, 'imbalance')):
		os.mkdir(os.path.join(donor, 'imbalance'))
	run_imbalance(donor, cluster)

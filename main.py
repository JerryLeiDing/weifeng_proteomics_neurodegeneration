import numpy as np
import pandas as pd
import re, sys, time
sys.dont_write_bytecode = True  # Avoid caching problems


from scripts.anova import *
from scripts.stat_tests import *
from scripts.plots import *


# Example workflow for KO95 vs GFP comparison, using protein report 
# data = pd.read_csv('protein.csv')
# ko95_anova, _ = anova_general(data, ko93=False, ko95=True, interaction=False)
# ko95_norm = normalize_plexb(ko95_anova)
# ko95_anova.iloc[:,:9] = ko95_norm.iloc[:,:9]
# res = run_protein(ko95_anova, 'KO95', 'GFP')

class RunMarch:
    """ Class to encapsulate the March dataset and associated procedures """
    data_cols = ['CKF_A1','CKF_A2','CKF_A3','CKF_A4',
                 'P25F_A1','P25F_A2','P25F_A3','P25F_A4','P25F_A5','P25F_A6']
    groups = [                                                      
            ('CKF', data_cols[:4]),                 
            ('P25F', data_cols[4:])
    ]
    # List of (str, str) pairs which define all valid pairwise comparisons
    valid_comparisons = [('P25F', 'CKF')]
    protein_results = dict.fromkeys(valid_comparisons, None)
    phospho_results = dict.fromkeys(valid_comparisons, None)

    def __init__(self):
        self.protein = pd.read_csv('data/March_2017/protein.csv')
        self.phospho_nonorm = pd.read_csv('data/March_2017/phospho.csv')

        # Normalize phosphosite by accession_number
        # Merge data portions of both datasets
        m_cols = self.data_cols + ['accession_number']
        phos_prot_merge = pd.merge(
                self.phospho_nonorm[m_cols],
                self.protein[m_cols],
                how='left', on='accession_number')
        # Normalize phosphosites by protein
        # Need to fill non-matched phosphosites with 0
        norm_values = (phos_prot_merge.iloc[:,:10].values -
                phos_prot_merge.iloc[:,11:].fillna(0).values)
        # Recenter norm_values around 0
        norm_values -= np.median(norm_values, axis=0, keepdims=True)
        self.phospho_norm = self.phospho_nonorm.copy()
        self.phospho_norm.iloc[:,:10] = norm_values

        # TODO: save quality figures to plot
        # Using clean_data + ri2py


    def plot_quality(self):
        """ Plots correlation for protein, phospho, and normalized phospho 
            Saved in 'plots/raw_quality/March_TYPE_channel_reps_GROUP.png"""
        compare_channel_replicates(
                self.protein, group=False,
                title='Mar_protein', col_groups=self.groups)
        compare_channel_replicates(
                self.phospho_nonorm, group=False,
                title='Mar_phospho_nonorm', col_groups=self.groups)
        compare_channel_replicates(
                self.phospho_norm, group=False,
                title='Mar_phospho_norm', col_groups=self.groups)


    def do_stat_tests(self, comparison=valid_comparisons[0]):
        """ Runs statistical tests 
            Results are stored in self.protein_results and self.phospho_results """
        if comparison not in self.valid_comparisons:
            raise ValueError("Invalid comparison: see self.valid_comparison")

        prot_res = run_protein(self.protein, comparison[0], comparison[1])
        phos_res = run_protein(self.phospho_norm, comparison[0], comparison[1])

        self.protein_results[comparison] = prot_res
        self.phospho_results[comparison] = phos_res
        
        return (prot_res, phos_res) 

    
    def plot_results(self):
        for (k,v) in self.protein_results.iteritems():
            if v is not None:
                generate_figures(v, LABEL_COLS_PROT,
                        filename='Mar_total_prot_%s_x_%s' % k)
            else:
                print "%s vs %s skipped in March figure gen" % comparison
        for (k,v) in self.phospho_results.iteritems():
            if v is not None:
                generate_figures(v, LABEL_COLS_PROT,
                        filename='Mar_phospho_%s_x_%s' % k)
            else:
                print "%s vs %s skipped in March figure gen" % comparison

        

### ALL DEPRECATED ###
### USE NORMALIZATION PROCEDURE ANOVA ###

class OldFunctions:
    def rerun_proteins(data_file='protein.csv'):
        data = pd.read_csv(data_file)
        res_93_gfp = run_protein(data, 'KO93', 'GFP')
        res_95_gfp = run_protein(data, 'KO95', 'GFP')
        res_dko_gfp = run_protein(data, 'DKO', 'GFP')

        res_93_gfp.to_csv('results/ns_KO93_vs_GFP.csv', index=False)
        res_95_gfp.to_csv('results/ns_KO95_vs_GFP.csv', index=False)
        res_dko_gfp.to_csv('results/ns_DKO_vs_GFP.csv', index=False)

        generate_figures(res_93_gfp, L2, 'ns_KO93_vs_GFP')
        generate_figures(res_95_gfp, L2, 'ns_KO95_vs_GFP')
        generate_figures(res_dko_gfp, L2, 'ns_DKO_vs_GFP')


    def rerun_proteins_onplex(data_file='protein.csv'):
        data = pd.read_csv(data_file)
        res_95_gfp_A = run_protein(data, 'KO95', 'GFP', plex='A')
        res_95_gfp_B = run_protein(data, 'KO95', 'GFP', plex='B')

        combined = pd.merge(res_95_gfp_A, res_95_gfp_B, on='accession_number')
        return combined

    def rerun_psm(data_file='psm.csv'):
        data = pd.read_csv(data_file)
        res_93_gfp = run_psm(data, 'KO93', 'GFP')
        res_95_gfp = run_psm(data, 'KO95', 'GFP')
        res_dko_gfp = run_psm(data, 'DKO', 'GFP')

        res_93_gfp.to_csv('results/psm_KO93_vs_GFP.csv', index=False)
        res_95_gfp.to_csv('results/psm_KO95_vs_GFP.csv', index=False)
        res_dko_gfp.to_csv('results/psm_DKO_vs_GFP.csv', index=False)

        generate_figures(res_93_gfp, 'psm_KO93_vs_GFP')
        generate_figures(res_95_gfp, 'psm_KO95_vs_GFP')
        generate_figures(res_dko_gfp, 'psm_DKO_vs_GFP')


    def replot_psm():
        res_93_gfp = pd.read_csv('results/psm_KO93_vs_GFP.csv')
        res_95_gfp = pd.read_csv('results/psm_KO95_vs_GFP.csv')
        res_dko_gfp = pd.read_csv('results/psm_DKO_vs_GFP.csv')

        generate_figures(res_93_gfp, 'psm_KO93_vs_GFP')
        generate_figures(res_95_gfp, 'psm_KO95_vs_GFP')
        generate_figures(res_dko_gfp, 'psm_DKO_vs_GFP')

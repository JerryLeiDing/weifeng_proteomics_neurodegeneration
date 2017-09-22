import numpy as np
import pandas as pd
import re, sys, time
sys.dont_write_bytecode = True  # Avoid caching problems

from abc import ABCMeta, abstractmethod
from scripts.anova import *
from scripts.stat_tests import *
from scripts.plots import *


# Example workflow for KO95 vs GFP comparison, using protein report
# data = pd.read_csv('protein.csv')
# ko95_anova, _ = anova_general(data, ko93=False, ko95=True, interaction=False)
# ko95_norm = normalize_plexb(ko95_anova)
# ko95_anova.iloc[:,:9] = ko95_norm.iloc[:,:9]
# res = run_protein(ko95_anova, 'KO95', 'GFP')

# TODO create parent class for all these running classes
#   Will allow better docstring
#   Make abstract class: abstract methods include init, stat_tests, plot_quality???
#   Filename string and docstring

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
                self.protein, group=False, cross=True,
                title='Mar/raw_quality/Mar_protein', col_groups=self.groups)
        compare_channel_replicates(
                self.phospho_nonorm, group=False, cross=True,
                title='Mar/raw_quality/Mar_phospho_nonorm', col_groups=self.groups)
        compare_channel_replicates(
                self.phospho_norm, group=False, cross=True,
                title='Mar/raw_quality/Mar_phospho_norm', col_groups=self.groups)

    def do_stat_test(self, comparison=valid_comparisons[0]):
        """ Runs statistical tests
            Results are stored in self.protein_results and self.phospho_results """
        if comparison not in self.valid_comparisons:
            raise ValueError("Invalid comparison: see self.valid_comparison")

        prot_res = run_protein(self.protein, comparison[0], comparison[1])
        phos_res = run_protein(self.phospho_norm, comparison[0], comparison[1])

        self.protein_results[comparison] = prot_res
        self.phospho_results[comparison] = phos_res

        return (prot_res, phos_res)

    def do_all_stat_tests(self):
        for comp in self.valid_comparison:
            self.do_stat_test(comp)

    def plot_results(self):
        for (k,v) in self.protein_results.iteritems():
            if v is not None:
                generate_figures(v, LABEL_COLS_PROT,
                        filename='Mar/Mar_total_prot_%s_x_%s' % k)
            else:
                print "%s vs %s skipped in March figure gen" % comparison
        for (k,v) in self.phospho_results.iteritems():
            if v is not None:
                generate_figures(v, LABEL_COLS_PROT,
                        filename='Mar/Mar_phospho_%s_x_%s' % k)
            else:
                print "%s vs %s skipped in March figure gen" % comparison


    def write_results(self):
        for (k,v) in self.protein_results.iteritems():
            if v is not None:
                v.to_csv('./results/Mar/Mar_protein_%s_x_%s.csv', index=False)
            else:
                print "%s vs %s skipped in March results saving" % comparison
        for (k,v) in self.phospho_results.iteritems():
            if v is not None:
                v.to_csv('./results/Mar/Mar_phospho_%s_x_%s.csv', index=False)
            else:
                print "%s vs %s skipped in March results saving" % comparison


class RunAug:
    """ Class to encapsulate the August dataset and associated procedures"""
    data_cols = [u'GFP_A1', u'GFP_A2', u'GFP_B1', u'GFP_B2',
                 u'KO95_A1', u'KO95_A2', u'KO95_A3', u'KO95_B1', u'KO95_B2', 
                 u'KO93_A1', u'KO93_A2', u'KO93_B1', u'KO93_B2', u'KO93_B3', 
                 u'DKO_A1', u'DKO_A2', u'DKO_B1', u'DKO_B2'       
                ]      
    # TODO finish this class


class RunTyr:
    """ Class to encapsulate the Tyrosine dataset and associated procedures """
    data_cols = ['KO93_A1', 'KO93_A2', 'KO93_A3',
                'KO95_A1', 'KO95_A2', 'KO95_A3',
                'KOSAP_A1', 'KOSAP_A2',
                'CT_A1', 'CT_A2']
    groups= [
            ('KO93', ['KO93_A1','KO93_A2','KO93_A3']),
            ('KO95', ['KO95_A1','KO95_A2','KO95_A3']),
            ('KOSAP', ['KOSAP_A1','KOSAP_A2']),
            ('CT', ['CT_A1','CT_A2']),
            ]
    # List of (str, str) pairs which define all valid pairwise comparisons
    valid_comparisons = [('KO93', 'CT'), ('KO95', 'CT'), ('KOSAP', 'CT')]
    protein_results = dict.fromkeys(valid_comparisons, None)
    phospho_results = dict.fromkeys(valid_comparisons, None)

    def __init__(self):
        self.phospho_norm = pd.read_csv('data/Tyr/phospho.csv')

    def plot_quality(self):
        """ Plots correlation for protein, phospho, and normalized phospho
            Saved in 'plots/Tyr/raw_quality/Tyr_TYPE_channel_reps_GROUP.png"""
        compare_channel_replicates(
                self.phospho_norm, group=False, cross=True,
                title='Tyr/raw_quality/Tyr_phospho_norm', col_groups=self.groups)

    def do_all_stat_tests(self):
        for comp in self.valid_comparison:
            self.do_stat_test(comp)

    def do_stat_test(self, comparison=valid_comparisons[0]):
        """ Runs statistical tests
            Results are stored in self.protein_results and self.phospho_results """
        if comparison not in self.valid_comparisons:
            raise ValueError("Invalid comparison: see self.valid_comparison")

        phos_res = run_protein(self.phospho_norm, comparison[0], comparison[1])
        self.phospho_results[comparison] = phos_res

        return phos_res

    def do_anova(self, bayes=True):
        """ Runs 1-way ANOVA """

        # Create design matrix
        design = pd.DataFrame.from_items([
            ('Intercept', np.ones(len(self.data_cols))),
            ('KO93', [1,1,1,0,0,0,0,0,0,0]),
            ('KO95', [0,0,0,1,1,1,0,0,0,0]),
            ('KOSAP', [0,0,0,0,0,0,1,1,0,0]),
            ])

        if bayes:
            self.phospho_anova,_ = anova_modt(
                    self.phospho_norm, self.data_cols, design)
        else:
            self.phospho_anova,_ = anova_noreg(
                    self.phospho_norm, self.data_cols, design)

    def plot_results(self):
        # Not used right now: if we get total protein data keep this
        """
        for (k,v) in self.protein_results.iteritems():
            if v is not None:
                generate_figures(v, LABEL_COLS_PROT,
                        filename='Tyr/Tyr_total_prot_%s_x_%s' % k)
            else:
                print "%s vs %s skipped in Tyr figure gen" % comparison
        """
        for (k,v) in self.phospho_results.iteritems():
            if v is not None:
                generate_figures(v, LABEL_COLS_PROT,
                        filename='Tyr/Tyr_phospho_%s_x_%s' % k)
            else:
                print "%s vs %s skipped in Tyr figure gen" % comparison

    def write_results(self):
        # Not used right now: if we get total protein data keep this
        """
        for (k,v) in self.protein_results.iteritems():
            if v is not None:
                v.to_csv('./results/Tyr/Tyr_protein_%s_x_%s.csv', index=False)
            else:
                print "%s vs %s skipped in Tyr results saving" % comparison
        """
        for (k,v) in self.phospho_results.iteritems():
            if v is not None:
                v.to_csv('./results/Tyr/Tyr_phospho_%s_x_%s.csv', index=False)
            else:
                print "%s vs %s skipped in Tyr results saving" % comparison

        if hasattr(self, 'phospho_anova'):
            self.phospho_anova.to_csv(
                    './results/Tyr/Tyr_phospho_ANOVA.csv', 
                    index=False)


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

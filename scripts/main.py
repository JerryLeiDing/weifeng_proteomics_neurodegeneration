import numpy as np
import pandas as pd
import statsmodels.api as sm
import re, sys, time
sys.dont_write_bytecode = True  # Avoid caching problems

from enum import Enum
from statsmodels.formula.api import ols
from statsmodels.sandbox.stats.multicomp import multipletests

from stat_tests import *
from plots import *



def validate_comp_subset_data(data, comp1, comp2):
    """ Validates that comp1 and comp2 define a valid comparison

    data: pandas df with columns GFP/KO93/KO95/DKO as necessary for comparison
          comp1, comp2: one of 'GFP', 'KO93', 'KO95', 'DKO'
    """
    valid = ['GFP', 'KO93', 'KO95', 'DKO']
    if comp1 not in valid or comp2 not in valid:
        raise ValueError('Invalid specification for comparison')
    columns = list(data.columns)

    comp1_pattern = re.compile('^' + comp1 + '_[A-B][0-9]$')
    comp1_cols = [c for c in columns if comp1_pattern.match(c) ]
    comp2_pattern = re.compile('^' + comp2 + '_[A-B][0-9]$')
    comp2_cols = [c for c in columns if comp2_pattern.match(c) ]

    c1 = data[comp1_cols]
    c2 = data[comp2_cols]

    return (c1, c2)


def run_protein(data, comp1, comp2, plex='both'):
    """
    Args:
        data: dataframe to test
        trial = (cond1, cond2), where cond1, cond2 are one of
                'GFP', 'KO93', 'KO95', 'DKO'
    """
    c1, c2 = validate_comp_subset_data(data, comp1, comp2)

    if plex == 'A' or plex == 'B':
        # Filter to corresponding plex
        c1 = c1[[col for col in c1.columns if col[-2] == plex]]
        c2 = c2[[col for col in c2.columns if col[-2] == plex]]
    elif plex == 'both':
        # Do nothing
        pass
    else:
        raise ValueError('Invalid specification of plex')

    pvals = do_stat_tests(c1, c2, True)
    # Delete pvals which are all NaN, i.e. skipped
    # Otherwise adjust pvals
    for c in pvals.columns:
        if pvals[c].isnull().all():
            del pvals[c]
        elif c == u'fold_change':
            continue  # Don't adjust the pval for fold change
        else:
            _, adj, _, _ = multipletests(pvals[c], alpha=0.05, method='fdr_bh')
            pvals[c+'_adj'] = adj

    keep_cols = [
            'numPepsUnique',
            'accession_number',
            'accession_numbers',
            'geneSymbol',
            'entry_name',
        ]
    return pd.concat((c1, c2, data[keep_cols], pvals), axis=1)


### TEMPORARY LOCATION FOR THIS CODE
### TODO: Compare the ANOVA results from eBayes to normal ANOVA without regularization
### TODO: Compare additive vs interaction ANOVA for experimental term
### TODO: Plots for evidence of interaction term? But how to do this in limit of many genes? I don't think the average response will be useful
### Let's just try fitting the ANOVA model
### MOST IMPORTANT: Is there an interaction term here?
### TODO: residual plots, hopefully they do show a positive correlation now
### We'll see!


# NOTE: don't change me unless you also change _make_design_matrix(!)
col_order = [
       u'GFP_A1', u'GFP_A2', u'GFP_B1', u'GFP_B2', u'KO95_A1', u'KO95_A2',
       u'KO95_A3', u'KO95_B1', u'KO95_B2', u'KO93_A1', u'KO93_A2', u'KO93_B1',
       u'KO93_B2', u'KO93_B3', u'DKO_A1', u'DKO_A2', u'DKO_B1', u'DKO_B2'
       ]
def _make_design_matrix(ko93, ko95, interact, samples):
    """
    Key for ko93 and ko95
    0=exclude, 1=no interaction term, 2=use interaction term

    Always returns an 18xCOND dataframe
    Drop rows with all zeros to yield a usable 
    """
    # if ko93 not in [0,1,2] or ko95 not in [0,1,2]:
    #     raise ValueError("ko93 and ko95 must be one of 0,1,2 in anova_general")
    # 
    # if ko93 + ko95 == 0:
    #     raise ValueError("At least one of ko93 and ko95 must be nonzero")

    design_cols = []
    coef_labels = []
    n = len(samples)
   
    # Determine correct indices for each ko marker
    ko93_col_idx = [i for i,col in enumerate(samples) if
            ('KO93' in col) or ('DKO' in col)]
    ko95_col_idx = [i for i,col in enumerate(samples) if
            ('KO95' in col) or ('DKO' in col)]
    plexB_col_idx = [i for i,col in enumerate(samples) if col[-2] == 'B']

    # Intercept term
    coef_labels.append('Intercept')
    design_cols.append(np.ones(n, dtype=int))
    # PlexB term
    coef_labels.append('PlexB')
    plexB = np.zeros(n, dtype=int)
    plexB[plexB_col_idx] = 1
    design_cols.append(plexB)

    if ko93:
        coef_labels.append('KO93')
        ko93_col = np.zeros(n, dtype=int)
        ko93_col[ko93_col_idx] = 1
        design_cols.append(ko93_col)
    if ko95:
        coef_labels.append('KO95')
        ko95_col = np.zeros(n, dtype=int)
        ko95_col[ko95_col_idx] = 1
        design_cols.append(ko95_col)

    if interact and ko93 and ko95:
        coef_labels.append('Interaction')
        interact = np.zeros(n, dtype=int)
        interact[[x for x in ko93_col_idx if x in ko95_col_idx]] = 1
        design_cols.append(interact)
    elif interact:
        print "Warning: interact has no effect if both ko93 and ko95 are not true"

    return pd.DataFrame.from_items(zip(coef_labels, design_cols))


def anova_general(df, ko93=True, ko95=True, interact=True):
    # TODO documentation
    # Mention that blocking term is additive
    
    # Select relevant columns for the comparison
    cols = [0,1,2,3]
    if ko95:
        cols += [4,5,6,7,8]
    if ko93:
        cols += [9,10,11,12,13]
    if ko95 and ko93:
        cols += [14,15,16,17]

    # Set up design, coefficients, and subset of data
    columns = np.array(col_order)[cols]
    data = df[columns]
    design = _make_design_matrix(ko93, ko95, interact, columns)
    coefs = list(design.columns)

    # Note that result is an R object
    # We can't do much with it directly except call topTable
    result = r['moderated.t'](data, design=design)

    # Now obtain best estimates for each coefficient
    res_coef = pandas2ri.ri2py(r['topTable'](
            result,
            number=data.shape[0],
            sort_by='none')).iloc[:,:len(coefs)]
    # F-test for significance for terms OTHER than intercept and PlexB
    coefs.remove('Intercept')
    coefs.remove('PlexB')
    res_f = pandas2ri.ri2py(r['topTable'](
            result,
            coef=np.array(coefs),
            number=data.shape[0],
            sort_by='none'))[['F', 'P.Value', 'adj.P.Val']]
    # Now bind together everything into one df
    aux_data = (
            df[['accession_number', 'geneSymbol', 'entry_name', 'numPepsUnique']].
            reset_index(drop=True))
    data.reset_index(drop=True, inplace=True)
    res_f.reset_index(drop=True, inplace=True)
    res_coef.reset_index(drop=True, inplace=True)
    res_df = pd.concat((data, res_coef, res_f, aux_data), axis=1)
    return res_df, result


# NOTE: deprecated in favor of anova_general
def anova_2way_interaction(df, use='KO93'):
    if use == 'KO93':
        ko93 = True
        KO = 'KO93'
    elif use == 'KO95':
        ko93 = False
        KO = 'KO95'
    else:
        raise ValueError("Invalid specification for use. Must be 'KO93' or 'KO95'")
    # Order data columns correctly
    if ko93:
        data = df[np.array(col_order)[[0,1,2,3,9,10,11,12,13]]]
    else:
        data = df[np.array(col_order)[[0,1,2,3,4,5,6,7,8]]]

    # Design matrix has 4 terms
    # Intercept, KO95, PlexB, both
    design_arr = np.zeros((data.shape[1], 4), dtype=int)
    design_arr[:,0] = 1
    design_arr[:,1] = np.array([0,0,0,0,1,1,1,1,1])  # KO term
    if ko93:
        design_arr[:,2] = np.array([0,0,1,1,0,0,1,1,1])  # PlexB term
        design_arr[:,3] = np.array([0,0,0,0,0,0,1,1,1])  # Interaction term
    else:
        design_arr[:,2] = np.array([0,0,1,1,0,0,0,1,1])  # PlexB term
        design_arr[:,3] = np.array([0,0,0,0,0,0,0,1,1])  # Interaction term
    design = pd.DataFrame(design_arr)
    # Interaction term coefficient is 'KO9[x]xPlexB'
    design.columns = ["Intercept", KO, "PlexB", KO+'xPlexB']

    # Note that result is an R object
    # We can't do much with it directly except call topTable
    result = r['moderated.t'](data, design=design)

    # Now obtain best estimates for each coefficient
    res_coef = pandas2ri.ri2py(r['topTable'](
            result,
            number=data.shape[0],
            sort_by='none')).iloc[:,:4]
    # F-test for significance for OTHER three terms
    res_f = pandas2ri.ri2py(r['topTable'](
            result,
            coef=np.array([KO, 'PlexB', KO+'xPlexB']),
            number=data.shape[0],
            sort_by='none'))[['F', 'P.Value', 'adj.P.Val']]
    # Now bind together everything into one df
    aux_data = (
            df[['accession_number', 'geneSymbol', 'entry_name', 'numPepsUnique']].
            reset_index(drop=True))
    data.reset_index(drop=True, inplace=True)
    res_f.reset_index(drop=True, inplace=True)
    res_coef.reset_index(drop=True, inplace=True)
    res_df = pd.concat((data, res_coef, res_f, aux_data), axis=1)
    return res_df, result


# NOTE: deprecated in favor of anova_general
def anova_2way_nointeraction(df, use='KO93'):
    if use == 'KO93':
        ko93 = True
        KO = 'KO93'
    elif use == 'KO95':
        ko93 = False
        KO = 'KO95'
    else:
        raise ValueError("Invalid specification for use. Must be 'KO93' or 'KO95'")
    # Order data columns correctly, usingn only GFP/KO95 cols
    if ko93:
        data = df[np.array(col_order)[[0,1,2,3,9,10,11,12,13]]]
    else:
        data = df[np.array(col_order)[[0,1,2,3,4,5,6,7,8]]]

    # Design matrix has 4 terms
    # Intercept, KO95, PlexB, both
    design_arr = np.zeros((data.shape[1], 3), dtype=int)
    design_arr[:,0] = 1
    design_arr[:,1] = np.array([0,0,0,0,1,1,1,1,1])  # KO term
    if ko93:
        design_arr[:,2] = np.array([0,0,1,1,0,0,1,1,1])  # PlexB term
    else:
        design_arr[:,2] = np.array([0,0,1,1,0,0,0,1,1])  # PlexB term
    design = pd.DataFrame(design_arr)
    design.columns = ["Intercept", KO, "PlexB"]

    # Note that result is an R object
    # We can't do much with it directly except call topTable
    result = r['moderated.t'](data, design=design)

    # Now obtain best estimates for each coefficient
    res_coef = pandas2ri.ri2py(r['topTable'](
            result,
            number=data.shape[0],
            sort_by='none')).iloc[:,:3]
    # F-test for significance for OTHER three terms
    res_f = pandas2ri.ri2py(r['topTable'](
            result,
            coef=np.array([KO, 'PlexB']),
            number=data.shape[0],
            sort_by='none'))[['F', 'P.Value', 'adj.P.Val']]
    # Now bind together everything into one df
    aux_data = (
            df[['accession_number', 'geneSymbol', 'entry_name', 'numPepsUnique']].
            reset_index(drop=True))
    data.reset_index(drop=True, inplace=True)
    res_f.reset_index(drop=True, inplace=True)
    res_coef.reset_index(drop=True, inplace=True)
    res_df = pd.concat((data, res_coef, res_f, aux_data), axis=1)
    return res_df, result


def normalize_plexB(anova_res_df):
    """ Requires columns PlexB, accession_number
    Can optionally include interaction terms [KO95xPlexB, KO93xPlexB, KO93xKO95xPlexB]

    This should subtract out the effects from the PlexB and interaction terms
    E[y] = INT + C_(ko95)X_(ko95) + C_(plxb)X_(plxb) + C_(ko95*plxb)X_(plxb)X_(ko95)
    E[y] - C_(plxb)X_(plxb) = INT + (C_(ko95) + C_(ko95*plxb)X_(plxb))X_ko95
    # For KO95, X_(ko95) === 1 always
    # Therefore, we can write
    E[y] - C_(plxb)X_(plxb) - C_(ko95*plxb)X_(plxb) = INT + C_(ko95)
    # So these should be correlated. Likewise for GFP, but no interaction term
    """
    normed_data = pd.DataFrame(index=np.arange(anova_res_df.shape[0]))
    # normed_data = anova_res_df.copy()
    # TODO: modify array in-place (preserving other columns?)

    if 'GFP_A1' in anova_res_df.columns:
        # Normalize GFP PlexB channels
        normed_data['GFP_A1'] = anova_res_df.GFP_A1
        normed_data['GFP_A2'] = anova_res_df.GFP_A2
        normed_data['GFP_B1'] = anova_res_df.GFP_B1 - anova_res_df.PlexB
        normed_data['GFP_B2'] = anova_res_df.GFP_B2 - anova_res_df.PlexB

    if 'KO95_A1' in anova_res_df.columns:
        # Normalize KO95 PlexB channels
        normed_data['KO95_A1'] = anova_res_df.KO95_A1
        normed_data['KO95_A2'] = anova_res_df.KO95_A2
        normed_data['KO95_A3'] = anova_res_df.KO95_A3
        normed_data['KO95_B1'] = anova_res_df.KO95_B1 - anova_res_df.PlexB
        normed_data['KO95_B2'] = anova_res_df.KO95_B2 - anova_res_df.PlexB
        if 'KO95xPlexB' in anova_res_df.columns:
            # If present, use interaction term between PlexB and KO95
            normed_data['KO95_B1'] -= anova_res_df.KO95xPlexB
            normed_data['KO95_B2'] -= anova_res_df.KO95xPlexB

    if 'KO93_A1' in anova_res_df.columns:
        # Normalize KO93 PlexB channels
        normed_data['KO93_A1'] = anova_res_df.KO93_A1
        normed_data['KO93_A2'] = anova_res_df.KO93_A2
        normed_data['KO93_B1'] = anova_res_df.KO93_B1 - anova_res_df.PlexB
        normed_data['KO93_B2'] = anova_res_df.KO93_B2 - anova_res_df.PlexB
        normed_data['KO93_B3'] = anova_res_df.KO93_B3 - anova_res_df.PlexB
        if 'KO93xPlexB' in anova_res_df.columns:
            # If present, use interaction term between PlexB and KO93
            normed_data['KO93_B1'] -= anova_res_df.KO93xPlexB
            normed_data['KO93_B2'] -= anova_res_df.KO93xPlexB
            normed_data['KO93_B3'] -= anova_res_df.KO93xPlexB

    if 'DKO_A1' in anova_res_df.columns:
        # Normalize DKO PlexB channels
        normed_data['DKO_A1'] = anova_res_df.DKO_A1
        normed_data['DKO_A2'] = anova_res_df.DKO_A2
        normed_data['DKO_B1'] = anova_res_df.DKO_B1 - anova_res_df.PlexB
        normed_data['DKO_B2'] = anova_res_df.DKO_B2 - anova_res_df.PlexB
        # If present, use double interaction terms
        # 93-Plex, 95-PlexB, 93-95-PlexB
        if 'KO95xPlexB' in anova_res_df.columns:
            # If present, use interaction term between PlexB and KO95
            normed_data['DKO_B1'] -= anova_res_df.KO95xPlexB
            normed_data['DKO_B2'] -= anova_res_df.KO95xPlexB
        if 'KO93xPlexB' in anova_res_df.columns:
            # If present, use interaction term between PlexB and KO93
            normed_data['DKO_B1'] -= anova_res_df.KO93xPlexB
            normed_data['DKO_B2'] -= anova_res_df.KO93xPlexB
        if 'KO93xKO95xPlexB' in anova_res_df.columns:
            # If present, use triple interaction term
            normed_data['DKO_B1'] -= anova_res_df.KO93xKO95xPlexB
            normed_data['DKO_B2'] -= anova_res_df.KO93xKO95xPlexB

    normed_data['accession_number'] = anova_res_df.accession_number
    return normed_data


### END TEMPORARY CODE

def run_psm(data, comp1, comp2, plex='both'):
    """
    Use PSM and roll up to the level of unique accession_number
    """
    start = time.time()
    c1, c2 = validate_comp_subset_data(data, comp1, comp2)

    if plex == 'A' or plex == 'B':
        # Filter to corresponding plex
        c1 = c1[[col for col in c1.columns if col[-2] == plex]]
        c2 = c2[[col for col in c2.columns if col[-2] == plex]]
    elif plex == 'both':
        # Do nothing
        pass
    else:
        raise ValueError('Invalid specification of plex')

    pvals = do_stat_tests_protein(c1, c2, data.accession_number)
    # Delete pvals which are all NaN, i.e. skipped
    # Otherwise adjust pvals
    for c in pvals.columns:
        if pvals[c].isnull().all():
            del pvals[c]
        elif c == u'fold_change_med':
            continue  # Don't adjust the pval for fold change
        elif pvals[c].dtype == np.number:
            # Mask NaN values so we don't bias the adjusted test
            pv = pvals[c]
            mask = np.isfinite(pv)
            pv_corr = np.full(pv.shape, np.nan)
            pv_corr[mask] = multipletests(
                    pv[mask], alpha=0.05, method='fdr_bh')[1]
            pvals[c+'_adj'] = pv_corr

    pvals.rename(columns={'protein_id': 'accession_number',
                          'fold_change_med': 'fold_change'},
                 inplace=True)
    # Make auxiliary info
    aux_info = data[[
        'accession_number',
        'geneSymbol',
    ]]
    aux_info.drop_duplicates(inplace=True)

    def get_group_counts(x):
        return pd.Series({'n_pep': len(x),
                          'n_valid': np.sum(1 - np.isnan(x).values)})
    tmp = (pd.concat((c1, c2), axis=1)
            .groupby(data.accession_number)
            .apply(get_group_counts)
            .reset_index())

    out = pd.merge(pvals, aux_info, on='accession_number')
    out = pd.merge(out, tmp, on='accession_number')

    print time.time() - start
    return out


def run_peptide(data, comp1, comp2):
    # TODO
    pass


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

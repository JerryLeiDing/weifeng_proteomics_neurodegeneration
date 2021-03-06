import numpy as np
import pandas as pd
import statsmodels.api as sm
import re, sys, time
sys.dont_write_bytecode = True  # Avoid caching problems

from statsmodels.formula.api import ols
from statsmodels.sandbox.stats.multicomp import multipletests

from stat_tests import *
from plots import *


ALPHA = 0.05  # Significance threshold for FDR

def validate_comp_subset_data(data, comp1, comp2):
    """ Validates that comp1 and comp2 define a valid comparison

    data: pandas df with columns GFP/KO93/KO95/DKO as necessary for comparison
          comp1, comp2: one of 'GFP', 'KO93', 'KO95', 'DKO'
    """
    # Old validation logic
    # valid = ['GFP', 'KO93', 'KO95', 'DKO']
    # if comp1 not in valid or comp2 not in valid:
    #     raise ValueError('Invalid specification for comparison')
    columns = list(data.columns)

    comp1_pattern = re.compile('^' + comp1 + '_[A-B][0-9]$')
    comp1_cols = [c for c in columns if comp1_pattern.match(c) ]
    comp2_pattern = re.compile('^' + comp2 + '_[A-B][0-9]$')
    comp2_cols = [c for c in columns if comp2_pattern.match(c) ]

    if len(comp1_cols) == 0:
        raise ValueError('No columns found for condition %s' % comp1)
    if len(comp2_cols) == 0:
        raise ValueError('No columns found for condition %s' % comp2)

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
            _, adj, _, _ = multipletests(pvals[c], alpha=ALPHA, method='fdr_bh')
            pvals[c+'_adj'] = adj

    # TODO check for optional columns
    # Should really make this a global option in a settings file somewhere
    # OR just keep everything that isn't explicitly removed?
    # Now bind together everything into one df
    aux_data = data.drop(list(c1.columns) + list(c2.columns), axis=1)
    return pd.concat((c1, c2, aux_data, pvals), axis=1)


def anova_modt(df, columns, design):
    """ Runs ANOVA on the subset of df defined by columns with the specified design matrix, using the moderated T linear model.

    Args:
        df (Pandas DataFrame): DataFrame with one column per measurement, rows=proteins
        columns (list(columns)): list of column names in df which have data
        design (Pandas DataFrame)): design matrix (see limma documentation for details)
            Note that the columns names of design are the returned coefficient names

    Returns:
        (res_df, result)
        res_df (Pandas DataFrame): DataFrame with one row per protein  
            - Same order as df
            - One column for each best fit coefficient
            - Has columns 'F_<COEF>' and 'PVal_<COEF>' for each coefficient
        result (R object)
    """
    coefs = list(design.columns)
    data = df[columns]
    # Note that result is an R object
    # We can't do much with it directly except call topTable
    result = r['moderated.t'](data, design=design)

    # Now obtain best estimates for each coefficient
    res_coef = pandas2ri.ri2py(r['topTable'](
            result,
            number=data.shape[0],
            sort_by='none')).iloc[:,:len(coefs)+3]
    # Adjust overall p-value
    res_coef['P.Value.Adj'] = multipletests(
            res_coef['P.Value'], alpha=ALPHA, method='fdr_bh')[1]
    # F-test for significance for terms OTHER than intercept and PlexB
    # Do this iteratively and obtain a p-value and F-value for every coefficient
    coefs.remove('Intercept')
    # Create mapping of coefficients to F and PVal columns
    coef_col_map = {c: ['F_%s'%c, 'PVal_%s'%c, 'PVal_%s_Adj'%c] for c in coefs}
    result_colnames = [col for cols in coef_col_map.values() for col in cols]
    # Create empty pvalue df
    res_f = pd.DataFrame(
            index=np.arange(data.shape[0]),
            columns=result_colnames,
            dtype=float)
    for c in coef_col_map.keys():
        # Find F and pvals/adj_pvals for each coefficient
        F_pv = pandas2ri.ri2py(r['topTable'](
                result,
                coef=c,
                number=data.shape[0],
                sort_by='none'))[['t', 'P.Value']]
        _, pv_adj, _, _ = multipletests(F_pv['P.Value'], alpha=ALPHA, method='fdr_bh')
        res_f[coef_col_map[c]] = np.concatenate(
                (F_pv.values, pv_adj[:,np.newaxis]), axis=1)

    # Now bind together everything into one df
    aux_data = df.drop(columns, axis=1).reset_index(drop=True)
    data.reset_index(drop=True, inplace=True)
    res_f.reset_index(drop=True, inplace=True)
    res_coef.reset_index(drop=True, inplace=True)
    res_df = pd.concat((data, res_coef, res_f, aux_data), axis=1)
    return res_df, result


def anova_noreg(df, columns, design):
    """ Runs ANOVA on the subset of df defined by columns with the specified design matrix

    Args:
        df (Pandas DataFrame): DataFrame with one column per measurement, rows=proteins
        columns (list(columns)): list of column names in df which have data
        design (Pandas DataFrame)): design matrix (see limma documentation for details)
            Note that the columns names of design are the returned coefficient names

    Returns:
        (res_df, result)
        res_df (Pandas DataFrame): DataFrame with one row per protein  
            - Same order as df
            - One column for each best fit coefficient
            - Has columns 'F_<COEF>' and 'PVal_<COEF>' for each coefficient
        result (R object)
    """
    # Set up design, coefficients, and subset of data
    data = df[columns]
    coefs = list(design.columns)
    # Make copy of coefs and remove intercept
    coefs_no_intercept = list(coefs)
    coefs_no_intercept.remove('Intercept')

    # Create mapping of coefficients to F and PVal columns
    coef_col_map = {c: ['F_'+c, 'PVal_'+c] for c in coefs_no_intercept}
    result_colnames = coefs + [y for x in coef_col_map.values() for y in x]
    # Create empty result df
    result = pd.DataFrame(
            index=np.arange(data.shape[0]),
            columns=result_colnames,
            dtype=float)
    # Make copy of design dataframe and add a column
    for_anova = design.copy()
    # Now for each row in data, run ANOVA separately
    for i in xrange(data.shape[0]):
        if i % 1000 == 0:
            print i
        for_anova['Y'] = data.iloc[i].values.astype(float)
        fit = ols('Y ~ ' + ' + '.join(coefs_no_intercept), for_anova).fit()
        # Set coefficients
        result.iloc[i][coefs] = fit.params.values
        # Extract ANOVA p-values and F-values
        anova_res = sm.stats.anova_lm(fit, typ=2)
        # Map to correct result columns
        for c in anova_res.index:
            if c in coef_col_map:
                result.iloc[i][coef_col_map[c]]= anova_res.loc[c][-2:].values
    # Adjust pvals
    pv_adj = pd.DataFrame({
        'PVal_%s_Adj' % coef: multipletests(
            result['PVal_%s'%coef], alpha=ALPHA, method='fdr_bh')[1]
         for coef in coefs_no_intercept
        })
    # Now bind together everything into one df
    aux_data = df.drop(columns, axis=1).reset_index(drop=True)
    data.reset_index(drop=True, inplace=True)
    res_df = pd.concat((data, result, pv_adj, aux_data), axis=1)
    return res_df

# The following code is specifically for the 2x2 ANOVA in the August dataset
# NOTE: don't change me unless you also change _make_design_matrix(!)
col_order = [
       u'GFP_A1', u'GFP_A2', u'GFP_B1', u'GFP_B2', u'KO95_A1', u'KO95_A2',
       u'KO95_A3', u'KO95_B1', u'KO95_B2', u'KO93_A1', u'KO93_A2', u'KO93_B1',
       u'KO93_B2', u'KO93_B3', u'DKO_A1', u'DKO_A2', u'DKO_B1', u'DKO_B2'
       ]
def _make_design_matrix(samples, ko93, ko95, interact, plexb_interact):
    """
    Key for ko93 and ko95
    0=exclude, 1=no interaction term, 2=use interaction term

    Always returns an 18xCOND dataframe
    Drop rows with all zeros to yield a usable 
    """

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

    if plexb_interact and ko93:
        coef_labels.append('KO93xPleXB')
        interact = np.zeros(n, dtype=int)
        interact[[x for x in ko93_col_idx if x in plexB_col_idx]] = 1
        design_cols.append(interact)
    if plexb_interact and ko95:
        coef_labels.append('KO95xPleXB')
        interact = np.zeros(n, dtype=int)
        interact[[x for x in ko95_col_idx if x in plexB_col_idx]] = 1
        design_cols.append(interact)

    return pd.DataFrame.from_items(zip(coef_labels, design_cols))


def anova_aug(df, ko93=True, ko95=True, interact=True, plexb_interact=False):
    # TODO documentation
    # Mention that blocking term is additive
    
    # Select relevant columns for the comparison
    cols = [0,1,2,3]
    if ko95:
        cols += [4,5,6,7,8]
    if ko93:
        cols += [9,10,11,12,13]
    if ko95 and ko93 and interact:  # Drop DKO unless interact
        cols += [14,15,16,17]

    # Set up design, coefficients, and subset of data
    columns = np.array(col_order)[cols]
    data = df[columns]
    design = _make_design_matrix(columns, ko93, ko95, interact, plexb_interact)
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
    # Do this iteratively and obtain a p-value and F-value for every coefficient
    coefs.remove('Intercept')
    # Create mapping of coefficients to F and PVal columns
    coef_col_map = {c: ['F_'+c, 'PVal_'+c] for c in coefs}
    result_colnames = [y for x in coef_col_map.values() for y in x]
    # Create empty pvalue df
    res_f = pd.DataFrame(
            index=np.arange(data.shape[0]),
            columns=result_colnames,
            dtype=float)
    for c in coef_col_map.keys():
        res_f[coef_col_map[c]] = pandas2ri.ri2py(r['topTable'](
                result,
                coef=c,
                number=data.shape[0],
                sort_by='none'))[['t', 'P.Value']]
    # Now bind together everything into one df
    aux_data = (
            df[['accession_number', 'geneSymbol', 'entry_name', 'numPepsUnique']].
            reset_index(drop=True))
    data.reset_index(drop=True, inplace=True)
    res_f.reset_index(drop=True, inplace=True)
    res_coef.reset_index(drop=True, inplace=True)
    res_df = pd.concat((data, res_coef, res_f, aux_data), axis=1)
    return res_df, result


def anova_noreg_aug(df, ko93=True, ko95=True, interact=True, plexb_interact=True):
    # TODO documentation
    # Mention that blocking term is additive
    
    # Select relevant columns for the comparison
    cols = [0,1,2,3]
    if ko95:
        cols += [4,5,6,7,8]
    if ko93:
        cols += [9,10,11,12,13]
    if ko95 and ko93 and interact:  # Drop DKO unless interact
        cols += [14,15,16,17]

    # Set up design, coefficients, and subset of data
    columns = np.array(col_order)[cols]
    data = df[columns]
    design = _make_design_matrix(ko93, ko95, interact, columns)
    coefs = list(design.columns)
    # Make copy of coefs and remove intercept
    coefs_no_intercept = list(coefs)
    coefs_no_intercept.remove('Intercept')

    # Create mapping of coefficients to F and PVal columns
    coef_col_map = {c: ['F_'+c, 'PVal_'+c] for c in coefs_no_intercept}
    result_colnames = coefs + [y for x in coef_col_map.values() for y in x]
    # Create empty result df
    result = pd.DataFrame(
            index=np.arange(data.shape[0]),
            columns=result_colnames,
            dtype=float)
    # Make copy of design dataframe and add a column
    for_anova = design.copy()
    # Now for each row in data, run ANOVA separately
    for i in xrange(data.shape[0]):
        if i % 100 == 0:
            print i
        for_anova['Y'] = data.iloc[i].values.astype(float)
        fit = ols('Y ~ ' + ' + '.join(coefs_no_intercept), for_anova).fit()
        # Set coefficients
        result.iloc[i][coefs] = fit.params.values
        # Extract ANOVA p-values and F-values
        anova_res = sm.stats.anova_lm(fit, typ=2)
        # Map to correct result columns
        for c in anova_res.index:
            if c in coef_col_map:
                result.iloc[i][coef_col_map[c]]= anova_res.loc[c][-2:].values
    # Now bind together everything into one df
    aux_data = (
            df[['accession_number', 'geneSymbol', 'entry_name', 'numPepsUnique']].
            reset_index(drop=True))
    data.reset_index(drop=True, inplace=True)
    res_df = pd.concat((data, result, aux_data), axis=1)
    return res_df
        


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
        print "Normalizing GFP"
        normed_data['GFP_A1'] = anova_res_df.GFP_A1
        normed_data['GFP_A2'] = anova_res_df.GFP_A2
        normed_data['GFP_B1'] = anova_res_df.GFP_B1 - anova_res_df.PlexB
        normed_data['GFP_B2'] = anova_res_df.GFP_B2 - anova_res_df.PlexB

    if 'KO95_A1' in anova_res_df.columns:
        # Normalize KO95 PlexB channels
        print "Normalizing KO95"
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
        print "Normalizing KO93"
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
        print "Normalizing DKO"
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


def compare_quality_before_after_anova(df, *args):

    # 93_interaction = anova_2way_interaction(df, use='KO93')[0]
    # 95_interaction = anova_2way_interaction(df, use='KO95')[0]
    # 93_nointeraction = anova_2way_nointeraction(df, use='KO93')[0]
    # 95_nointeraction = anova_2way_nointeraction(df, use='KO95')[0]
    full_nointeraction, _ = anova_general(df, *args)
    normalized = normalize_plexB(full_nointeraction)
    orig_norm_col = normalized.columns
    normalized.columns = [col + '_NORM' for col in normalized.columns]

    normalized.reset_index(drop=True, inplace=True)
    df_subset = df[[c for c in df.columns if c in orig_norm_col]].reset_index(drop=True)

    final_df = pd.concat((df_subset, normalized), axis=1)
    if 'numPepsUnique' in df.columns:
        final_df['numPepsUnique'] = df.numPepsUnique

    def plot_comparison(cols=[0,1,2,3]):
        label_cols = final_df.columns[cols]
        n = len(label_cols)
        f = compare_measures(final_df, label_cols, count=False, equal=True)
        label_cols = [col + '_NORM' for col in label_cols] 
        label_cols2 = label_cols
        axarr = np.reshape(f.axes, (n, n))

        if 'numPepsUnique' in final_df.columns: 
            alpha = np.log2(final_df.numPepsUnique)/np.max(np.log2(final_df.numPepsUnique)) 
        elif 'n_pep' in final_df.columns:
            alpha = np.log2(final_df.n_pep) / np.max(np.log2(final_df.n_pep))
        else:                                                                       
            alpha = np.full((final_df.shape[0],), 1.0)

        for i in xrange(n):
            for j in xrange(n):
                if i <= j:
                    continue
                mask = (np.isfinite(final_df[label_cols[i]]) &
                    np.isfinite(final_df[label_cols2[j]]))                           
                val = np.sum(mask)
                # Set alpha correctly                                               
                colors = matplotlib.cm.Reds(np.arange(0, 1, 1./val))                
                # colors = np.zeros((n, 4)) 
                # colors[:,0] = 0.9  # Make red 
                # colors[:,1] = 0.7 # Make red
                # colors[:,2] = 0.6 # Make red                                      
                colors[:,-1] = alpha[mask]                                          
                axarr[i][j].scatter(   # Note that the x-y order is MIRRORED
                        final_df[label_cols[i]].values[mask],
                        final_df[label_cols2[j]].values[mask],
                        color=colors, lw=0)     
        # Set axes for better visualization 
        f.axes[0].set_ylim(-5, 5)
        f.axes[0].set_xlim(-5, 5)
        return f
    return plot_comparison


def find_uncorrelated_gfp_proteins(normed):
    # Here, we treat all proteins with mean GFP.PlexB of >= 1 as uncorrelated
    # TODO: we can also try getting multiple categories from each channel

    is_uncorr = (np.mean(normed[['GFP_B1', 'GFP_B2']], axis=1) >= 1.0)
    return normed[is_uncorr]


### END ANOVA CODE

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
    raise ValueError("Not implemented!")

# Example workflow for KO95 vs GFP comparison, using protein report 
# data = pd.read_csv('protein.csv')
# ko95_anova, _ = anova_general(data, ko93=False, ko95=True, interaction=False)
# ko95_norm = normalize_plexb(ko95_anova)
# ko95_anova.iloc[:,:9] = ko95_norm.iloc[:,:9]
# res = run_protein(ko95_anova, 'KO95', 'GFP')


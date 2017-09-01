import math
import numbers

import numpy as np
import pandas as pd
import rpy2
import scipy as sp
import statsmodels.api as sm

from rpy2.robjects import r
from rpy2.robjects import pandas2ri
from scipy import stats
from statsmodels.tools import add_constant


#################
#     SETUP     #
#################

# Setup R environment
pandas2ri.activate()
r['source']('modT.r')
r['source']('bayesreg.R')
# INSTRUCTIONS FOR USE
# Call functions r['modT_test'], r['bayesT'], r['find_protein_medians']
# r['modT_test'](data, "placeholder", dframe=TRUE, data.col=[DATA COLS])
# r['bayesT'](data, numC, numE)
# r['find_protein_medians'](pepdf)
# r['proteinBayesT'](data, numC, numE, pool_intensity=True)


########################
#  STATISTICAL TESTS   #
########################

def modT(ctrl, exp):
    """
    NOTE: only defined if ncol(ctrl) == 1 or ncol(ctrl) = ncol(exp)
    Be sure ctrl and exp are presented in the same order
    both ctrl and exp should be log2 ratios
    """
    if ctrl.shape[0] != exp.shape[0]:
        raise ValueError('Length of exp and ctrl data frames not identical')

    if ctrl.shape[1] == 1 or ctrl.shape[1] == exp.shape[1]:
        data = pd.DataFrame(exp.values - ctrl.values)
    else:
        raise ValueError('Not valid number of ctrl columns, see documentation')
    # Also create id col for convenience
    data.columns = ['Log2_Ratio_%d' % i for i in xrange(data.shape[1])]
    data_cols = data.columns
    data['id'] = np.arange(len(data))

    res = r['modT_test'](data, "placeholder", id_col='id', data_col=data_cols, dframe=True)
    res = pandas2ri.ri2py(res)
    res['std_err'] = (res['CI.R'] - res['CI.L'])/3.92

    return res


def modT_2sample(ctrl, exp):
    """
    Runs moderated T with 2 sample t test
    """
    if ctrl.shape[0] != exp.shape[0]:
        raise ValueError('Length of exp and ctrl data frames not identical')

    data = pd.DataFrame(np.concatenate((ctrl.values, exp.values), axis=1))
    data.columns = (['C%d' % (i+1) for i in xrange(ctrl.shape[1])] +
            ['E%d' % (i+1) for i in xrange(exp.shape[1])])
    data_cols= data.columns
    data['id'] = np.arange(len(data))

    design = np.array(([-1] * ctrl.shape[1]) + ([1] * exp.shape[1]), dtype=int)

    res = r['modT_test'](data, "placeholder", id_col='id', data_col=data_cols, dframe=True, design=design)
    res = pandas2ri.ri2py(res)
    res['std_err'] = (res['CI.R'] - res['CI.L'])/3.92

    return res


def cyberT(ctrl, exp, **kwargs):
    if ctrl.shape[0] != exp.shape[0]:
        raise ValueError('Length of exp and ctrl data frames not identical')
    ctrl.reset_index(drop=True, inplace=True)
    exp.reset_index(drop=True, inplace=True)

    df = pd.concat([ctrl, exp], axis=1, ignore_index=True)
    df.columns = (['C%d' % (i+1) for i in xrange(ctrl.shape[1])] +
            ['E%d' % (i+1) for i in xrange(exp.shape[1])])
    res = r['bayesT'](df, numC = ctrl.shape[1], numE = exp.shape[1], ppde=False, doMulttest=True, **kwargs)
    res = pandas2ri.ri2py(res)

    return res


def t_test(ctrl, exp, ratio_test = False):
    """
    Ratio test flag indicates if t test should be performed as 1-sample t test
    """

    if ctrl.shape[0] != exp.shape[0]:
        raise ValueError('Length of exp and ctrl data frames not identical')

    if ratio_test:
        if ctrl.shape[1] == 1 or ctrl.shape[1] == exp.shape[1]:
            data = np.array(exp.values - ctrl.values)
            _, pvals = stats.ttest_1samp(data, 0, axis=1, nan_policy='omit')
        else:
            raise ValueError('Not valid number of ctrl columns, see documentation')
    else:
        _, pvals = stats.ttest_ind(ctrl, exp, axis=1, nan_policy='omit')

    return pvals


def protein_wls_test(ctrl, exp, protein, onesample=True, use_bayes=False, reg=10, adj_var=False):
    # Do cyberT
    cyberT_peps = cyberT(ctrl, exp)
    cyberT_peps['protein'] = protein
    return protein_wls_test_cyberT(
            cyberT_peps, onesample=onesample, use_bayes=use_bayes, reg=reg, adj_var=adj_var)


def protein_wls_test_cyberT(cyberT_peps, onesample=True, use_bayes=False, reg=10, adj_var=False):
    """ TODO documentation
    Especially what columns are required from cyberT_peps
    and what columns are in output

    reg = number >= 1: Default 10, prior weight on regularization std
                Set to None for no regularization
    """
    # Group by protein
    grouped = cyberT_peps.groupby('protein', sort=False)
    res = pd.DataFrame(
            columns=['protein_id','fold_change','raw_std_err', 'pval_no_reg', 'n'],
            index=np.arange(len(grouped)))
    res['protein_id'] = res['protein_id'].astype(str)
    res['fold_change'] = res['fold_change'].astype(float)
    res['raw_std_err'] = res['raw_std_err'].astype(float)
    res['pval_no_reg'] = res['pval_no_reg'].astype(float)
    res['n'] = 0
    for j, (name, group) in enumerate(grouped):
        # Do WLS
        if onesample:
            fold_change, std_err, pval = _wls_onesample(
                    group, use_bayes=use_bayes, adj_var=adj_var)
        else:
            fold_change, std_err, pval = _wls(
                    group, use_bayes=use_bayes)

        nC, nE = cyberT_peps['nC'].astype(float), cyberT_peps['nE'].astype(float)
        if use_bayes:
            stdC, stdE = cyberT_peps['bayesSDC'], cyberT_peps['bayesSDE']
        else:
            stdC, stdE = cyberT_peps['stdC'], cyberT_peps['stdE']
        res.ix[j] = [name, fold_change, std_err, pval, group.shape[0]]
        # If we regularize, adjust the std
    
    # Add regularization
    # Smooth the protein level variances by moving toward mean
    if reg is not None:
        # Rolling average of the protein level standard deviations
        WINSIZE = 101  # Window size: must be odd positive integer
        m = res.shape[0]
        if m <= WINSIZE:
            p_sd = np.ones(m, dtype=float)*np.mean(res['raw_std_err'].values)
        else:
            # Rolling average sorted in order of fold change
            sorted_std_err = res['raw_std_err'].reindex(
                    np.argsort(res['fold_change'].values))
            p_sd = sorted_std_err.rolling(window=WINSIZE, center=True).mean()
            # Pad ends to avoid NaNs
            pad_size = (WINSIZE-1) / 2
            p_sd.iloc[:pad_size] = p_sd.iloc[pad_size]
            p_sd.iloc[-pad_size:] = p_sd.iloc[-pad_size-1]
            # Now make sure we are in protein id order again
            p_sd.sort_index(inplace=True)

        std_err = ((reg*(p_sd**2) + (res['n']-1)*(res['raw_std_err']**2)) / 
                   (reg+res['n']-1))**0.5
        # Two-tailed test
        p_val = 2*stats.t.cdf(-abs(res['fold_change']) / std_err, df=res['n']-1)
        res['std_err'] = std_err
        res['pval'] = p_val

        # Copy over degenerate pvalues if n=1
        res.loc[res['n'] == 1,'pval'] = res.loc[res['n'] == 1,'pval_no_reg']
        res.loc[res['n'] == 1,'std_err'] = res.loc[res['n'] == 1,'raw_std_err']
    else:
        res['std_err'] = res['raw_std_err']
        res['pval'] = res['pval_no_reg']

    return res


def _wls_degenerate(data, use_bayes):
    """ Weighted least squares with just one peptide

    Peptide mean directly estimates protein mean
    Valid for both one and two sample approach
    data should be df with one row
    """
    data = data.ix[0]
    nC, nE = data['nC'], data['nE']
    if use_bayes:
        stdC, stdE = data['bayesSDC'], data['bayesSDE']
    else:
        stdC, stdE = data['stdC'], data['stdE']
    return (data['meanE'] - data['meanC'],
            (((nC-1)*(stdC**2) + (nE-1)*(stdE**2)) /
                (nC+nE-2)) * (float(nC+nE) / (nC*nE)),
            data['pVal'])

def _wls(data, use_bayes=False):
    """ Weighted least squares for peptides in protein.
    DEPRECATED, do not use
    Args:
        data - df with columns including
            meanC, meanE, stdC, stdE, bayesSDC, bayesSDE
        use_bayes - bool, def False. If True, use bayesSDC/bayesSDE
        reg - (mean, n_obs) mean, pseudocount of of prior observation
            n_obs > 0
    Returns:
        (beta, std_err, pval, n)
    """
    # Degenerate case, only one peptide
    # Regularization doesn't affect answer here
    if data.shape[0] == 1:
        return _wls_degenerate(data, use_bayes)

    y = data['meanC'].append(data['meanE']).values
    if use_bayes:
        w = data['bayesSDC'].append(data['bayesSDE']).values**2
    else:
        w = data['stdC'].append(data['stdE']).values**2
    x = np.ones(data.shape[0]*2)
    x[:data.shape[0]] = 0

    mod_wls = sm.WLS(y, add_constant(x, prepend=False), weights=1./w)
    res_wls = mod_wls.fit()
    return (res_wls.params[0], res_wls.bse[0], res_wls.pvalues[0], data.shape[0])

def _wls_onesample(data, use_bayes=False, adj_var=True):
    """ Weighted least squares for one-sample peptides in protein
    Args:
        data - df with columns including
            meanC, meanE, stdC, stdE, bayesSDC, bayesSDE
        use_bayes - bool, def False. If True, use bayesSDC/bayesSDE
        reg - (std, n_obs) standard dev, pseudocount of of prior on std dev
            std, n_obs > 0
    Returns:
        (beta, std_err, pval)
    """
    n = data.shape[0]
    # Degenerate case, only one peptide
    # Regularization doesn't affect answer here
    if n == 1:
        return _wls_degenerate(data, use_bayes)

    y = data['meanE'].values - data['meanC'].values
    x = np.ones(data.shape[0])
    nC, nE = data['nC'].values.astype(float), data['nE'].values.astype(float)
    if use_bayes:
        stdC, stdE = data['bayesSDC'].values, data['bayesSDE'].values
    else:
        stdC, stdE = data['stdC'].values, data['stdE'].values
    # Variance is additive for control and experimental conditions
    w = 1. / (stdC**2 + stdE**2)
    # Weighted least squares with only constant term
    mod_wls = sm.WLS(y, x, weights=w)
    res_wls = mod_wls.fit()
    beta_hat, std_err, p_val = (res_wls.params[0], res_wls.bse[0], res_wls.pvalues[0])
    if adj_var:
        # Adjust the standard error for peptide level uncertainty
        # = reciprocal of sum of reciprocals of peptide variances
        std_err += 1. / sum(w)
        p_val = 2*stats.t.cdf(-abs(beta_hat) / std_err, df=(n-1))

    return (beta_hat, std_err, p_val) 


###############
#    MISC     #
###############

def protein_rollup(protein_df):
    """
    Runs cyberT on peptides and rolls up to protein rollup using R script

    Args:
        ctrl, exp should be dataframes of PEPTIDE level measurements
        protein is vector which designates the protein for each row of ctrl, exp

    Returns (WIP):

    """
    protein_df = r['find_protein_medians'](protein_df, use_isoform=True)
    out = pandas2ri.ri2py(protein_df)

    return out


def do_stat_tests(ctrl, exp, run_modT_2sample=False):
    """
    Runs modT, cyberT, t-test, fold_change analysis

    Returns modT_pvals,
            cyberT_pvals,
            ttest_pvals,
            ttest_ratio_pvals,
            fold_change (+),
            [modT_2sample]
    Any of these may be None if number of channels is not suitable
    (+) = Proper measure is inverted: i.e. x* = max(x) - x
    """

    do_ratio = (ctrl.shape[1] == 1 or ctrl.shape[1] == exp.shape[1]) and exp.shape[1] > 1
    do_t = (ctrl.shape[1] > 1)

    if ctrl.shape[1] == exp.shape[1]:
        modT_pvals = modT(ctrl, exp)['P.Value'].values
        ttest_ratio_pvals = t_test(ctrl, exp, ratio_test=True)
    else:
        min_dim = min(ctrl.shape[1], exp.shape[1])
        modT_pvals = modT(
                ctrl.iloc[:,:min_dim],
                exp.iloc[:,:min_dim]
        )['P.Value'].values
        ttest_ratio_pvals = t_test(
                ctrl.iloc[:,:min_dim], 
                exp.iloc[:,:min_dim], ratio_test=True)
    # TODO troublesome
    """
    if do_ratio:
        modT_res = modT(ctrl, exp)
        modT_pvals = modT_res['P.Value'].values
        print "Ran moderated T test"

        ttest_ratio_pvals = t_test(ctrl, exp, ratio_test=True)
        print "Ran one sample t test"
    else:
        modT_pvals = np.nan
        print "Skipped moderated T test, dimensions not suitable"

        ttest_ratio_pvals = np.nan
        print "Skipped one sample t test, dimensions not suitable"
    """

    if do_t:
        ttest_pvals = t_test(ctrl, exp)
        print "Ran two sample t test"
    else:
        ttest_pvals = np.nan
        print "Skipped two sample t test, too few channels"

    cyberT_res = cyberT(ctrl, exp)
    cyberT_pvals = cyberT_res['pVal'].values
    print "Ran cyberT test"

    fold_change = np.mean(ctrl.values, axis=1) - np.mean(exp.values, axis=1)

    res = pd.DataFrame({'modT_2samp': modT_pvals,
        'cyberT':cyberT_pvals,
        't_test_2samp': ttest_pvals,
        't_test_1samp': ttest_ratio_pvals,
        'fold_change': fold_change})
    if run_modT_2sample:
        modT_2sample_pvals = modT_2sample(ctrl, exp)['P.Value'].values
        res['modT_2samp'] = modT_2sample_pvals

    return res


def do_stat_tests_protein(ctrl, exp, protein):
    """
    Runs modT, cyberT, t-test, fold_change on protein level

    Returns:
        (modT_pvals (median)
         cyberT_pvals (median)
         fold_change (median)
         ), df

    pvals are all series, sorted in increasing order of protein id
    df has columns for pvals. Can pass directly into plot visualization
    """
    # TODO make sure all protein indices are sorted consistently
    # TODO handle nans gracefully
    if ctrl.shape[0] != exp.shape[0]:
        raise ValueError('Length of exp and ctrl data frames not identical')

    ctrl.columns = ['C%d' % (i+1) for i in xrange(ctrl.shape[1])]
    exp.columns = ['E%d' % (i+1) for i in xrange(exp.shape[1])]

    ## Do stat tests individually
    # (modT_pvals,
    #  cyberT_pvals,
    #  ttest_pvals,
    #  ttest_ratio_pvals,
    #  fold_change) = do_stat_tests(ctrl, exp)
    # Revert fold change inversion
    # fold_change = np.max(fold_change) + 0.01 - fold_change
    ctrl_mean = ctrl.groupby(protein).median()
    exp_mean = exp.groupby(protein).median()

    # TODO this is a little sketchy
    # Drop arbitrary column to make ctrl and exp numbers line up
    if ctrl_mean.shape[1] == exp_mean.shape[1]:
        modT_pvals = modT(ctrl_mean, exp_mean)['P.Value']
    else:
        min_dim = min(ctrl_mean.shape[1], exp_mean.shape[1])
        modT_pvals = modT(
                ctrl_mean.iloc[:,:min_dim],
                exp_mean.iloc[:,:min_dim]
        )['P.Value']

    cyberT_pvals = cyberT(ctrl_mean, exp_mean)['pVal']
    ttest_pvals = t_test(ctrl_mean, exp_mean)
    fold_change = (np.mean(np.ma.masked_invalid(ctrl_mean.values), axis=1) -
                  np.mean(np.ma.masked_invalid(exp_mean.values), axis=1))

    pval_df = pd.DataFrame({
        'protein_id': sorted(pd.unique(protein)),  # groupby sorts keys in asc
        'fold_change_med': fold_change,
        'modT_PVal_med': modT_pvals,
        'cyberT_PVal_med': cyberT_pvals,
        'ttest_PVal_med': ttest_pvals,
        })
    out = pval_df

    # Do naive aggregates of peptide level tests
    # Sorts in increasing order of id
    # out = pval_df.groupby('protein_id').median()

    # ctrl.reset_index(drop=True, inplace=True)
    # exp.reset_index(drop=True, inplace=True)
    # meanC, meanE = cyberT_res['meanC'].values, cyberT_res['meanE'].values
    meanC = np.mean(np.ma.masked_invalid(ctrl.values), axis=1)
    meanE = np.mean(np.ma.masked_invalid(exp.values), axis=1)

    # Note that all these are returned in ascending order of id
    # CyberT by peptide
    protein_cyberT_bypep = r['bayesT.pair'](
            pd.DataFrame.from_items([
                ('ratio', meanE - meanC),
                ('index', meanE + meanC),
            ]),
            1,
            aggregate_by=protein)
    protein_cyberT_bypep = pandas2ri.ri2py(protein_cyberT_bypep)
    protein_ttest_bypep = r['bayesT.pair'](
            pd.DataFrame.from_items([
                ('ratio', meanE - meanC),
                ('index', meanE + meanC),
            ]),
            1,
            bayes=False,
            aggregate_by=protein)
    protein_ttest_bypep = pandas2ri.ri2py(protein_ttest_bypep)

    cyberT_res = cyberT(ctrl, exp)
    cyberT_res['protein'] = protein
    # Testing different regularization coefficients
    # for reg in [2, 4, 8, 12, 16]:
    #     out['wls_pval_reg_%02d' % reg] = protein_wls_test_cyberT(
    #             cyberT_res, onesample=True, reg=reg)['pval'].values


    # out['wls'] = protein_wls_test_cyberT(
    #        cyberT_res, use_bayes=True, onesample=True)['pval'].values

    # out['wls_pval'] = wls_pval_reg['pval']
    out['cyberT_bypep'] = protein_cyberT_bypep['pVal'].values
    out['ttest_bypep'] = protein_ttest_bypep['pVal'].values
    out.reset_index(inplace=True)
    return out


def do_stat_tests_phospho(ctrl, exp, protein_labels, protein_df):
    """ Do statistical tests for phosphoproteins

    protein_df should have columns 'protein_id', 'fold_change', 'std_err'
        If 'std_err' not provided, will default to 0
        Optional: 'fold_change_nwls'
    """
    if 'std_err'in protein_df:
        std_err = protein_df['std_err'][protein_labels].values
    else:
        std_err = 0
        print "Warning! No standard errors were included in protein_df for phosphoprotein statistical tests, defaulting to zero"

    no_norm = cyberT(ctrl, exp)
    norm_exp = exp.subtract(
            protein_df['fold_change'][protein_labels].reset_index(drop=True),
            axis=0)
    mean_norm = cyberT(ctrl, norm_exp)
    norm_with_err = cyberT(ctrl, norm_exp, base_vars=std_err)

    res = pd.DataFrame({
        'no_adj': no_norm['pVal'],
        'mean_adj': mean_norm['pVal'],
        'var_adj': norm_with_err['pVal'],
        })

    if 'fold_change_nwls' in protein_df:
        # TODO what about non-wls mean???
        norm_non_wls_exp = exp.subtract(
                protein_df['fold_change_nwls'][protein_labels].reset_index(drop=True),
                axis=0)
        mean_norm_nwls = cyberT(ctrl, norm_non_wls_exp)
        res['norm_nwls_adj'] = mean_norm_nwls['pVal']

    return res


## TODO this is VERY JANKY
## Convenience function for protein pval labels
def protein_pval_labels():
    ctrl, exp, is_changed, protein = sample_proteins(500, 0, 2, 2)

    tmp = do_stat_tests_protein(ctrl, exp, protein)
    return tmp.columns

## Convenience function for peptide pval labels
def peptide_pval_labels(run_modT_2sample=True):
    ctrl, exp, is_changed = sample_no_ctrl_uniform(500, 0, 2)

    tmp = do_stat_tests(ctrl, exp, run_modT_2sample)
    return tmp.columns

def phospho_pval_labels(nowls=True):
    ctrl, exp, is_changed, protein = sample_proteins(500, 0, 2, 2)
    protein_df = pd.DataFrame({'protein_id': np.arange(500),
                               'fold_change': np.ones(500)})
    ctrl_p, exp_p, is_changed, mapping = sample_phospho(500, 0, 0, protein_df)

    tmp = do_stat_tests_phospho(ctrl_p, exp_p, mapping, protein_df)
    if nowls:
        return list(tmp.columns) + [u'mean_adj_nowls']
    else:
        return tmp.columns

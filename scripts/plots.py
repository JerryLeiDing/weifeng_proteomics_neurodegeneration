import matplotlib
import matplotlib.pyplot as plt

import numpy as np
import pandas as pd

from scipy.stats.stats import pearsonr

def plot_pvalue_dist(pvals, axes=None):
    """
    Pvals should be a Pandas dataframe
    """
    valid_labels = sorted(list(pvals.columns))
    m = len(valid_labels)
    if axes is None:
        f, axes = plt.subplots(1, m, sharey='row', squeeze=True)
    else:
        f = None

    # f, (hist_axs, log_axs) = plt.subplots(2, m, sharey='row', squeeze=False)

    axes[0].set_ylabel('Density')
    # log_axs[0].set_ylabel('Observed p-value')

    for i, name in enumerate(valid_labels):
        ax = axes[i]

        _, _, rects = ax.hist(pvals[name].dropna(),
                bins=20, range=(0,1), normed=True, alpha=0.5)
        ax.plot([0, 1], [1, 1], color='grey', lw=1, linestyle='--')
        # ax.set_title(name if name not in LABEL_MAPPING else LABEL_MAPPING[name])
        ax.set_title(name)
        ax.set_xlabel('p-value')

    return f


def volcano_plots(pval_df, fold_change, alpha=0.05, axes=None):
    valid_labels = sorted(list(pval_df.columns))
    m = len(valid_labels)
    if axes is None:
        f, axes = plt.subplots(1, m, sharey='row', squeeze=True)
    else:
        f = None

    axes[0].set_ylabel('$-\log_{10}$(p-value)')

    ALPHA = - np.log10(alpha)
    # Warn if fold change is nan
    if np.any(np.isnan(fold_change)):
        fold_change = fold_change.copy()
        print "[Volcano plot] WARNING: some fold changes were NaN"
        print "%d NaN were set to 0" % np.sum(np.isnan(fold_change))
        fold_change[np.isnan(fold_change)] = 0

    # Volcano plots
    for i, name in enumerate(valid_labels):
        ax = axes[i]
        pv = pval_df[name].values
        mask = np.isfinite(pv)
        thresh = (pv[mask] <= alpha)

        ax.scatter(fold_change[mask], -np.log10(pv[mask]),
                   c=thresh, alpha=0.15, cmap=matplotlib.cm.coolwarm, lw=0)
        (l, r) = ax.get_xlim()
        ax.plot([l, r], [ALPHA, ALPHA], color='grey', lw=1, linestyle='--')
        # ax.set_title(name if name not in LABEL_MAPPING else LABEL_MAPPING[name])
        ax.set_title(name + ('\n%d Significant' % np.sum(thresh)))
        ax.set_ylim(bottom=0)
        ax.set_xlim(l, r)
        ax.set_xlabel('$\log_2$(fold change)')

    return f


# Pvals to use for psm
LABEL_COLS = [
        'cyberT_PVal_med',
        'cyberT_bypep',
        'modT_PVal_med',
        'ttest_PVal_med',
        'ttest_bypep'
]
LABEL_COLS_PROT = ['cyberT', 'modT_2samp', 't_test_1samp', 't_test_2samp']
LABEL_COLS_ADJ = [l + '_adj' for l in LABEL_COLS]
LABEL_COLS_PROT_ADJ = [l + '_adj' for l in LABEL_COLS_PROT]
def compare_measures(res, label_cols, label_cols2=None, title='', 
                     count=True, equal=True, sharey=True, sharex=True, corr=False):
    """ Generate scatterplots comparing pvals with each other
    """
    # TODO: for convenience, should check if integer locations
    if label_cols2 is None:
        label_cols2 = label_cols
    m, n = (len(label_cols), len(label_cols2))
    f, axarr = plt.subplots(m, n, sharex=sharex, sharey=sharey)
    f.suptitle(title)

    if 'numPepsUnique' in res.columns:
        alpha = np.log2(res.numPepsUnique) / np.max(np.log2(res.numPepsUnique))
    elif 'n_pep' in res.columns:
        alpha = np.log2(res.n_pep) / np.max(np.log2(res.n_pep))
    else:
        alpha = np.full((res.shape[0],), 0.25)

    for i in xrange(m):
        for j in xrange(n):
            ax = axarr[i][j]
            if i == 0:
                ax.set_title(label_cols2[j])
            if j == 0:
                ax.set_ylabel(label_cols[i])
            if i > j:
                # ax.axis('off')
                continue
            if i == j:
                # TODO fix this
                # Plot histogram of data distribution instead of scatterplot
                # ax.hist(res[label_cols[i]])
                pass
            x = res[label_cols2[j]].values
            y = res[label_cols[i]].values
            mask = (np.isfinite(x) & np.isfinite(y))
            val = np.sum(mask)
            # Set alpha correctly
            colors = matplotlib.cm.Reds(np.arange(0, val, dtype=float) / val)
            # colors = np.zeros((n, 4))
            # colors[:,0] = 0.9  # Make red
            # colors[:,1] = 0.7 # Make red
            # colors[:,2] = 0.6 # Make red
            colors[:,-1] = alpha[mask]
            ax.scatter(x[mask],
                       y[mask],
                       color=colors, lw=0)
            if count:
                # Count number of significant
                nsig = np.sum((x[mask] <= 0.05) &
                              (y[mask] <= 0.05))
                ax.set_xlabel('%d sig' % nsig)
            if corr:
                # Find correlation
                r,_ = pearsonr(x[mask], y[mask])
                ax.set_xlabel(ax.xaxis.get_label_text() + (' r=%.3f' % r))
            if equal:
                # Equalize axis scales
                ax.axis('equal')

    f.tight_layout()
    return f



def generate_figures(res, label_cols=LABEL_COLS, filename='psm'):
    """
    Generate QA figures and related for results df
    """
    # TODO
    # Plot psm vs protein level p-vals
    # Plot distribution of peptides per protein in psm
    # Find correlation (if any) between number of peptides and p-value
    # Determine variance characteristics of data
        # Peptide level? Protein level???
        # This is awkward because there are only 1-2 samples per condition(!)
        # at peptide level, but at protein level it doesn't make sense
        # to directly compare all the peptide channels

    LABEL_COLS = label_cols
    LABEL_COLS_ADJ = [l + '_adj' for l in LABEL_COLS]
 
    # Plot raw p-values (histogram, volcano plot, fc distribution)
    f, axarr = plt.subplots(2, len(LABEL_COLS), figsize=(16, 10))
    plot_pvalue_dist(res[LABEL_COLS], axes=axarr[0])
    volcano_plots(res[LABEL_COLS], res.fold_change, axes=axarr[1])
    f.tight_layout(rect=(0, 0, 1, 0.95))
    f.suptitle('%s raw p-values' % filename)
    f.savefig('figures/%s_pvals_volcano_raw.png' % filename, dpi=100)

    # Plot adj p-values (histogram, volcano plot, fc distribution)
    f, axarr = plt.subplots(2, len(LABEL_COLS), figsize=(16, 10))
    plot_pvalue_dist(res[LABEL_COLS_ADJ], axes=axarr[0])
    volcano_plots(res[LABEL_COLS_ADJ], res.fold_change, axes=axarr[1])
    f.tight_layout(rect=(0, 0, 1, 0.95))
    f.suptitle('%s adjusted p-values' % filename)
    f.savefig('figures/%s_pvals_volcano_adj.png' % filename, dpi=100)

    # Histogram fc distribution
    f, ax = plt.subplots(figsize=(10, 10))
    ax.hist(res.fold_change.dropna(), bins=100)
    f.tight_layout()
    f.suptitle('%s log2 fold changes' % filename)
    f.savefig('figures/%s_fc_hist.png' % filename, dpi=100)

    # Compare pval measures
    f = compare_measures(res, LABEL_COLS, title='%s raw p-values' % filename)
    f.set_size_inches(10, 10)
    f.tight_layout(rect=(0, 0, 1, 0.95))
    f.savefig('figures/%s_compare_raw.png' % filename, dpi=100)

    # Compare pval measures adjusted
    f = compare_measures(res, LABEL_COLS_ADJ, title='%s adjusted p-values' % filename)
    f.set_size_inches(10, 10)
    f.tight_layout(rect=(0, 0, 1, 0.95))
    f.savefig('figures/%s_compare_adj.png' % filename, dpi=100)


def compare_channel_replicates(data, group=True, title='', col_groups=None, cross=False): 
    """ 
    Plot (ncols x ncols) grid of scatterplots comparing the measures in each column
    data must have 'accession_number' column to aggregate by if group is True
    col_groups is a list of name-value pairs to split the plot into several
        saved figures, where name is the label of the group and value is a 
        list of columns
        OR one of 'Mar', 'Aug', 'Sep', 'Tyr'
    cross [bool] - designates whether to include one panel with cross of all channels
        Note that cross will not with with group

    """
    # Default column groups for each dataset
    if col_groups is None:
        col_groups = [('all', np.xrange(data.shape[1]))]
    elif isinstance(col_groups, basestring):
        if col_groups == 'Aug':
            col_groups = [
            ('Control', ['GFP_A1', 'GFP_A2', 'GFP_B1', 'GFP_B2']),
            ('KO93',    ['KO93_A1','KO93_A2','KO93_B1','KO93_B2','KO93_B3']),
            ('KO95',    ['KO95_A1','KO93_A2','KO95_A3','KO95_B1','KO95_B2']),
            ('DKO',     ['DKO_A1', 'DKO_A2', 'DKO_B1', 'DKO_B2']),
            ]
        elif col_groups == 'Sep':
            col_groups = [
            ('P25_EE', ['P25_EE_A1','P25_EE_A2','P25_EE_A3']),
            ('EE', ['CT_EE_A1','CT_EE_A2','CT_EE_A3']),
            ('P25', ['P25_HC_A1','P25_HC_A2']),
            ('Control', ['CT_HC_A1','CT_HC_A2']),
            ]
    
    # Obtain all named columns
    all_cols = sum([g[1] for g in col_groups] ,[])
    if group:
        count = data.accession_number.groupby(data.accession_number).count()
        aggregated = data[all_cols].groupby(data.accession_number).mean()
        aggregated['n_pep'] = count
    else:
        aggregated = data

    for name, cols in col_groups:
        f = compare_measures(aggregated, cols, title=title,
                corr=True, count=False)
        f.set_size_inches(10, 10)
        f.tight_layout(rect=(0, 0, 1, 0.95))
        f.savefig('figures/%s_channel_reps_%s.png' % (title, name), dpi=100)
    
    # TODO separate aggregation code so cross_groups works with agg
    if cross:
        f = compare_measures(aggregated, all_cols, title="ALL",
                corr=True, count=False)
        f.set_size_inches(2*len(all_cols), 2*len(all_cols))
        f.tight_layout(rect=(0, 0, 1, 0.95))
        f.savefig('figures/%s_channel_reps_ALL.png' % (title,), dpi=100)


def intensity_psd_nonpsd(res, psd, ax=None):
    # Only plots average of GFP columns
    psd_accnums = set(psd.Accession)
    in_psd = np.array([1 if ac in psd_accnums else 0 
                       for ac in res.accession_number])

    res['in_psd'] = in_psd
    
    if ax is None:
        f, ax = plt.subplots()
    else:
        f = None

    # Plot non-psd elements
    gfp_cols = ['GFP_A1', 'GFP_A2', 'GFP_B1', 'GFP_B2']
    rng = (-8, 8)
    nb = 40
    a = 0.5

    ax.hist(np.mean(res[gfp_cols][res.in_psd == 1], axis=1),
            bins=nb, alpha=a, normed=True, range=rng, label='PSD proteins')
    ax.hist(np.mean(res[gfp_cols][res.in_psd == 0], axis=1),
            bins=nb, alpha=a, normed=True, range=rng, label='Other proteins')
    ax.legend()
    ax.set_title('Intensity of PSD proteins vs other proteins')
    ax.set_xlabel('$log_2$ Intensity (normalized to reference channel')
    ax.set_ylabel('Density')

    return f, ax


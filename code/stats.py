"""
MPACT
Copyright 2022, Robert M. Samples, Sara P. Puckett, and Marcy J. Balunas
Modifications Copyright 2026, Robert M. Samples
"""

import shutil
from pathlib import Path

import pandas as pd
import numpy as np
from scipy import stats
from filter import listfilter

def tstat(m1, sd1, n1, m2, sd2, n2):
    """
    Calculates the t-statistic for two independent samples with unequal variances.
    
    Args:
        m1 (float): Sample mean of the first sample.
        sd1 (float): Sample standard deviation of the first sample.
        n1 (int): Sample size of the first sample.
        m2 (float): Sample mean of the second sample.
        sd2 (float): Sample standard deviation of the second sample.
        n2 (int): Sample size of the second sample.
    
    Returns:
        float: The calculated t-statistic.
    """
    numerator = m1 - m2
    den1 = (sd1**2/n1).replace([np.inf, -np.inf, np.nan], 0)
    den2 = (sd2**2/n2).replace([np.inf, -np.inf, np.nan], 0)
    denominator = (den1 + den2)**0.5
    return(np.absolute(numerator/denominator))

def ws(sd1, n1, sd2, n2):
    """
    Calculates the effective sample size using the Welch-Satterthwaite equation.

    Args:
        sd1 (float): Sample standard deviation of the first sample.
        n1 (int): Sample size of the first sample.
        sd2 (float): Sample standard deviation of the second sample.
        n2 (int): Sample size of the second sample.

    Returns:
        float: The calculated effective sample size.
    """
    s1 = sd1/(n1**.5)
    s2 = sd2/(n2**.5)
    numerator = (s1**2/n1 + s2**2/n2)**2
    den1 = s1**4/((n1**2)*(n1-1))
    den2 = s2**4/((n2**2)*(n2-1))
    denominator = den1.replace([np.inf, -np.inf, np.nan], 0) + den2.replace([np.inf, -np.inf, np.nan], 0)
    return(numerator/denominator)

def groupave(analysis_params):
    """
    Averages groups and calculates technical and biological RSDs.

    Args:
        analysis_params (object): A custom object that contains various analysis parameters.

    Returns:
        None
    """
    print('Averaging groups')

    # Works on the wide table (features x injections) and averages across
    # columns. Stacking it into one row per value, as before, used several GB
    # for a large table (40,000 features x 400 injections is 17 million rows
    # with a six-level index) and got the process killed on small servers.
    msdata = pd.read_csv(
        analysis_params.outputdir / (analysis_params.filename.stem + '_formatted.csv'),
        sep=',', header=[0, 1, 2], index_col=[0, 1, 2])
    msdata.index.names = ['Compound', 'm/z', 'Retention time']
    msdata.columns.names = ['Group', 'Sample', 'Injection']
    values = msdata.T   # one row per injection column

    # Mean per injection (repeated columns for the same injection are averaged)
    injection_mean = values.groupby(level=[0, 1, 2], sort=True).mean()
    del values, msdata

    # Technical RSDs and counts: over the injections of each sample
    by_sample = injection_mean.groupby(level=[0, 1], sort=True)
    tech_mean = by_sample.mean()
    tech_stddev = by_sample.std().fillna(0)
    tech_rsd = (tech_stddev / tech_mean).fillna(0)
    techn = by_sample.size()   # number of injections per sample
    del injection_mean, tech_stddev

    # Biological RSDs and counts: over the samples in each group
    by_group = tech_mean.groupby(level=0, sort=True)
    biol_mean = by_group.mean()
    biol_stddev = by_group.std().fillna(0)
    biol_rsd = (biol_stddev / biol_mean).fillna(0)
    bioln = by_group.size()   # number of samples per group
    group_tech_rsd = tech_rsd.groupby(level=0, sort=True).mean()
    group_techn = techn.groupby(level=0, sort=True).mean()

    def long(frame):
        """Group x feature frame -> Series indexed (Compound, m/z, RT, Group)."""
        try:
            return frame.T.stack(future_stack=True)
        except TypeError:
            return frame.T.stack(dropna=False)

    average = long(biol_mean)
    groups = average.index.get_level_values('Group')
    msdata_errprop_mean = pd.DataFrame({
        'average': average,
        'biolRSD': long(biol_rsd),
        'bioln': bioln.reindex(groups).values,
        'techRSD': long(group_tech_rsd),
        'techn': group_techn.reindex(groups).values,
    }, index=average.index).sort_index()

    # Save results
    msdata_errprop_mean.to_csv(
        analysis_params.outputdir / (analysis_params.filename.stem + '_summarydata.csv'),
        header=True, index=True
    )
    msdata_errprop_grpav = msdata_errprop_mean[['average']]
    msdata_errprop_grpav.to_csv(
        analysis_params.outputdir / (analysis_params.filename.stem + '_groupaverages.csv'),
        header=True, index=True
    )

# samplenormalize warns about a sample whose signal is more than this many
# times above or below the median sample.
EXTREME_FACTOR = 5.0


def normalization_factors(msdata, method, skip_groups=()):
    """
    Per-sample normalization factors for a wide peak table.

    Args:
        msdata (DataFrame): features x injections, columns (group, sample, injection).
        method (str): 'sum' (total signal) or 'pqn' (probabilistic quotient,
            Dieterle et al. 2006, Anal. Chem. 78, 4281).
        skip_groups: groups left as they are (factor 1), such as blanks.

    Returns:
        DataFrame indexed (group, sample) with the number of injections, the
        sample's signal relative to the median sample, and the factor that
        divides it out (1 / relative signal).
    """
    # One profile per sample: the mean of its technical replicate injections,
    # so every injection of a sample gets the same factor and technical CVs
    # are unchanged.
    profiles = msdata.T.groupby(level=[0, 1], sort=True).mean().T
    injections = msdata.T.groupby(level=[0, 1], sort=True).size()
    used = ~profiles.columns.get_level_values(0).isin(list(skip_groups))
    relative = pd.Series(np.nan, index=profiles.columns)
    if used.any():
        sample_profiles = profiles.loc[:, used]
        if method == 'sum':
            totals = sample_profiles.sum(axis=0)
            relative[used] = totals / totals[totals > 0].median()
        elif method == 'pqn':
            reference = sample_profiles.median(axis=1)
            ratios = sample_profiles.div(reference, axis=0)
            ratios = ratios.where((sample_profiles > 0) & (reference > 0).values[:, None])
            quotients = ratios.median(axis=0)
            relative[used] = quotients / quotients[quotients > 0].median()
        else:
            raise ValueError('unknown normalization method: ' + str(method))
    # A sample with no usable signal is left as it is.
    relative = relative.where(relative > 0)
    factors = pd.DataFrame({'injections': injections,
                            'relative signal': relative,
                            'factor': (1 / relative).fillna(1.0)})
    factors.index.names = ['Group', 'Sample']
    return factors


def _scale_table(path, factors):
    """Multiply a peak table on disk (three header rows, then Compound, m/z,
    RT and one column per injection) by its sample's factor, in place."""
    table = pd.read_csv(path, sep=',', header=None, index_col=None, low_memory=False)
    header, body = table.iloc[:3], table.iloc[3:]
    lookup = factors['factor'].to_dict()
    column_factors = [lookup.get((str(header.iat[0, i]), str(header.iat[1, i])), 1.0)
                      for i in range(3, table.shape[1])]
    values = body.iloc[:, 3:].astype(float) * np.array(column_factors)
    body = pd.concat([body.iloc[:, :3], values], axis=1)
    pd.concat([header, body], ignore_index=True).to_csv(path, index=False, header=False)


def samplenormalize(analysis_params):
    """
    Beta: normalize each biological sample to the same overall signal.

    The factor is worked out from the features that passed every filter, so
    noise, solvent and blank peaks do not set it, and it is applied to every
    injection of the sample alike, so technical CVs are unchanged. The
    formatted and filtered tables are rescaled in place; the unscaled filtered
    table is kept as <name>_filtered_raw.csv and the factors are written to
    <name>_normfactors.csv. The blank group, and any group named in
    samplenorm_skip (a media blank, say), is left unscaled and is not used to
    work out the factors of the others.

    Args:
        analysis_params (object): needs samplenorm ('sum' or 'pqn'), outputdir,
            filename and blnkgrp; samplenorm_skip is optional.

    Returns:
        DataFrame: the factors (see normalization_factors).
    """
    method = getattr(analysis_params, 'samplenorm', 'none')
    print(f'Normalizing samples (beta: {method})')
    stem = analysis_params.outputdir / analysis_params.filename.stem
    filtered_path = Path(str(stem) + '_filtered.csv')
    formatted_path = Path(str(stem) + '_formatted.csv')

    filtered = pd.read_csv(filtered_path, sep=',', header=[0, 1, 2], index_col=[0, 1, 2])
    blank = getattr(analysis_params, 'blnkgrp', '') or ''
    skip = [str(g) for g in getattr(analysis_params, 'samplenorm_skip', None) or []]
    if blank and blank not in skip:
        skip.insert(0, blank)
    if skip:
        print('  left unscaled: ' + ', '.join(skip))
    factors = normalization_factors(filtered, method, skip_groups=skip)
    factors.to_csv(Path(str(stem) + '_normfactors.csv'), header=True, index=True)
    names = factors.index.get_level_values(1)
    unused = factors['relative signal'].isna() & ~factors.index.get_level_values(0).isin(skip)
    if unused.any():
        print('  warning: no usable signal in ' + ', '.join(names[unused]) + '; left unscaled')
    # A factor this far from 1 is more likely a group with different contents
    # (a media blank, a strain that makes more) than a dilution.
    relative = factors['relative signal']
    extreme = (relative < 1 / EXTREME_FACTOR) | (relative > EXTREME_FACTOR)
    if extreme.any():
        listed = ', '.join(f'{n} ({r:.3g}x)' for n, r in zip(names[extreme], relative[extreme]))
        print(f'  warning: signal more than {EXTREME_FACTOR:g}x from the median sample in '
              f'{listed}; check these are dilutions, or leave their groups out')

    shutil.copyfile(filtered_path, Path(str(stem) + '_filtered_raw.csv'))
    for path in (formatted_path, filtered_path):
        _scale_table(path, factors)
    return factors


def properr(analysis_params):
    """
    Propagates error and calculates the combined relative standard deviation,
    absolute combined standard deviation, and effective sample size.

    Args:
        analysis_params (object): A custom object that contains various analysis parameters.

    Returns:
        None
    """
    print('Propagating error')
    msdata_errprop = pd.read_csv(
        analysis_params.outputdir / (analysis_params.filename.stem + '_summarydata.csv'),
        sep=',', header=0, index_col=[0, 1, 2, 3]
    )

    # Calculate combined RSD and absolute combined SD
    msdata_errprop['combRSD'] = np.sqrt(msdata_errprop['techRSD'] ** 2 + msdata_errprop['biolRSD'] ** 2)
    msdata_errprop['combASD'] = msdata_errprop['combRSD'] * msdata_errprop['average']

    # Set NaN values to 0
    msdata_errprop[['combASD', 'biolRSD', 'techRSD']] = msdata_errprop[['combASD', 'biolRSD', 'techRSD']].fillna(0)

    # Calculate effective sample size and replace infinite and NaN values with 0
    msdata_errprop['neff'] = (ws(
        msdata_errprop['techRSD'], msdata_errprop['techn'],
        msdata_errprop['biolRSD'], msdata_errprop['bioln']
    ) + 1).replace([np.inf, -np.inf, np.nan], 0)

    msdata_errprop = msdata_errprop.unstack()
    msdata_errprop.to_csv(
        analysis_params.outputdir / (analysis_params.filename.stem + '_summarydata.csv'),
        sep=',', index=True, header=True
    )


    
def runfc(analysis_params, statstgrps):
    """
    Runs fold change analysis and saves formatted backup for later use.

    Args:
        analysis_params (object): A custom object that contains various analysis parameters.
        statstgrps (list): A list of length two containing the names of the two sample groups to compare.

    Returns:
        None
    """
    iondict = pd.read_csv(analysis_params.outputdir / 'iondict.csv', sep = ',', header = [0], index_col = [0])
    maxval=100
    minval=.01
    msdata_errprop = pd.read_csv(analysis_params.outputdir / (analysis_params.filename.stem + '_groupaverages.csv'), sep = ',', header = [0], index_col = [0, 1, 2, 3]).unstack()
    msdata_errprop.columns = msdata_errprop.columns.droplevel()
    maxvals = msdata_errprop.max(axis=1) #the max stuff is used to base point opacity based on abundance, have to put in a better location than this
    msdata_errprop = msdata_errprop.loc[:, msdata_errprop.columns.intersection(statstgrps)]
    msdata_errprop['FC']=msdata_errprop[statstgrps[0]]/msdata_errprop[statstgrps[1]]
    # Use .loc (single-step assignment) instead of chained assignment, which is
    # deprecated and will stop updating the frame in pandas 3.0.
    msdata_errprop.loc[msdata_errprop['FC'] >= maxval, 'FC'] = maxval
    msdata_errprop.loc[msdata_errprop['FC'] <= minval, 'FC'] = minval
    msdata_errprop['max'] = maxvals # this has to be done before resetting the index or else tuple/string expectation error thrown
    msdata_errprop = msdata_errprop.reset_index([1,2])
    iondict['fc'] = msdata_errprop['FC']
    iondict['max'] = msdata_errprop['max']
    iondict['logmax'] = np.log10(iondict['max']) 
    #iondict['logfc'] = np.log2(iondict['fc']) maybe include this later, can replace relevant line in fc3d
    iondict.to_csv(analysis_params.outputdir / 'iondict.csv', header = True, index = True) #saves formatted backup for later use
    
    
def runttest(analysis_params, statstgrps, groupsets):
    """
    Runs a two-sample t-test for each ion and calculates the q-value.

    Args:
        analysis_params (object): A custom object that contains various analysis parameters.
        statstgrps (list): A list of length two containing the names of the two sample groups to compare.
        groupsets (dict): A dictionary containing ion sets for each comparison group.

    Returns:
        None
    """
    # Load iondict and data
    iondict = pd.read_csv(analysis_params.outputdir / 'iondict.csv', sep=',', header=[0], index_col=[0])
    msdata_errprop = pd.read_csv(analysis_params.outputdir / (analysis_params.filename.stem + '_summarydata.csv'), sep=',', header=[0, 1], index_col=[0, 1, 2])

    # Set minimum p-value
    minval = 0.00001
        
    # Calculate t-test statistics
    msdata_teststats = msdata_errprop
    msdata_teststats['T'] = tstat(msdata_teststats[('average', statstgrps[0])], msdata_teststats[('combASD', statstgrps[0])], msdata_teststats[('neff', statstgrps[0])], msdata_teststats[('average', statstgrps[1])], msdata_teststats[('combASD', statstgrps[1])], msdata_teststats[('neff', statstgrps[1])])
    msdata_teststats['deg'] = msdata_teststats[('neff', statstgrps[0])] + msdata_teststats[('neff', statstgrps[1])] - 2
    msdata_teststats['p'] = (1 - stats.t.cdf(msdata_teststats['T'], msdata_teststats['deg'])) * 2
    msdata_teststats.loc[msdata_teststats['p'] <= minval, 'p'] = minval
    msdata_teststats['logp'] = np.log10(msdata_teststats['p']) 
    
    # Save msdata_teststats. Previously written to the current working
    # directory as 'msdata_teststats_test.csv' (a debug-named file that
    # littered code/ and was never read back); now written into the run's
    # output directory under a descriptive name alongside the other outputs.
    msdata_teststats = msdata_teststats.reset_index([1, 2])
    msdata_teststats.to_csv(
        analysis_params.outputdir / (analysis_params.filename.stem + '_teststats.csv'),
        header=True, index=True)
    
    # Update iondict with -logp
    iondict['-logp'] = -msdata_teststats['logp']
    iondict.to_csv(analysis_params.outputdir / 'iondict.csv', header=True, index=True)
    
    # Load iondict and msdata
    iondict = pd.read_csv(analysis_params.outputdir / 'iondict.csv', sep=',', header=[0], index_col=[0])
    msdata = pd.read_csv(analysis_params.outputdir / (analysis_params.filename.stem + '_filtered.csv'), sep=',', header=[2], index_col=0).index.tolist()
    
    # Filter iondict based on msdata
    iondict = iondict.loc[msdata].reset_index()
    
    # Create list of ions to plot and test
    plottedions = []
    for elem in analysis_params.querylist:
        plottedions = list(set(plottedions + groupsets[elem].ionlist))
    
    # Filter iondict based on plottedions
    iondict = listfilter(iondict, plottedions, True)
    
    # Compute p-values for each ion
    iondict['p'] = 10 ** (-iondict['-logp'])
    
    # Sort ions by p-value
    iondict = iondict.sort_values(by=['p'])
    
    # Compute q-values for each ion
    num_ions = len(iondict[iondict['p'] <= 1])
    iondict['qval'] = iondict.reset_index().index + 1
    iondict['qval'] = iondict['p'] * num_ions / iondict['qval']
    iondict = iondict.sort_values(by=['p'], ascending=False)
    
    # Compute minimum q-value
    min_qval = 1
    qvals = iondict['qval'].tolist()
    for pos in range(len(qvals)):
        if qvals[pos] < min_qval:
            min_qval = qvals[pos]
        else:
            qvals[pos] = min_qval
    
    # Compute -logq for each ion
    iondict['qval'] = qvals
    iondict['-logq'] = -np.log10(iondict['qval'])
    
    # Save results to CSV files
    iondict = iondict.set_index('Compound')
    # Written into the run's output directory rather than the current working
    # directory (was a bare 'qdata.csv' that landed in code/ and was never
    # read back); the canonical -logq still goes into iondict.csv below.
    iondict.to_csv(
        analysis_params.outputdir / (analysis_params.filename.stem + '_qvalues.csv'),
        index=True, header=True)
    iondict2 = pd.read_csv(analysis_params.outputdir / 'iondict.csv', sep=',', header=[0], index_col=[0])
    iondict2['-logq'] = np.nan
    iondict2.loc[iondict.index.tolist(), '-logq'] = iondict['-logq']
    iondict2.to_csv(analysis_params.outputdir / 'iondict.csv', header=True, index=True)
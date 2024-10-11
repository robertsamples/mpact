"""
MPACT
Copyright 2022, Robert M. Samples, Sara P. Puckett, and Marcy J. Balunas
"""

import pandas as pd
import numpy as np
from scipy import stats
from filter import listfilter
from tqdm import tqdm

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

    # Calculate total number of rows and set chunk size
    total_rows = sum(1 for _ in open(analysis_params.outputdir / (analysis_params.filename.stem + '_formatted.csv'))) - 1
    header = pd.read_csv(analysis_params.outputdir / (analysis_params.filename.stem + '_formatted.csv'), nrows=0, header=[0, 1, 2])
    total_columns = header.shape[1]
    cells_per_chunk = 100_000
    chunk_size = max(1, cells_per_chunk // total_columns)  # Rows per chunk

    # Initialize lists to collect results from each chunk
    sum_values_list = []
    sum_squares_list = []
    counts_list = []

    # Process data in chunks
    with tqdm(total=(total_rows // chunk_size) + 1, desc='Processing chunks') as pbar:
        for chunk in pd.read_csv(
            analysis_params.outputdir / (analysis_params.filename.stem + '_formatted.csv'),
            sep=',', header=[0, 1, 2], index_col=[0, 1, 2], chunksize=chunk_size
        ):
            # Stack the DataFrame to move column levels into the index
            try:
                chunk_stacked = chunk.stack(level=[0, 1, 2], future_stack=True)
            except TypeError:
                chunk_stacked = chunk.stack(level=[0, 1, 2])

            # Set index names according to your data structure
            chunk_stacked.index.names = ['Compound', 'm/z', 'Retention time', 'Group', 'Sample', 'Injection']

            # Compute sum, sum of squares, and counts per group
            group_levels = ['Compound', 'm/z', 'Retention time', 'Group', 'Sample', 'Injection']
            sum_values_chunk = chunk_stacked.groupby(level=group_levels).sum()
            sum_squares_chunk = (chunk_stacked ** 2).groupby(level=group_levels).sum()
            count_chunk = chunk_stacked.groupby(level=group_levels).count()

            # Append results to lists
            sum_values_list.append(sum_values_chunk)
            sum_squares_list.append(sum_squares_chunk)
            counts_list.append(count_chunk)

            pbar.update(1)

    # Concatenate all results
    all_sum_values = pd.concat(sum_values_list)
    all_sum_squares = pd.concat(sum_squares_list)
    all_counts = pd.concat(counts_list)

    # Aggregate over the entire dataset
    sum_values_df = all_sum_values.groupby(level=group_levels).sum()
    sum_squares_df = all_sum_squares.groupby(level=group_levels).sum()
    counts_df = all_counts.groupby(level=group_levels).sum()

    # Calculate mean and variance per injection
    mean_values = sum_values_df / counts_df
    variance_values = (sum_squares_df / counts_df) - (mean_values ** 2)
    stddev_values = variance_values ** 0.5

    # Calculate technical RSDs and counts
    # Group over technical replicates within each sample
    tech_group_levels = ['Compound', 'm/z', 'Retention time', 'Group', 'Sample']
    tech_mean = mean_values.groupby(level=tech_group_levels).mean()
    tech_stddev = mean_values.groupby(level=tech_group_levels).std().fillna(0)
    tech_rsd = (tech_stddev / tech_mean).fillna(0)
    techn = counts_df.groupby(level=tech_group_levels).size()  # Number of injections per sample

    # Calculate biological RSDs and counts
    # Group over samples within each group
    biol_group_levels = ['Compound', 'm/z', 'Retention time', 'Group']
    biol_mean = tech_mean.groupby(level=biol_group_levels).mean()
    biol_stddev = tech_mean.groupby(level=biol_group_levels).std().fillna(0)
    biol_rsd = (biol_stddev / biol_mean).fillna(0)
    bioln = tech_mean.groupby(level=biol_group_levels).size()  # Number of samples per group

    # Prepare final DataFrame
    msdata_errprop_mean = pd.DataFrame({
        'average': biol_mean,
        'biolRSD': biol_rsd,
        'bioln': bioln,
        'techRSD': tech_rsd.groupby(level=biol_group_levels).mean(),
        'techn': techn.groupby(level=biol_group_levels).mean()
    })

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
    msdata_errprop['FC'][msdata_errprop['FC'] >= maxval] = maxval
    msdata_errprop['FC'][msdata_errprop['FC'] <= minval] = minval
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
    msdata_teststats['p'][msdata_teststats['p'] <= minval] = minval
    msdata_teststats['logp'] = np.log10(msdata_teststats['p']) 
    
    # Save msdata_teststats
    msdata_teststats = msdata_teststats.reset_index([1, 2])
    msdata_teststats.to_csv('msdata_teststats_test.csv', header=True, index=True)
    
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
    iondict.to_csv('qdata.csv', index=True, header=True)
    iondict2 = pd.read_csv(analysis_params.outputdir / 'iondict.csv', sep=',', header=[0], index_col=[0])
    iondict2['-logq'] = np.nan
    iondict2.loc[iondict.index.tolist(), '-logq'] = iondict['-logq']
    iondict2.to_csv(analysis_params.outputdir / 'iondict.csv', header=True, index=True)
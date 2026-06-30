"""
MPACT
Copyright 2022, Robert M. Samples, Sara P. Puckett, and Marcy J. Balunas
"""

import pandas as pd
import numpy as np
import filter
import stats
from groupsets import normalize_graphfilters
from datetime import datetime
import time

#---Classes---


class groupset:  # parses the ions to be plotted as a specific color based on thier presence in user input sample groups
    """
    Parses the ions to be plotted as a specific color based on their presence in user input sample groups.

    Args:
    - name (str): The name of the group set.
    - query (query_parameters): A query_parameters object.
    - iondict (DataFrame): The ion dictionary, indexed by ``Compound``, with
      a ``groups`` column. Callers already have this loaded in memory (it's
      built immediately before every groupset is constructed) -- pass that
      directly rather than a path, so each groupset doesn't re-read the same
      file from disk (previously: once per Plot Feature Set configured).

    Attributes:
    - legendname (str): The legend name of the group set.
    - excl (list of str): The exclusion groups in the group set.
    - incl (list of str): The inclusion groups in the group set.
    - plotcol (str): The plot color of the group set.
    - ionlist (list of str): The ions to be plotted in the group set.
    """
    def __init__(self, name, query, iondict):
        self.legendname = query.name
        self.excl = query.excl
        self.incl = query.incl
        self.plotcol = query.colour

        iondict = iondict.loc[iondict['groups'].notna(), 'groups'].to_frame()
        exclgrps = ' ' + '| '.join(self.excl)
        if len(exclgrps)>3: #only runs below line if there is an exclusion group
            iondict = iondict.loc[~iondict['groups'].str.contains(exclgrps), 'groups'].to_frame()#this does not work right if there are no exclusion groups so an artificial one was added
        for group in self.incl:
            iondict = iondict.loc[iondict['groups'].str.contains(' ' + str(group)), 'groups'].to_frame()
        self.ionlist = iondict.index.to_list()
       
class analysis_parameters:
    """
    An empty class that is used as a placeholder to hold analysis parameters.

    Attributes:
    - init (str): An initial string.
    """
    def __init__(self):
        self.init = ''


class AnalysisResult:
    """The output of ``run_MSFaST()``: everything it computes, bundled into
    one plain object instead of written onto whichever object was passed in.

    Attributes:
    - ionfilters (dict): Per-filter-type ``ionfilter`` results (cv/relfil/insource).
    - groupionlists (dict): Per-biological-group lists of failing ion IDs.
    - groupsets (dict): Per-Plot-Feature-Set ``groupset`` objects.
    - filtereddfs (dict): Per-Plot-Feature-Set filtered peak tables.
    """
    def __init__(self, ionfilters, groupionlists, groupsets, filtereddfs):
        self.ionfilters = ionfilters
        self.groupionlists = groupionlists
        self.groupsets = groupsets
        self.filtereddfs = filtereddfs


#---Methods---

def start_time(): 
    """Function to start calculating runtime"""
    global initial
    initial = time.time()
    return initial

def stop_time(): 
    """Function to stop calculating runtime and return the elapsed time"""
    final = time.time()
    return(final - initial)

def importdata():
    """Function to import data from files and format it"""
    print('Loading files')

    # Import sample/extract metadata
    extractmetadata = pd.read_csv(analysis_params.extractmetadatafilename,
                                   sep=',', header=0, index_col=None)

    # Import instrument sample list
    samplelist = pd.read_csv(analysis_params.samplelistfilename,
                              sep=',', header=0, index_col=None)

    # Join extract metadata and sample list by the sample code
    combinedmetadata = (extractmetadata.set_index('Sample_Code')
                        .join(samplelist.set_index('Sample_Code'))
                        .reset_index().set_index('Injection'))


    # Import feature list
    msdata = pd.read_csv(analysis_params.filename, sep=',', header=None,
                          index_col=[0, 1, 2], low_memory=False)

    # Iterate over header to format
    for position, elem in enumerate(msdata.iloc[1]):
        msdata.iloc[1, position] = combinedmetadata.loc[
            msdata.iloc[2, position], 'Sample_Code']
        msdata.iloc[0, position] = combinedmetadata.loc[
            msdata.iloc[2, position], 'Biological_Group']

    # Write formatted peak table
    print('Writing formatted peak table')
    msdata.to_csv(analysis_params.outputdir / (analysis_params.filename.stem + '_formatted.csv'),
                   header=False, index=True)

    # Import data from formatted peak table. This re-read is a genuine
    # reshape, not a redundant one: writing with header=False/index=True
    # then reading back with header=[0, 1, 2]/index_col=[0] promotes the
    # first 3 (already-edited) data rows into a 3-level column header and
    # collapses the 3-level row index down to one level -- replicating that
    # in memory would mean re-deriving exactly what pd.read_csv's header/
    # index parsing does (dtype inference, NaN handling, etc.) by hand, the
    # same category of subtle mismatch that caused the LossySetitemError
    # fixed earlier. Not worth the risk for one disk round-trip.
    msdata = pd.read_csv(analysis_params.outputdir /
                          (analysis_params.filename.stem + '_formatted.csv'),
                          sep=',', header=[0, 1, 2], index_col=[0])

    # Calculate mean of each row and drop rows with mean of 0
    msdataave = msdata.iloc[:, 2:].astype(float)
    msdataave['mean'] = msdataave.mean(axis=1, skipna=True)
    msdata = msdata[msdataave['mean'] != 0]

    # Save formatted peak table and original feature list
    msdata.to_csv(analysis_params.outputdir / (analysis_params.filename.stem + '_formatted.csv'),
                   header=True, index=True)
    msdata.to_csv(analysis_params.outputdir / (analysis_params.filename.name),
                   header=True, index=True)

    # Build the feature dictionary (with KMD) from msdata already in memory
    # instead of re-reading the file just written above with yet another
    # header shape (header=[2], flattening the 3-level column header to its
    # bottom level only). Equivalent to msdata.columns.get_level_values(2)
    # plus resetting the index back to a plain Compound column -- verified
    # against the previous disk-read version with real example data
    # (identical values and dtypes) before relying on it here.
    iondict = msdata.copy()
    iondict.columns = iondict.columns.get_level_values(2)
    iondict = iondict.reset_index().iloc[:, :3]
    iondict.columns = ['Compound', 'm/z', 'Retention time (min)']
    iondict['kmd'] = iondict['m/z'] - np.floor(iondict['m/z'])
    iondict.to_csv(analysis_params.outputdir / 'iondict.csv', header=True, index=False)


def run_MSFaST(params):
    """
    Runs the MSFaST analysis pipeline for the given analysis parameters.

    Args:
    params (analysis_parameters): The analysis parameters/inputs to run.

    Returns:
    AnalysisResult: everything the pipeline computed (ionfilters,
    groupionlists, groupsets, filtereddfs) -- nothing is written onto
    ``params`` or any other caller-supplied object.
    """
    start_time()

    # importdata()/filter.py/stats.py read and mutate this same object via
    # the `global analysis_params` they declare -- unchanged from before,
    # just no longer aliased through a caller-supplied `self.analysis_paramsgui`.
    global analysis_params
    analysis_params = params
    importdata()

    # Filtering and error propagation
    print('Filtering data')
    ionfilters = {}
    # Initialise here (not only inside `if analysis_params.grpave:`) so the
    # unconditional groupionlists[...] writes further down (and the blank
    # filter, which reads it) can't raise NameError if grpave is ever off.
    # The GUI currently forces grpave=True, but loaded sessions/tests need not.
    groupionlists = {}
    if analysis_params.relfil:
        ionfilters = filter.relationalfilter(analysis_params, ionfilters)
        if analysis_params.merge:
            filter.mergeions(analysis_params, ionfilters)
    if analysis_params.grpave:
        stats.groupave(analysis_params)
        print('Parsing ion lists')
        groupionlists = filter.parsionlists(analysis_params)
    if analysis_params.blnkfltr:
        msdata = pd.read_csv(analysis_params.outputdir / (analysis_params.filename.stem + '_formatted.csv'), sep=',', header=[0, 1, 2], index_col=None)
        msdata = filter.listfilter(msdata, groupionlists[analysis_params.blnkgrp], False)
        msdata = msdata.drop(analysis_params.blnkgrp, axis=1, level=0)
        msdata.to_csv(analysis_params.outputdir / (analysis_params.filename.stem + '_formatted.csv'), header=True, index=False)
        iondict = pd.read_csv(analysis_params.outputdir / 'iondict.csv', sep=',', header=[0], index_col=0)
        iondict['pass_blnkfil'] = ~iondict.index.isin(groupionlists[analysis_params.blnkgrp])
        iondict.to_csv(analysis_params.outputdir / 'iondict.csv', header=True, index=True)
    if analysis_params.CVfil:
        ionfilters = filter.cvfilter(analysis_params, ionfilters, analysis_params.cvthresh)
    if analysis_params.decon:
        ionfilters = filter.decon(analysis_params, ionfilters)
    filter.applyfilters(analysis_params, ionfilters)
    if analysis_params.prperr:
        stats.properr(analysis_params)

    # Parse ion lists and add filter lists
    groupionlists['cv'] = ionfilters['cv'].ions if analysis_params.CVfil else []
    groupionlists['relfil'] = ionfilters['relfil'].ions if analysis_params.relfil else []
    groupionlists['insource'] = ionfilters['insource'].ions if analysis_params.decon else []

    # Add groups column to iondict csv
    iondict = pd.read_csv(analysis_params.outputdir / 'iondict.csv', sep=',', header=0, index_col=None)
    iondict['groups'] = ''
    for group in groupionlists:
        iondict.loc[iondict['Compound'].isin(groupionlists[group]), 'groups'] += (' ' + str(group))
    iondict.to_csv(analysis_params.outputdir / 'iondict.csv', header=True, index=None)

    # Add default filters to querylist. Reuse the iondict already in memory
    # (just written above, unchanged since) instead of re-reading it from
    # disk -- and build the Compound-indexed view groupset needs once here,
    # rather than each groupset() call below re-reading the file from disk
    # itself (previously: once per Plot Feature Set configured).
    iondict_by_compound = iondict.set_index('Compound')
    groupsets, filtereddfs = {}, {}
    for elem in analysis_params.querylist:
        groupsets[elem] = groupset(elem, analysis_params.querydict[elem], iondict_by_compound)
        filtereddfs[elem] = filter.listfilter(iondict, groupsets[elem].ionlist, True)


    #block creates user specified plots some of these can be changed to eliminate a few arguments with data pulled from analysis_params
    if analysis_params.FC:
        stats.runfc(analysis_params, analysis_params.statstgrps)
    if analysis_params.Ttest:
        stats.runttest(analysis_params, analysis_params.statstgrps, groupsets)

    #---Analysis info file writing---
    runtime = stop_time()
    msdata_filtered = pd.read_csv(analysis_params.outputdir / (analysis_params.filename.stem + '_filtered.csv'), sep = ',', header = [0, 1, 2], index_col = [0])
    msdata_header = pd.read_csv(analysis_params.outputdir / (analysis_params.filename.stem + '_filtered.csv'), sep = ',', header = None, index_col = [0,1,2]).iloc[:3,:].transpose()
    msdata_header.columns = ['Biolgroup', 'Sample', 'Injection']

    iondict = pd.read_csv(analysis_params.outputdir / 'iondict.csv', sep = ',', header = [0], index_col = [0])

    msdata_unformatted = pd.read_csv(analysis_params.filename, sep = ',', header = [0, 1, 2], index_col = [0, 1, 2]) #imports feature list
    msdata = pd.read_csv(analysis_params.outputdir / (analysis_params.filename.name), sep = ',', header = [0, 1, 2], index_col = [0, 1, 2])
    msdata_filtered = pd.read_csv(analysis_params.outputdir / (analysis_params.filename.stem + '_filtered.csv'), sep = ',', header = [0, 1, 2], index_col = [0, 1, 2])
    analysisrec = open(analysis_params.outputdir / 'analysisinfo.txt',"w")
    analysisrec.writelines(['Analysis Date: ' + str(datetime.now()) + '\n',
                            'Runtime: ' + str(round(runtime, 2)) + ' seconds\n',
                            'Input file: ' + str(analysis_params.filename) + '\n',
                            'Sample list: ' + str(analysis_params.samplelistfilename) + '\n',
                            'Extract metadata file: ' + str(analysis_params.extractmetadatafilename) + '\n',
                            'Number of injections: ' + str(msdata_header['Injection'].nunique()) + '\n',
                            'Number of samples: ' + str(msdata_header['Sample'].nunique())  + '\n',
                            'Number of biological groups: ' + str(msdata_header['Biolgroup'].nunique())  + '\n',
                            '\n',
                            '\n',
                            '---Filtering---\n',
                            '-CV filtering-\n',
                            'CV threshold: ' + str(analysis_params.cvthresh) + '\n',
                            '\n',
                            '-Relative filtering-\n',
                            'Ringing window: ' + str(analysis_params.ringingwin) + '\n',
                            'Isotopic peak window: ' + str(analysis_params.isopeakwin) + '\n',
                            'Dimer filter window: ' + str(analysis_params.dimerpeakwin) + '\n',
                            'RT window: ' + str(analysis_params.RTwin) + '\n',
                            'max isotopic peak shift: ' + str(analysis_params.maxisowin) + '\n',
                            '\n',
                            'Total features: ' + str(len(msdata.index)) + '\n'])


    text = ''
    if analysis_params.relfil:
        text += 'Features failing peak correction filtering: ' + str(len(ionfilters['relfil'].ions)) + '/' + str(len(msdata_unformatted.index)) + ' ' + str(round(100 * len(ionfilters['relfil'].ions) / len(msdata_unformatted.index), 2)) + '%\n'
    if analysis_params.blnkfltr:
        text += 'Features failing blank filtering: ' + str(len(groupionlists[analysis_params.blnkgrp])) + '/' + str(len(msdata_unformatted.index)) + ' ' + str(round(100 * len(groupionlists[analysis_params.blnkgrp]) / len(msdata_unformatted.index), 2)) + '%\n'
    if analysis_params.decon:
        text += 'Features failing in-source/deconvolution filtering: ' + str(len(ionfilters['insource'].ions)) + '/' + str(len(msdata_unformatted.index)) + ' ' + str(round(100 * len(ionfilters['insource'].ions) / len(msdata_unformatted.index), 2)) + '%\n'
    if analysis_params.CVfil:
        text += 'Features failing CV filtering: ' + str(len(ionfilters['cv'].ions)) + '/' + str(len(msdata_unformatted.index)) + ' ' + str(round(100 * len(ionfilters['cv'].ions) / len(msdata_unformatted.index), 2)) + '%\n'
    text += 'Features failing any filters: ' + str(len(msdata_unformatted.index) - len(msdata_filtered.index)) + '/' + str(len(msdata_unformatted.index)) + ' ' + str(round(100 * (len(msdata_unformatted.index) - len(msdata_filtered.index)) / len(msdata_unformatted.index), 2)) + '%\n'
    text += 'Features passing all filters: ' + str(len(msdata_filtered.index)) + '/' + str(len(msdata_unformatted.index)) + ' ' + str(round(100 * len(msdata_filtered.index) / len(msdata_unformatted.index), 2)) + '%\n'

    analysisrec.writelines([text,
                            '\n',
                            '\n',
                            '---Graphing Parameters---\n',
                            'Filters: \n'])

    for elem in normalize_graphfilters(analysis_params.graphfilters):
        analysisrec.write(elem + '\n')

    analysisrec.writelines(['\n',
                            '-Groups-\n'])

    for elem in analysis_params.querylist:
        analysisrec.write(elem + '\n')

    analysisrec.writelines(['\n',
                            '-Plots generated-\n',
                            'RT/mz: ' + str(analysis_params.MZRTplt) + '\n',
                            'RT/mz/FC: ' + str(analysis_params.FC3Dplt) + ' ' + str(analysis_params.statstgrps) + '\n',
                            'KMD/mz ' + str(analysis_params.KMD) + '\n',
                            #'KMD/mz/RT ' + str(analysis_params.___) + '\n',
                            'PCA unfiltered: ' + str(analysis_params.PCA) + '\n',
                            'PCA filtered: ' + str(analysis_params.PCA) + '\n',
                            'Dendrogram (ward) unfiltered: ' + str(analysis_params.Dendrogram) + '\n',
                            'Dendrogram (ward) Filtered: ' + str(analysis_params.Dendrogram) + '\n',
                            'Volcano plot: ' + str(analysis_params.Volcanoplt) + ' ' + str(analysis_params.statstgrps) + '\n',])
    analysisrec.close()

    return AnalysisResult(
        ionfilters=ionfilters,
        groupionlists=groupionlists,
        groupsets=groupsets,
        filtereddfs=filtereddfs,
    )

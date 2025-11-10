"""
MPACT
Copyright 2022, Robert M. Samples, Sara P. Puckett, and Marcy J. Balunas
"""


import pandas as pd
import numpy as np
from pathlib import Path


def format_check(parent):
    """
    Check if the file has the right format, and if not, reformat it.
    """
    try:
        if parent.filename.suffix == '.txt':
            reformat_msdial(Path(parent.filename))
            parent.filename = Path(str(parent.filename)[:-4] + '.csv')

            
        if parent.filename.suffix == '.csv':
            msdata = pd.read_csv(parent.filename, sep = ',', header = None, index_col = None) #imports feature list
            if msdata.iloc[0, 0] == 'row ID':
                reformat_mzmine(Path(parent.filename))
            if msdata.iloc[0, 0] == 'Bucket label':
                reformat_metaboscape(Path(parent.filename))
                
            msdata = pd.read_csv(Path(parent.filename), sep = ',', header = [0,1,2], index_col = [0]) #imports data
            if len(msdata[~msdata.index.duplicated()].index) > 0:
                rename_duplicates(Path(parent.filename))
            
    except Exception:
        pass
        return()


def reformat_metaboscape(file):
    """Reformat a file from the Metaboscape software.
    
    Args:
    - file (Path): the path of the file to reformat
    """
    data = pd.read_csv(file, sep = ',', header = None, index_col = None)
    data.iloc[0:3, 0] = ('','','Compound')
    data.iloc[0:3, 1] = ('','','Retention time (min)')
    data.iloc[0:3, 2] = ('','','m/z')
    
    mz = list(data.iloc[:,2])
    data.iloc[:,2] = data.iloc[:,1]
    data.iloc[:,1] = mz
    
    injections = list(data.iloc[0,3:])
    groups = list(data.iloc[1,3:])
    samples = list(data.iloc[2,3:])
    
    data.iloc[0,3:] = groups
    data.iloc[1,3:] = samples
    data.iloc[2,3:] = injections
    
    data=data.drop(data.columns[3:5], axis=1)
    data.to_csv(file, header = False, index = False) #saves formatted backup for later use This line might break on import of new files...
    data = pd.read_csv(file, sep = ',', header = [0,1,2], index_col = [0]) #imports data
    
    data.index = data.index + data.groupby(level=0).cumcount().astype(str).replace('0','')
    data.to_csv(file, header = True, index = True) #saves formatted backup for later use


def reformat_mzmine(file):
    """Reformat a file from the MZmine software (updated for newer versions)."""
    data = pd.read_csv(file, sep=',', header=None, index_col=None)
    
    # Find where sample columns start (look for file extensions or "Peak area")
    sample_start_col = None
    for idx, val in enumerate(data.iloc[0, :]):
        if pd.notna(val) and ('.mzML' in str(val) or '.mzXML' in str(val) or 
                               '.raw' in str(val) or '.cdf' in str(val) or
                               'Peak area' in str(val)):
            sample_start_col = idx
            break
    
    if sample_start_col is None:
        # Find last metadata column
        metadata_cols = ['row ID', 'row m/z', 'row retention time', 'row ion mobility', 
                        'row ion mobility unit', 'row CCS', 'correlation', 'annotation']
        last_metadata_idx = 0
        for idx, val in enumerate(data.iloc[0, :]):
            if pd.notna(val) and any(meta in str(val) for meta in metadata_cols):
                last_metadata_idx = max(last_metadata_idx, idx)
        sample_start_col = last_metadata_idx + 1
    
    # Find m/z and RT columns
    mz_col = rt_col = id_col = None
    for idx, val in enumerate(data.iloc[0, :sample_start_col]):
        if pd.notna(val):
            val_str = str(val).lower()
            if 'row m/z' in val_str:
                mz_col = idx
            elif 'row retention time' in val_str:
                rt_col = idx
            elif 'row id' in val_str:
                id_col = idx
    
    # Defaults
    id_col = id_col if id_col is not None else 0
    mz_col = mz_col if mz_col is not None else 1
    rt_col = rt_col if rt_col is not None else 2
    
    # Find where valid sample columns end (exclude empty columns)
    sample_end_col = len(data.columns)
    for col_idx in range(len(data.columns) - 1, sample_start_col - 1, -1):
        # Check if column header is empty or column has all NaN/zero values
        if pd.isna(data.iloc[0, col_idx]) or str(data.iloc[0, col_idx]).strip() == '':
            sample_end_col = col_idx
        else:
            break  # Found a valid column, stop looking
    
    # Build new data structure
    new_data = []
    
    # Two empty header rows
    num_cols = 3 + (sample_end_col - sample_start_col)
    new_data.append([np.nan] * num_cols)
    new_data.append([np.nan] * num_cols)
    
    # Header row
    header_row = ['Compound', 'm/z', 'Retention time (min)']
    for col_idx in range(sample_start_col, sample_end_col):
        sample_name = str(data.iloc[0, col_idx])
        if pd.notna(data.iloc[0, col_idx]):
            # Clean sample name
            for ext in ['.mzML', '.mzXML', '.raw', '.cdf', '.mzData']:
                sample_name = sample_name.replace(ext, '')
            sample_name = sample_name.replace(' Peak area', '').strip()
            header_row.append(sample_name)
    
    new_data.append(header_row)
    
    # Data rows
    for row_idx in range(1, len(data)):
        mz_val = data.iloc[row_idx, mz_col]
        rt_val = data.iloc[row_idx, rt_col]
        
        # Create compound ID as RT_mz
        if pd.notna(rt_val) and pd.notna(mz_val):
            compound_id = f"{rt_val}_{mz_val}"
        else:
            compound_id = str(data.iloc[row_idx, id_col] if pd.notna(data.iloc[row_idx, id_col]) else row_idx)
        
        new_row = [compound_id, mz_val if pd.notna(mz_val) else '', rt_val if pd.notna(rt_val) else '']
        
        # Add sample data (only valid columns)
        for col_idx in range(sample_start_col, sample_end_col):
            new_row.append(data.iloc[row_idx, col_idx])
        
        new_data.append(new_row)
    
    # Save
    result_df = pd.DataFrame(new_data)
    result_df.to_csv(file, header=False, index=False)


def reformat_msdial(file):
    """Reformat a file from the MS-DIAL software.
    
    Args:
    - file (Path): the path of the file to reformat
    """
    database = pd.read_csv(file, sep = '\t', header = None, index_col = None) #imports data
    database.iloc[:,3] = database.iloc[:,2]
    database.iloc[:,2] = database.iloc[:,1]
    database.iloc[:,1] = database.iloc[:,3]
    headser = database.iloc[3,:] == 'Stdev'
    cutlen = len(headser[headser == True]) * 2 #cut cols for stdev and average
    database2 = database.iloc[:,0:3].join(database.iloc[:,35: -cutlen])
    database2.iloc[2,:] = database2.iloc[1,:]
    database2.iloc[3,:] = database2.iloc[1,:]
    database2 = database2.iloc[2:,:]
    database2.iloc[2,:3] = ['Compound', 'm/z', 'Retention time (min)']
    database2.to_csv(file.parent / (str(file.stem) + '.csv'), header = False, index = False) #saves formatted backup for later use This line might break on import of new files...
    database = pd.read_csv(file.parent / (str(file.stem) + '.csv'), sep = ',', header = [0,1,2], index_col = None) #imports data
    database.iloc[:,0] = database.iloc[:,2].astype(str)+"_"+database.iloc[:,1].astype(str)
    database.to_csv(file.parent / (str(file.stem) + '.csv'), header = True, index = False) #saves formatted backup for later use This line might break on import of new files...

def rename_duplicates(file):
    """Rename duplicate entries in the file's index to make them unique.

    Args:
    - file (Path): the path of the file to reformat
    """
    data = pd.read_csv(file, sep = ',', header = [0,1,2], index_col = [0]) #imports data
    
    data.index = data.index + data.groupby(level=0).cumcount().astype(str).replace('0','')
    data.to_csv(file, header = True, index = True) #saves formatted backup for later use
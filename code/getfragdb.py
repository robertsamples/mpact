"""
MPACT
Copyright 2022, Robert M. Samples, Sara P. Puckett, and Marcy J. Balunas
Modifications Copyright 2026, Robert M. Samples
"""


import numpy as np

class ion():
    """
    A class representing a ion with associated fragmentation data.
    
    Attributes:
    -----------
    fragparams: dict
        A dictionary of parameters for the ion.
    pattern: np.ndarray
        A numpy array representing the ion's pattern.
    """
    
    def __init__(self, fragparams, pattern):
        self.fragparams = fragparams
        self.pattern = pattern

class fragmentation_db():
    """
    A class representing a database of fragmentation ions.
    
    Attributes:
    -----------
    regex: list
        A list of regular expressions for the database.
    ions: dict
        A dictionary of ions in the database.
    """
    
    def __init__(self, regex):
        self.regex = regex
        self.ions = {}

def importfrag_v1(fragfile):
    """
    Imports a fragmentation database from a Progenesis-style MSP file.
    
    Parameters:
    -----------
    fragfile: str
        The path to the MSP file.
    
    Returns:
    --------
    fragmentation_db
        A database of fragmentation ions.
    """
    
    fragmsp = open(fragfile, 'r')
    regex = []
    while True:
        line = fragmsp.readline()
        if (not line) or (':' not in line):
            break
        regex.append(line.split(':')[0])
    fragmsp.close()
    
    fragmsp = open(fragfile, 'r')
    namelist = []
    while True:
        line = fragmsp.readline()
        if not line:
            break
        if 'Name:' in line:
            namelist.append(line.split('(')[1].split(')')[0])
    fragmsp.close()
    
    
    fragmsp = open(fragfile, 'r')
    pairlist = []
    fragparams = {}
    pair = []
    fragdb = fragmentation_db(regex)
    while True:
        line = fragmsp.readline()
        if (not line):
            break
        if 'Name:' in line:
            fragparams = {}
            if line.split('(')[1].split(')')[0] in namelist:
                name = line.split('(')[1].split(')')[0]
                pairlist = []
                while len(line)>5:
                    if ':' in line:
                        fragparams[line.split(':')[0]] = line.split(':')[1]
                    else:
                        pair = line.split('\n')[0].split(' ') #cuts line break an splits at space to seperate mass and abundance
                        pairlist.append([float(pair[0]), float(pair[1])]) #appeds float list
                    line = fragmsp.readline()
                outputarr = np.array(pairlist)
                fragdb.ions[name] = ion(fragparams, outputarr)
    fragmsp.close()
    
    return(fragdb)

def importfrag_v2(fragfile):
    """
    Imports a fragmentation database from an MS-DIAL-style MSP file.
    
    Parameters:
    -----------
    fragfile: str
        The path to the MSP file.
    
    Returns:
    --------
    fragmentation_db
        A database of fragmentation ions.
    """
    fragdb = fragmentation_db([])
    with open(fragfile, 'r') as fragmsp:
        content = fragmsp.read()
    
    # Split the file content into entries for each ion
    entries = content.strip().split('\n\n')
    for entry in entries:
        lines = entry.strip().split('\n')
        if not lines:
            continue
        fragparams = {}
        pairlist = []
        for line in lines:
            if ':' in line:
                parts = line.split(':', 1)
                key = parts[0].strip().upper()  # Use upper case for consistency
                value = parts[1].strip()
                fragparams[key] = value
            elif line.strip():
                pair = line.strip().split()
                if len(pair) >= 2:
                    try:
                        mz = float(pair[0])
                        intensity = float(pair[1])
                        pairlist.append([mz, intensity])
                    except ValueError:
                        # Handle cases where conversion to float fails
                        continue
        if 'RETENTIONTIME' in fragparams and 'PRECURSORMZ' in fragparams:
            try:
                retention_time = float(fragparams['RETENTIONTIME'])
                precursor_mz = fragparams['PRECURSORMZ']
                name = f"{round(retention_time, 3)}_{precursor_mz}"
                outputarr = np.array(pairlist)
                fragdb.ions[name] = ion(fragparams, outputarr)
            except ValueError:
                # Handle cases where conversion to float fails
                continue
    return fragdb


def importfrag_mgf(fragfile):
    """
    Imports a fragmentation database from an MS-DIAL-style MGF file.

    Entries are keyed the same way as importfrag_v2 keys the matching MSP
    export (rounded retention time and precursor m/z), so either file works
    with the same peak table. Entries without a retention time are skipped.

    Parameters:
    -----------
    fragfile: str
        The path to the MGF file.

    Returns:
    --------
    fragmentation_db
        A database of fragmentation ions.
    """
    # MGF keys and the MSP keys importfrag_v2 stores for the same values.
    msp_keys = {'TITLE': 'NAME', 'PEPMASS': 'PRECURSORMZ', 'ION': 'PRECURSORTYPE',
                'RTINMINUTES': 'RETENTIONTIME'}
    fragdb = fragmentation_db([])
    fragparams, pairlist = None, []
    with open(fragfile, 'r') as fragmgf:
        for line in fragmgf:
            stripped = line.strip()
            upper = stripped.upper()
            if upper == 'BEGIN IONS':
                fragparams, pairlist = {}, []
            elif upper == 'END IONS':
                if fragparams is not None:
                    _add_mgf_entry(fragdb, fragparams, pairlist)
                fragparams = None
            elif fragparams is None or not stripped:
                continue
            elif '=' in stripped:
                key, value = stripped.split('=', 1)
                key = key.strip().upper()
                value = value.strip()
                if key == 'PEPMASS':
                    value = value.split()[0] if value else value
                elif key == 'RTINSECONDS':
                    try:
                        key, value = 'RTINMINUTES', str(float(value) / 60)
                    except ValueError:
                        continue
                fragparams[msp_keys.get(key, key)] = value
            else:
                pair = stripped.split()
                if len(pair) >= 2:
                    try:
                        pairlist.append([float(pair[0]), float(pair[1])])
                    except ValueError:
                        continue
    return fragdb


def _add_mgf_entry(fragdb, fragparams, pairlist):
    """Store one parsed MGF entry under the importfrag_v2 name."""
    try:
        retention_time = float(fragparams['RETENTIONTIME'])
        precursor_mz = fragparams['PRECURSORMZ']
    except (KeyError, ValueError):
        return
    name = f"{round(retention_time, 3)}_{precursor_mz}"
    fragdb.ions[name] = ion(fragparams, np.array(pairlist))


def importfrag(fragfile):
    """
    Imports a fragmentation database from an MSP or MGF file, detecting the file type automatically.

    Parameters:
    -----------
    fragfile: str
        The path to the MSP or MGF file.

    Returns:
    --------
    fragmentation_db
        A database of fragmentation ions.
    """
    has_parentheses = False
    has_pipe = False

    with open(fragfile, 'r') as fragmsp:
        for line in fragmsp:
            if line.strip().upper() == 'BEGIN IONS':
                print('MGF file Detected')
                return importfrag_mgf(fragfile)
            if line.upper().startswith('NAME:'):
                if '|' in line:
                    has_pipe = True
                if '(' in line and ')' in line:
                    has_parentheses = True
                # Break after checking the first NAME: line
                break

    if has_parentheses:
        print('Progenesis MSP file Detected')
        return importfrag_v1(fragfile)
    elif has_pipe:
        print('MS-DIAL MSP file Detected')
        return importfrag_v2(fragfile)
    else:
        print('Unknown MSP file format. Attempting MS-DIAL parsing by default.')
        return importfrag_v2(fragfile)

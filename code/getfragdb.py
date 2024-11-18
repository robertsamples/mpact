
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


def importfrag(fragfile):
    """
    Imports a fragmentation database from an MSP file, detecting the file type automatically.

    Parameters:
    -----------
    fragfile: str
        The path to the MSP file.

    Returns:
    --------
    fragmentation_db
        A database of fragmentation ions.
    """
    has_parentheses = False
    has_pipe = False

    with open(fragfile, 'r') as fragmsp:
        for line in fragmsp:
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

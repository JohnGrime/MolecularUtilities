#!/usr/bin/python

#
#   Author: John Grime, The University of Chicago.
#

import sys, math, PDB


def expand_range( range_string, range_sep ):
    """
    Convert a range string into a list of integers

    Args:
        range_string (string): text representation of integer range
        range_sep (string): range separator character

    Returns:
        list of integers from start index to stop index INCLUSIVE, or None if
        failed to convert input string.

    Example:

        >>> expand_range( '1-10', '-' )
        [1,2,3,4,5,6,7,8,9,10]

        >>> expand_range( '5:9', ':' )
        [5,6,7,8,9]
    """
    tokens = range_string.split( range_sep )
    if len(tokens) < 2: return None

    try:
        start = int(tokens[0])
        stop = int(tokens[1])
    except:
        return None
    
    indices = range( start, stop+1 )
    return indices


#
# Calls expand_range() on a set of range strings separated as specified.
#
def expand_ranges( range_strings, entry_sep, range_sep ):
    """
    Convert a list of range strings into a combined list of integers

    Args:
        range_strings (string): text representation of set of integer ranges
        entry_sep (string): separator for the individual range strings
        range_sep (string): range separator character

    Returns:
        combined list of integers from all defined ranges, or None if
        failure to parse/convert the range string.

    Example:

        >>> expand_range( '1-4:21-25:111-113', ':', '-' )
        [1,2,3,4,21,22,23,24,25,111,112,113]
    """
    indices = []
    tokens = range_strings.split( entry_sep )
    for t in tokens:
        ind = expand_range( t, range_sep )
        if ind == None: return None
        indices.append( ind )
    return indices


def renumber( molecules, delta ):
    """
    Add specified delta to the resSeq entry of all specified PDB molecules

    Args:
        molecules (list of PDB molecules): mmolecular data to modify
        delta (integer): offset to add to resSeq in each molecule

    Returns:
        list of updated PDB molecules
    """
    new_molecules = []
    for mol in molecules:
        new_mol = []
        for a in mol:
            new_a = dict( a )
            new_a['resSeq'] += delta
            new_mol.append( new_a )
        new_molecules.append( new_mol )
    return new_molecules


#
# Rename chains sequentially, according to PDB spec
#
def fix_chains( molecules ):
    """
    Rename chains in PDB molecules to be sequential as per PDB specification

    Args:
        molecules (list of PDB molecules): mmolecular data to modify

    Returns:
        list of updated PDB molecules
    """
    chains = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789'
    new_molecules = []
    for mol in molecules:
        new_mol = []
        chain = chains[ len(new_molecules) % len(chains) ]
        for a in mol:
            new_a = dict( a )
            new_a['chainID'] = chain
            new_mol.append( new_a )
        new_molecules.append( new_mol )
    return new_molecules


def extract_chains( molecules, chain_IDs, preserve_order = False ):
    """
    Extract specified chains from PDB molecules where

    Args:
        molecules (list of PDB molecules): mmolecular data to modify
        chain_IDs (list of strings): chain IDs to extract
        preserve_order (boolean): if True, maintain order of chains in the file. Otherwise, use order of appearance from chain_IDs

    Returns:
        list of extracted chains as PDB molecules
    """
    new_molecules = []

    # preserve ordering from PDB file
    if preserve_order == True:
        for mol in molecules:
            if (len(mol)<1) or (mol[0]['chainID'] not in chain_IDs): continue
            new_mol = [ dict(a) for a in mol ]
            new_molecules.append( new_mol )
    # output in order specified by chain_IDs
    else:
        chain_to_index = {}
        for i in range( 0, len(molecules) ):
            if (len(molecules[i])<1): continue
            chain_to_index[ molecules[i][0]['chainID'] ] = i
        for cID in chain_IDs:
            if cID not in chain_to_index: continue
            new_mol = [ dict(a) for a in molecules[chain_to_index[cID]] ]
            new_molecules.append( new_mol )

    return new_molecules


def extract_mols( molecules, set_size, set_indices ):
    """
    Extract specified sets of contiguous PDB molecules, with a specified set size

    Args:
        molecules (list of PDB molecules): mmolecular data
        set_size (integer): number of consecutive atoms to view as a distinct set
        set_indices (list of integers): indices of sets to extract (UNIT BASED!)

    Returns:
        list of extracted atom sets, with each set as a separate PDB molecule
    """
    new_molecules = []
    for set_i in set_indices:
        i = ( (set_i-1)*set_size ) # convert to UNIT BASED!
        j = ( set_i*set_size )
        for mol_i in range(i,j):
            new_mol = [ dict(a) for a in molecules[mol_i] ]
            new_molecules.append( new_mol )
    return new_molecules


def filter_molecules( molecules, filter_strings ):
    """
    Filter molecules to produce a new set of molecules containing only atoms that passed the filter(s)

    Args:
        molecules (list of PDB molecules): mmolecular data to filter
        filter_strings (list of strings): filters to apply, using PDB ATOM-style fields

    Returns:
        list of molecules containing filtered data.

    Notes:

        - filter strings are in the form `key=value,key=value,...`, where numerical ranges are INCLUSIVE and
          indicated using a dash to separate the start and end indices, e.g. `name=CA,resSeq=2-14`
        - if no atoms in a molecule pass the filter, an empty molecule results
    """
    filters = {}
    for f in filter_strings:
        PDB.UpdateFilters( f, filters, '=', ',', '-' )
    
    new_molecules = []
    for mol in molecules:
        new_mol = PDB.FilterAtoms( mol, filters )
        new_molecules.append( new_mol )
    
    return new_molecules


def scale_molecules( molecules, scale ):
    """
    Scale all coordinates in molecules

    Args:
        molecules (list of PDB molecules): mmolecular data to scale
        scale (float): scale factor to apply

    Returns:
        list of molecules containing scaled data.
    """
    new_molecules = []
    for mol in molecules:
        new_mol = []
        for a in mol:
            new_a = dict( a )
            new_a['x'] *= scale
            new_a['y'] *= scale
            new_a['z'] *= scale
            new_mol.append( new_a )
        new_molecules.append( new_mol )
    return new_molecules


def recentre( molecules ):
    """
    Translate molecules such that their combined centre-of-geometry is at the origin

    Args:
        molecules (list of PDB molecules): mmolecular data to translate

    Returns:
        list of translated molecules
    """
    new_molecules = []
    x, y, z, N = 0.0, 0.0, 0.0, 0

    for mol in molecules:
        for a in mol:
            x += a['x']
            y += a['y']
            z += a['z']
            N += 1
    x /= N
    y /= N
    z /= N

    for mol in molecules:
        new_mol = []
        for a in mol:
            new_a = dict( a )
            new_a['x'] -= x
            new_a['y'] -= y
            new_a['z'] -= z
            new_mol.append( new_a )
        new_molecules.append( new_mol )
    return new_molecules


def print_usage( prog ):
    print ''
    print 'Usage: cat whatever.pdb | %s [renumber resSeq_offset] | [fix_chains] | [recentre] | [scale X] | [extract_chains A B C ...] | [extract_mols set_size set_i set_j set_k ...] | [filter key=val,val,... key=val,val,... ]' % ( prog )
    print ''
    print 'Where:'
    print '  - renumber        : add resSeq_offset to all residue sequence numbers'
    print '  - fix_chains      : sequentially rename chains to A B C ...'
    print '  - recentre        : move centre of geometry to the origin'
    print '  - scale           : scale coordinates by specified factor'
    print '  - extract_chains  : extract the specified chains (IN THE ORDER SPECIFIED!)'
    print '  - extract_mols    : extract the UNIT_BASED sets of molecules with sets of specified size from the PDB file (molecule = a TER- or MODEL-separated entry)'
    print '  - filter          : filter input PDB with specified filters. Range separator char for resSeq is "-"'
    print ''
    print 'In all cases, the resultant PDB file is printed to stdout.'
    print ''
    sys.exit( -1 )


#
# Do we have sufficient command line arguments?
#
if len(sys.argv) < 2:
    print_usage( sys.argv[0] )

#
# Get params
#
action = sys.argv[1]
params = sys.argv[2:]

#
# Get PDB data
#
if sys.stdin.isatty() == True:
    print_usage( sys.argv[0] )

lines = sys.stdin.readlines()
molecules = PDB.GetPDBMolecules( lines )

#
# Figure out what to do, and generate new_molecules
# 
if action == 'renumber':
    if len(params) < 1:
        print_usage( sys.argv[0] )
    try:
        delta = int( params[0] )
    except:
        print 'Unable to convert "%s" into an integer' % ( params[0] )
        sys.exit( -1 )
    new_molecules = renumber( molecules, delta )
    
elif action == 'fix_chains':
    new_molecules = fix_chains( molecules )
    
elif action == 'recentre':
    new_molecules = recentre( molecules )

elif action == 'extract_chains':
    # split param tokens into individual characters, if needed.
    chains = []
    for tok in params:
        chains += [ c for c in tok if c not in ' \t\n\r' ]
    new_molecules = extract_chains( molecules, chains )

elif action == 'extract_mols':
	if len(params) < 2:
		print_usage( sys.argv[0] )
	new_molecules = extract_mols( molecules, int(params[0]), [ int(mol_i) for mol_i in params[1:] ] )

elif action == 'filter':
    new_molecules = filter_molecules( molecules, params )

elif action == 'scale':
    new_molecules = scale_molecules( molecules, float(params[0]) )

#
# Print new_molecules
#
for mol in new_molecules:
    for a in mol:
        outline = PDB.MakePDBAtomLine( a )
        print outline
    print 'TER   '

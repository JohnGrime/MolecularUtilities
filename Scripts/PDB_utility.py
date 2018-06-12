#
#   Author: John Grime, The University of Chicago.
#

import sys, math, PDBTools

#
# Convert string to the type 'i<range_sep>j' into a range of numbers from i to j inclusive.
#
def expand_range( range_string, range_sep ):
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
    indices = []
    tokens = range_strings.split( entry_sep )
    for t in tokens:
        ind = expand_range( t, range_sep )
        if ind == None: return None
        indices.append( ind )
    return indices

#
# Add delta to all resSeq
#
def renumber( molecules, delta ):
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


#
# Extract specified chains
#
def extract_chains( molecules, chain_IDs, preserve_order = False ):
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

#
# Extract specified molecules (set_indices are UNIT BASED!)
#
def extract_mols( molecules, set_size, set_indices ):
	new_molecules = []
	for set_i in set_indices:
		i = ( (set_i-1)*set_size ) # convert to UNIT BASED!
		j = ( set_i*set_size )
		for mol_i in range(i,j):
			new_mol = [ dict(a) for a in molecules[mol_i] ]
			new_molecules.append( new_mol )
	return new_molecules

#
# Filter atoms
#
def filter_molecules( molecules, filter_strings ):
    filters = {}
    for f in filter_strings:
        PDBTools.UpdateFilters( f, filters, '=', ',', '-' )
    
    new_molecules = []
    for mol in molecules:
        new_mol = PDBTools.FilterAtoms( mol, filters )
        new_molecules.append( new_mol )
    
    return new_molecules

#
# Scale all coords
#
def scale_molecules( molecules, scale ):
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

#
# Reset centre of geometry
#
def recentre( molecules ):
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
    print ''
    print '\t renumber        : add resSeq_offset to all residue sequence numbers'
    print '\t fix_chains      : sequentially rename chains to A B C ...'
    print '\t recentre        : move centre of geometry to the origin'
    print '\t scale           : scale coordinates by specified factor'
    print '\t extract_chains  : extract the specified chains (IN THE ORDER SPECIFIED!)'
    print '\t extract_mols    : extract the UNIT_BASED sets of molecules with sets of specified size from the PDB file (molecule = a TER- or MODEL-separated entry)'
    print '\t filter          : filter input PDB with specified filters. Range separator char for resSeq is "-"'
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
#pdb_path = sys.argv[1]
action = sys.argv[1]
params = sys.argv[2:]

#
# Get PDB data
#
if sys.stdin.isatty() == True:
    print_usage( sys.argv[0] )

#lines = open( pdb_path ).readlines()
lines = sys.stdin.readlines()
molecules = PDBTools.GetPDBMolecules( lines )

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
        outline = PDBTools.MakePDBAtomLine( a )
        print outline
    print 'TER   '

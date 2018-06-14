#!/usr/bin/python

#
#   Author: John Grime, The University of Chicago.
#

import sys, math

#
# This script identifies and extracts structures from a PDB file continaing an HIV CA assembly.
# Depending on what is included and missing in the source file, you may need to modify the atoms
# that the script uses to determine what's connected to what.
#
# The script outputs the following files:
#
# dimers.pdb : CTD/CTD dimer pairs
# rings.N.pdb : a set of files containing cyclic NTD/NTD structures, with N denoting the number of monomers in the cyclic structure.
# trimer.pdb : trimer-of-dimer structures
#
# IMPORTANT: the script assumes monomers are composed of a certain number of consecutive "subunits" in the PDB file, with a subunit
# delineated by a TER entry. I typically split the NTD and CTD of the CA proteins into separate subunits, which are arranged as
# NTD then CTD for each monomer in the source files. These subunit delineations are preserved in the output files.
#



#############################################################
# First section of this file is just PDB file IO stuff.
#############################################################



#
# The following dictionary is used to parse PDB atoms; 'start' and 'end' are
# the UNIT BASED text column indices (as per the 'ATOM' definition in the RCSB),
# and 'converter' is a function to convert the raw text string into something
# more natural for that specific data.
#
PDB_atom_line_info = {
    'type':       { 'start': 1, 'stop': 6, 'converter':str },
    'serial':     { 'start': 7, 'stop':11, 'converter':int },
    'serial':     { 'start': 7, 'stop':11, 'converter':lambda x: int(x,16) }, # for weird big PDBs with hex serial: use base 16.
    'name':       { 'start':13, 'stop':16, 'converter':str },
    'altLoc':     { 'start':17, 'stop':17, 'converter':str },
    'resName':    { 'start':18, 'stop':20, 'converter':str },
    'chainID':    { 'start':22, 'stop':22, 'converter':str },
    'resSeq':     { 'start':23, 'stop':26, 'converter':int },
    'iCode':      { 'start':27, 'stop':27, 'converter':str },
    'x':          { 'start':31, 'stop':38, 'converter':float },
    'y':          { 'start':39, 'stop':46, 'converter':float },
    'z':          { 'start':47, 'stop':54, 'converter':float },
    'occupancy':  { 'start':55, 'stop':60, 'converter':float },
    'tempFactor': { 'start':61, 'stop':66, 'converter':float },
    'element':    { 'start':77, 'stop':78, 'converter':str },
    'charge':     { 'start':79, 'stop':80, 'converter':str }
}

#
# Text line => PDB atom dictionary, assumes valid ATOM line passed in.
#
def parse_PDB_line( line, line_info ):
    parsed = {}
    for key in line_info.keys():
        li = line_info[key]
        converter, i, j = li['converter'], li['start']-1, li['stop']

        if i >= len(line): continue
        if j > len(line): j = len(line)

        # special case; big PDB files from VMD can have serial=='*****' where > 5 digits
        if (key=='serial') and (line[i] == '*'): parsed[key] = 1
        else: parsed[key] = converter( line[i:j].strip() ) # remove leading/trailing whitespace before conversion
    return parsed

#
# Revert atom data => PDB compatible text string
#
def make_PDB_line( atom ):
    PDB_atom_line_format = '%-6.6s%5.5s %4.4s%1.1s%3.3s %1.1s%4.4s%1.1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2.2s%2.2s'
    defaults = { 'occupancy':1.0, 'tempFactor':0.0, 'element':'', 'charge':'' }

    # fill in common missing values with defaults, if needed
    checked_values = {}
    for key in defaults.keys():
        checked_values[key] = defaults[key]
        if key in atom:
            checked_values[key] = atom[key]

    return PDB_atom_line_format % (
        atom['type'],
        str(atom['serial']),
        atom['name'],
        atom['altLoc'],
        atom['resName'],
        atom['chainID'],
        str(atom['resSeq']),
        atom['iCode'],
        atom['x'],
        atom['y'],
        atom['z'],
        checked_values['occupancy'],
        checked_values['tempFactor'],
        checked_values['element'],
        checked_values['charge'] )



##############################################################################
# Utility class to load PDB data, with different ways to reference the atoms.
# The PDBData class also allows you to build distance tables, and look for 
# connectivity between sets of atoms.
##############################################################################



#
# Class to simplify mapping between atoms/subunits/monomers,
# with some functionality to create filtered sets of atoms and distance tables.
# This treatment allows several "subunits" to be considered as a single monomer,
# so we can handle PDB files with eg the NTD and CTD as separate entries separated
# by a TER line.
#
class PDBData:
    
    def __init__( self, pdb_path, subunits_per_monomer ):
        # map of global atom id => [subunit_i,subunit_atom_i], and vice versa.
        # l2g is actually a 2D map, used as l2g[subunit][atom] => global atom id
        self.g2l = {}
        self.l2g = {}

        # each monomer is simply a list of subunit indices; "monomers" is
        # therefore a list of subunit index lists.
        self.monomers = []

        # global atom id => monomer, and subunit index => monomer
        self.g2monomer = {}
        self.s2monomer = {}

        #
        # Parse line-by-line rather than loading entire file, better for large PDB files.
        #
        last_resSeq = None
        self.subunits = []
        current_subunit = []
        valid_atom_prefixes = [ 'ATOM' ]
        f = open( pdb_path, 'r' )
        for line in f:

            if line[0:4] in valid_atom_prefixes:
                line_data = parse_PDB_line( line, PDB_atom_line_info )

                # skip non-CA
                if line_data['name'] != 'CA': continue

                # check for wrap ...
                resSeq = line_data['resSeq']
                if (last_resSeq != None) and (resSeq < last_resSeq):
                    if len(current_subunit) > 0:
                        self.subunits.append( current_subunit )
                    current_subunit = []
#                    print 'Read %d subunits so far ...' % ( len(self.subunits) )
                last_resSeq = resSeq

                current_subunit.append( line_data )
            elif line[0:3] == 'TER':
                if len(current_subunit) > 0:
                    self.subunits.append( current_subunit )
                current_subunit = []
#                if len(self.subunits) % 10 == 0: print 'Read %d subunits so far ...' % ( len(self.subunits) )

        #
        # Include trailing atom data, if no final TER present in file
        #
        if len(current_subunit) > 0:
            self.subunits.append( current_subunit )
        f.close()

        #
        # Add subunit indices to monomers
        #
        for i in range( 0, len(self.subunits), subunits_per_monomer ):
            monomer = range( i, i+subunits_per_monomer )
            self.monomers.append( monomer )

        #
        # Generate maps
        #
        global_i = 0
        for monomer_i in range( 0, len(self.monomers) ):
            for si in self.monomers[monomer_i]:
                for ai in range( 0, len(self.subunits[si]) ):

                    # global atom index => local [subunit,offset]
                    self.g2l[global_i] = [si,ai]

                    # reverse mapping for the above
                    if si not in self.l2g:
                        self.l2g[si] = {}
                    self.l2g[si][ai] = global_i

                    # mapping from global atom index into monomer,
                    # and from subunit into monomer.
                    self.g2monomer[global_i] = monomer_i
                    self.s2monomer[si] = monomer_i

                    # increment global atom identifier
                    global_i += 1

    #
    # Convert between global and local atom indices
    #
    def local_to_global( self, si, ai ): return self.l2g[si][ai]
    def global_to_local( self, gi ): return self.g2l[gi]

    #
    # Return atom data for global or local atom indices
    #
    def local_to_atom( self, si, ai ): return self.subunits[si][ai]
    def global_to_atom( self, gi ):
        si, ai = self.g2l[gi]
        return self.subunits[si][ai]

    #
    # Look up monomer from subunit, and vice versa
    #
    def subunit_to_monomer( self, si ): return self.s2monomer[si]
    def monomer_to_subunits( self, mi ): return self.monomers[mi]

    #
    # Return filtered subset of atoms as a list of global ids
    #
    def filter_atoms( self, filters ):
        filtered_atoms = []
        for gi in self.g2l:
            atom, add_me = self.global_to_atom( gi ), True
            for key in filters:
                if (key in atom) and (atom[key] not in filters[key]):
                    add_me = False
                    break
            if add_me == True: filtered_atoms.append( gi )
        return filtered_atoms

    #
    # Table coords [i][j] are monomer indices, but i_gids and j_gids are
    # lists of global atom indices.
    #
    def dr2_table_as_monomers( self, i_gids, j_gids ):
        dr2 = [ [None for j in self.monomers] for i in self.monomers ]
        for i in range( 0, len(i_gids) ):
            si, ai = self.global_to_local( i_gids[i] )
            mi = self.subunit_to_monomer( si )
            a = self.local_to_atom( si, ai )
            for j in range( 0, len(j_gids) ):
                sj, aj = self.global_to_local( j_gids[j] )
                mj = self.subunit_to_monomer( sj )
                if mi == mj:
                    continue # exclude self-monomers
                b = self.local_to_atom( sj, aj )
                dx = a['x'] - b['x']
                dy = a['y'] - b['y']
                dz = a['z'] - b['z']
                dr2[mi][mj] = dx*dx + dy*dy + dz*dz
        return dr2

##############################################################################
# HIV capsid protein structural identification routines.
##############################################################################

class HIVStructures:

    def __init__( self ):
        pass

    #
    # Return CLOSEST paired entry for each index. pairs[i] == None if nothing < rcut from i
    # NOTE - returns 2 copies of a pair if symmetrical (i->j and j->i), so be careful.
    #
    def get_closest_pairs( self, dr2_table, rcut ):
        rcut2 = rcut*rcut
        pairs = {}
        for i in range( 0, len(dr2_table) ):
            closest = None
            for j in range( 0, len(dr2_table[i]) ):
                dr2 = dr2_table[i][j]
                if (i==j) or (dr2==None) or (dr2>rcut2):
                    continue
                if (closest==None) or (dr2<dr2_table[i][closest]):
                    closest = j
            pairs[i] = closest
        return pairs

    #
    # Identify dimers via closest pairs in dr2_table.
    # Returns list of [monomer_i, monomer_j] entries.
    #
    def get_dimers( self, dr2_table, rcut, symmetrical ):
        dimers = []
        pairs = self.get_closest_pairs( dr2_table, rcut )
        for mi in pairs:
            mj = pairs[mi]
            if (mj==None) or ((symmetrical==True) and (mi>mj)):
                continue # don't double count the dimers!
            dimers.append( [mi,mj] )
        return dimers

    #
    # Recursive cyclic structure identifier; uses list not dictionary, so we preserve ring ordering.
    # Returns map of cyclic structure length => list of structures of that length.
    #
    def cyclic_recurse( self, pairs, upto, previous ):
        if upto==None:
            return False, previous
        if upto in previous:
            return True, previous
        previous.append( upto )
        # TEMP - can we actually find a connection from the current?
        if upto not in pairs:
            return False, previous
        return self.cyclic_recurse( pairs, pairs[upto], previous )

    def get_cyclic( self, dat, dr2_table, rcut ):
        rings = {}
        pairs = self.get_closest_pairs( dr2_table, rcut )
        while len(pairs) > 0:
            start = pairs.keys()[0]
            result, ring = self.cyclic_recurse( pairs, start, [] )
            if result == True:
                rkey = len(ring)
                if rkey not in rings:
                    rings[rkey] = []
                rings[rkey].append( ring )
            # delete any monomers we considered from the current info, regardless of whether it formed a ring
            for m in ring:
                if m in pairs:
                    del pairs[m]
        return rings

    #
    # Identify trimer-of-dimers, given NTD/NTD and CTD/CTD monomer pair tables
    # Returns map of unique_trimer_key => trimer structure
    #
    def get_trimers( self, ntd_dr2_table, ntd_rcut, ctd_dr2_table, ctd_rcut ):
        ntd_pairs = self.get_closest_pairs( ntd_dr2_table, ntd_rcut )
        ctd_pairs = self.get_closest_pairs( ctd_dr2_table, ctd_rcut )
        
        #
        # Note: as NTD/NTD dimer directional, we can't walk around a trimer "the wrong way" by starting
        # with following a CTD/CTD dimer interface even though the CTD/CTD interface is undirected: the
        # following NTD/NTD dimer test would fail, due to the directionality of that interface. We can
        # therefore be assured that the resultant trimer-of-dimers structures we detect will always be
        # superposable onto one another, as the monomer orderings will be compatible.
        #
        
        trimers = {}
        # m1 is start monomer for each trimer identification attempt
        for m1 in range( 0, len(ntd_dr2_table) ):
            # get next monomer in dimer via CTD/CTD interface (ie m1-m2 is a CTD dimer)
            m2 = ctd_pairs[m1]
            if m2 == None: continue

            # see if m2 is connected to another monomer m3 via NTD/NTD pair
            m3 = ntd_pairs[m2]
            if m3 == None: continue

            # follow monomer m3 along CTD pair (ie m3-m4 is a CTD dimer)
            m4 = ctd_pairs[m3]
            if m4 == None: continue

            # see if m4 is connected to another monomer via NTD/NTD pair
            m5 = ntd_pairs[m4]
            if m4 == None: continue

            # follow monomer m5 along CTD pair (ie m5-m6 is a CTD dimer)
            m6 = ctd_pairs[m5]
            if m6 == None: continue

            # if m6 is bound to m1 via a NTD/NTD connection, we have a trimer-of-dimers
            test = ntd_pairs[m6]
            if (test==None) or (test!=m1): continue

            # use sorted monomer indices for the trimer key, to ensure unique
            # lists of monomers, but KEEP THE TRIMER DATA IN THE ORDER WE FOUND IT!
            # this ensures consecutive monomers in the trimer data represent CTD/CTD dimers!
            trimer = [m1,m2,m3,m4,m5,m6]
            key = ' '.join( [ '%d'%(m) for m in sorted(trimer) ] )
            if key not in trimers:
                #
                # Trick: ensure first monomer in trimer has an NTD neighbour; due to the "direction dependence"
                # of the defined NTD/NTD interface, and the threefold symmetry, making sure the first monomer
                # listed has an NTD neighbour means we can have all the trimer listed in a comparable order.
                #
                # You can justr everse the monomer order for the trimer to move in the opposite direction!
                #
                if ntd_pairs[ trimer[0] ] == None:
                    trimer = list( reversed(trimer) )
                trimers[key] = trimer

        return [ trimers[key] for key in trimers ]


#
# Slightly different to the save atom routine in PDBData;
# each structure has same "internal" chainID values, and
# structures are separated into PDB 'MODEL' sections.
# once again, we assume structures are lists of monomer ids.
#
def save_structures( fpath, structures, dat ):
    chains = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz123456789'
    f = open( fpath, 'w' )
    serial = 1
    for structure in structures:
        chain_i=0
        serial = 1 # for if we wish to reset count each time (e.g. lots of structures as in the capsid!)
        print >>f, 'MODEL'
        for monomer_i in structure:
            for subunit_i in dat.monomers[monomer_i]:
                for a in dat.subunits[subunit_i]:
                    a2 = dict( a ) # modify a copy of the data, not the original
                    a2['chainID'] = chains[ chain_i%len(chains) ]

                    a2['serial'] = serial # change this to retain original serial?
                    serial += 1

                    print >>f, make_PDB_line( a2 )
                print >>f, 'TER'
            chain_i += 1
        print >>f, 'ENDMDL'
    f.close()


##########################################################
# Command line parameter wrangling routine & error method
##########################################################


#
# convert command line parameters into something more useful, in particular
# expand the resSeq data into a list of indices. Not currently used!
#
def get_params( argv ):
    params = {}
    for i in range( 0, len(argv) ):
        tokens = argv[i].split('=')
        key = tokens[0]
        params[key] = None

        if len(tokens) <= 1:
            continue

        subtokens = tokens[1].split(',')
        params[key] = subtokens
    return params


def print_usage( prog ):
    print 'Usage: %s input=source.pdb CTD=name,resSeq,rcut NTD=name1,resSeq1,name2,resSeq2,rcut output_prefix=x subunits_per_monomer=x' % ( prog )
    print 'Where:'
    print '\t input: PDB file with subunits separated by TER lines'
    print '\t output_prefix: prefix for all output files'
    print '\t subunits_per_monomer: number of consecutive PDB subunits which form a monomer'
    sys.exit( -1 )



###########################################
# Main program starts here
###########################################


structure_identifier = HIVStructures()

#
# Get command line parameters, check we have required info
#

params = get_params( sys.argv )

required_info = { 'input':1, 'output_prefix':1, 'subunits_per_monomer':1, 'CTD':3, 'NTD':5 }
for key in required_info.keys():
    if (key not in params) or (len(params[key])!=required_info[key]):
        print_usage( sys.argv[0] )

input_PDB = params['input'][0]
output_prefix = params['output_prefix'][0]
subunits_per_monomer = int( params['subunits_per_monomer'][0] )

temp = params['CTD']
ctd_name1, ctd_resSeq1 = temp[0], int(temp[1])
ctd_name2, ctd_resSeq2 = ctd_name1, ctd_resSeq1
ctd_rcut = float(temp[2])

temp = params['NTD']
ntd_name1, ntd_resSeq1 = temp[0], int(temp[1])
ntd_name2, ntd_resSeq2 = temp[2], int(temp[3])
ntd_rcut = float(temp[4])

#
# Load input PDB data
#

print 'Loading data from "%s" ...' % ( input_PDB )
dat = PDBData( input_PDB, subunits_per_monomer )
print '\t %d subunits, %d monomers' % ( len(dat.subunits), len(dat.monomers) )

#
# Build sets of filtered atoms for CTD/CTD and NTD/NTD interface detection
#
ctd_gi1  = dat.filter_atoms( {'name':[ctd_name1], 'resSeq':[ctd_resSeq1]} )
ctd_gi2  = dat.filter_atoms( {'name':[ctd_name2], 'resSeq':[ctd_resSeq2]} )

ntd_gi1 = dat.filter_atoms( {'name':[ntd_name1], 'resSeq':[ntd_resSeq1]} )
ntd_gi2 = dat.filter_atoms( {'name':[ntd_name2], 'resSeq':[ntd_resSeq2]} )

#
# CTD/CTD dimer structures as paired monomers
#

print 'Detecting CTD dimers ...'
ctd_dr2_table = dat.dr2_table_as_monomers( ctd_gi1, ctd_gi2 )
ctd_dimers = structure_identifier.get_dimers( ctd_dr2_table, ctd_rcut, symmetrical=True )
print '\t %d CTD dimer pairs found.' % ( len(ctd_dimers) )

#
# NTD/NTD dimer structures as paired monomers
#

print 'Detecting NTD dimers ...'
ntd_dr2_table = dat.dr2_table_as_monomers( ntd_gi1, ntd_gi2 )
ntd_dimers = structure_identifier.get_dimers( ntd_dr2_table, ntd_rcut, symmetrical=False )
print '\t %d NTD dimer pairs found.' % ( len(ntd_dimers) )

#
# Cyclic NTD/NTD rings as sets of monomers
#

print 'Detecting cyclic structures ...'
ntd_rings = structure_identifier.get_cyclic( dat, ntd_dr2_table, ntd_rcut )
for key in ntd_rings:
    print '\t %d rings of length %d' % ( len(ntd_rings[key]), key )

#
# Trimer-of-dimers as monomers - needs ntd_dr2 and ctd_dr2 tables from before.
#

print 'Detecting trimer-of-dimers ...'
trimers = structure_identifier.get_trimers( ntd_dr2_table, ntd_rcut, ctd_dr2_table, ctd_rcut )
print '\t %d trimer-of-dimers found.' % ( len(trimers) )

#
# Save structural mappings to file (subunit numbers denote zero-based indices into the original PDB file)
#

fpath = '%s.structures.txt'%(output_prefix)
f=open( fpath, 'w' )
print >>f, '#'
print >>f, '# Monomers are defined in terms of subunits, zero-based indices into the molecular units of the source file'
print >>f, '# Dimers, trimers, and rings are defined in terms of monomers'
print >>f, '#'
print >>f, '# Source file: %s' % ( params['input'][0] )
print >>f, '#'

print >>f, ''
print >>f, '#'
print >>f, '# Monomers : %d' % ( len(dat.monomers) )
print >>f, '#'
print >>f, ''
for monomer in dat.monomers:
    print >>f, 'monomer ', ' '.join( [ '%6d'%(subunit) for subunit in monomer ] )
	
print >>f, ''
print >>f, '#'
print >>f, '# CTD_dimers : %d' % ( len(ctd_dimers) )
print >>f, '#'
print >>f, ''
for structure in ctd_dimers:
    print >>f, 'ctd_dimer', ' '.join( [ '%6d'%(monomer) for monomer in structure ] )
	
print >>f, ''
print >>f, '#'
print >>f, '# NTD_dimers : %d' % ( len(ntd_dimers) )
print >>f, '#'
print >>f, ''
for structure in ntd_dimers:
    print >>f, 'ntd_dimer', ' '.join( [ '%6d'%(monomer) for monomer in structure ] )
	
print >>f, ''
print >>f, '#'
print >>f, '# Trimers : %d' % ( len(trimers) )
print >>f, '#'
print >>f, ''
for structure in trimers:
    print >>f, 'trimer', ' '.join( [ '%6d'%(monomer) for monomer in structure ] )
	
for key in ntd_rings:
    print >>f, ''
    print >>f, '#'
    print >>f, '# Ring %d : %d' % ( key, len(ntd_rings[key]) )
    print >>f, '#'
    print >>f, ''
    for structure in ntd_rings[key]:
        print >>f, 'cyc%d'%(key), ' '.join( [ '%6d'%(monomer) for monomer in structure ] )
f.close()

#
# Save structures as PDB files, if they exist
#

if len(ctd_dimers) > 0:
    fpath = '%s.ctd_dimers.pdb'%( output_prefix )
    save_structures( fpath, ctd_dimers, dat )
    print '\t => %s' % ( fpath )

if len(ntd_dimers) > 0:
    fpath = '%s.ntd_dimers.pdb'%( output_prefix )
    save_structures( fpath, ntd_dimers, dat )
    print '\t => %s' % ( fpath )

if len(trimers) > 0:
    fpath = '%s.trimers.pdb'%( output_prefix )
    save_structures( fpath, trimers, dat )
    print '\t => %s' % ( fpath )

for key in ntd_rings:
    if len( ntd_rings[key]) > 0:
        fpath = '%s.rings.%d.pdb' % (output_prefix,key)
        save_structures( fpath, ntd_rings[key], dat )
        print '\t => %s' % ( fpath )

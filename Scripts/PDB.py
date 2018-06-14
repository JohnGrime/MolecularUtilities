#
#	Author: John Grime, The University of Chicago.
#

import sys

#
# Default descriptions of the elements in a PDB ATOM/HETATM line. Maps a PDB element
# name onto the following information:
#   - start/stop : UNIT BASED column indices into a PDB line string specifying the element location
#   - converter : routine to convert a string into a suitable data type for the element.
# These are basically the standard PDB fields
#
PDB_atom_line_info = {
	'type':       { 'start':1,  'stop':6,  'converter':str },
	'serial':     { 'start':7,  'stop':11, 'converter':int },
#	'serial':     { 'start':7,  'stop':11, 'converter':lambda x: int(x,16) }, # for weird big PDBs with hex serial: use base 16.
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


def parse_line( line, line_info ):
	"""
	Convert a line of text from a PDB file into a parsed atom dictionary. This is not intended to be called
	directly by user code!

	Args:
	line (string): PDB file line (assumes ATOM or HETATM entry)
	line_info (dictionary): provides retrieval of start/end columns and conversion routines

	Returns:
	Dictionary of parsed atom data
	"""
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


def GetPDBMolecules( lines ):
	"""
	Convert a set of lines from a PDB file into a set of PDB molecules.

	Args:
	  lines (list of strings): a list of lines in PDB file format.

	Returns:
	  A list of molecules, where a molecule is a list of PDB atom dictionaries.
	"""
	molecules = []
	current_molecule = []
	
	for i in range( 0, len(lines) ):
		line = lines[i].strip()
		
		# catches both 'TER   ' and 'TER', although only the
		# former is technically well-formed for PDB entries.
		if line[0:3] == 'TER':
			if len(current_molecule) > 0:
				molecules.append( current_molecule )
				current_molecule = []
			continue
		
		# convert ATOM / HETATM line into a dictionary of appropriate values.
		# this looks a little odd, but it overcomes badly formatted PDB files,
		# such as where too many digits exist in the serial columns, so some
		# digits flood into the type columns etc.
		if line[0:4]=='ATOM':
			line_data = parse_line( line, PDB_atom_line_info )
			line_data['type'] = 'ATOM  ' # force correct !
			current_molecule.append( line_data )
		elif line[0:6]=='HETATM':
			line_data = parse_line( line, PDB_atom_line_info )
			line_data['type'] = 'HETATM' # force correct !
			current_molecule.append( line_data )
			
	# include any data we've not added yet, in case of missing a final TER etc.
	if len(current_molecule) > 0:
		molecules.append( current_molecule )

	return molecules


def MakeEmptyPDBAtom():
	"""
	Utility routine to return an empty PDB atom dictionary.

	Args:
	  None.

	Returns:
	  An 'empty' PDB atom dictionary.
	"""
	atom = {}
	atom['type'] = 'ATOM'
	atom['serial'] = 1
	atom['name'] = '?'
	atom['altLoc'] = ''
	atom['resName'] = '?'
	atom['chainID'] = 'A'
	atom['resSeq'] = 1
	atom['iCode'] = ''
	atom['x'] = 0
	atom['y'] = 0
	atom['z'] = 0
	atom['occupancy'] = 0
	atom['tempFactor'] = 0
	atom['element'] = ''
	atom['charge'] =''
	return atom


def MakePDBAtomLine( atom ):
	"""
	Generates a line of text suitable for writing to a PDB file from a PDB atom dictionary.

	Args:
	  atom (dictionary): A PDB atom dictionary

	Returns:
	  A string suitable for writing as an ATOM or HETATM line in a PDB file.
	"""
	PDB_atom_line_format = '%-6.6s%5.5s %4.4s%1.1s%3.3s %1.1s%4d%1.1s   %8.2f%8.2f%8.2f%6.6s%6.6s          %2.2s%2.2s'

	#check a few things that might be missing, and add dummy values if needed
	checks = { 'occupancy':'%6.2f', 'tempFactor':'%6.2f', 'element':'%2.2s', 'charge':'%2.2s' }
	checked_values = {}
	for key in checks.keys():
		if key in atom:
			checked_values[key] = checks[key] % ( atom[key] )
		else:
			checked_values[ key ] = ' '

	string = PDB_atom_line_format % (
		atom['type'],
		str(atom['serial']),
		atom['name'],
		atom['altLoc'],
		atom['resName'],
		atom['chainID'],
		atom['resSeq'],
		atom['iCode'],
		atom['x'],
		atom['y'],
		atom['z'],
		checked_values['occupancy'],
		checked_values['tempFactor'],
		checked_values['element'],
		checked_values['charge']
	)
	
	
	return string

	
def UpdateFilters( filter_string, filters, keyval_sep = '=', val_sep = ',', resSeq_sep = ':' ):
	"""
	Update the information stored in a set of PDB atom filters. Provides progressive generation of sets of filters.

	Args:
	  filter_string (string): new filter information to add
	  filters (dictionary): maps a PDB element key onto a list of acceptable values for that element
	  keyval_sep (string): separator character for key/values strings (default: '=', as in 'key=val')
	  val_sep (string): delimiter character for lists of values (default: ',', as in 'val1,val2,val3')
	  resSeq_Sep (string): delimiter for INCLUSIVE integer ranges specified as a resSeq filter (default is ':', as in 'resSeq=1:10')

	Returns:
	  An updated filter dictionary.

	Notes:
		Filters are simple dictionaries: key => [ value1, value2, ... ]
		where the key is one of the entries found in the PDB atoms,
		and the value list is a set of acceptable values for that key.
	"""
	tokens = filter_string.split()
	for token in tokens:
		#
		# Expecting key=value1,value2,...
		#
		subtoks = token.split( keyval_sep )
		if len(subtoks) < 2: continue
		key = subtoks[0]
		values = subtoks[1].split( val_sep )

		#
		# Key not defined in the PDB info, so ignore this.
		#
		if key not in PDB_atom_line_info: continue

		#
		# If key not defined in filters previously, add empty list to start
		#
		if key not in filters: filters[key] = []

		#
		# Special case: resSeq values can be defined as i:j
		#
		if key == 'resSeq':
			for v in values:
				rtoks = v.split( resSeq_sep )
				i = int( rtoks[0] )
				j = i

				if len(rtoks) > 1:
					j = int( rtoks[1] )

				if j < i:
					print( 'AddFilters: Bad range %d => %d, from "%s"' % (i,j,v) )
					sys.exit( -1 )

				for resSeq in range( i, j+1 ):
					filters[key].append( resSeq )
		#
		# Otherwise, use the converter method defined in the PDBTools module for this key
		#
		else:
			converter = PDB_atom_line_info[key]['converter']
			for v in values:
				filters[key].append( converter(v) )

	return filters


def FilterAtoms( in_atoms, filters ):
	"""
	Return a list of PDB atom dictionaries that pass the specified filers.

	Args:
	  in_atoms (list of PDB atom dictionaries): input atoms to filter
	  filters (dictionary): map filter element keys to a list of acceptable value for that key

	Returns:
	  A list of PDB atom dictionaries that pass the specified filters.
	"""
	out_atoms = []
	for a in in_atoms:
		#
		# See whether "a" passes the filters
		#
		add_atom = True
		for key in filters:
			if (key in a) and (a[key] not in filters[key]):
				add_atom = False
				break
		#
		# Add "a" to "atoms", if appropriate
		#
		if add_atom == True:
			out_atoms.append( dict(a) )
	return out_atoms


def FilterAtomsByIndex( in_atoms, in_indices, filters ):
	"""
	Return a list of INDICES of the PDB atom dictionaries that pass the specified filters.

	Args:
	  in_atoms (list of PDB atom dictionaries): input atom data for filtering
	  in_indices (list of integer): indices into in_atoms list to apply filters to
	  filters (dictionary): map filter element keys to a list of acceptable value for that key

	Returns:
	  A list of INDICES into in_atoms that passed the specified filters.
	"""
    out_indices = []
    for index in in_indices:
        a = in_atoms[index]
        #
        # See whether "a" passes the filters
        #
        add_atom = True
        for key in filters:
            if (key in a) and (a[key] not in filters[key]):
                add_atom = False
                break
        #
        # Add "a" to "atoms", if appropriate
        #
        if add_atom == True:
            out_indices.append( index )
    return out_indices

#
#	Author: John Grime, The University of Chicago.
#

import sys
import xml.etree.ElementTree as ET
import PDB

#
# Convert Simone Mattei's geometrical HIV data from the XML output into PDB data.
# See e.g. http://science.sciencemag.org/content/354/6318/1434.long
#

def tripletise( string, delim ):
	"""
	Convert string of tokens into floating point [x,y,z] triplets.

	Args:
		string (string) : the text string to convert
		delim (string) : the delimiters to use in tokenization

	Returns:
		list of [x,y,z] entries, where elements are floats

	Example:
	
		>>> tripletise( '1 2 3, 4 5 6', ',' )
		[ [1.0,2.0,3.0], [4.0,5.0,6.0] ]
	"""
	tokens = string.strip().split(delim)
	stuff = []
	for token in tokens:
		stuff.append( [ float(t) for t in token.split() ] )
	return stuff


def print_PDB( f, atoms, bonds, unit_length = None ):
	"""
	Print geometrical information to PDB file format for visualization.

	Args:
		f (file) : output destination
		atoms (list of PDB-style atom) : vertices
		bonds (list of integer pairs) : vertex connections
		unit_length (integer) : number of consecutive vertices in a notional geometrical unit

	Returns:
		Nothing
	"""
	conect_format = 'CONECT%5d%5d'
	
	if unit_length == None:
		unit_length = len(atoms)
	
	counter = 0
	for a in atoms:
		if (counter>0) and (counter%unit_length==0):
			print( 'TER   ', file=f )
		line = PDB.MakePDBAtomLine( a )
		print( line, file=f )
		counter += 1

	print( 'TER   ', file=f )

	for b in bonds:
		line = conect_format % ( b[0], b[1] )
		print( line, file=f )


def parse_arguments( args ):
	"""
	Tokenizes a list of strings in the format `key=value` and generates a dictionary.

	Args:
		args (list of strings) : input strings
		atoms (list of PDB-style atom) : vertices
		bonds (list of integer pairs) : vertex connections
		unit_length (integer) : number of consecutive vertices in a notional geometrical unit

	Returns:
		Dictionary of parameters, mapping string => string
	"""
	results = {}
	for arg in args:
		tokens = arg.strip().split('=')
		key, value = tokens[0], None
		if len( tokens ) > 1: value = tokens[1]
		results[key] = value
	return results


def print_usage( prog ):
	print( '' )
	print( 'Usage: %s input=file.x3d [edge_name=X] [edge_resSeq=X] [scale=X] [output_prefix=X]' % (prog) )
	print( '' )
	print( 'Where:' )
	print( '' )
	print( '- input : x3d file generated by Chimera' )
	print( '- edge_name : OPTIONAL atom name for edge points (default: CA)' )
	print( '- edge_resSeq : OPTIONAL residue number for edge points (default: 1)' )
	print( '- scale : OPTIONAL scaling of coords (default: 1.0)' )
	print( '- output_prefix : OPTIONAL output prefix for files (dfault: output)' )
	print( '' )
	sys.exit( -1 )

#
# Check command line args
#
output_prefix = 'blah'

if len(sys.argv) < 2:
	print_usage( sys.argv[0] );

if len(sys.argv) > 2:
	output_prefix = sys.argv[2]

args = parse_arguments( sys.argv[1:] )

if 'input' not in args:
	print_usage( sys.argv[0] )

input_x3d = args['input']
edge_name, edge_resSeq, scale, output_prefix = 'CA', 1, 1.0, 'output'

if 'edge_name'     in args: edge_name     = args['edge_name']
if 'edge_resSeq'   in args: edge_resSeq   = int( args['edge_resSeq'] )
if 'scale'         in args: scale         = float( args['scale'] )
if 'output_prefix' in args: output_prefix = args['output_prefix']


#
# Parse speified input
#
try:
	tree = ET.parse( input_x3d )
	root  = tree.getroot()
except:
	print( 'Unable to parse input x3d file "%s"!'%(input_x3d) )
	sys.exit( -1 )

#
# Strip out the color and coordinate data from the file.
# We're assuming sequential matching of those data.
#
color_vec, coord_vec = [], []
for node in tree.iter():
	if node.tag == 'Color':
		data_string = node.attrib['color']
		color_vec.append( tripletise(data_string,',') )
	if node.tag == 'Coordinate':
		data_string = node.attrib['point']
		coord_vec.append( tripletise(data_string,',') )

#
# Remove the trailing (0,0,0) triplet that seems to be
# present for every element in the color and coord vecs
#
color_vec = [ x[0:len(x)-1] for x in color_vec ]
coord_vec = [ x[0:len(x)-1] for x in coord_vec ]

#
# Scale coordinates
#
for i in range( 0, len(coord_vec) ):
	capsomer_coords = coord_vec[i]
	for j in range( 0, len(capsomer_coords) ):
		x,y,z = coord_vec[i][j]
		coord_vec[i][j] = [ x*scale, y*scale, z*scale ]

#
# Print a PDB file, connecting entries as appropriate
#
capsomer_atom_names = { 5:'P',   6:'H' }
capsomer_res_names  = { 5:'PEN', 6:'HEX' }

atoms, bonds, serials, resSeqs = {}, {}, {}, {} # key = length of capsomer structure
for key in capsomer_atom_names:
	atoms[key] = []
	bonds[key] = []
	serials[key] = 1
	resSeqs[key] = 1

atom = PDB.MakeEmptyPDBAtom()
for capsomer_index in range(0,len(color_vec)):
	colors = color_vec[capsomer_index]
	coords = coord_vec[capsomer_index]
	capsomer_length = len(coords)
	
	if capsomer_length not in capsomer_atom_names:
		print( 'Unknown capsomer length %d!'%(capsomer_length) )
		sys.exit(-1)

	serial = serials[capsomer_length]
	resSeq = resSeqs[capsomer_length]

	#
	# Generate an ATOM entry for every vertex in the
	# capsomer, with sequential bonds for edges.
	#
	for j in range(0,capsomer_length):
		x,y,z = coords[j]
		atom['name'] = capsomer_atom_names[capsomer_length]
		atom['serial'] = serial
		atom['resName'] = capsomer_res_names[capsomer_length]
		atom['resSeq'] = resSeq
		atom['x'], atom['y'], atom['z'] = x, y, z

		atoms[capsomer_length].append( dict(atom) )
		if( j > 0 ): bonds[capsomer_length].append( [serial-1,serial] )

		serial += 1

	resSeq += 1

	#
	# Add trailing bond to connect last to first.
	# Note ascending order of serials!
	#
	bonds[capsomer_length].append( [serial-capsomer_length,serial-1] )

	serials[capsomer_length] = serial
	resSeqs[capsomer_length] = resSeq

#
# Save separate capsomer files
#
for key in atoms:
	f = open( '%s.%d.pdb'%(output_prefix,key), 'w' )
	print_PDB( f, atoms[key], bonds[key], key )
	f.close()

#
# Generate midpoints between vertices: useful for superpositions of tops of H9
#
serial = 1
new_stuff = {}
for key in bonds:
	new_stuff[key] = []
	for b in bonds[key]:
		ai = atoms[key][ b[0]-1 ]
		xi, yi, zi = [ ai[axis] for axis in ['x','y','z'] ]

		aj = atoms[key][ b[1]-1 ]
		xj, yj, zj = [ aj[axis] for axis in ['x','y','z'] ]
		
		x,y,z = (xi+xj)/2, (yi+yj)/2, (zi+zj)/2

		a = dict(ai)
		a['name'] = edge_name
		a['resSeq'] = edge_resSeq
		a['serial'] = serial
		a['x'], a['y'], a['z'] = x, y, z
		new_stuff[key].append( dict(a) )

		serial += 1

for key in new_stuff:
	f = open( '%s.%d.midpoints.pdb'%(output_prefix,key), 'w' )
#	print_PDB( f, new_stuff[key], [], key )
	print_PDB( f, new_stuff[key], [], 1 )
	f.close()


sys.exit()
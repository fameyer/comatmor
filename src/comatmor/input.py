# Routines to control data input to python/pymor
# by Falk Meyer, 09.02.2015
# Import methods to read sparse .mat files from matlab
from scipy import sparse, io
# column sparse matrix
from scipy.sparse import csc_matrix
# Import parameters (FIND NICER WAYS TO DO SO)
import parameter
# MAYBE MAKE A WHOLE INPUT CLASS INCLUDING A DECOMPOSER?

def getMatrix(type='disc', row=None, col=None, data=None):
	"""
	DOC ME
	"""
	# Read sparse matrices from harddisc into dictionary
	if type == 'disc':
		dic = {}
		for key in parameter.matfile:
		print 'Reading '+key+'...'
		dic[key] = io.loadmat(parameter.matfile[key][0])[parameter.matfile[key][1]]
		# DECOMPOSER here?
		# dic stores name and state of a matrix
		return dic
	# Generate sparse column matrix
	if type == 'direct':
		dic = {}
		print 'Mmh?'
		#print data
		A = csc_matrix((data[1:],(row[1:],col[1:])),shape=(len(data[1:]),len(data[1:])))
		return A
	else:
		pass

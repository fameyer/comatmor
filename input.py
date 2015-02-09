# Routines to control data input to python/pymor
# by Falk Meyer, 09.02.2015

# Import methods to read sparse .mat files from matlab
from scipy import sparse, io

# Import parameters (FIND NICER WAYS TO DO SO)
import parameter

# MAYBE MAKE A WHOLE INPUT CLASS INCLUDING A DECOMPOSER?
def getMatrix(type='disc', matrixDict = None):
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

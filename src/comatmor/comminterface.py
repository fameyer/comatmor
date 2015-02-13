# Interface class to control data exchange between outer source (e.g. matlab)
# and pymor.
# by Falk Meyer, 10.02.2015

import numpy as np

# Import methods to read sparse .mat files from matlab
from scipy import sparse, io
# column sparse matrix
from scipy.sparse import csc_matrix

# pymor includes
from pymor.parameters.functionals import GenericParameterFunctional
from pymor.parameters.base import Parameter

from pymor.la.numpyvectorarray import NumpyVectorArray

import parameter

class comminterface(object):
	"""
	DOC ME
	"""
	def __init__(self, name="Interface for comatmor communications", type = 'direct'):
		"""
		DOC ME
		"""
		self._name = name
		self._matDict = {}
		self._type = type

	def pushMat(self, name=None, row=None, col=None, data=None, paramName=None, paramValue=None, paramSpace=None):
		"""
		DOC ME
		"""
		# Direct argument call
		if self._type == 'direct':
			
			# check for correct data types
			assert isinstance(row, np.ndarray)
			assert isinstance(col, np.ndarray)
			assert isinstance(data, np.ndarray)

			# get dimension of system
			dim = row.max()+1 if row.max() > col.max() else col.max()+1

			# transform to scipymatrix
			matrix = csc_matrix((data,(row,col)),shape=(dim,dim))

			# create parameterfunctional
			paramType = Parameter({ paramName: np.array(paramValue)})
			paramFunc = GenericParameterFunctional(lambda mu: mu[paramName], parameter_type = paramType)

			# decompose matrix comes together with above line later on
			self.decompose(matrix)

			# save matrix to dic
			self._matDict[name]=(matrix, paramFunc, paramSpace, paramType)
		
		# Calls by harddisc access
		if self._type == 'disc':
	                for key in parameter.matfile:
				# Obtain information from the parameter file
				paramName = parameter.matfile[key][1]
				paramValue = parameter.matfile[key][2]
				paramType = Parameter({ paramName: np.array(paramValue)})
				paramSpace = parameter.matfile[key][3]
				paramFunc = GenericParameterFunctional(lambda mu: mu[paramName], parameter_type = paramType)
                        	print 'Reading '+key+'...'

	                        self._matDict[key] = (io.loadmat(parameter.matfile[key][0])[key],paramFunc,paramSpace,paramType)

	def getMat(self):
		"""
		DOC ME
		"""
		return self._matDict
	
	def decompose(self,mat):
		"""
		DOC ME
		"""
		pass

	def pushRhs(self):
		"""
		DOC ME
		"""
		# do that better without for-loop
		for key in parameter.rhsfile:
			print 'Reading rhs...'
			return io.loadmat(parameter.rhsfile[key])[key]

	def writeSolutions(self, u):
		"""
		Write given u to disc
		"""		

		assert self._type == 'disc'
		assert isinstance(u,dict) 
	
		io.savemat('RBsolutions',u)

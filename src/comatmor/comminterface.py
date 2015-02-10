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

	def pushMat(self, name, row, col, data, paramName, paramValue, paramSpace, training_set):
		"""
		DOC ME
		"""
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
			paramType = Parameter({ paramName: np.array([paramValue])})
			paramFunc = GenericParameterFunctional(lambda mu: mu[paramName], parameter_type = paramType)

			# decompose matrix comes together with above line later on
			self.decompose(matrix)

			# save matrix to dic
			self._matDict[name]=(matrix,paramFunc,paramSpace,training_set, paramType)

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

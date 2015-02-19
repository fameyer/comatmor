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
from pymor.parameters.base import ParameterType

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
		self._matDict = ({}, ParameterType({}), {})
		self._type = type

	def pushMat(self, name=None, row=None, col=None, data=None, paramName=None, paramShape=None, paramRange=None):
		"""
		DOC ME
		"""
		# Direct argument call
		if self._type == 'direct':
			
			# check for correct data types
			assert isinstance(row, np.ndarray)
			assert isinstance(col, np.ndarray)
			assert isinstance(data, np.ndarray)

			# check for all input arguments
			assert not isinstance(paramName,None)
			assert not isinstance(paramShape,None)
			assert not isinstance(paramSpace,None)

			# get dimension of system
			dim = row.max()+1 if row.max() > col.max() else col.max()+1

			# transform to scipymatrix
			matrix = csc_matrix((data,(row,col)),shape=(dim,dim))

			# create parameterfunctional
			paramType = ParameterType({ paramName: paramShape})
			paramFunc = GenericParameterFunctional(lambda mu: mu[paramName], parameter_type = paramType)
			# save matrix to dic
			self._matDict[0][name]=(matrix, paramFunc)
				
			# Update of parametertypes and ranges
			self._matDict[1][paramName]=paramType[paramName]
			self._matDict[2][paramName]=paramRange

		# Call by harddisc access
		if self._type == 'disc':
		
	                for key in parameter.matfile:

				# Obtain information from the parameter file
				paramName = parameter.matfile[key][1]
				paramShape = parameter.matfile[key][2]
				paramRange = parameter.matfile[key][3]
				
				# assert linear, scalar parameterdependence
				assert len(paramName)== 1
			
				# If just one scalar parameter given, we have to correct due to sparse scipy matrices
				paramType = ParameterType({ paramName[0]: 0})
					
				paramFunc = GenericParameterFunctional(lambda mu: mu[paramName[0]], parameter_type = paramType)		
			
                        	print 'Reading '+key+'...'

	                        self._matDict[0][key] = (io.loadmat(parameter.matfile[key][0])[key],paramFunc)
				# add parametertypes and parameterRanges
				self._matDict[1][paramName[0]] = paramType[paramName[0]]
				self._matDict[2][paramName[0]] = paramRange[0]

	def getMat(self):
		"""
		DOC ME
		"""
		return self._matDict
	
	def pushRhs(self):
		"""
		DOC ME
		"""
		# do that better without for-loop
		for key in parameter.rhsfile:
			print 'Reading rhs...'
			return io.loadmat(parameter.rhsfile[key])[key]

	def writeSolutions(self, u, file=None):
		"""
		Write given u to disc
		"""		

		assert self._type == 'disc'
		assert isinstance(u,dict) 
	
		# check if user-given filename is available, otherwise use default
		if file == None:
			io.savemat('RBsolutions',u)
		else:
			io.savemat(file,u)

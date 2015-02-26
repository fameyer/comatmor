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
from pymor.parameters.base import Parameter

from pymor.la.numpyvectorarray import NumpyVectorArray

from pymor.operators.numpy import NumpyMatrixOperator
from pymor.operators.constructions import LincombOperator

# make that NICER
#from comatmor.elliptic import parameterHeateq as parameter
from comatmor.heat import parameterHeateq as parameter


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
					
	                 	print 'Reading '+key+'...'
				# add parametertypes and parameterRanges

				self._matDict[1][paramName[0]] = paramType[paramName[0]]
				self._matDict[2][paramName[0]] = paramRange[0]

			for key in parameter.matfile:
				paramName = parameter.matfile[key][1]
				print 'Reading parameter: '+paramName[0]
				name = paramName[0]
				paramFunc = GenericParameterFunctional(lambda mu, name=name: mu[name], parameter_type = self._matDict[1])
				self._matDict[0][key] = (io.loadmat(parameter.matfile[key][0], mat_dtype=True)[key],paramFunc) 

        def assembleOperators(self):
                """
                Assemble operators to be prepared for RB calculations 
                """
                matDict = self._matDict

                # Create lincomboperator for all involved matrices
                stiffOps, rhsOps, stiffPops, rhsPops = [], [], [], []
                
                # get stiffness matrices and rhs matrices
                for key in parameter.stiffNames:
                #for key in matDict[0]:
                        stiffOps.append(NumpyMatrixOperator(matDict[0][key][0]))
                        stiffPops.append(matDict[0][key][1])
                for key in parameter.rhsNames:
                        rhsOps.append(NumpyMatrixOperator(matDict[0][key][0].T))
                        rhsPops.append(matDict[0][key][1])

                stiffOp = LincombOperator(stiffOps, coefficients=stiffPops)
                rhsOp = LincombOperator(rhsOps,coefficients=rhsPops)
                return stiffOp, rhsOp

	def getMat(self):
		"""
		DOC ME
		"""
		return self._matDict
	
	def pushRhs(self):
		"""
		DOC ME
		"""
		if self._type == 'disc':
			# do that better without for-loop
			for key in parameter.rhsfile:
				print 'Reading rhs...'
				return io.loadmat(parameter.rhsfile[key])[key]
		else:
			pass
	
	def readMass(self):
		"""
		Read mass/damping matrix
		"""
		if self._type == 'disc':
			for key in parameter.massfile:
				print 'Reading mass...'
				return io.loadmat(parameter.massfile[key])[key]
		else:
			pass

	def readU0(self):
		"""
		Read initial solution for time-dependent problems
		"""
		if self._type == 'disc':
			for key in parameter.u0file:
				print 'Reading initial solution...'
				return io.loadmat(parameter.u0file[key])[key]
		else:
			pass

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

	def getTrainingSet(self):
		"""
		Read training-set from disc
		"""
		if self._type == 'disc':
			for key in parameter.trainingSetfile:
				print 'Obtain training_set...'
				# transform to correct format
				training_set = io.loadmat(parameter.trainingSetfile[key], mat_dtype=True)[key]
				return [tuple(training_set[i]) for i in range(0,len(training_set))]
		else:
			pass

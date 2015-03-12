# Interface class to control data exchange between outer source (e.g. matlab)
# and pymor.
# by Falk Meyer, 10.02.2015

import numpy as np
import csv

# Import methods to read sparse .mat files from matlab
from scipy import sparse, io
# column sparse matrix
from scipy.sparse import csc_matrix

# pymor includes
from pymor.parameters.functionals import GenericParameterFunctional
from pymor.parameters.functionals import ExpressionParameterFunctional
from pymor.parameters.base import ParameterType
from pymor.parameters.base import Parameter

from pymor.la.numpyvectorarray import NumpyVectorArray

from pymor.operators.numpy import NumpyMatrixOperator
from pymor.operators.constructions import LincombOperator, VectorFunctional

# make that NICER
#from comatmor.elliptic import parameterHeateq as parameter
from comatmor.IRT import parameterIRT as parameter


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
				# For just one parameter so far
				paramType = ParameterType({ paramName[0]: 0})
					
	                 	print 'Reading '+key+'...'
				# add parametertypes and parameterRanges

				self._matDict[1][paramName[0]] = paramType[paramName[0]]
				self._matDict[2][paramName[0]] = paramRange[0]

			for key in parameter.matfile:
				paramName = parameter.matfile[key][1]
				print 'Reading parameter: '+paramName[0]
				name = paramName[0]
				#paramFunc = GenericParameterFunctional(lambda mu, name=name: mu[name], parameter_type = self._matDict[1])
				paramFunc = ExpressionParameterFunctional(name, parameter_type = self._matDict[1])

				self._matDict[0][key] = (io.loadmat(parameter.matfile[key][0], mat_dtype=True)[key],paramFunc) 

        def assembleOperators(self):
                """
                Assemble operators to be prepared for RB calculations 
                """
                matDict = self._matDict

                # Create lincomboperator for all involved matrices
                stiffOps, rhsOps, massOps, stiffPops, rhsPops, massPops = [], [], [], [], [], []
                # get stiffness matrices and rhs matrices
		for key in parameter.stiffNames:
                        stiffOps.append(NumpyMatrixOperator(matDict[0][key][0]))
                        stiffPops.append(matDict[0][key][1])
                for key in parameter.rhsNames:
                        rhsOps.append(NumpyMatrixOperator(matDict[0][key][0].T))
                        rhsPops.append(matDict[0][key][1])
		for key in parameter.massname:
			massOps.append(NumpyMatrixOperator(matDict[0][key][0]))
			massPops.append(matDict[0][key][1])

		# Correct for Dirichlet data in right-hand side
		# Load given COMSOL-array to find dirichlet indices
                L = self.getDirichletIndex()
		lenl = len(L)

		# read values of dirichlet function
   		values = self.getDirichletValues()			
      
		# Add _t to parametertype
		localtype = matDict[1].copy()
		localtype['_t'] = 0
		
		diriRhs = NumpyMatrixOperator(np.zeros((1,rhsOps[0]._matrix.shape[1])))
	        for j in range(0,lenl):
	        	diriRhs._matrix[0,L[j]-1] = 1.0

		rhsOps.append(diriRhs)
		rhsPops.append(ExpressionParameterFunctional('k*_t', parameter_type=localtype))

		# Shift RHS for dirichlet
                U_d = NumpyVectorArray(self.readU0().T.copy())
                for i in range(len(U_d.data[0])):
	                if i+1 in L:
        		        U_d.data[0][i] = 1.0
	                else:
        		        U_d.data[0][i] = 0.0 
	
		# Subtract appropr. mass matrices from RHS	
                for key in parameter.massname:
        	       rhsOps.append(VectorFunctional(NumpyMatrixOperator(matDict[0][key][0]).apply(U_d)*(-1)))
		       paramName = parameter.matfile[key][1]
                       name = paramName[0]
		       rhsPops.append(ExpressionParameterFunctional(name+'*_t', parameter_type=localtype))

		# Subtract appropr. stiffness matrices from RHS
                for key in parameter.stiffNames:
                       rhsOps.append(VectorFunctional(NumpyMatrixOperator(matDict[0][key][0]).apply(U_d)*(-1)))
		       paramName = parameter.matfile[key][1]
                       name = paramName[0]
                       rhsPops.append(ExpressionParameterFunctional(name+'*_t', parameter_type=localtype))
		stiffOp = LincombOperator(stiffOps, coefficients=stiffPops)
                rhsOp = LincombOperator(rhsOps,coefficients=rhsPops)
		massOp = LincombOperator(massOps, coefficients=massPops)
                return stiffOp, rhsOp, massOp

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

	def getParameterSet(self):
		"""
		Read parameter-set from disc
		"""
		if self._type == 'disc':
			for key in parameter.parameterSetfile:
				print 'Obtain parameter_set...'
				# transform to correct format
				parameter_set = io.loadmat(parameter.parameterSetfile[key], mat_dtype=True)[key]
				return [tuple(parameter_set[i]) for i in range(0,len(parameter_set))]
		else:
			pass

	def getDirichletIndex(self, file = 'dirichletIndex.mat'):
		"""
		Read indices of dirichlet values from matlab saved file
		"""
		return (io.loadmat(file)['index'])[0]

	def getDirichletValues(self, file='dirichlet.csv'):
		"""
		Read time-dependent values from comsol saved csv file
		(First column time, second column value)
		"""
		f = open(file)
		values = []
		for row in csv.reader(f,delimiter=','):
			values.append(row)
		return values

	def getSignature(self, num_samples, steps, T):
		"""
		Construct signature for given RB object
		"""
		print 'Generating signature...' 
		# create signature
		sig = 'IRT'
		dim = 0
		for key in self._matDict[0]:
			if dim == 0:
				dim = self._matDict[0][key][0].shape[0]
			sig = sig+'_'+key
		sig = sig+'_'+str(num_samples)
		sig = sig+'_'+str(steps)
		sig = sig+'_'+str(T)
		sig = sig+'_'+str(dim)
		return sig  

	def saveSignature(self, file, signature):
		"""
		Save signature for underlying object
		"""
		print 'Writing new signature to file...'
		e = open(file,'a')
		e.write(signature+'\n')
		e.close()

	def checkSignature(self, file, signature):
		"""
	 	Check if given signature has already been saved
		"""
		try:
			print 'Checking signature...'
			e = open(file,'r')
			contents = e.readlines()
			# check if signature already given
			for i in contents:
				if i == signature+'\n':
					e.close()
					return True
			# No signature found
			e.close()
			return False
		except:
			print 'No signature file found...'
			return False

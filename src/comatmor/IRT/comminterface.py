#    Copyright (C) 2015 Falk Meyer
#
#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License along
#    with this program; if not, write to the Free Software Foundation, Inc.,
#    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

#===========================================================================

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

#from pymor.la.numpyvectorarray import NumpyVectorArray
from pymor.vectorarrays.numpy import NumpyVectorArray
from pymor.operators.numpy import NumpyMatrixOperator
from pymor.operators.constructions import LincombOperator, VectorFunctional

# make that NICER
#from comatmor.elliptic import parameterHeateq as parameter
from comatmor.IRT import parameterIRT as parameter


class comminterface(object):
	"""
	Interface class providing methods for matrix/vector input	
	"""
	def __init__(self, name="Interface for comatmor communications"):
		"""
		Constructor
		"""
		self._name = name
		self._matDict = ({}, ParameterType({}), {})

	def pushMat(self):
		"""
		Insert matrices from harddisc
		"""
		# Call by harddisc access
	
	
		for key in parameter.matfile:

			# Obtain information from the parameter file
			paramName = parameter.matfile[key][1]
			paramRange = parameter.matfile[key][2]
			
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
		for key in parameter.massNames:
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
                for key in parameter.massNames:
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
		Return matrix dictionary		
		"""
		return self._matDict
	
	def pushRhs(self):
		"""
		Insert Right-Hand side vector from harddisk
		"""
		for key in parameter.rhsfile:
			print 'Reading rhs...'
			return io.loadmat(parameter.rhsfile[key])[key]

	def readU0(self):
		"""
		Read initial solution for time-dependent problems
		"""
		for key in parameter.u0file:
			print 'Reading initial solution...'
			return io.loadmat(parameter.u0file[key])[key]

	def writeSolutions(self, u, file=None):
		"""
		Write given u to disc
		"""		
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
		for key in parameter.parameterSetfile:
			print 'Obtain parameter_set...'
			# transform to correct format
			parameter_set = io.loadmat(parameter.parameterSetfile[key], mat_dtype=True)[key]
			return [tuple(parameter_set[i]) for i in range(0,len(parameter_set))]
		
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

	def getSignature(self, num_samples, steps, T, max_extensions):
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
		sig = sig+'_'+str(dim)		
		sig = sig+'_'+str(num_samples)
		sig = sig+'_'+str(steps)
		sig = sig+'_'+str(T)
		sig = sig+'_'+str(max_extensions)
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

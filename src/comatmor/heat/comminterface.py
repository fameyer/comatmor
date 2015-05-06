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
		Returns inserted matrices
		"""
		return self._matDict
	
	def pushRhs(self):
		"""
		Inserts right-hand side
		"""
		for key in parameter.rhsfile:
			print 'Reading rhs...'
			return io.loadmat(parameter.rhsfile[key])[key]

	def readMass(self):
		"""
		Read mass/damping \matrix
		"""
		for key in parameter.massfile:
			print 'Reading mass...'
			return io.loadmat(parameter.massfile[key])[key]

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

	def getSignature(self, num_samples, steps, T):
		"""
		Construct signature for given RB object
		"""
		print 'Generating signature...' 
		# create signature
		sig = 'heatequation'
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

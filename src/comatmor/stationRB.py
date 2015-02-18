# Class to apply pymor RB on a stationary PDE
# by Falk Meyer, 10.02.2015

import numpy as np

from functools import partial

# for debugging
#import pdb

# pymor includes
from pymor.operators.constructions import LincombOperator
from pymor.operators.numpy import NumpyMatrixOperator

from pymor.la.numpyvectorarray import NumpyVectorArray

from pymor.discretizations.basic import StationaryDiscretization

from pymor.reductors.linear import reduce_stationary_affine_linear

from pymor.parameters.spaces import CubicParameterSpace

from pymor.algorithms.greedy import greedy
from pymor.algorithms.basisextension import trivial_basis_extension

# local imports
from comminterface import comminterface as CI

class stationRB(object):
	"""
	DOC ME
	"""
	def __init__(self, name = 'RB on stationary problem', inputmethod = 'direct'):
		"""
		DOC ME
		"""
		self._name = name
		self._type = inputmethod
		# communication interface
		self._CI = CI(type = self._type)
		# mathematical objects
		self._rhs = None
		self._matDict = None
		# RB-basis
		self._rb = None
		# RB discretization
		self._rd = None
		# RB reconstructor
		self._rc = None

	def addMatrix(self, name=None, row=None, col=None, data=None, paramName=None, paramValue=None, paramSpace=None):
		"""
		paramSpace = [1,2] for intervall so to speak, even discrete or continuous
		"""
		# just call comminterface
		self._CI.pushMat(name, row, col, data, paramName, paramValue, paramSpace)
		self._matDict = self._CI.getMat()

	def getMatrix(self):
		"""
		Namen umbenennen vermoege besseren Wissens!
		"""
		return self._CI.getMat()

	def addRhs(self, rhs=None):
		"""
		DOC ME
		"""
		# also assert right dimension?
		if rhs != None:
			assert isinstance(rhs, np.ndarray)
			self._rhs = NumpyMatrixOperator(rhs)
		else:
			# try to get rhs from matDictionary(has to be transposed!)
			self._rhs = NumpyMatrixOperator(self._CI.pushRhs().T)

	def getRhs(self):
		"""
		DOC ME
		"""
		return self._rhs

	def constructRB(self):
		"""
		TO DO: add possibility to change basis extension for instance
		"""
		print 'Constructing reduced basis...' 

		matDict = self._matDict

		# Create lincomboperator for all involved matrices
		ops, pops, paramSpaces, paramTypes = [], [], [], []
		# maybe improve method below - put it somewhere else prob. !!!
		for key in matDict:
			ops.append(NumpyMatrixOperator(matDict[key][0]))
			pops.append(matDict[key][1])
			paramSpaces.append(matDict[key][2])
			paramTypes.append(matDict[key][3])

		op = LincombOperator(ops, coefficients=pops)
		print op	
		print self._rhs
		# create discretization
		dis = StationaryDiscretization(operator=op, rhs=self._rhs)
	
		# create parameterSpace TO DO IN GENERAL - NOW JUST FOR ONE PARAM
		print paramSpaces[0][0]
		print paramSpaces[0][1]
		paramSpace = CubicParameterSpace(paramTypes[0], minimum=paramSpaces[0][0], maximum=paramSpaces[0][1])
		print paramSpace	
		# create reductor
	 	reductor = partial(reduce_stationary_affine_linear, error_product = None)

		# greedy search to construct RB 		
		self._rb = greedy(dis, reductor, paramSpace.sample_uniformly(10), use_estimator=True, extension_algorithm=trivial_basis_extension, target_error=1e-10, max_extensions = 10) 

		# get the reduced discretization and the reconstructor
		self._rd, self._rc = self._rb['reduced_discretization'], self._rb['reconstructor']

	def compute(self, training_set=None, error=False):
		"""
		Think about error estimators and parameters etc etc
		Compute rb solutions for training set - with errors?
		"""
		# assert right set structure
		assert isinstance(training_set,np.ndarray) or isinstance(training_set, list)
		if self._type == 'direct': 
			# just supports return of one solution so far!!! Due to matlab restrictions
			for mu in training_set:
				u = self._rd.solve(mu)
				ur = self._rc.reconstruct(u)
		
				print ur
				return ur	
		if self._type == 'disc':
			
			solutions = {}
			i = 0
			# save solutions for all parameters	
			for mu in training_set:
				i=i+1
				u = self._rd.solve(mu)
				# Use data function to transform NumpyVectorArray to standard NumpyArray
				# Have to define valid matlab variable names
				# 'mu'+str(int(mu*100))
				solutions['mu'+str(i)]=(self._rc.reconstruct(u)).data
				
			# save solutions to disk
			self._CI.writeSolutions(solutions)
		
		
	def getRB(self):
		"""
		DOC ME
		"""
		return self._rb
	
	def saveObj(self):
		"""
		provide opportunity to save object
		"""
		pass

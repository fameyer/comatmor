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
		self._CI = CI()
		# mathematical objects
		self._rhs = None
		self._matDict = None
		# RB-basis
		self._rb = None
		# RB discretization
		self._rd = None
		# RB reconstructor
		self._rc = None

	def addMatrix(self, name, row, col, data, paramName=None, paramValue=None, paramSpace=None, training_set=None):
		"""
		paramSpace = [1,2] for intervall so to speak, even discrete or continuous
		"""
		# just call comminterface
		self._CI.pushMat(name, row, col, data, paramName, paramValue, paramSpace, training_set)
		self._matDict = self._CI.getMat()
	def getMatrix(self):
		"""
		Namen umbenennen vermoege besseren Wissens!
		"""
		return self._CI.getMat()

	def addRhs(self, rhs):
		"""
		DOC ME
		"""
		# also assert right dimension?
		assert isinstance(rhs, np.ndarray)
		self._rhs = NumpyMatrixOperator(rhs)

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
		ops, pops, paramSpaces, trainingSets, paramTypes = [], [], [], [], []
		# maybe improve method below
		for key in matDict:
			ops.append(NumpyMatrixOperator(matDict[key][0]))
			pops.append(matDict[key][1])
			paramSpaces.append(matDict[key][2])
			trainingSets.append(matDict[key][3])
			paramTypes.append(matDict[key][4])

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

		# just supports return of one solution so far!!! Due to matlab restrictions
		for mu in training_set:
			u = self._rd.solve(mu)
			ur = self._rc.reconstruct(u)

		print ur
		return ur	

	def getRB(self):
		"""
		DOC ME
		"""
		return self._rb

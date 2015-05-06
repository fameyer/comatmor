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
# Class to apply pymor RB on the simplified IRT model
# by Falk Meyer,12.03.2015

import numpy as np
from scipy import sparse, io

from functools import partial

# pymor includes
from pymor.operators.constructions import LincombOperator
from pymor.operators.numpy import NumpyMatrixOperator

from pymor.la.numpyvectorarray import NumpyVectorArray

from pymor.discretizations.basic import InstationaryDiscretization

from pymor.reductors.basic import reduce_generic_rb

from pymor.parameters.spaces import CubicParameterSpace

from pymor.algorithms.greedy import greedy
from pymor.algorithms.basisextension import trivial_basis_extension
from pymor.algorithms.basisextension import gram_schmidt_basis_extension
from pymor.algorithms.basisextension import pod_basis_extension
#from pymor.algorithms.timestepping import ImplicitEulerTimeStepper

# for Output
from pymor.core import logger
logger.MAX_HIERACHY_LEVEL = 2
logger.set_log_levels({'pymor.algorithms': 'INFO',
                       'pymor.discretizations': 'INFO',
                       'pymor.la': 'INFO',
                       'pymor.reductors': 'INFO'})

from timestepping import ImplicitEulerTimeStepper

from pymor.core.pickle import dump, load

# local imports
#from ..comminterface import comminterface as CI
from comminterface import comminterface as CI

class IRTRB(object):
	"""
        This is the class realizing input of saved matrices and communcation with 
        pyMOR to solve the simplified IRT model.
	"""
	def __init__(self, name = 'RB on stationary problem'):
		"""
		Constructor
		"""
		self._name = name
		# communication interface
		self._CI = CI()
		# mathematical objects
		self._rhs = None
		self._matDict = None
		self._u0 = None
		self._mass = None
		# RB-basis
		self._rb = None
		# RB discretization
		self._rd = None
		# RB reconstructor
		self._rc = None
		# Desired endtime of time-stepping
		self._T = 0
		# Number of timesteps
		self._steps = 0	
		# for inputmethod = disc call all necessary get-functions
		self.addMatrix()
		self.getU0()
		# Bool to check whether saving is desired
		self._save = True
		# Given signature file
		self._signFile = 'sign.txt'
		# Get dirichlet indices
		self._dirichletIndex = self._CI.getDirichletIndex()
		# Get dirichlet values over time
		self._dirichletValues = self._CI.getDirichletValues()

	def addMatrix(self):
		"""
		Add matrices to given instationHeatRB object - direct or by disc access
		"""
		# just call comminterface
		self._CI.pushMat()
		self._matDict = self._CI.getMat()

	def getMatrix(self):
		"""
		Insert matrix from harddisc
		"""
		return self._CI.getMat()

	def getRhs(self):
		"""
		DOC ME
		"""
		return self._rhs

	def getU0(self):
		"""
		Read initial solution for time-stepping
		"""
		self._u0 = NumpyMatrixOperator(self._CI.readU0())
		
	def constructRB(self, num_samples = 10, T=0, steps=0, max_extensions = 30, target_error =1e-10):
		"""
		Construct reduced basis with end-time T and number of steps steps
		"""
		print('Constructing reduced basis...')
	
		# Save parameters for compute
		self._num_samples = num_samples
		self._T = T
		self._steps = steps

		# check if reduced basis was already constructed before, then load it before contructing a new one
		if self._save:
			signature = self._CI.getSignature(num_samples, steps, T, max_extensions)
			if self._CI.checkSignature(self._signFile,signature):
				with open(signature,'r') as f:
					self._rd, self._rc = load(f)
					print('Reduced basis and reconstructor already computed before, loading...')
					return

		# call assembleOperators and get right operators
		stiffOp, rhsOp, massOp = self._CI.assembleOperators()
		
		# get parametertypes and ranges
		paramTypes  = self._matDict[1]
		paramRanges = self._matDict[2]
		
		# create timestepper
		time_stepper = ImplicitEulerTimeStepper(steps)

		# create discretization with induced norm
		ones =tuple([1 for i in range(len(paramTypes))])
	
		# Get dirichlet indices
    	        L = self._dirichletIndex
		lenl = len(L)

		# Extract initial solution without pure Dirichlet solution
		U_d = NumpyVectorArray(self._u0._matrix.T.copy())
		for i in range(len(U_d.data[0])):
  			if i+1 in L:
                       		U_d.data[0][i] = self._u0._matrix[i]
               		else:
                       		U_d.data[0][i] = 0.0 
		UNull = NumpyMatrixOperator((self._u0._matrix.T - U_d.data[0]).T)

		# Get matrix to be used for error_norm in greedy algorithm
		Normmatrix = NumpyMatrixOperator(io.loadmat('KNorm.mat', mat_dtype=True)['K']) 

		dis = InstationaryDiscretization(operator=stiffOp, rhs=rhsOp, initial_data=UNull, T=T, time_stepper=time_stepper, mass=massOp, products={'l2': Normmatrix})

		# create parameterSpace
		paramSpace = CubicParameterSpace(parameter_type = paramTypes, ranges = paramRanges)
		print('Given parameterspace: '+str(paramSpace))

		# create reductor
		def reductor(discretization, rb, extends = None):
			return reduce_generic_rb(dis, rb, extends=extends)
		
		print('Do greedy search...')
		
		# greedy search to construct RB 		
		self._rb = greedy(dis, reductor, paramSpace.sample_uniformly(num_samples), use_estimator=False, extension_algorithm=pod_basis_extension, target_error=target_error, max_extensions = max_extensions, error_norm=lambda U: np.max(dis.l2_norm(U))) 
		# get the reduced discretization and the reconstructor
		self._rd, self._rc = self._rb['reduced_discretization'], self._rb['reconstructor']
		self._bas = self._rb['basis']

		print('Greedy search successfull! Reduced basis has dimension: '+str(len(self._rb['basis'])))
		# If saving desired, save the reduced basis and the reconstructor to the disc
		if self._save:
			print('Saving reduced basis and reconstructor...')
			self._CI.saveSignature(self._signFile, signature) 
			with open(signature, 'w') as f:
				dump((self._rd, self._rc),f)		

	def compute(self, parameter_set=None, error=False, file=None):
		"""
		Compute rb solutions for parameter_set
		"""
		# assert right set structure
		assert isinstance(parameter_set,np.ndarray) or isinstance(parameter_set, list) or parameter_set == None
		assert not self._rd is None and not self._rc is None	
		# Get parameter_set from disc
		assert parameter_set == None	
		print('Computing solutions for given parameter set in respect to reduced basis...')
		parameter_set = self._CI.getParameterSet()
		solutions = {}
		i = 0

		# Open files to get right dirichlet values
		values = self._dirichletValues

		L = self._dirichletIndex
		lenl = len(L)

		# save solutions for all parameters	
		for mu in parameter_set:
			i=i+1
			u = self._rd.solve(mu)
			# Use data function to transform NumpyVectorArray to standard NumpyArra             			# Have to define valid matlab variable names
			solutions['mu'+str(i)]=(self._rc.reconstruct(u)).data
			
			# Add default dirichlet values to the solution
			t = 0.0
			dt = float(self._T)/float(self._steps);
			# Imply Dirichlet - will be linearly interpolated when time values do not match
			for j in range(0,len(solutions['mu'+str(i)])):
				for k in range(1,len(values)):
					if t >= float(values[k-1][0]) and t <= float(values[k][0]):
						diriValue = float(values[k-1][1]) + (float(values[k][1]) - float(values[k-1][1]))/(float(values[k][0])-float(values[k-1][0]))*(t-float(values[k-1][0]))
						break

				for l in range(len(solutions['mu'+str(i)][j])):
					if l+1 in L:
						solutions['mu'+str(i)][j][l] = diriValue	
				t += dt

		# save solutions to disk
		self._CI.writeSolutions(solutions,file)
			
	def getRB(self):
		"""
		Returns computed reduced basis
		"""
		assert self._rb 
		return self._rb
	
	def __str__(self):
		"""
		Information about the object
		"""
		return 'Reduced basis method with greedy search basis generation method for the IRT-model'

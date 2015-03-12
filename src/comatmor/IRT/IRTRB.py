# Class to apply pymor RB on the simplified IRT model
# by Falk Meyer, 10.02.2015

import pickle
import time

import numpy as np
from scipy import sparse, io
import copy

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
from pymor.algorithms.timestepping import ImplicitEulerTimeStepper

from pymor.core.pickle import dump, load

# local imports
#from ..comminterface import comminterface as CI
from comminterface import comminterface as CI
from problems import deaffinize_discretization

class IRTRB(object):
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
		self._u0 = None
		self._mass = None
		# RB-basis
		self._rb = None
		# RB discretization
		self._rd = None
		# RB reconstructor
		self._rc = None
		# Default
		self._save = False

		# for inputmethod = disc call all necessary get-functions
		if self._type == 'disc':
			self.addMatrix()
			self.getU0()
			# Bool to check whether saving is desired
			self._save = True
			# Given signature file
			self._signFile = 'sign.txt'

	def addMatrix(self, name=None, row=None, col=None, data=None, paramName=None, paramShape=None, paramRange=None):
		"""
		Add matrices to given instationHeatRB object - direct or by disc access
		"""
		# just call comminterface
		self._CI.pushMat(name, row, col, data, paramName, paramShape, paramRange)
		self._matDict = self._CI.getMat()

	def getMatrix(self):
		"""
		Namen umbenennen vermoege besseren Wissens!
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
		if self._type == 'disc':
			self._u0 = NumpyMatrixOperator(self._CI.readU0())
		else:
			pass

	def constructRB(self, num_samples = 10, T=0, steps=0):
		"""
		Construct reduced basis with end-time T and number of steps n
		"""
		print 'Constructing reduced basis...' 
	
		# check if reduced basis was already constructed before, then load it before contructing a new one
		if self._save:
			signature = self._CI.getSignature(num_samples, steps, T)
			if self._CI.checkSignature(self._signFile,signature):
				with open(signature,'r') as f:
					self._rd, self._rc = load(f)
					print 'Reduced basis and reconstructor already computed before, loading...'
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

    	        L = io.loadmat('/home/310191226/pymorDir/comatmor/src/comatmor/IRT/dirichletIndex.mat')['index']
    		L = L[0]
		lenl = len(L)

		#Extract initial solution without Dirichlet solution
		U_d = NumpyVectorArray(self._u0._matrix.T.copy())
		for i in range(len(U_d.data[0])):
  			if i+1 in L:
                       		U_d.data[0][i] = 50.0
               		else:
                       		U_d.data[0][i] = 0.0 
		UNull = NumpyMatrixOperator((self._u0._matrix.T - U_d.data[0]).T)

		Normmatrix = NumpyMatrixOperator(io.loadmat('KcNorm.mat', mat_dtype=True)['Kc']) 

		dis = InstationaryDiscretization(operator=stiffOp, rhs=rhsOp, initial_data=UNull, T=T, time_stepper=time_stepper, mass=massOp, products={'l2': Normmatrix})#stiffOp.assemble(ones)})
		#L = np.linalg.cholesky(stiffOp.assemble(ones)._matrix.todense())
		#print('is positive definite!')

		#loese = dis.solve((1.0,1.0,1.0))
		#io.savemat('loesung',{'loese':loese.data})
		#exit()

		#R = dis.solve((1.0,40.0))
		#self._CI.writeSolutions({'R': R.data},file='Test11')
		#exit()
		# create parameterSpace
		paramSpace = CubicParameterSpace(parameter_type = paramTypes, ranges = paramRanges)
		print 'Given parameterspace: '+str(paramSpace)
		print next(paramSpace.sample_uniformly(2))
		# create reductor
		def reductor(discretization, rb, extends = None):
			return reduce_generic_rb(dis, rb, extends=extends)
		
		print 'Do greedy search...'
		
		# greedy search to construct RB 		
		self._rb = greedy(dis, reductor, paramSpace.sample_uniformly(num_samples), use_estimator=False, extension_algorithm=pod_basis_extension, target_error=1e-10, max_extensions = 10, error_norm=lambda U: np.max(dis.l2_norm(U)))  
		# get the reduced discretization and the reconstructor
		self._rd, self._rc = self._rb['reduced_discretization'], self._rb['reconstructor']
		self._bas = self._rb['basis']
		#sjkdg =adas

		print 'Greedy search successfull! Reduced basis has dimension: '+str(len(self._rb['basis']))
		# If saving desired, save the reduced basis and the reconstructor to the disc
		if self._save:
			print 'Saving reduced basis and reconstructor...'
			self._CI.saveSignature(self._signFile, signature) 
			with open(signature, 'w') as f:
				dump((self._rd, self._rc),f)		

	def compute(self, parameter_set=None, error=False, file=None):
		"""
		Think about error estimators and parameters etc etc
		Compute rb solutions for parameter set - with errors?
		"""
		# assert right set structure
		assert isinstance(parameter_set,np.ndarray) or isinstance(parameter_set, list) or parameter_set == None
		if self._type == 'direct': 
			# just supports return of one solution so far!!! Due to matlab restrictions - or glue them all together in the end and decompose them in matlab
			for mu in parameter_set:
				u = self._rd.solve(mu)
				ur = self._rc.reconstruct(u)
		
				print ur
				return ur	
		if self._type == 'disc':
			# Get parameter_set from disc
		 	assert parameter_set == None	
			print 'Computing solutions for given parameter set in respect to reduced basis...'
			parameter_set = self._CI.getParameterSet()
			solutions = {}
			i = 0
			# save solutions for all parameters	
			for mu in parameter_set:
				i=i+1
				u = self._rd.solve(mu)
				# Use data function to transform NumpyVectorArray to standard NumpyArray

				# Add dirichlet values				

				# Have to define valid matlab variable names
				# 'mu'+str(int(mu*100))
				solutions['mu'+str(i)]=(self._rc.reconstruct(u)).data
				#skajgs =sdfkasf
			# save solutions to disk
			self._CI.writeSolutions(solutions,file)
			
	def getRB(self):
		"""
		DOC ME
		"""
		return self._rb
	
	def __str__(self):
		"""
		Give detailed information about the class(which equation can be solved? etc.)
		"""
		pass

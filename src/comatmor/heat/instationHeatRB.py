# Class to apply pymor RB on a stationary PDE
# by Falk Meyer, 10.02.2015

import pickle
import time

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
from pymor.algorithms.timestepping import ImplicitEulerTimeStepper


# local imports
#from ..comminterface import comminterface as CI
from comminterface import comminterface as CI
import parameterHeateq as parameter

class instationHeatRB(object):
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

	def addMatrix(self, name=None, row=None, col=None, data=None, paramName=None, paramShape=None, paramRange=None):
		"""
		paramRange = [1,2] for intervall so to speak, even discrete or continuous
		"""
		# just call comminterface
		self._CI.pushMat(name, row, col, data, paramName, paramShape, paramRange)
		self._matDict = self._CI.getMat()

	def getMatrix(self):
		"""
		Namen umbenennen vermoege besseren Wissens!
		"""
		return self._CI.getMat()

	# DEPRECATED
	#def addRhs(self, rhs=None):
	#	"""
	#	DOC ME
	#	"""
	#	# also assert right dimension?
	#	if rhs != None:
	#		assert isinstance(rhs, np.ndarray)
	#		self._rhs = NumpyMatrixOperator(rhs)
	#	else:
	#		# try to get rhs from matDictionary(has to be transposed!)
	#		self._rhs = NumpyMatrixOperator(self._CI.pushRhs().T)

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

	def getMass(self):
		"""
		Read mass matrix
		"""
		if self._type == 'disc':
			self._mass = NumpyMatrixOperator(self._CI.readMass())
		else:
			pass

	def assembleOperators(self):
		"""
		Assemble operators to be prepared for RB calculations 
		"""
		matDict = self._matDict

		# Create lincomboperator for all involved matrices
		stiffOps, rhsOps, stiffPops, rhsPops = [], [], [], []
		# maybe improve method below - put it somewhere else prob. !!!
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

	def constructRB(self, T=0, n=0):
		"""
		Construct reduced basis with end-time T
		"""
		print 'Constructing reduced basis...' 
		# call assembleOperators and get right operators
		stiffOp, rhsOp = self.assembleOperators()
		#print stiffOp	
		#print self._rhs._matrix
		#raw_input()		
		
		# get parametertypes and ranges
		paramTypes  = self._matDict[1]
		paramRanges = self._matDict[2]
		
		# create timestepper
		time_stepper = ImplicitEulerTimeStepper(n)

		# create mass matrix
		#dimension = op.source.dim
		#mass = NumpyMatrixOperator(np.eye(dimension))

		#print stiffOp.assemble((4,2))._matrix
		#raw_input()
		#print rhsOp.assemble((4,2))._matrix
		#io.savemat('Rhs',{'r': rhsOp.assemble((4,2))._matrix})
		#raw_input()

		# create discretization
		dis = InstationaryDiscretization(operator=stiffOp, rhs=rhsOp, initial_data=self._u0, T=T, time_stepper=time_stepper, mass=self._mass)

		io.savemat('mass',{'mass': self._mass._matrix})	
		exit()
			
		R = dis.solve((1,1))
		self._CI.writeSolutions({'R': R.data},file='Test11')
		exit()
		# create parameterSpace
		paramSpace = CubicParameterSpace(parameter_type = paramTypes, ranges = paramRanges)
		print 'Given parameterspace: '+str(paramSpace)

		# create reductor
		def reductor(discretization, rb, extends = None):
			return reduce_generic_rb(dis, rb, extends=extends)
		
		print 'Do greedy search...'
		
		# greedy search to construct RB 		
		self._rb = greedy(dis, reductor, paramSpace.sample_uniformly(10), use_estimator=False, extension_algorithm=pod_basis_extension, target_error=1e-10, max_extensions = 10) 
		# get the reduced discretization and the reconstructor
		self._rd, self._rc = self._rb['reduced_discretization'], self._rb['reconstructor']

		print 'Greedy search successfull! Reduced basis has dimension: '+str(len(self._rb['basis']))


	def compute(self, training_set=None, error=False, file=None):
		"""
		Think about error estimators and parameters etc etc
		Compute rb solutions for training set - with errors?
		"""
		# assert right set structure
		assert isinstance(training_set,np.ndarray) or isinstance(training_set, list)
		if self._type == 'direct': 
			# just supports return of one solution so far!!! Due to matlab restrictions - or glue them all together in the end and decompose them in matlab
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
			self._CI.writeSolutions(solutions,file)
			
	def getRB(self):
		"""
		DOC ME
		"""
		return self._rb
	
	def save(self, file=None):
		"""
		provide opportunity to save current object with pickle
		perhaps think of deleting self._matDict first, 'cause reduced basis essential.
		"""
		# if no filename given, take standard one
		if file==None:
			# give fileName a local time stamp depending on the current day and time
			t = time.localtime()
			file = str(t[2])+'_'+str(t[1])+'_'+str(t[0])+'_'+str(t[3])+str(t[4])+'_StationRB.save'
		# Save current object
		e = open(file,'w')
		pickle.dump(self,e)
	
	@staticmethod
	def load(path):
		"""
		CLASS METHOD (no instance of class needed to call it)
		Load existing stationRB object located in path
		"""
		e = open(path,'r')
		return pickle.load(e)
	
	def __str__(self):
		"""
		Give detailed information about the class(which equation can be solved? etc.)
		"""
		pass

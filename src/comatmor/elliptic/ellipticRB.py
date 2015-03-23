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

# Class to apply pymor RB on an elliptic PDE
# by Falk Meyer, 10.02.2015

import pickle
import time

import numpy as np

from scipy import sparse,io
from functools import partial

# pymor includes
from pymor.operators.constructions import LincombOperator
from pymor.operators.numpy import NumpyMatrixOperator

from pymor.la.numpyvectorarray import NumpyVectorArray

from pymor.discretizations.basic import StationaryDiscretization

from pymor.reductors.linear import reduce_stationary_affine_linear

from pymor.parameters.spaces import CubicParameterSpace
from pymor.parameters.base import Parameter

from pymor.algorithms.greedy import greedy
from pymor.algorithms.basisextension import trivial_basis_extension
from pymor.algorithms.basisextension import gram_schmidt_basis_extension

from pymor.core.pickle import dump, load

# comatmor imports
#from ..comminterface import comminterface as CI
from comminterface import comminterface as CI


class ellipticRB(object):
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
		# Default
		self._save = False

		# for inputmethod = disc call all necessary get-functions
		if self._type == 'disc':
			self.addMatrix()
                        # Bool to check whether saving is desired
                        self._save = True
                        # Given signature file
                        self._signFile = 'sign.txt'	

	def addMatrix(self, name=None, row=None, col=None, data=None, paramName=None, paramShape=None, paramRange=None):
		"""
		Add matrices to given ellipticRB object - direct or by disc access
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

	# moved to comminterface
        #def assembleOperators(self):
        #        """
        #        Assemble operators to be prepared for RB calculations 
        #        """
        #        matDict = self._matDict
	#
        #        # Create lincomboperator for all involved matrices
        #        stiffOps, rhsOps, stiffPops, rhsPops = [], [], [], []
        #        # maybe improve method below - put it somewhere else prob. !!!
        #        # get stiffness matrices and rhs matrices
        #        for key in parameter.stiffNames:
        #        #for key in matDict[0]:
        #                stiffOps.append(NumpyMatrixOperator(matDict[0][key][0]))
        #                stiffPops.append(matDict[0][key][1])
        #        for key in parameter.rhsNames:
	#		print(matDict[0][key][0].dtype)
        #                rhsOps.append(NumpyMatrixOperator(matDict[0][key][0].T))
        #                rhsPops.append(matDict[0][key][1])
        #        stiffOp = LincombOperator(stiffOps, coefficients=stiffPops)
        #        rhsOp = LincombOperator(rhsOps,coefficients=rhsPops)
        #        return stiffOp, rhsOp


	def constructRB(self, num_samples = 10):
		"""
		TO DO: add possibility to change basis extension for instance
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
                stiffOp, rhsOp = self._CI.assembleOperators()

		paramTypes  = self._matDict[1]
		paramRanges = self._matDict[2]
	
		#print self._rhs._matrix
		#io.savemat('pyrhs',{'pyrhs': rhsOp.assemble((4,2))._matrix})	
		#io.savemat('pyop',{'pyop': stiffOp.assemble((4,2))._matrix})	
		
		# create discretization
		dis = StationaryDiscretization(operator=stiffOp, rhs=rhsOp)

		# create parameterSpace
		paramSpace = CubicParameterSpace(parameter_type = paramTypes, ranges = paramRanges)
		print 'Given parameterspace: '+str(paramSpace)

		# create reductor
	 	reductor = partial(reduce_stationary_affine_linear, error_product = None)
		print 'Do greedy search...'
		
		# greedy search to construct RB 		
		self._rb = greedy(dis, reductor, paramSpace.sample_uniformly(num_samples), use_estimator=False, extension_algorithm=gram_schmidt_basis_extension, target_error=1e-10, max_extensions = 10) 
		# get the reduced discretization and the reconstructor
		self._rd, self._rc = self._rb['reduced_discretization'], self._rb['reconstructor']

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
			# Get parameter set from disc
                        assert parameter_set == None
			print 'Computing solutions for given parameter set in respect to reduced basis...'

			parameter_set = self._CI.getTrainingSet()
			solutions = {}
			i = 0
			# save solutions for all parameters	
			for mu in parameter_set:
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
	
	def __str__(self):
		"""
		Give detailed information about the class(which equation can be solved? etc.)
		"""
	

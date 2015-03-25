# This file was part of the pyMOR project (http://www.pymor.org).
# Copyright Holders: Rene Milk, Stephan Rave, Felix Schindler
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
# modified by Falk Meyer for the comatmor model 

""" This module provides generic time-stepping algorithms for the solution of
instationary problems.

The algorithms are generic in the sense that each algorithms operates exclusively
on |Operators| and |VectorArrays|. In particular, the algorithms
can also be used to turn an arbitrary stationary |Discretization| provided
by an external library into an instationary |Discretization|.

Currently, implementations of :func:`explicit_euler` and :func:`implicit_euler`
time-stepping are provided. The :class:`TimeStepperInterface` defines a
common interface that has to be fulfilled by the time-steppers that are used
by |InstationaryDiscretization|. The classes :class:`ExplicitEulerTimeStepper`
and :class:`ImplicitEulerTimeStepper` encapsulate :func:`explicit_euler` and
:func:`implicit_euler` to provide this interface.
"""

from __future__ import absolute_import, division, print_function

from pymor.core.interfaces import ImmutableInterface, abstractmethod
from pymor.la.interfaces import VectorArrayInterface
from pymor.operators.interfaces import OperatorInterface
from pymor.algorithms.timestepping import TimeStepperInterface

# local and additional includes
from scipy import io
import csv
import copy


class ImplicitEulerTimeStepper(TimeStepperInterface):
    """Implict-Euler time-stepper.

    Solves equations of the form ::

        M * d_t u + A(u, mu, t) = F(mu, t).

    Parameters
    ----------
    nt
        The number of time-steps the time-stepper will perform.
    invert_options
        The :attr:`~pymor.operators.interfaces.OperatorInterface.invert_options` used
        to invert `M + dt*A`.
    """

    def __init__(self, nt, invert_options=None):
        self.nt = nt
        self.invert_options = invert_options

    def solve(self, initial_time, end_time, initial_data, operator, rhs=None, mass=None, mu=None, num_values=None):
        return implicit_euler(operator, rhs, mass, initial_data, initial_time, end_time, self.nt, mu,
                              self.invert_options, num_values)

def implicit_euler(A, F, M, U0, t0, t1, nt, mu=None, invert_options=None, num_values=None):
    assert isinstance(A, OperatorInterface)
    assert isinstance(F, (OperatorInterface, VectorArrayInterface))
    assert isinstance(M, OperatorInterface)
    # Changed assertion for M to not non parametric
    assert A.source == A.range == M.source == M.range
    num_values = num_values or nt + 1
    dt = (t1 - t0) / nt
    DT = (t1 - t0) / (num_values - 1)

    if isinstance(F, OperatorInterface):
        assert F.range.dim == 1
        assert F.source == A.range
        F_time_dep = F.parametric and '_t' in F.parameter_type
        if not F_time_dep:
            dt_F = F.as_vector(mu) * dt
    else:
        assert len(F) == 1
        assert F in A.range
        F_time_dep = False
        dt_F = F * dt

    assert U0 in A.source
    assert len(U0) == 1

    A_time_dep = A.parametric and '_t' in A.parameter_type

    R = A.source.empty(reserve=nt+1)
    R.append(U0)

    M_dt_A = M + A * dt
    if not A_time_dep:
        M_dt_A = M_dt_A.assemble(mu)

    t = t0
    U = U0.copy()

    # read values of dirichlet function
    f=open('dirichlet.csv')
    values = [];
    for row in csv.reader(f,delimiter=','):
	values.append(row)	

    for n in xrange(nt):
        t += dt

	# Imply Dirichlet - will be linearly interpolated when time values do not match
	for i in range(1,len(values)):
		if t >= float(values[i-1][0]) and t <= float(values[i][0]):
			mu['_t'] = float(values[i-1][1]) + (float(values[i][1]) - float(values[i-1][1]))/(float(values[i][0])-float(values[i-1][0]))*(t-float(values[i-1][0]))
			break

	#mu['_t'] = [float(values[i][1]) for i in range(0,len(values)) if t==float(values[i][0])][0]
        if F_time_dep:
            dt_F = F.as_vector(mu) * dt
	
	U = M_dt_A.apply_inverse(M.apply(U,mu=mu)+dt_F, mu=mu,options=invert_options)	
        while t - t0 + (min(dt, DT) * 0.5) >= len(R) * DT:
            R.append(U)

    return R

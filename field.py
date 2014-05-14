from stdlib import *
from numpy import savez
from os.path import isfile
from copy import copy
from scipy.ndimage.interpolation import map_coordinates

# cartesian products, variable dimension of common coordinate space, grids of lower rank than coordinate space dimension (rank, dimension are identifiers)
# TODO extrapolation inserts nans, support counts nans as zero
# Field as a subclass of ndarray

# Glossary
#
# grid: a set of grid coordinates
# points: a set of common coordinates

_albls = ['_s', '_x', '_y', '_z']

class Grid:
	"""
All operations on a single grid refer to that grid's coordinates.
"""

		# integrating a field can return a lower-dimensional field
		
	def __str__(self):
		return "<grid " + " in R^" + str(self.dim()) + \
			" shape " + str(self.shape) + \
			" over " + str([tuple(b) for b in self.bounds()]) + ">"
			
	def __repr__(self):
		return str(self)
		
	#
	# construction
	#
		
	@classmethod
	def from_axes(cls, *axes):
		# FIXME assert result.axes() are close to axes
		axes = [array(q).flatten() for q in axes]
		N = array([q.size for q in axes])
		origin = [q[0] for q in axes]
		return cls(
			shape=N,
			h=[ptp(q)/(n-1) for n, q in zip(N, axes)],
			U=eye(len(N)),
			p=origin,
			o=origin)
	
	@classmethod
	def default(cls):
		ffile = load('fields.npz')
		return cls.from_axes(*[ffile[q] for q in _albls])

	def __init__(self, shape, o, p, h, U):
		# instance variables explained in geometry.tex
		self.U = array(U, dtype=float)
		assert allclose(dot(self.U.T, self.U), eye(self.rank()))
		self.shape = tuple(int(n) for n in shape)
		assert len(self.shape) == self.rank()
		assert all(n>1 for n in self.shape)
		self.h = array(h, dtype=float).reshape((self.rank(),))
		assert (self.h>0).all()
		self.p = array(p, dtype=float).reshape((self.rank(),))
		self.o = array(o, dtype=float).reshape((self.dim(),))
			
	def __mul__(self, other):
		"cartesian product of grids on cartesian product of common coordinates"
		U = zeros((self.dim()+other.dim(), self.rank()+other.rank()))
		U[:self.dim(),:self.rank()] = self.U
		U[self.dim():,self.rank():] = other.U
		return Grid(
			shape = self.shape + other.shape,
			o = concatenate((self.o, other.o)),
			p = concatenate((self.p, other.p)),
			h = concatenate((self.h, other.h)),
			U = U)
	
	#
	# basic properties
	#
	
	def rank(self):
		"return the rank of arrays of samples"
		return self.U.shape[1]
		
	def dim(self):
		"return the dimension of the common coordinate space"
		return self.U.shape[0]
		
	def axes(self):
		return [p+h*arange(n) for p, h, n in zip(self.p, self.h, self.shape)]
		
	def bounds(self):
		return array([[q.min(), q.max()] for q in self.axes()])
			
	#
	# indices and coordinates
	# these assume the first axis of the coordinate arrays is the components
	# if we don't have full rank, project othogonally
	#
		
	def r(self, i=None, R=None):
		assert i is None or R is None
		if i is not None:
			i = array(i)
			h, p, o = self._vectors(i)
			return p + h*i
		if R is not None:
			R = array(R)
			h, p, o = self._vectors(R)
			return p + tensordot(self.U.T, R-o, 1)		
		else:
			return Field(self.r(i=self.i()), self)
		
	def R(self, i=None, r=None):
		"""convert indices or grid coordinates to common coordinates.

when called with no arguments, return a field of the common coordinates for each grid point.  otherwise, return the common coordinates of the points with indices i or grid coordinates r: the first dimension is coordinate components.
"""
		assert i is None or r is None
		if i is not None:
			i = array(i)
			h, p, o = self._vectors(i)
			return o + tensordot(self.U, h*i, 1)
		if r is not None:
			r = array(r)
			h, p, o = self._vectors(r)
			return o + tensordot(self.U, r-p, 1)
		else:
			return Field(self.R(i=self.i()), self)
		
	def i(self, r=None, R=None):
		assert r is None or R is None
		if r is not None:
			r = array(r)
			h, p, o = self._vectors(r)
			return (r-p)/h
		if R is not None:
			R = array(R)
			h, p, o = self._vectors(R)
			return tensordot(self.U.T, R-o, 1)/h
		else:
			return Field(indices(self.shape), self)
		
	def rr(self, r=None):
		"inertia tensors about the grid origin for a unit mass at each point of r, or the grid by default.  unlear if this is useful for rank greater than three."
		if r is None:
			r = self.r()
		else:
			r = ascontiguousarray(r)
		P = -r[newaxis,::]*r[:,newaxis,::]
		# hack strides to address the diagonal plane of P
		Pd = ndarray(buffer=P, dtype=P.dtype, \
			shape=r.shape, strides=((3+1)*P.strides[1],)+ P.strides[2:])
		Pd += (r**2).sum(0)
		return P
		
	#
	# utility
	# 
	
	def _vectors(self, A):
		"return h, p,and  o, reshaped to broadcast over A"
		tail = (1,)*(A.ndim-1)
		return self.h.reshape((self.rank(),)+tail), \
			self.p.reshape((self.rank(),)+tail), \
			self.o.reshape((self.dim(),)+tail)
			
	def _clone(self, **args):
		"return a Grid like self, but with some things changed"
		d = dict([(x, getattr(self, x)) for x in ['shape', 'h', 'o', 'p', 'U']])
		d.update(args)
		return Grid(**d)
	
	#
	# transformation
	#
	
	def __getitem__(self, subidx):
		"for now, we only handle slices.  could have integers return a lower rank grid."
		def bound(i, n, m):
			if i is None:
				i = m
			if i<0:
				i += n
			assert 0 <= i and i <= n
			return i
			
		assert all([type(x) is slice and x.step is None for x in subidx])
		lowcorner = array([bound(x.start, n, 0) for x, n in zip(subidx, self.shape)])
		highcorner = array([bound(x.stop, n, n) for x, n in zip(subidx, self.shape)])
		return self._clone(shape=highcorner - lowcorner,
			p = self.p + self.h*lowcorner,
			o = self.o + self.h*lowcorner)
	
	def shifted(self, new_origin):
		"translate grid coordinates, while leaving the grid fixed in common coordinates"
		return Grid(p=self.p-new_origin, o=self.o, h=self.h, shape=self.shape, U=self.U)
		
	def translated(self, dr):
		"Shift the grid points, instead of the coordinates"
		assert False
		
	def rotated(self, V, centre=None):
		"return the same grid, representing points transfomed by V about centre (centre is in grid coordinates)"
		if centre is None:
			centre = zeros(self.rank())
		assert V.shape == 2*(self.dim(),)
		Rc = self.R(r=centre)
		return self._clone(o=Rc+dot(V, self.o-Rc), U=dot(V,self.U))
	
	# loading and saving to file
		
	def be_default(self):
		ftab = dict(load('fields.npz'))
		ftab.update(dict(zip(_albls, self.axes)))
		savez('fields.npz', **ftab)
		
	def blank(self):
		return empty(self.shape)
		
	#
	# comparision
	#
		
	def __eq__(self, other):
		assert False	# testing floats for equality is a sin
		
	def __neq__(self, other):
		return not self.__eq__(other)
		
	def close_points(self, other):
		"test if the points of self and other have similar common coordinates"
		return self.shape == other.shape and \
			all([allclose(getattr(self, x), getattr(other, x))
				for x in ['h', 'o', 'U']])
	
	def close_grid(self, other):
		"test if the points of self and other have similar grid coordinates"
		return self.shape == other.shape and \
			all([allclose(getattr(self, x), getattr(other, x))
				for x in ['h', 'p']])
		
	def close(self, other):
		return self.close_points(other) and self.close_grid(other)
		
	#
	# Fourierology
	#
		
	#
	# integration
	#
		
	def S(self, ordinates):
		# integrate: should allow axes to be picked
		# the method ndarray.sum returns an ndarray, not a Field,
		# so this doesn't trigger shape checks.
		return prod(self.h)*ordinates.sum(tuple(range(-self.rank(),0)))
		
		
	
# For now, load and store read and write the whole file every time a
# saved field changes. The other way would be a singleton object that
# cached the file as a dict, and automatically saved the fields on exit; doing that quickly is beyond my Python
# ability.

def load_field(lbl):
	ffile = load('fields.npz')
	return Field(ffile[lbl], label=lbl)
	

class Field(ndarray):
	
	# FIXME binary ufuncs should check the grids are the same
	
	def __new__(cls, ordinates, abscissae=Grid.default(), label=None):
		obj = asarray(ordinates).view(cls)
		obj.abscissae = abscissae
		assert obj.shape[-obj.abscissae.rank():] == obj.abscissae.shape
		if label is not None:
			obj.label = label
		return obj
				
	def __array_finalize__(self, obj):
		if obj is None: return
		self.abscissae = getattr(obj, 'abscissae', None)
		# catch reductions
		if self.abscissae:
			assert self.shape[-self.abscissae.rank():] == self.abscissae.shape
			
	def __str__(self):
		return '<You lose>'
		
	def __repr__(self):
		return '<You lose>'
	
	#
	# methods delegated to abscissae
	#
	
	def r(self):
		return self.abscissae.r()
	def R(self):
		return self.abscissae.R()
	def i(self):
		return self.abscissae.i()
	def rr(self):
		return self.abscissae.rr()
	
	def S(self):
		"should allow dimensions to be specified"
		return self.abscissae.S(self)
	
	#
	# loading and saving
	#
	
	def save(self, label=None):
		if label:
			self._label = label
		assert self.abscissae == Grid.default()
		ftab = dict(load('fields.npz'))
		ftab[self._label] = self.ordinates
		savez('fields.npz', **ftab)
	
	#
	# interpolation
	#
	
	def sampled(self, abscissae):
		# shortcut: interpolation on the same grid
		if abscissae.close_points(self.abscissae):
			return Field(self, abscissae)
		# shortcut: time slicing n.y.i.
		else:
			return Field(
				map_coordinates(self, self.abscissae.i(R=abscissae.R()), cval=nan),
				abscissae)
	
	#
	# misc
	#

	def support(self, cut=0):
		"return a subgrid of abscissae, on which ordinates exceed cutoff"
		A = array(self)
		A[isnan(A)] = -inf
		bools = array(nonzero(A>cut))
		return self.abscissae[[slice(m,n) for m, n in zip(bools.min(axis=1), bools.max(axis=1) + 1)]]
	
	def central_frame(self):
		"return a Grid with origin at my centre of mass, and axes aligned to my principle axes.  only sensible for scalar fields."
		assert False

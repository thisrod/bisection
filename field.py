from stdlib import *
from numpy import savez
from os.path import isfile
from copy import copy
from scipy.ndimage.interpolation import map_coordinates

# cartesian products, variable dimension of common coordinate space, grids of lower rank than coordinate space dimension (rank, dimension are identifiers)
# TODO extrapolation inserts nans, support counts nans as zero
# Field as a subclass of ndarray

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
			return array(meshgrid(*self.axes(), indexing='ij'))
		
	def R(self, i=None, r=None):
		assert i is None or r is None
		if i is None and r is None:
			i = indices(self.shape)
		if i is not None:
			i = array(i)
			h, p, o = self._vectors(i)
			return o + tensordot(self.U, h*i, 1)
		if r is not None:
			r = array(r)
			h, p, o = self._vectors(r)
			return o + tensordot(self.U, r-p, 1)	
		
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
			return indices(self.shape)
		
	#
	# utility
	# 
	
	def _vectors(self, A):
		"return h, p,and  o, reshaped to broadcast over A"
		tail = (1,)*(A.ndim-1)
		return self.h.reshape((self.rank(),)+tail), \
			self.p.reshape((self.rank(),)+tail), \
			self.o.reshape((self.dim(),)+tail)
	
	#
	# transformation
	#
	
	def shifted(self, new_origin):
		# index origin not changed in common coordinates,
		# moved by -new_origin in grid coordinates.
		S = copy(self)
		new_origin = concatenate(([0], new_origin))
		S._configure(p=self.p - new_origin)
		return S
		
	def translated(self, dr):
		"Shift the grid points, instead of the coordinates"
		assert False
		
	def rotated(self, V, centre=None):
		"return grid transfomed by V about centre (centre is in grid coordinates).  unlike shifted, this returns a different set of points"
		if centre is None:
			centre = self.origin()
		Rc = self.o + dot(self.U, centre-self.p)
		S = copy(self)
		S._configure(U=V)
		S._configure(o=Rc + dot(S.U,self.o-Rc))
		S._configure(U=dot(V, self.U3))
		return S
	
	# loading and saving to file
		
	def be_default(self):
		ftab = dict(load('fields.npz'))
		ftab.update(dict(zip(_albls, self.axes)))
		savez('fields.npz', **ftab)
		
	def blank(self):
		return empty(self.shape)
		
	def __eq__(self, other):
		assert False
		return all([(q == p).all() for q, p in zip(self.axes, other.axes)])
		
	def __neq__(self, other):
		return not self.__eq__(other)
		
	def rr(self):
		"Return a 3*3*1*shape[1:] array of inertia tensors about the grid origin for a unit mass at each point."
		r = self.r()[1:,0:1,::]
		P = -r[newaxis,::]*r[:,newaxis,::]
		# hack strides to address the diagonal plane of P
		Pd = ndarray(buffer=P, dtype=P.dtype, \
			shape=r.shape, strides=((3+1)*P.strides[1],)+ P.strides[2:])
		Pd += (r**2).sum(0)
		return P
		
	def S(self, ordinates):
		# integrate: should allow axes to be picked
		# the method ndarray.sum returns an ndarray, not a Field,
		# so this doesn't trigger shape checks.
		return prod(self.h)*ordinates.sum(tuple(range(-self.rank(),0)))
		
	def sample(self, fld):
		return Field(fld.samples(self), self)
	
	def subgrid(self, corner, shape):
		return Grid(*[q[c:c+s] for q, c, s in zip(self.axes(), corner, shape)],
			orientation=self.U3, origin=self.origin())
		
		
	
# For now, load and store read and write the whole file every time a
# saved field changes. The other way would be a singleton object that
# cached the file as a dict, and automatically saved the fields on exit; doing that quickly is beyond my Python
# ability.

def load_field(lbl):
	ffile = load('fields.npz')
	return Field(ffile[lbl], label=lbl)
	

class Field(ndarray):
	
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
			
	def S(self):
		"should allow dimensions to be specified"
		return self.abscissae.S(self)
		
	def save(self, label=None):
		if label:
			self._label = label
		assert self.abscissae == Grid.default()
		ftab = dict(load('fields.npz'))
		ftab[self._label] = self.ordinates
		savez('fields.npz', **ftab)
	
	def samples(self, points):
		return map_coordinates(self.ordinates, self.abscissae.indices_for(points.R()), cval=nan)

	def support(self, cut=0):
		"return a subset of the sampling grid, where my values exceed cutoff"
		bools = array(nonzero(self.ordinates>cut))
		return self.abscissae.subgrid(bools.min(axis=1), bools.ptp(axis=1) + 1)
	
	def central_frame(self):
		"return a Grid with origin at my centre of mass, and axes aligned to my principle axes.  only sensible for scalar fields."
		assert False

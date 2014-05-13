from stdlib import *
from numpy import savez
from os.path import isfile
from copy import copy
from scipy.ndimage.interpolation import map_coordinates

# TODO: extrapolation inserts nans, support counts nans as zero

_albls = ['_s', '_x', '_y', '_z']
_vector = (4,1,1,1,1)		# shape in which Grid stores its defining vectors

class Grid:
	"""
All operations on a single grid refer to that grid's coordinates. To manage the bookeeping for interpolation, these grids are related to a common set of coordinates by the relations r = p+Hi and R = o+HUi, where r is a set of grid coordinates, i is an index in the grid, and R is a set of common coordinates.  We will refer to p as the index origin, and o as the grid origin.  U is a unitary matrix, that doesn't rotate the time axis.
"""
	

	# See geometry.tex for the geometric definitions of the instance variables.
	# U is a unitary 3*3 matrix, the rest are vectors.  H is the matrix with the vector h along its diagonal.  The vectors are stored with shape (4,1,1,1,1), for broadcasting over grids.
	# the vectors come in o3, p3, h3 and U3 versions, with only the space dimensions.  
	# there are also ov, pv and h3, which are constant vector fields
	
	def __init__(self, s, x, y, z, origin=zeros(4), orientation=eye(3)):
		"""
The grid coordinates are defined by the axes s, x, y and z.  These should be ndarrays; the constructor assumes the values are evenly spaced, and deduces the spacing from the first and last coordinates.  The origin of the grid is implied by where the values in the axis arrays.

The parameters origin and orientation define where the grid sits in R^4.  Origin is the point in common coordinates at which s=x=y=z=0; orientation is the unitary matrix U from the coordinate relation.  Note that the coordinate origin, provided as a parameter, is not generally the same as the index origin, stored as self.o.
"""
		# TODO: implement sections by allowing U to be rectangular
		# integrating a field can return a lower-dimensional field
		
		# step 1 is used for singleton axes to make H invertible
		axes = [q.flatten() for q in [s, x, y, z]]
		N = array([q.size for q in axes])
		self._configure(
			shape=N,
			h=[1 if n==1 else ptp(q)/(n-1) for n, q in zip(N, axes)],
			U=orientation,
			p=[q[0] for q in axes],
			o=empty(4))
		self._configure(o=array(origin) + dot(self.U,self.p))

	def _configure(self, **args):
		"set up instance variables and assert sanity"
		
		assert set(args.keys()).issubset(set(['shape', 'o', 'p', 'h', 'U']))
		if 'shape' in args:
			self.shape = tuple(args['shape'])
			assert len(self.shape) == 4
			assert all(n>0 for n in self.shape)
		if 'U' in args:
			W = eye(4);  self.U = W
			self.U3 = W[1:,1:]
			self.U3[::] = args['U']
			assert allclose(dot(self.U3.T, self.U3), eye(3))
		if 'h' in args:
			self.h = array(args['h'])
			assert self.h.size == 4
			assert (self.h>0).all()
			self.h3 = self.h[1:]
			self.hv = self.h.reshape(_vector)
		if 'p' in args:
			self.p = array(args['p'])
			assert self.p.size == 4
			self.p3 = self.p[1:]
			self.pv = self.p.reshape(_vector)
		if 'o' in args:
			self.o = array(args['o'])
			assert self.o.size == 4
			self.o3 = self.o[1:]
			self.ov = self.o.reshape(_vector)
		
	def axes(self):
		return [p+h*arange(n) for p, h, n in zip(self.p, self.h, self.shape)]
		
	def origin(self):
		"return origin as would be passed to constructor"
		return self.o - dot(self.U, self.p)
		
	def r(self):
		"grid coordinates of points"
		return array(meshgrid(*self.axes(), indexing='ij'))
		
	def R(self, i=None):
		"""common coordinates of points at given indices, by default the grid indices.
		"""
		assert i is None
		return self.ov + tensordot(self.U, self.hv*indices(self.shape), 1)
		
	def indices_for(self, R):
		"R has shape 4*whatever, as grid coordinates conventionally do"
		return tensordot(self.U.T, R-self.ov, 1)/self.hv
		
	def bounds(self):
		return array([[q.min(), q.max()] for q in self.axes()]).T
		
	# transformation
	
	def shifted(self, new_origin):
		# index origin not changed in common coordinates,
		# moved by -new_origin in grid coordinates.
		S = copy(self)
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
	
	@classmethod
	def default(cls):
		ffile = load('fields.npz')
		return cls(*[ffile[q] for q in _albls])
		
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
		# integrate over space
		return prod(self.h3)*ordinates.sum((-3, -2, -1))
		
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
	

class Field:
	"""A three plus 1 dimensional field.
	
	Field just packs up a grid and a set of samples.  plan: anything we can do with a Field, we can do with the Grid and the array of samples.
	"""
	
	def __init__(self, ordinates, abscissae=Grid.default(), label=None):
		assert abscissae.shape == ordinates.shape
		self.abscissae, self.ordinates = abscissae, ordinates
		if label:
			self._label = label
		
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
		
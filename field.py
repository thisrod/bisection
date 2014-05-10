from stdlib import *
from numpy import savez
from os.path import isfile

# these classes are minimal: we add things to them that are required in several scripts.

_albls = ['_s', '_x', '_y', '_z']

class Grid:
	"The coordinates (s,x,y,z) on the grid represent the point (s,O+U(x,y,z)) in 3+1 dimensional space.  This is just bookeeping for the resampling methods: all operations on a single grid refer to that grid's coordinates, and don't depend on the origin or orientation.  The grid must be even and rectangular for resampling to work."
	
	# the stored attributes are:
	#	self.shape
	#	self.o
	#	self.h
	#	self.U
	# U is a unitary matrix, the rest are vectors.  shape is the size of the grid; a point with index vector i is located at r=o+h*(U.i)
	
	def __init__(self, s, x, y, z, origin=zeros(4), orientation=eye(3)):
		# FIXME: origin is redundant, because the axes can start away from zero.  it should denote the point in R^n that stays fixed when the grid is rotated.
		axes = [q.flatten() for q in [s, x, y, z]]
		self.shape = tuple(q.size for q in axes)
		self.o = origin.reshape((4,1,1,1,1))
		assert allclose(dot(orientation.T, orientation), eye(3))
		self.U = zeros((4,4));  self.U[1:,1:] = orientation;  self.U[0,0] = 1
		self.h = array([ptp(q) for q in axes])/array(self.shape)
		self.h = self.h.reshape((4,1,1,1,1))
		
	def axes(self):
		return [a*arange(n) for a, n in zip(self.h, self.shape)]
		
	def r(self):
		"Coordinates of points wrt grid"
		return array(meshgrid(*self.axes(), indexing='ij'))
		
	def R(self):
		"Coordinates of points in R^n"
		return self.o + self.h*dot(self.U, self.r())
		
	def indices_for(self, R):
		"this assumes R has shape ...*4*?*?*?*?, as grid coordinates conventionally do"
		return dot(self.U.T, (R-self.o)/self.h)
		
	def bounds(self):
		return array([[q.min(), q.max()] for q in self.axes])
	
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
		return all([(q == p).all() for q, p in zip(self.axes, other.axes)])
		
	def __neq__(self, other):
		return not self.__eq__(other)
		
	def rr(self):
		"Return a 3*3*1*shape[1:] array of inertia tensors about the origin for a unit mass at each point."
		r = self.r[1:,0:1,::]
		P = -r[newaxis,::]*r[:,newaxis,::]
		# hack strides to address the diagonal plane of P
		Pd = ndarray(buffer=P, dtype=P.dtype, \
			shape=r.shape, strides=((3+1)*P.strides[1],)+ P.strides[2:])
		Pd += (r**2).sum(0)
		return P
		
	def S(self, ordinates):
		# integrate over space - for now, this is just a sum, dV assumed to cancel out.
		return(ordinates.sum((-3, -2, -1)))
		
		
	
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
		return map_coordinates(self.ordinates, self.abscissae.indices_for(points.R))
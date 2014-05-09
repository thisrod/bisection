from stdlib import *
from numpy import savez
from os.path import isfile

# these classes are minimal: we add things to them that are required in several scripts.

class Grid:
	"A 3D grid"
	
	def __init__(self, s, x, y, z):
		self.ss = s
		self.axes = [s, x, y, z]
		self.gaxes = [x.reshape((-1,1,1)), y.reshape((1,-1,1)), z.reshape((1,1,-1))]
		self.shape = tuple(a.size for a in self.axes)
		
	@classmethod
	def default(cls):
		ffile = load('fields.npz')
		return cls(ffile['_s'], ffile['_x'], ffile['_y'], ffile['_z'])
		
	def __eq__(self, other):
		return all([(a == b).all() for a, b in zip(self.axes, other.axes)])
	
# For now, load and store read and write the whole file every time a
# saved field changes. The other way would be a singleton object that
# cached the file as a dict, and automatically saved the fields on exit; doing that quickly is beyond my Python
# ability.

def load_field(lbl):
	ffile = load('fields.npz')
	return Field(ffile[lbl], label=lbl)
	

class Field:
	"""A three plus 1 dimensional field.
	
	Fields sampled at one time are constant.	
	Idea: Grid knows how to do many things with arrays of samples, Field just packs up a grid and a set of samples.
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
	
	def be_sampled_on(self, abscissae):
		"""g is a grid: the stored field is sampled on it."""
		pass

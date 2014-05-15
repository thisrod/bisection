Field library
===

by Rodney E. S. Polkinghorne

In general, a field is a scalar or vector quantity that is a function of position and time.  For the purposes of this Python library, a `Field` is an `ndarray` containing samples of such a quantity, that remembers the points in R<sup>n</sup> at which it was sampled.  Knowing this, an array can integrate, differentiate and Fourier transform itself, generate Gaussian noise with unit variance, and other such things.  This spares the user from inordinate bookeeping, prevents entire classes of bugs from possibly occuring, and makes the code focussed enough on the physical problem that the remaining errors can be found by physical nous.

The mechanics of this are done by a class `Grid`.  This is a abstract representation of a rectangular grid of points in R<sup>n</sup>, along with a system of grid coordinates.  The points must form a rectangular grid, but this does not have to be aligned with the usual axes, start at the origin, or have rank n.  We can represent a skew plane in space, and things like that.  

The class `Field` is thin: it basically wraps up an array of samples with the grid they were sampled on, and dispatches some operations to the appropriate grid.  The know-how is contained in |Grid|.  This stores an abstract specification of a grid.  Many of its methods take an array of samples as an argument, that is interpreted as a set of samples on that grid.  The methods do things such as integrate the samples, interpolate them on another grid, etc.  There is also a method |r| to return the coordinates of the grid points, and |rr| to give their moments of inertia about the grid origin.

Arrays of samples are assumed to have the shape components*grid_shape.  This has two advantages.  Firstly, multiplying on the left using numpy.dot has the effect of a matrix product in the components, evaluated at each grid point.  Secondly, the arrays returned by |r| and |rr| can be broadcast over the components of a field, and three-dimensional arrays of samples in space can be broadcast over times.

Grids can be shifted or rotated.  In the common frame, the specified transformation is applied to the grid, while fields and other grids stay fixed.

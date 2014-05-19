from field import *

x, y, z, s = [Grid.from_axes(q) for q in
	[arange(2), 0.3*arange(3), pi+arange(5), [0, 17]]]
R = x*y*z
T = s*R
f = Field(R.blank(), R);  f[::] = 1

# check index and coordinate functions

e = array([[1,1,1],[1,3,2]]).T
assert all(type(x) is Field for x in [f.i(), f.r(), f.R(), f.rr()])
assert all(type(x) is ndarray for x in [R.i(r=e), R.r(i=e), R.R(i=3), R.rr(r=e)])
assert allclose(R.r(), R.r(i=indices(R.shape)))
assert allclose(R.r(R=e), e)

assert allclose(e, R.i(R=R.R(i=e)))
assert allclose(e, R.R(i=R.i(R=e)))
# fill the round-trip tests out

# test rotation
U = array([[cos(pi/6), -sin(pi/6), 0], [sin(pi/6), cos(pi/6), 0], [0, 0, 1]])
origin = R.R(r=(0,0,0)).reshape((3,1,1,1))
assert allclose(R.rotated(U).R() - origin, tensordot(U, R.R() - origin, 1))

# test subgrids
Rs = R[0:2, 1:2, 1:3]
assert allclose(Rs.r(), array(R.r())[:, 0:2, 1:2, 1:3])
assert allclose(Rs.R(), array(R.R())[:, 0:2, 1:2, 1:3])
assert Rs.close(R[0:2, 1:2, 1:-2])

# test low rank
assert R.close(R[:,0,0]*R[0,:,0]*R[0,0,:])
# assert y.close(R[pi,:,exp(1)])

# possible bugs
# interpolation on low-rank grids that are close but unequal
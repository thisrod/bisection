from field import *

x = arange(2);  y = 0.3*arange(3);  z = pi+arange(5)
R = Grid.from_axes(x, y, z)
s = Grid.from_axes([0, 17])
T = s*R

e = array([[1,1,1],[1,3,2]]).T
assert allclose(R.r(), R.r(i=indices(R.shape)))
assert allclose(R.r(R=e), e)

assert allclose(e, R.i(R=R.R(i=e)))
assert allclose(e, R.R(i=R.i(R=e)))
# fill the round-trip tests out

# test rotation
T = Grid(*axs)
U = array([[cos(pi/6), -sin(pi/6), 0], [sin(pi/6), cos(pi/6), 0], [0, 0, 1]])
R = T.rotated(U).R() - T.origin().reshape((4,1,1,1,1))
S = T.R() - T.origin().reshape((4,1,1,1,1))
assert allclose(R[0,::], S[0,::])
assert allclose(R[1:,::], tensordot(U, S[1:,::], 1))
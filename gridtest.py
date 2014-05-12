from field import *

s = array([-5]);  x = arange(2);  y = 0.3*arange(3);  z = pi+arange(5)
axs = [s, x, y, z]
bds = array([[q.min() for q in axs], [q.max() for q in axs]])
T = Grid(*axs)

tax = [q.flatten() for q in T.axes()]
assert allclose(axs[0], tax[0])
assert allclose(axs[1], tax[1])
assert allclose(axs[2], tax[2])
assert allclose(axs[3], tax[3])
assert allclose(T.bounds(), bds)
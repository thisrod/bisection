"""
Rotate the potential to a trap-aligned frame, find the initial ground state by imaginary time propagation, and save the resulting fields.

The trap-aligned frame is defined by the centre of mass and principle axes of the final trap, after bisection.
"""

from field import *
from stdlib import *
from scipy.io import loadmat

# everything is in dispersion units, 10^-7 m
# T is the grid supplied by Vienna, shifted to the centre of mass
# S is the computation grid, in the trap-aligned frame, covering just the volume where atoms are to be found.

Kcut = 10		# potential at which to truncate grid V
N = 7e3		# number of atoms
s = 1e-3*arange(18)/1.368e-5; 		# times samples supplied

# load supplied grid
vfile = loadmat('potentials/RWA_X_3D_0.mat')
T = Grid(array(nan), *[10*vfile[q] for q in ['x', 'y', 'z']])

# load final potential, shift to avoid underflow
K = 0.1719*loadmat('potentials/RWA_X_3D_17.mat')['v']
K -= K.min()

# shift grid origin to centre of weight
wgt = exp(-K/(2*2.88**2))
r0 = T.S(T.r()*wgt)/T.S(wgt)
T = T.shifted(r0.flatten())

# find the principle axes and
# set U to the rotation from the trap-aligned frame to the original data frame
ew, ev = eig(T.S(wgt*T.rr())[:,:,0])
U = ev[:,(1,2,0)]	# order of increasing moments is z, x, y
U = dot(U, diagflat(sign(diag(U))))	# align senses
S = T.rotated(U)
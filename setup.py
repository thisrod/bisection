"""
Rotate the potential to a trap-aligned frame, find the initial ground state by imaginary time propagation, and save the resulting fields.

The trap-aligned frame is defined by the centre of mass and principle axes of the final trap, after bisection.
"""

from field import *
from scipy.io import loadmat

# everything is in dispersion units, 10^-7 m
# T is the grid supplied by Vienna, shifted to the centre of mass
# S is the computation grid, in the trap-aligned frame, covering just the volume where atoms are to be found.

Kcut = 10		# potential at which to truncate grid V
N = 7e3		# number of atoms
s = Grid.from_axes(1e-3*arange(18)/1.368e-5)	# times samples supplied

# load supplied grid
vfile = loadmat('potentials/RWA_X_3D_0.mat')
T = Grid.from_axes(*[10*vfile[q] for q in ['x', 'y', 'z']])

# load final potential, shift to avoid underflow
K = Field(0.1719*loadmat('potentials/RWA_X_3D_17.mat')['v'], T)
K -= K.min()

# shift grid origin to centre of weight
wgt = exp(-K/(2*2.88**2))
T = T.shifted((wgt.r()*wgt).S()/wgt.S())
wgt = wgt.sampled(T)

# find the principal axes and
# set U to the rotation from the trap-aligned frame to the original data frame
ew, ev = eig((wgt*wgt.rr()).S())
U = ev[:,(1,2,0)]	# order of increasing moments is z, x, y
U = dot(U, diagflat(sign(diag(U))))	# align senses
S = T.rotated(U)

# load inital potential, combine weights
K = Field(0.1719*loadmat('potentials/RWA_X_3D_0.mat')['v'], T)
K -= K.min()
wgt += exp(-K/(2*2.35**2))

# truncate S to bounds of atoms in initial and final potentials
S = wgt.sampled(S).support(cut=exp(-Kcut/(2*2.5**2)))
# trim z axis to avoid extrapolation
S = S[:,:,1:-1]

# load and rotate potentials
V = Ts.blank()
for i in range(s.size):
	vfile = loadmat('potentials/RWA_X_3D_' + str(i) + '.mat')
	V[i,:,:,:] = vfile['v']
K = Ss.sample(Field(0.1719*V,Ts))
assert not isnan(K.ordinates).any()

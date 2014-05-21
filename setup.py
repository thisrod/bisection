"""
Rotate the potential to a trap-aligned frame, find the initial ground state by imaginary time propagation, and save the resulting fields.

The trap-aligned frame is defined by the centre of mass and principle axes of the final trap, after bisection.
"""

from field import *
from scipy.io import loadmat
import sys

# everything is in dispersion units, 10^-7 m
# T is the grid supplied by Vienna, shifted to the centre of mass
# S is the computation grid, in the trap-aligned frame, covering just the volume where atoms are to be found.

Kcut = 10		# potential at which to truncate grid V
N = 7e3		# number of atoms
l = 2e-1		# timestep for imag. time integration
s = Grid.from_axes(1e-3*arange(18)/1.368e-5)	# times samples supplied

# load supplied grid
vfile = loadmat('potentials/RWA_X_3D_0.mat')
T = Grid.from_axes(*[10*vfile[q] for q in ['x', 'y', 'z']])

# load final potential, shift to avoid underflow
K = Field(0.1719*loadmat('potentials/RWA_X_3D_17.mat')['v'], T)
K -= K.min()

# shift grid origin to centre of weight
wgt = exp(-K/(2*2.88**2))
T = T.shifted((wgt.w()*wgt).S()/wgt.S())
wgt = wgt.sampled(T)

# find the principal axes and
# set U to the rotation from the trap-aligned frame to the original data frame
ew, ev = eig((wgt*wgt.ww()).S())
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
K = (s*T).blank()
for i in range(len(s)):
	vfile = loadmat('potentials/RWA_X_3D_' + str(i) + '.mat')
	K[i,:,:,:] = vfile['v']
K = 0.1719*K.sampled(s*S)
assert not isnan(K).any()

# find ground states
q = (s*S).blank()
for j in range(1):
	w = 1+0*K[j,:,:,:]; w *= sqrt(N/(w**2).S())
	print('\nSplit step')
	for i in range(100):
		print('#', end='')
		sys.stdout.flush()
		w = w.expdsq(-l)
		w *= exp(-l*(K[j,:,:,:]+1.330*abs(w)**2))
		w *= sqrt(N/(w**2).S())
	print()
	q[j,:,:,:] = w

# plot sections through centre of weight
x, y, z = S[;,0,0], S[0,:,0], S[0,0,:]
y0 = S
slice = x*Grid.from_axes([0])*z
figure()
# extent=(z[0], z[-1], x[0], x[-1])
imshow(w.sampled(slice), aspect='auto', interpolation='nearest').set_cmap('gray')
xlabel('z');  ylabel('x')
savefig('gsxz.pdf')

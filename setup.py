"""
Rotate the potential to a trap-aligned frame, find the initial ground state by imaginary time propagation, and save the resulting fields.

The trap-aligned frame is defined by the centre of mass and principle axes of the final trap, after bisection.
"""

from field import *
from scipy.io import loadmat
import sys

# everything is in dispersion units, 10^-7 m
# R is the computation grid, in the trap-aligned frame, covering just the volume where atoms are to be found.

Kcut = 10		# potential at which to truncate grid V
N = 7e3		# number of atoms
l = 2e-1		# timestep for imag. time integration

# T is the 3D grid on which the supplied potential was sampled
# 3D resampling on a 4D common space n.y.i.
vfile = loadmat('potentials/RWA_X_3D_0.mat')
T = Grid.from_axes(*[10*vfile[q] for q in ['x', 'y', 'z']])

# load final potential, shift to avoid underflow
K1 = SampledField(0.1719*loadmat('potentials/RWA_X_3D_17.mat')['v'], T)
K1 -= K1.min()

# shift T grid origin to thermal cloud centre of mass
wgt = exp(-K1/(2*2.88**2))
T = T.shifted((wgt.r()*wgt).S()/wgt.S())
wgt = wgt.sampled(T)

# find the principal axes and
# set U to the rotation from the trap-aligned frame to the original data frame
ew, ev = eig((wgt*wgt.rr()).S())
U = ev[:,(1,2,0)]	# order of increasing moments is z, x, y
U = dot(U, diagflat(sign(diag(U))))	# align senses

# load inital potential, combine weights
K0 = SampledField(0.1719*loadmat('potentials/RWA_X_3D_0.mat')['v'], T)
K0 -= K0.min()
wgt += exp(-K0/(2*2.35**2))

# set up trap frame S, and extend to 3+1D
S = T.rotated(U)
S = wgt.sampled(S).support(cut=exp(-Kcut/(2*2.5**2)))
# trim z axis to avoid extrapolation
S = S[:,:,1:-1]
S = Grid.from_axes(1e-3*arange(18)/1.368e-5)*S
s, x, y, z = S.axes

# load and rotate potentials
K = (Grid.from_axes(1e-3*arange(18)/1.368e-5)*T).blank()
for i in range(len(s)):
	vfile = loadmat('potentials/RWA_X_3D_' + str(i) + '.mat')
	K[i,:,:,:] = vfile['v']
K = 0.1719*K.sampled(S)
assert not isnan(K).any()

# find ground states
q = 0j*S.blank()
for j in range(1):
	w = 1+0*K[j,:,:,:]; w *= sqrt(N/(w**2).S())
	print('\nSplit step')
	for i in range(100):
		print('#', end=''); sys.stdout.flush()
		w = w.expdsq(-l)
		w *= exp(-l*(K[j,:,:,:]+1.330*abs(w)**2))
		w *= sqrt(N/(w**2).S())
	print()
	q[j,:,:,:] = w

# plot sections through centre of weight
# FIXME bodgy.  delta needs to be generalised
dfun = Grid.delta(0)
D = dfun*(dfun*dfun*dfun).rotated(U)
ds, dx, dy, dz = D.axes
slice = x*dy*z
figure()
(abs(w)**2).sampled(slice).section_positive()
xlabel('z');  ylabel('x')
savefig('gsxz.pdf')

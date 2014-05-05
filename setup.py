"""
Rotate the potential to a trap-aligned frame, find the initial ground state by imaginary time propagation, and save

The trap-aligned frame is defined as follows:
	* the origin is the center of mass of e^{-V/mu}, mu being around the energy quantum for the tight axis
	* coordinates are rotated around the y axis to align the z axis with the long principle axis of e^{-V/mu}
	* this makes the x axis approximate the splitting direction

the transform matrix is saved to align.npy
"""

from stdlib import *
import sys
from numpy import savez
from scipy.io import loadmat
from scipy.integrate import simps
from scipy.interpolate import interp1d, interp2d
from scipy.ndimage.interpolation import map_coordinates

# coordinates:
#
# everything is in dispersion units, 10^-7 m
# ux, uy, uz, uR are the grid supplied by Vienna
# x, y, z, R is the inscribed grid in the trap-aligned frame
# the columns of R are [x y z] for each point in the grid
# q is often used to interate over x, y, z

# times at which potential is sampled
t = 1e-3*arange(18)/2.737e-5
Kcut = 10		# potential at which to truncate grid
M = 7e3		# atoms

# load u grid and convert units
vfile = loadmat('potentials/RWA_X_3D_0.mat')
ushape = vfile['v'].shape
ux, uy, uz = [10*vfile[q].flatten() for q in ['x', 'y', 'z']]
uR = array([q.flatten() for q in meshgrid(ux, uy, uz, indexing='ij')])

# load initial and final potentials, shift to avoid underflow
uK0 = 0.1719*vfile['v']
vfile = loadmat('potentials/RWA_X_3D_17.mat')
uK1 = 0.1719*vfile['v']
Koff = uK0.min()
uK0 -= Koff;  uK1 -= Koff

# shift origin to centre of weight
wgt = exp(-uK0.flatten()/(2*2.35**2))
r0 = (uR*wgt).sum(1)/wgt.sum()
uR -= r0.reshape((3,1))

# set MI to the inertia tensor of the weights
P = -uR.reshape((3,1,-1))*uR.reshape((1,3,-1))
# hack strides to address the diagonal plane of P
Pd = ndarray(buffer=P, dtype=P.dtype, \
	shape=uR.shape, strides=((3+1)*P.strides[1], P.strides[2]))
Pd += (uR**2).sum(0)
MI = (wgt*P).sum(2)

# set U to the rotation from the trap-aligned frame to the original data frame
# hU is the equivalent transform on grid indices
ew, ev = eig(MI)
W = zeros((3,3));  W[1,0] = 1;  W[:,1] = ev[:,0];  W[0,2] = 1
U, R = qr(W)
U = U[:,argmax(abs(U), axis=1)]		# order of original axes
U = dot(U, diagflat(sign(U.sum(0))))		# sense of original axes
usteps = array(ushape)-1
h = (uR[:,-1]-uR[:,0])/usteps
hU = (1/h).reshape((3,1))*U*h.reshape((1,3))

# rotate initial and final potentials, find restricted interpolation grid
# changed to test convergence and grid
grid = indices(ushape).reshape((3,-1))
rotgrid = dot(hU, grid-usteps[:,newaxis]/2) + usteps[:,newaxis]/2
K0 = map_coordinates(uK0, rotgrid, cval=nan).reshape(ushape)
K1 = map_coordinates(uK1, rotgrid, cval=nan).reshape(ushape)
support = array(nonzero(logical_or(K0 < Kcut, K1 < Kcut)))
corner = support.min(axis=1)
extent = support.max(axis=1) - corner + 1

# trim one point from each end of z axis to avoid extrapolation
corner[2] += 1;  extent[2] -= 2;
grid = indices(extent).reshape((3,-1)) + corner[:,newaxis]
cent = (usteps - corner)[:,newaxis]/2
rotgrid = dot(hU, grid-cent) + cent

# trim axes
x, y, z = [(ux, uy, uz)[i][corner[i]:corner[i]+extent[i]] for i in range(3)]
r = array([Q.flatten() for Q in meshgrid(x, y, z, indexing='ij')])

# load and rotate potentials
K = empty((t.size,)+tuple(extent))
for i in range(t.size):
	vfile = loadmat('potentials/RWA_X_3D_' + str(i) + '.mat')
	K[i,:,:,:] = 0.1719*map_coordinates(vfile['v'], rotgrid, cval=nan).reshape(extent)
assert not isnan(K).any()
K -= Koff

# find centres of mass.  only to plot sections, so trapezoid rule suffices
wgt = exp(-K/(2*2.35**2)).reshape((t.size,1,-1))
r0 = (r[newaxis,:,:]*wgt).sum(2)/wgt.sum(2)


# Set up to find ground state
steps = extent - 1
h = (r[:,-1]-r[:,0])/steps	# grid steps
dV = prod(h)
l = 2e-1	# timestep
def kspace(q):
	return 2*pi*fftfreq(q.size, ptp(q)/(q.size-1))
kx, ky, kz = meshgrid(kspace(x), kspace(y), kspace(z), indexing='ij')
dvec = exp(-l*(kx**2+ky**2+kz**2))
def D(q):	return ifftn(dvec*fftn(q))
def N(q):	return -l*((JJ[0,:,:,:]+1.330*abs(q)**2-E)*q)
w = 1+0*K[0,:,:,:]; w *= sqrt(M/dV)/norm(w)
E = ((K[0,:,:,:]+1.330*abs(w)**2)*abs(w)**2).sum()/(abs(w)**2).sum()

# Propagate in imaginary time to find ground state
print('\nSplit step')
for i in range(100):
	print('#', end='')
	sys.stdout.flush()
	w = D(w)
	w *= exp(-l*(K[0,:,:,:]+1.330*abs(w)**2-E))
	mm = simps(simps(simps(abs(w)**2)))*dV
	w *= sqrt(M/mm)
print()

# Save potentials, centres of weight and ground state
savez('trap.npz', r0=r0, t=t, x=x, y=y, z=z, K=K, w0=w, N=M)


# plot sections through centre of weight
rcm = r0[0,:]
figure()
fy = interp1d(y-rcm[1], abs(w)**2, axis=1)
imshow(fy(0), extent=(z[0], z[-1], x[0], x[-1]), aspect='auto').set_cmap('gray')
xlabel('z');  ylabel('x')
savefig('gsxz.pdf')
figure()
fz = interp1d(z-rcm[2], abs(w)**2, axis=2)
imshow(fz(0), extent=(y[0], y[-1], x[0], x[-1]), aspect='auto').set_cmap('gray')
xlabel('y');  ylabel('x')
savefig('gsxy.pdf')

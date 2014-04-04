"""Bisect the Vienna condensate according to the GPE

Solve the scaled Gross-Pitaevskii Equation, 
	dq/ds = (-del^2 + K^2 + 8*pi*a*|q|^2)q
in imaginary time.
"""


# seem to have lost stability by truncating the trap-sampling grid

from stdlib import *
import sys
from scipy.integrate import simps
from scipy.interpolate import interp1d, interp2d

# Load trap data
tfile = load('trap.npz')
r0, w0, y, x, K, z, t, N = [tfile[s] for s in ['r0', 'w0', 'y', 'x', 'K', 'z', 't', 'N']]

r = array([Q.flatten() for Q in meshgrid(x, y, z, indexing='ij')])
h = (r[:,-1]-r[:,0])/(array(w0.shape)-1)
dV = prod(h)

def kspace(x):
	return 2*pi*fftfreq(x.size, ptp(x)/(x.size-1))

kx, ky, kz = meshgrid(kspace(x), kspace(y), kspace(z), indexing='ij')
def D(q):	return ifftn(dvec*fftn(q))
def N(q):	return -l*((K[0,:,:,:]+1.330*abs(q)**2-E)*q)

# solve GPE
l = 1		# timestep
w = w0
E = ((K[0,:,:,:]+1.330*abs(2)**2)*abs(2)**2).sum()/(abs(2)**2).sum()
nn = int(t[-1]/l)
print('Real time up to %.1f/%.1f step = %.2f' % (nn*l, t[-1], l))
dvec = exp(1j*l*(kx**2+ky**2+kz**2))
Kt = interp1d(t,K,axis=0,kind='linear')
for i in range(nn):
	print('+', end='')
	sys.stdout.flush()
	w = D(w)
	w *= exp(1j*l*(Kt(l*(i+0.5))+1.330*abs(w)**2-E))
print()

rcm = interp1d(t, r0, axis=0)(nn*l)

# plot sections
figure()
fy = interp1d(y-rcm[1], abs(w)**2, axis=1)
imshow(fy(0), extent=(z[0], z[-1], x[0], x[-1]), aspect='auto').set_cmap('gray')
xlabel('z');  ylabel('x')
title('Real time up to %.1f/%.1f step = %.2f' % (nn*l, t[-1], l))
savefig('gpexz.pdf')
figure()
fz = interp1d(z-rcm[2], abs(w)**2, axis=2)
imshow(fz(0), extent=(y[0], y[-1], x[0], x[-1]), aspect='auto').set_cmap('gray')
xlabel('y');  ylabel('x')
savefig('gpexy.pdf')



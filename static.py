"""In a constant potential, of the initial trap, check that the atomic density in the Wigner approximation remains the ground-state density.
"""

from stdlib import *
import sys
from scipy.integrate import simps
from scipy.interpolate import interp1d, interp2d

samples = 1000

# Load trap data
tfile = load('trap.npz')
r0, w0, y, x, K, z, t, M = [tfile[s] for s in ['r0', 'w0', 'y', 'x', 'K', 'z', 't', 'N']]

r = array([Q.flatten() for Q in meshgrid(x, y, z, indexing='ij')])
h = (r[:,-1]-r[:,0])/(array(w0.shape)-1)
dV = prod(h)

tsteps = concatenate(([nan], logspace(-0.5,-2.5,5)))	# nan is initial state

def kspace(q):
	return 2*pi*fftfreq(q.size, ptp(q)/(q.size-1))

kx, ky, kz = meshgrid(kspace(x), kspace(y), kspace(z), indexing='ij')
def D(q):	return ifftn(dvec*fftn(q))
def N(q):	return -l*((K[0,:,:,:]+1.330*abs(q)**2-E)*q)
# E = ((K[0,:,:,:]+1.330*abs(w0)**2)*abs(w0)**2).sum()/(abs(w0)**2).sum()
E=0

ss = logspace(0.5,4,20).astype(int)
ss = ss[ss<samples]
dcrps = [[] for s in tsteps]
dacc = zeros((tsteps.size,)+w0.shape)

for j in range(samples):

	if j in ss:
		print('Sampling %d' % j)
		dens = dacc/j - 0.5/dV
		for k in range(tsteps.size):
			dcrps[k].append(dV*norm(dens[k,:,:,:]-abs(w0)**2))

	# Set up Wigner initial state
	noise = dot([1,1j], randn(2,w0.size)) / sqrt(2)
	w1 = w0 + sqrt(1/(2*dV))*noise.reshape(w0.shape)
	dacc[0,:,:,:] += abs(w1)**2
	
	# solve GPE for each timestep
	for k in range(1,tsteps.size):
		l = tsteps[k]
		dvec = exp(1j*l*(kx**2+ky**2+kz**2))
		nn = int(30/l)
		w = w1.copy()
		for i in range(nn):
			w = D(w)
			w *= exp(1j*l*(K[0,:,:,:]+1.330*abs(w)**2-E))
		dacc[k,:,:,:] += abs(w)**2

# Plot convergence
figure()
loglog(ss, dcrps[0], '+k', ss, dcrps[1], 'vw', ss, dcrps[2], '^w', ss, dcrps[3], 'sw', \
	ss, dcrps[4], 'pw', ss, dcrps[5], '*k', ss, 100*ss**-0.5, '--k')
xlabel('samples');  ylabel('norm of density discrepency')
title('Convergence of TW statics.  Interval %.1f, %d atoms' % (nn*l, M))
legend(['initial'] + ['%.3f' % h for h in tsteps[1:]], loc='lower left')
savefig('dcrps.pdf')

rcm = r0[0,:]

# plot sections
#figure()
#fy = interp1d(y-rcm[1], dens, axis=1)
#imshow(fy(0), extent=(z[0], z[-1], x[0], x[-1]), aspect='auto', interpolation='nearest').set_cmap('gray')
#xlabel('z');  ylabel('x')
#title('Real time up to %.1f/%.1f step = %.2f' % (nn*l, t[-1], l))
#savefig('rhoxz.pdf')
#figure()
#fz = interp1d(z-rcm[2], dens, axis=2)
#imshow(fz(0), extent=(y[0], y[-1], x[0], x[-1]), aspect='auto', interpolation='nearest').set_cmap('gray')
#xlabel('y');  ylabel('x')
#savefig('rhoxy.pdf')



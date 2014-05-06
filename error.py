"""In a constant potential, of the initial trap, check that the atomic density in the Wigner approximation remains the ground-state density.
"""

from stdlib import *
import sys
from scipy.integrate import simps
from scipy.interpolate import interp1d

# samples = 150
samples = 10

# Load trap data
tfile = load('trap.npz')
r0, w0, y, x, K, z, t, M = [tfile[s] for s in ['r0', 'w0', 'y', 'x', 'K', 'z', 't', 'N']]
w0 = w0[0,:,:,:]

r = array([Q.flatten() for Q in meshgrid(x, y, z, indexing='ij')])
h = (r[:,-1]-r[:,0])/(array(w0.shape)-1)
dV = prod(h)

ls = array([0.09, 0.03])
# ss = concatenate((logspace(-1,log10(620/49),10), linspace(0,620,50)[2:]))
ss = concatenate((logspace(-1,log10(10/4),3), linspace(0,10,5)[2:]))
nss = (ss/ls[0]).astype(int)

# for each sample:
#	integrate statics with timesteps 0.03, 0.09
#	accumulate rho, rho^2, |psi1-psi2| (take small ts for rho) at log short times and lin times
#	integrate dynamics, collect same moments
# find moments and plot

# set up for split-step spectral differentiation
def kspace(q):
	return 2*pi*fftfreq(q.size, ptp(q)/(q.size-1))
kx, ky, kz = meshgrid(kspace(x), kspace(y), kspace(z), indexing='ij')
dvecs = [exp(1j*l*(kx**2+ky**2+kz**2)) for l in ls]
def D(q):	return ifftn(dvec*fftn(q))
def N(q):	return -l*((K[0,:,:,:]+1.330*abs(q)**2-E)*q)
# E = ((K[0,:,:,:]+1.330*abs(w0)**2)*abs(w0)**2).sum()/(abs(w0)**2).sum()
E=0
# nn = int(t[-1]/ls[0])
nn = int(10.1/ls[0])

# tables for data
# Can we store moments instead of fields?  Not if we compare density to initial density
# (static dynamic) (r^2 r r_longstep) ss x y z
acc = zeros((2, 3, (nss<nn).sum()) + w0.shape)

for j in range(samples):
	
	print('Sample %d' % (j+1))
	
	# Set up Wigner initial state
	noise = dot([1,1j], randn(2,w0.size)) / sqrt(2)
	w1 = w0 + sqrt(1/(2*dV))*noise.reshape(w0.shape)

	# solve static GPE for both timesteps
	w2 = w1.copy();  w3 = w1.copy(); rec = 0
	for i in range(nn):
		dvec = dvecs[0]
		w2 = D(w2)
		w2 *= exp(1j*ls[0]*(K[0,:,:,:]+1.330*abs(w2)**2-E))
		dvec = dvecs[1]
		for p in range(3):
			w3 = D(w3)
			w3 *= exp(1j*ls[1]*(K[0,:,:,:]+1.330*abs(w3)**2-E))
		if i in nss:
			acc[0,:,rec,:,:,:] = array([abs(w3)**4, abs(w3)**2, abs(w2)**2])
			rec += 1
			
	# solve dynamic GPE for both timesteps
	w2 = w1.copy();  w3 = w1.copy(); rec = 0
	for i in range(nn):
		dvec = dvecs[0]
		w2 = D(w2)
		Kt = interp1d(t, K, axis=0)
		w2 *= exp(1j*ls[0]*(Kt(ls[0]*(i+0.5))+1.330*abs(w2)**2-E))
		dvec = dvecs[1]
		for p in range(3):
			w3 = D(w3)
			w3 *= exp(1j*ls[1]*(Kt(ls[1]*(3*i+p+0.5))+1.330*abs(w3)**2-E))
		if i in nss:
			acc[1,:,rec,:,:,:] = array([abs(w3)**4, abs(w3)**2, abs(w2)**2])
			rec += 1
			

# adjust moments
# (static dynamic) (r^2 r r_longstep) ss x y z
def moments(xs):
	return sqrt(dV*xs.reshape((xs.shape[0], -1)).sum(1))
def distances(x,y):
	return moments((x-y)**2)
def fieldstds(xsq,x):
	return moments(xsq-x**2)
ss = ss[nss<nn]
acc = acc/samples
dens0 = (abs(w0)**2)[newaxis,::]
dens = abs(acc[:,1:3,::])**2 - 0.5/dV

# Plot static convergence
figure()
semilogy(ss, distances(dens[0,0,::],dens0), '+k', ss, distances(dens[0,0,::],dens[0,1,::]), 'vk', \
	ss, fieldstds(acc[0,1,::],acc[0,0,::]), '^k')
legend(['vs. GPE ground state', 'time discretisation', 'sampling'], 'right')
xlabel('s');  ylabel('misplaced atoms (7000 total)')
savefig('serror.pdf')

# Plot dynamic convergence
figure()
semilogy(ss, distances(dens[1,0,::],dens[1,1,::]), 'vk', ss, fieldstds(acc[1,1,::],acc[1,0,::]), '^k')
xlabel('s');  ylabel('misplaced atoms (7000 total)')
legend(['time discretisation', 'sampling'], 'right')
savefig('derror.pdf')

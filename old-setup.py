
# find centres of mass.  only to plot sections, so trapezoid rule suffices
wgt = exp(-K/(2*2.35**2)).reshape((t.size,1,-1))
r0 = (r[newaxis,:,:]*wgt).sum(2)/wgt.sum(2)


# Set up to find ground state at each sampled time
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

w0 = 1j*nan*K

for j in range(t.size):
	w = 1+0*K[j,:,:,:]; w *= sqrt(M/dV)/norm(w)
	E = ((K[j,:,:,:]+1.330*abs(w)**2)*abs(w)**2).sum()/(abs(w)**2).sum()
	
	# Propagate in imaginary time to find ground state
	print('\nSplit step, t = %.1f' % t[j])
	for i in range(100):
		print('#', end='')
		sys.stdout.flush()
		w = D(w)
		w *= exp(-l*(K[j,:,:,:]+1.330*abs(w)**2-E))
		mm = simps(simps(simps(abs(w)**2)))*dV
		w *= sqrt(M/mm)
	print()
	w0[j,:,:,:] = w

# Save potentials, centres of weight and ground state
savez('trap.npz', r0=r0, t=t, x=x, y=y, z=z, K=K, w0=w0, N=M)

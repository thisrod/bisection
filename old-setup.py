
# hU is the equivalent transform on grid indices

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

# Find principal axes of the trap potential

N = 7e3;		# nominal number of atoms

# load parallel cigar potential, and convert to dispersion units
load potentials/RWA_X_3D_17.mat
x = 10*x;  y = 10*y;  z = 10*z;
K = 0.1719*v;
K = K - min(min(min(K)));

# set up grid
[X0,Y0,Z0] = ndgrid(x,y,z);

# find nominal density, shift origin to centre of mass
rho = exp(-K/(2*2.88^2));
rho = rho(:);
r0 = [X0(:) Y0(:) Z0(:)]'*rho/sum(rho);
xa = x-r0(1);  ya = y - r0(2);  za = z-r0(3);
[X,Y,Z] = ndgrid(xa,ya,za);
X = X(:);  Y = Y(:);  Z = Z(:);

# plot section through centre of mass
# figure
# zmin = round(interp1(za,1:length(za),0));
# contour(xa,ya,K(:,:,zmin),linspace(0,sqrt(2),10).^2)
# xlabel x;  ylabel y;  title('Transverse section near bottom of trap, za = ' + string(za(zmin)))

# inertia matrix
I = diag(([Y X X].^2+[Z Z Y].^2)'*rho);
I([2 3 6]) = -[X.*Y  X.*Z Y.*Z]'*rho;
I([4 7 8]) = I([2 3 6]);

# principal axes, in x y z order, projected into xz plane.
[U, ~] = eig(I);
[~, k] = max(abs(U'));
U = U(:,k);  U([2 4 6 8]) = 0;
for i = 1:3
	U(:,i) = sign(U(i,i))*U(:,i)/norm(U(:,i));
end

function [x, y, z, K] = loadk(t, U, r0)
	load(['potentials/RWA_X_3D_' int2str(t) '.mat'])
	x = 10*x;  y = 10*y;  z = 10*z;
	xa = x-r0(1);  ya = y - r0(2);  za = z-r0(3);
	K = 0.1719*v;
	K = K - min(min(min(K)));
	[X, Y, Z] = ndgrid(xa, ya, za);
	R = U*[X(:) Y(:) Z(:)]';
	# interp3 takes meshgrid order, reshape takes ndgrid order
	K = reshape(interp3(ya, xa, za, K, R(2,:), R(1,:), R(3,:)), size(K));
	# trim extrapolated points from rotated grid
	x = x(4:end-3);  y = y;  z = z(2:end-1);
	K = K(4:end-3, :, 2:end-1);
endfunction
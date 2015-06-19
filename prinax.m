# Find principal axes of the trap potential

N = 7e3;		# nominal number of atoms

function [x, y, z, K] = loadk(t, U, r0)
	% load potential in dispersion um, rotate to prinaxes.
	load(['potentials/RWA_X_3D_' int2str(t) '.mat'])
	K = 17.19*v;  K = K - min(K(:));
	if nargin == 1, return, endif
	x1 = x-r0(1);  y1 = y - r0(2);  z1 = z-r0(3);
	[X, Y, Z] = ndgrid(x1, y1, z1);
	R = U*[X(:) Y(:) Z(:)]';
	# interp3 takes meshgrid order, reshape takes ndgrid order
	K = interp3(y1, x1, z1, K, R(2,:), R(1,:), R(3,:));
	K = reshape(K, size(X));
endfunction

function [r, rho] = tcent(x, y, z, K)
	% centre of mass of a thermal cloud
	[X,Y,Z] = ndgrid(x,y,z);
	rho = exp(-K*0.25^2);  rho = rho(:);
	r = [X(:) Y(:) Z(:)]'*rho/sum(rho);
endfunction

[x, y, z, K] = loadk(17);
[r0, rho] = tcent(x,y,z,K);
x1 = x-r0(1);  y1 = y - r0(2);  z1 = z-r0(3);
[X,Y,Z] = ndgrid(x1,y1,z1);
X = X(:);  Y = Y(:);  Z = Z(:);

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

# plot section through centre of mass
[~,i0] = min(abs(z1));  [X,Y] = ndgrid(x,y);
clf, hold on, axis([-4 4 100 106.5]), axis equal
contour(x, y, K(:,:,i0).', (1:9).^2/0.288.^2), colormap([0 0 0])
plot(X, Y, '+k', 'markersize', 2, 'color', [0.7 0.7 0.7])
plot(r0(1), r0(2), 'ok', 'markersize', 4)
text(0.05, 103.35, 'r_0', 'fontsize', 14)
plot(1.35+103.15i+0.288*exp(2*pi*1i*(0:0.01:1)), '-r')
text(1.4, 102.7, 'x_0', 'color', 'r', 'fontsize', 14)
xlabel x,  ylabel y
title(['Section of split trap, t = 17 ms, z = ' num2str(z(i0), 3)], 'fontsize', 14)
text(1.5, 105.9, 'K^2=1/x_0^2, 4/x_0^2, 9/x_0^2, ...', 'fontsize', 14)

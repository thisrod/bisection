function [K, x, y, z] = loadk(t, U, r0)
	% load potential in dispersion um, rotate to prinaxes.
	load(['potentials/RWA_X_3D_' int2str(t) '.mat'])
	K = 17.19*v;  K = K - min(K(:));
	if nargin == 1, return, endif
	x1 = x-r0(1);  y1 = y - r0(2);  z1 = z-r0(3);
	[X, Y, Z] = ndgrid(x1, y1, z1);
	R = U*[X(:) Y(:) Z(:)]';
	% interp3 takes meshgrid order, reshape takes ndgrid order
	K = interp3(y1, x1, z1, K, R(2,:), R(1,:), R(3,:));
	K = reshape(K, size(X));
end

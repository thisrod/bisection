% initialise principal axes of the trap potential
% run prinaxfigs to plot cross section

N = 7e3;		# nominal number of atoms

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

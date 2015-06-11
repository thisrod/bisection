% fit quadratic/quartic potentials to the Vienna trap data, for interpolation around the bottom of the trap.

% load raw data, not interpolants

prinax, [x, y, z, K] = loadk(0, U, r0);
x = x(26:end-25);  y = y(31:44);  K = K(26:end-25, 31:44, :);
[X,Y,Z] = ndgrid(x,y,z);  r = [X(:) Y(:) Z(:)]-repmat(r0', numel(X), 1);

[i,j] = ndgrid(1:3, 1:3);		% should we force U to be symmetric?
wgt = exp(-K(:)/(2*2.88^2));
Ua = ([r(:,i(:)).*r(:,j(:)) r].*repmat(wgt(:),1,12)) \ (K(:).*wgt(:));
U = reshape(Ua(1:9),3,3);
a = Ua(10:12);  r1 = -(U\a)/2;



% fit quadratic/quartic potentials to the Vienna trap data, for interpolation around the bottom of the trap.
% run qrtfitfigs to print the figures

set(0, 'defaultaxesfontsize', 14, 'defaulttextfontsize', 14)
prinax

cfs = repmat(nan, 5, 18);		% cols hold coefficients [x^4 1 x^2 y^2 z^2]

tic
for t = 0:17
	[x, y, z, K] = loadk(t);
	[r0, wgt] = tcent(x,y,z,K);
	x1 = x - r0(1);  y1 = y - r0(2);  z1 = z - r0(3);
	
	% fit symmetric quartic to potential well below second excited state
	% fit doesn't justify interpolating from nearest grid point to exact minimum
	%	or rotating grid
	[~,i1] = min(abs(x1));  [~,i2] = min(abs(y1));  [~,i3] = min(abs(z1));
	K1 = squeeze(K(:,i2,i3));
	xf = x(K1 <= 50)';
	cfs(1:3, t+1) = [xf.^4 ones(size(xf)) xf.^2] \ K1(K1 <= 50);
	
	% r expands sample points over principal axes
	% TODO filter out rows with tiny weights
	[X,Y,Z] = ndgrid(x1,y1,z1);  r = [X(:) Y(:) Z(:)]*U;
	rsdl = K(:) - [r(:,1).^4 ones(size(r(:,1))) r(:,1).^2]*cfs(1:3, t+1);
	cfs(4:5, t+1) = ([r(:,2).^2 r(:,3).^2].*[wgt wgt]) \ (rsdl.*wgt);
	
	if t == 13
		Ks = K;  Ks1 = K1;  rr = r0;  j0 = i3;
		Kf = [r(:,1).^4 ones(size(r(:,1))) r.^2]*cfs(:,t+1);
	endif
endfor
toc

% set x^4 coefficient constant
cfs(1,:) = mean(cfs(1,10:end));

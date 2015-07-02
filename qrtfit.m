% fit quadratic/quartic potentials to the Vienna trap data, for interpolation around the bottom of the trap.
% run qrtfitfigs.m to print the figures.  currently needs to be done in Matlab.

prinax

cfs = repmat(nan, 5, 18);		% cols hold coefficients [x^4 1 x^2 y^2 z^2]

tplot = [0 13 17];
for t = 0:17
	[x, y, z, K] = loadk(t);
	if t == 0
		Ks = zeros([length(tplot) length(x) length(y) 2]);
		Kf = zeros([length(tplot) length(x) length(y) 2]);
		Ks1 = zeros([length(tplot) length(x)]);
		zsec = zeros(size(tplot));
	endif
	[r0, wgt] = tcent(x,y,z,K);  win = wgt > 0.01;  wgt = wgt(win);
	x1 = x - r0(1);  y1 = y - r0(2);  z1 = z - r0(3);
	
	% fit symmetric quartic to potential well below second excited state
	% fit doesn't justify interpolating from nearest grid point to exact minimum
	%	or rotating grid
	[~,i1] = min(abs(x1));  [~,i2] = min(abs(y1));  [~,i3] = min(abs(z1));
	K1 = squeeze(K(:,i2,i3));
	xf = x(K1 <= 50)';
	cfs(1:3, t+1) = [xf.^4 ones(size(xf)) xf.^2] \ K1(K1 <= 50);
	
	% r expands sample points over principal axes
	[X,Y,Z] = ndgrid(x1,y1,z1);  r = [X(:) Y(:) Z(:)]*U;  w = r(win,:);
	rsdl = K(:)(win) - [w(:,1).^4 ones(size(w(:,1))) w(:,1).^2]*cfs(1:3, t+1);
	cfs(4:5, t+1) = ([w(:,2).^2 w(:,3).^2].*[wgt wgt]) \ (rsdl.*wgt);
	
	[~,i] = ismember(t, [0 13 17]);
	if i
		Ks(i,:,:,:) = K(:,:,[i3 70]);  Ks1(i,:) = K1;
		Kfit = [r(:,1).^4 ones(size(r(:,1))) r.^2]*cfs(:,t+1);
		Kfit = reshape(Kfit, size(K));
		Kf(i,:,:,:) = Kfit(:,:,[i3 70]);
		zsec(i) = z(i3);
	endif
endfor

% set x^4 coefficient constant
cfs(1,:) = mean(cfs(1,10:end));

save -mat qrtfitdat.mat x y zsec tplot cfs Ks Ks1 Kf

% fit quadratic/quartic potentials to the Vienna trap data, for interpolation around the bottom of the trap.  Plot contours of original and interpolated potentials.

prinax

% The ith column of cfs holds the interpolation coefficients at time
% (i-1)ms, in the order [x^4 1 x^2 y^2 z^2]
cfs = nan(5,18);

tplot = [0 15 17];

[~, x, y, z] = loadk(0);
Ks = zeros([length(tplot) length(x) length(y) 2]);
Ks1 = zeros([length(tplot) length(x)]);
zsec = zeros([length(tplot) 2]);

for t = 0:17
	K = loadk(t);
	[r0, wgt] = tcent(x,y,z,K);  
	
	% shift origin to trap centre
	x = x - r0(1);  y = y - r0(2);  z = z - r0(3);
	[~,i0] = min(abs(x));  [~,j0] = min(abs(y));  [~,k0] = min(abs(z));
	
	% Fit a quartic potential along a line that isn't quite the
	% x-axis, but is close enough for the resulting error to
	% be small compared to the residual.
	%
	% Away from the origin, the potential is actually quadratic
	% not quartic.  We fit the quartic to the potential
	% around the wells, where the gas is.
	Kx = squeeze(K(:,j0,k0));
	ix = Kx <= 50;  wix = wgt > 0.01;  wgt = wgt(wix);
	basis = [x.^4; ones(size(x)); x.^2]';
	cfs(1:3, t+1) = basis(ix,:) \ Kx(ix);
	
	% Each row of r is the coordinates of a sample point along
	% the principal axes of the trap
	[X,Y,Z] = ndgrid(x,y,z);  r = [X(:) Y(:) Z(:)]*U;  w = r(wix,:);
	
	% Fit a quadratic in y and z, using the weights returned by tcent
	rsdl = K(wix) - [w(:,1).^4 ones(size(w(:,1))) w(:,1).^2]*cfs(1:3, t+1);
	cfs(4:5, t+1) = ([w(:,2).^2 w(:,3).^2].*[wgt wgt]) \ (rsdl.*wgt);
	
	[~,i] = ismember(t, tplot);
	if i
		Ks(i,:,:,:) = K(:,:,[k0 70]);  Ks1(i,:) = Kx;
		% Kfit = [r(:,1).^4 ones(size(r(:,1))) r.^2]*cfs(:,t+1);
		% Kfit = reshape(Kfit, size(K));
		% Kf(i,:,:,:) = Kfit(:,:,[k0 70]);
		zsec(i,:) = z([k0 70]);
	end
end

xoff = [0 0 z(70)]*U;  x = [x; x+xoff(1)];

% set x^4 coefficient constant - juggle this to get plots that fit
cfs(1,:) = mean(cfs(1,15:end));

% recompute fit with constant coeff of x^4
L = 1;  N = 16;  h = 2*L/N;  
x2 = h*(1:3*N)-3*L;  y2 = h*(1:N)-L; 
Kf = zeros([length(tplot) length(x2) length(y2) 2]);
for i = 1:length(tplot)
	t = tplot(i);
	[X,Y,Z] = ndgrid(x2,y2,zsec(i,:));
	r = [X(:) Y(:) Z(:)]*U;
	Kfit = [r(:,1).^4 ones(size(r(:,1))) r.^2]*cfs(:,t+1);
	Kf(i,:,:,:) = reshape(Kfit, size(X));
end

save cfs.mat cfs

% save qrtfitdat.mat x y x2 y2 xoff zsec tplot cfs Ks Ks1 Kf

set(0, 'defaultaxesfontsize', 14, 'defaulttextfontsize', 14)

xx = linspace(-2.2, 2.2, 100)';
ctrs = (1:10).^2/0.25.^2;

for i = 1:length(tplot)

t = tplot(i);

figure, plot(x(1,:), Ks1(i,:), '.k')
axis([-2.2 2.2 0 200]), hold on
plot(xx, [xx.^4 ones(size(xx)) xx.^2]*cfs(1:3, t+1), '-k')
title(['Fit along splitting axis at t = ' int2str(t) ' ms'])
s = sprintf('K^2 = %.2f x^4 %+.2f x^2 %+.2f', cfs(1, t+1), cfs(3, t+1), cfs(2, t+1));
text(-1.5, 125, s)
xlabel x, ylabel K^2

% TODO get Matlab to put a box around the figure

for j = 1:2
	figure
	contour(x(j,:), y, squeeze(Ks(i,:,:,j))', ctrs, '-k')
	hold on
	contour(x2, y2, squeeze(Kf(i,:,:,j))', ctrs, ':k')
	% xlabel x, ylabel y
	axis equal, axis([-4 4 -3 3]), axis off
	% title(sprintf('Fit at t = %d, z = %.1f', t, zsec(i,j)))
	if i == 1
		patch([-3.8 -3.8 -2.8 -2.8], [-1 -1.15 -1.16 -1], 0.5*[1 1 1])
		text(-3.5,-0.8,'1 \mum')
	end
	if j == 1, saveTightFigure(['kmode' int2str(t) '.pdf']), end
end

end	% i = 1:length(tplot)

figure
plot(0:17, cfs(1:4,:)), hold on
line(0:17, 1e4*cfs(5,:))
title 'coefficients of quartic fit', xlabel t
legend('x^4', '1', 'x^2', 'y^2', '10^4*z^2')


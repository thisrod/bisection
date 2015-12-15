function BoseOne()                        %%name of main function

% This function computes the observables, and saves them to a file
% BoseOut_N.mat, where N is a serial number for the run of the
% simuation. These results can be plotted with xgraph BoseOut_N.mat.

% This version computes the total number observables only. The code to
% compute the linear densities is in the tag prebranch. It might need to
% be updated to the current version of XPDE.

% The code calls the trap axis x in order to plot 1D graphs along
% it. The bisection axis is called y. The graph labels match my thesis.

icond.xlabels = {'t', 'z', 'x'};
icond.klabels = {'\omega', 'k_z', 'k_x'};

load cfs		% trap potential polynomials fitted by qrtfit.m

tmax = 5;

icond.name =		'find ground state order parameter for the initial trap';
icond.dimension =	3;
icond.fields =		1; 
icond.ranges =		[tmax 200 6];
icond.c.N =		7000;	% number of trapped atoms
icond.c.y0 = 		0.235;
icond.c.x0 = 		3.3;
icond.cfs =		cfs([1 2 3 5], :);  % column i has coefs of [y^4 1 y^2 x^2] at t/ms = i-1
% points(3) is even to ensure that all particles lie on a definite side of the trap with none on the centreline.
icond.points =		[49 70 28];
icond.steps =		20*tmax;
icond.step =		@nrmstp;
icond.initial =		@(w,r) ones(size(r.x));
icond.da =		@(a,~,r) -1i*dA(a,r,r.cfs(:, 1));
icond.linear =		@(D,r) 0.5*(D.x.^2 + D.y.^2); 
icond.ensembles =	[2 1 1];

icond.file =		'BoseOne.mat';
icond.images =		{2 0 0 0 0 0 0 0 0};
icond.pdimension =	{3 1 1 1 1 1 1 1 1};

icond.observe{1} =	@(a,r)  abs(a).^2;
icond.olabels{1} =		'<|\psi|^2>';

icond.olabels{2} = 		'N';
icond.compare{2} = 		@(t,in) 7000;
icond.observe{2} = 		@(a,r) xint(abs(a).^2, r.dx, r);

icond.observe{3} =	@(a,r)  xint(abs(a).^2.*potential(r.cfs(:, 1), r), r.dx, r)/r.c.N;
icond.compare{3} = @(t, in) 9;
icond.olabels{3} =		'<K^2>/N';

icond.observe{4} =	@(a,r)  xint(0.2255*abs(a).^4, r.dx, r)/r.c.N;
icond.olabels{4} =		'<V>/N';

icond.transforms{5} = [false true true];
icond.observe{5} =	@(a,r)  xint(abs(a).^2.*(r.kx.^2 + r.ky.^2), r.dk, r)/r.c.N;
icond.compare{5} = @(t, in) 1/(2*in.c.x0.^2) + 1/(2*in.c.y0.^2);
icond.olabels{5} =		'<T>/N';

icond.olabels{6} = 		'\Delta x';
icond.observe{6} = 		@(a,r) sqrt(xint(r.y.^2.*abs(a).^2, r.dx, r)/r.c.N);
icond.compare{6} = 		@(t,in) in.c.y0/sqrt(2);

icond.olabels{7} = 		'\Delta z';
icond.observe{7} = 		@(a,r) sqrt(xint(r.x.^2.*abs(a).^2, r.dx, r)/r.c.N);
icond.compare{7} = 		@(t,in) in.c.x0/sqrt(2);

icond.olabels{8} = 		'\Delta k_x';
icond.observe{8} = 		@(a,r) sqrt(xint(r.ky.^2.*abs(a).^2, r.dk, r)/r.c.N);
icond.compare{8} = 		@(t,in) 1/(sqrt(2)*in.c.y0);

icond.olabels{9} = 		'\Delta k_z';
icond.observe{9} = 		@(a,r) sqrt(xint(r.kx.^2.*abs(a).^2, r.dk, r)/r.c.N);
icond.compare{9} = 		@(t,in) 1/(sqrt(2)*in.c.x0);

% monte = rmfield(icond, {'observe', 'pdimension'});
monte.name =		'simulate Vienna trap splitting';
monte.seed =		29101977;
monte.randoms =	[2 0];	% Re and Im
monte.transfer =	@(w,r,a0,r0) a0 + [1 1i]*w/2;
% 1.368 converts t from dispersion units to interpolate into 0:17 ms
monte.da =		@(a,~,r) dA(a, r, interp1(0:17, r.cfs', 1.368*r.t)');
monte.linear =		@(D,r) 0.5*1i*(D.x.^2 + D.y.^2); 
monte.step =		@xMP;
monte.ranges =	[17/1.368 200 6];
monte.steps =		10;

% Total number observables.

monte.observe{1} =	@(a,r)  xint((r.y<0).*(abs(a).^2-1/(2*r.dV)), r.dx, r);
monte.olabels{1} =		{'N_l'};

monte.observe{2} =	@(a,r)  xint((r.y>0).*(abs(a).^2-1/(2*r.dV)), r.dx, r);
monte.olabels{2} =		{'N_r'};

monte.observe{3} = @(a,r) ...
	( xint(sign(r.y).*(abs(a).^2-1/(2*r.dV)), r.dx, r).^2 - prod(r.points(2:end))/4 ) ...
	./ xave(xint((r.y~=0).*(abs(a).^2-1/(2*r.dV)), r.dx, r)).';

monte.olabels{3} =		{'<(N_r-N_l)^2>/<N_l+N_r>'};

monte.observe{4} =	@(a,r)  xint((r.y~=0).*(abs(a).^2-1/(2*r.dV)), r.dx, r);
monte.olabels{4} =		{'N'};

monte.pdimension =	{1 1 1 1};

xsim(icond)

end	% function BoseOne


function b = nrmstp(a,xi,dt,r)
	% normalise the order parameter for imaginary time integration
	b = xMP(a,xi,dt,r);
	s = xint(abs(b).^2, r.dx, r);
	b = sqrt(r.c.N./s).*b;
end

function da = dA(a,r,c)
	da = -0.5*1i*(potential(c,r) + 0.2255*abs(a).^2).*a;
end

function K = potential(c,r)
	% sample quartic trap potential on grid
	K = [r.y(:).^4 ones(size(r.x(:))) r.y(:).^2 r.x(:).^2]*c;
	K = reshape(K, size(r.x));
	K = min(K, 100);		% truncate unphysical part of quartic fit
end
function [e] = BoseOne()                        %%name of main function

% This script calls the trap axis x in order to plot 1D graphs along
% it. The bisection axis is called y.  The graph labels match my thesis.

gr1.xlabels = {'t', 'z', 'x'};

global cfs,  load cfs
cfs = cfs([1 2 3 5], :);  % column i has coefs of [y^4 1 y^2 x^2] at t/ms = i-1

tmax = 5;

icond.name =		'two dimensional Vienna trap, imaginary time initial state';
icond.dimension =	3;
icond.fields =		1; 
icond.ranges =		[tmax,200,6];
icond.points =		[49 70 27];
icond.order =		0;
icond.steps =		10*tmax;
icond.step =		@nrmstp;
% icond.graphs =		[0 0 0 0];
icond.initial =		@(w,r) ones(size(r.x));
icond.da =		@Da;
icond.linear =		@(D,r) 0.5*(D.x.^2 + D.y.^2); 
% icond.observe =	@(a,~,~) abs(a).^2;
icond.ensembles =	[1 1 1];

gr1.images =		[0];
gr1.olabels =		{'<|\psi|^2>'};

monte = icond;
monte.name =		'two dimensional Vienna trap splitting';
monte.seed =		29101977;
monte.transfer =		@(w,r,a0,r0) a0 + [1 1i]*w/sqrt(2);
monte.da =		@Db;
monte.linear =		@(D,r) 0.5*1i*(D.x.^2 + D.y.^2); 
monte.step =		@xMP;
monte.ranges =	[17/1.368,200,6];
monte.noises =		[2 0];	% Re and Im
monte.steps =		50;

gr2 = gr1;
gr2.images =		[5];
gr2.pdimension =	[3 2];

monte.observe{1} =	@(a,r) abs(a).^2 - 1/(2*r.dV);
gr2.olabels{1} =		{'n'};

monte.observe{2} =	@(a,r) angle(aone(a,r,+1));
gr2.olabels{2} =		{'\phi_R'};

monte.observe{3} =	@(a,r) angle(aone(a,r,-1));
gr2.olabels{3} =		{'\phi_L'};

monte.observe{4} =	@(a,r) angle(conj(aone(a,r,+1)).*aone(a,r,-1));
gr2.olabels{4} =		{'\phi_-'};

% monte.observe = @foo;
% gr2.olabels =		{'\phi_-'};

e  =  xspde({icond, monte}, {gr1, gr2});

end	% function BoseOne

function b = nrmstp(a,r,xi,dt)
	b = xMP(a,r,xi,dt);
	s = xint(abs(b).^2, r, r.dx);
	b = sqrt(7000./s).*b;
end

function f = aone(a,r,dn)
	m = r.dx;
	m(2) = 0;	% IWRT y
	f = xint(a.*(dn*r.y > 0), r, m);
end

function c = crln(o,r,in)	% 1D autocorrelation
	o = o(r.y == r.y(1));	% take a line from the plane
	n = length(o);
	c = toeplitz(o, [o(1) zeros(1, n-1)])' * o(:) * in.dx(2);
	% FIXME: repeat with correct shape
end

% FIXME D.R.Y.

function da  =  Da(a,r,~)
	global cfs
	K = [r.y(:).^4 ones(size(r.x(:))) r.y(:).^2 r.x(:).^2]*cfs(:, 1);
	K = reshape(K, size(r.x));
	K = min(K, 100);		% trim unphysical part from quartic fit
	da = -0.5*(K + 0.2255*abs(a).^2).*a;
end

function da  =  Db(a,r,~)
	global cfs
	ts = (0:17)/1.368;
	c = interp1(ts, cfs', r.t)';
	K = [r.y(:).^4 ones(size(r.x(:))) r.y(:).^2 r.x(:).^2]*c;
	K = reshape(K, size(r.x));
	K = min(K, 100);		% trim unphysical part from quartic fit
	da = -0.5*1i*(K + 0.2255*abs(a).^2).*a;
end

function [e] = BoseOne()                        %%name of main function

global cfs
load cfs		% cols hold coefficients [x^4 1 x^2 y^2 z^2]

tmax = 5;

icond.name =		'two dimensional Vienna trap, imaginary time initial state';
icond.dimension =	3;
icond.fields =		1; 
icond.ranges =		[tmax,6,200];
icond.order =		0;
icond.steps =		10*tmax;
icond.step =		@nrmstp;
icond.graphs =		[1 0 0 0];
icond.initial =		@(r,w,in) ones(size(r.x));
icond.da =		@Da;
icond.linear =		@(D,in) 0.5*(D.x.^2 + D.y.^2); 
icond.observe =	@(a,~,~) abs(a).^2;
icond.ensembles =	[1 1 1];

gr1.images =		[0];
gr1.olabels =		{'<|\psi|^2>'};

monte = icond;
monte.name =		'two dimensional Vienna trap splitting';
monte.seed =		29101977;
monte.initials =		@(r,w,in,r0,a0,in0) a0 + [1 1i]*w/sqrt(2);
monte.da =		@Db;
monte.linear =		@(D,in) 0.5*1i*(D.x.^2 + D.y.^2); 
% monte.observe =	@(a,r,in) abs(a).^2 - 1/(2*in.dV);
monte.observe =	@(a,r,in) phase(a,r,in,+1);
monte.step =		@xMP;
monte.ranges =	[17/1.368,6,200];
monte.noises =		[2 0];	% Re and Im
monte.steps =		50;

gr2 = gr1;
gr2.images =		[5];
gr2.olabels =		{'\phi_R'};

e  =  xspde({icond, monte}, {gr1, gr2});

end	% function BoseOne

function b = nrmstp(a,r,dt,dw,ft,in)
	b = xMP(a,r,dt,dw,ft,in);
	s = xint(abs(b).^2, in.dx, in);
	b = sqrt(7000./s).*b;
end

function f = phase(a,r,in,dn)
	m = in.dx;
	m(3) = 0;
	f = xint(a.*(dn*r.x > 0), m, in);
	f = angle(f);
end

function da  =  Da(a,r,dt,dw)
	global cfs
	K = [r.x(:).^4 ones(size(r.x(:))) r.x(:).^2 r.y(:).^2]*cfs([1 2 3 5], 1);
	K = reshape(K, size(r.x));
	K = min(K, 100);		% trim unphysical part from quartic fit
	da = -0.5*(K + 0.2255*abs(a).^2).*a*dt;
end

function da  =  Db(a,r,dt,dw)
	global cfs
	ts = (0:17)/1.368;
	c = interp1(ts, cfs([1 2 3 5], :)', r.t)';
	K = [r.x(:).^4 ones(size(r.x(:))) r.x(:).^2 r.y(:).^2]*c;
	K = reshape(K, size(r.x));
	K = min(K, 100);		% trim unphysical part from quartic fit
	da = -0.5*1i*(K + 0.2255*abs(a).^2).*a*dt;
end

function in = vienna(s)
% XMDS input structures and standard observables for 2D Vienna trap simulations.

load cfs		% trap potential polynomials fitted by qrtfit.m

switch s, case 'trap'

in.name = 		'Vienna trap configuration';
in.dimension = 		3;
in.xlabels = 		{'t', 'z', 'x'};
in.ranges = 		[NaN 200 6];
in.c.N = 			7000;
in.c.y0 = 			0.235;
in.c.x0 = 			3.3;
in.c.rpsn = 		0.2255;
in.c.tfs = 			0:17;
in.c.cfs = 			cfs([1 2 3 5], :);
in.c.K =			@potential;
in.fields =			1;
in.da = 			@dA;
in.linear = 		@(D,r) 0.5*1i*(D.x.^2 + D.y.^2);
in.klabels = 		{'\omega', 'k_z', 'k_x'};
in.points =		[49 70 28];

case 'initial'

in = vienna('trap');
in.name =			'Find ground state order parameter for initial trap';
in.ranges(1) =		5;
in.steps =			100;
in.step =			@nrmstp;
in.initial =			@(w,r) ones(size(r.x));
in.c.tfs =			in.c.tfs([1 end]);
in.c.cfs =			in.c.cfs(:, [1 1]);
in.da =			@(a,w,r) -1i*dA(a,w,r);
in.linear =			@(D,r) 0.5*(D.x.^2 + D.y.^2);

case 'observables'

obs = { ...
	'flag', 'olabels', 'observe', 'transforms', 'pdimension', 'description'; ...
	'K', '<K^2>/N', @(a,r) xint(abs(a).^2.*r.c.K(r), r.dx, r)/r.c.N, [], {1}, 'trap potential energy'; ...
	'V', '<V>/N', @(a,r)  xint(0.2255*abs(a).^4, r.dx, r)/r.c.N, [], 1, 'repulsion energy' ...
};
xinstrument(cell2table(obs(2:end,2:end), 'VariableNames', obs(1,2:end), 'RowNames', obs(2:end,1)));

end, end	% function vienna


function b = nrmstp(a,xi,dt,r)
	% normalise the order parameter for imaginary time integration
	b = xMP(a,xi,dt,r);
	s = xint(abs(b).^2, r.dx, r);
	b = sqrt(r.c.N./s).*b;
end

function da = dA(a,~,r)
	da = -0.5*1i*(r.c.K(r) + r.c.rpsn*abs(a).^2).*a;
end

function K = potential(r)
	% sample quartic trap potential on grid
	c = interp1(r.c.tfs, r.c.cfs', 1.368*r.t)';
	K = [r.y(:).^4 ones(size(r.x(:))) r.y(:).^2 r.x(:).^2]*c;
	K = reshape(K, size(r.x));
	K = min(K, 100);		% truncate unphysical part of quartic fit
end
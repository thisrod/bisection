function in = vienna(s, varargin)
% XMDS input structures and standard observables for 2D Vienna trap simulations.

load cfs		% trap potential polynomials fitted by qrtfit.m

switch s, case 'trap'

in.name = 		'Vienna split trap';
in.ranges = 		[NaN 200 6];
in.xlabels = 		{'t', 'z', 'x'};
in.c.N = 			7000;
in.c.y0 = 			0.235;
in.c.x0 = 			3.3;
in.a.g = 			0.2255;
in.c.tfs = 			0:17;
in.c.cfs = 			cfs([1 2 3 5], :);
in.a.K =			@potential;
in.fields =			1;
% in.da = 			@(a,~,r) -0.5*1i*(r.c.K(r) + r.c.rpsn*abs(a).^2).*a;
% in.linear = 		@(D,r) 0.5*1i*(D.x.^2 + D.y.^2);
in.klabels = 		{'\omega', 'k_z', 'k_x'};
in.points =		[49 70 28];
in.steps =			30;
in = trap(in);

case 'stub'

in = rmfield(vienna('trap'), {'da' 'linear' 'steps'});
in.points(1) = 2;  in.ranges(1) = 1;

case 'initial'

in = vienna('trap');
in = order(in, in.c.N, 0);
in.ranges(1) =		5;
in.steps =			100;

case 'observables'

obs = { ...
	'flag', 'olabels', 'observe', 'transforms', 'pdimension', 'description'; ...
	'K', '<K^2>/N', @(a,r) xint(abs(a).^2.*r.a.K(r), r.dx, r)/r.a.N, [], {1}, 'trap potential energy'; ...
	'Kx', 'K^2', @(a,r) r.a.K(r), [], [], 'trap potential function'; ...
	'V', '<V>/N', @(a,r)  xint(r.a.g*abs(a).^4, r.dx, r)/r.a.N, [], 1, 'repulsion energy'; ...
	'rshl', '\xi^{-2}', @(a,r)  xave(r.a.g*abs(a).^2), [], [], 'reciprocal square healing length'; ...
	'g2tw', '<(\psi^+)^2\psi^2>', @(a,r) abs(a).^4 - 2*abs(a).^2/r.dV + 1/(2*r.dV^2), [], 2, 'pair correlations'; ...
	'ntw', 'n', @(a,r) abs(a).^2-1/(2*r.dV), [], [], 'atom number density'
};
xinstrument(cell2table(obs(2:end,2:end), 'VariableNames', obs(1,2:end), 'RowNames', obs(2:end,1)));

end, end	% function vienna


function b = nrmstp(a,xi,r)
	% normalise the order parameter for imaginary time integration
	b = xMP(a,xi,r);
	s = xint(abs(b).^2, r.dx, r);
	b = sqrt(r.a.N./s).*b;
end

function K = potential(r)
	% sample quartic trap potential on grid
	c = interp1(r.c.tfs, r.c.cfs', 1.368*r.t)';
	K = [r.y(:).^4 ones(size(r.x(:))) r.y(:).^2 r.x(:).^2]*c;
	K = reshape(K, size(r.x));
	K = min(K, 100);		% truncate unphysical part of quartic fit
end
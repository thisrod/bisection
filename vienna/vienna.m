function in = vienna()
% XMDS input structure for 2D Vienna trap simulations.

load cfs		% trap potential polynomials fitted by qrtfit.m

in.name = 		'Vienna split trap';
in.ranges = 		[NaN 200 6];
in.xlabels = 		{'t', 'z', 'x'};
in.a.N = 			7000;
in.a.y0 = 			0.235;
in.a.x0 = 			3.3;
in.a.g = 			0.2255;
in.c.tfs = 			0:17;
in.c.cfs = 			cfs([1 2 3 5], :);
in.a.K =			@potential;
in.fields =			1;
in.klabels = 		{'\omega', 'k_z', 'k_x'};
in.points =		[49 70 28];
in.steps =			30;
in = trap(in);

end	% function vienna

function K = potential(r)
	% sample quartic trap potential on grid
	c = interp1(r.c.tfs, r.c.cfs', 1.368*r.t)';
	K = [r.y(:).^4 ones(size(r.x(:))) r.y(:).^2 r.x(:).^2]*c;
	K = reshape(K, size(r.x));
	K = min(K, 100);		% truncate unphysical part of quartic fit
end
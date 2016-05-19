function in = xinstrument(in, varargin)
%XINSTRUMENT     Common XPDE observables.
%   XINSTRUMENT prints a table with names and descriptions of standard
%   observables.
%   
%   out = XINSTRUMENT(in, O1, O2) returns the input structure with Os appended
%   to its observables.
%   
%   out = XINSTRUMENT(in, O, N) sets the corresponding element of images
%   to the integer N.
%   
%   out = XINSTRUMENT(in, O, @F(t,in)) sets a comparison function.
%   
%   out = XINSTRUMENT(T) appends the table T to the standard observables
%   stored.
%
%   This clobbers any existing observables and plotting functions.

%obs = { ...
%	'flag', 'olabels', 'observe', 'transforms', 'pdimension', 'description'; ...
%	'n', 'n', @(a,~) abs(a).^2, [], [], 'square modulus of field'; ...
%	'nk', 'nk', @(a,~) abs(a).^2, [false true], [], 'square modulus of reciprocal field'; ...
%	'Re', 'Re', @(a,~) real(a), [], [], 'real part of atomic field'; ...
%	'Im', 'Im', @(a,~) imag(a), [], [], 'imaginary part of atomic field'; ...
%	'N', 'N', @(a,r) xint(abs(a).^2, r.dx, r), [], 1, 'total atom number'; ...
%	'Teqm', 'Teqm', @(a,r) xint(r.kx.^2.*abs(a).^2, r.dk, r), [false true], 1, 'kinetic energy'; ...
%	'Keqm', 'Keqm', @(a,r) xint(r.a.K(r).*abs(a).^2, r.dx, r), [], 1, 'external potential energy'; ...
%	'Reqm', 'Reqm', @(a,r) xint(r.a.g.*abs(a).^4, r.dx, r), [], 1, 'repulsion energy'; ...
%	'healing', 'healing', @(a,r) 1/sqrt(r.a.g)./abs(a), [], [], 'healing length'; ...
%	'ntw', 'ntw', @(a,r) abs(a).^2-1/(2*r.dV), [], [], 'atom number density'; ...
%	'Ntw', 'Ntw', @(a,r) xint(abs(a).^2-1/(2*r.dV), r.dx, r), [], 1, 'total atom number'; ...
%	'g2tw', 'g2tw', @(a,r) abs(a).^4 - 2*abs(a).^2/r.dV + 1/(2*r.dV^2), [], 2, 'pair correlations'; ...
%	'x', '<x>', @(a,r) xint(r.x.*abs(a).^2, r.dx, r)./xint(abs(a).^2, r.dx, r), [], 1, 'x component of centre of mass'; ...
%	'kx', '<kx>', @(a,r) xint(r.kx.*abs(a).^2, r.dx, r)./xint(abs(a).^2, r.dx, r), [false true true], 1, 'mean kx'; ...
%	'Vx', 'V_x', @(a,r) xint(r.x.^2.*abs(a).^2, r.dx, r)./xint(abs(a).^2, r.dx, r), [], 1, 'variance along x axis'; ...
%	'Vy', 'V_y', @(a,r) xint(r.y.^2.*abs(a).^2, r.dx, r)./xint(abs(a).^2, r.dx, r), [], 1, 'variance along y axis' ...
%};
%obs = cell2table(obs(2:end,2:end), 'VariableNames', obs(1,2:end), 'RowNames', obs(2:end,1));

persistent plotfs

if nargin == 0
	disp(plotfs(:,{'description'}))
	return
end

if isempty(plotfs)
	plotfs = { ...
		'flag', 'olabel', 'plotfun', 'moments', 'pdimension', 'description'; ...
		'n', 'n', @(d,~) size(d), [], [], 'density'; ...
	};
	plotfs = cell2table(plotfs(2:end,2:end), 'VariableNames', plotfs(1,2:end), 'RowNames', plotfs(2:end,1));
end

in.observe = {@(a,~) abs(a).^2};
in.transforms = {[]};

n = 0;
for o = varargin
	o = o{1};
	if ischar(o)
		n = n + 1;
		ls = 0;
		if isfield(in, 'xlabels'), ls = numel(in.xlabels); end
		s = plotfs.olabel(o);
		if ls >= 2, s = strrep(s, 'x', in.xlabels{2}); end
		if ls >= 3, s = strrep(s, 'y', in.xlabels{3}); end
		if ls >= 4, s = strrep(s, 'z', in.xlabels{4}); end
		in.olabels(n) = s;	
		in.functions(n) = plotfs.plotfun(o);
		in.pdimension(n) = plotfs.pdimension(o);
	elseif isnumeric(o)
		in.images{n} = o;
	elseif isa(o, 'function_handle')
		in.compare{n} = o;
	end
end

end

function o = foo(a,r) 
keyboard
o = false;
end
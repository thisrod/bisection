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
%   out = XINSTRUMENT(in, O, @F) sets a comparison function.
%   
%   out = XINSTRUMENT(T) appends the table T to the standard observables
%   stored.


% TODO: substitute axes from in.xlabels

persistent obs

if isempty(obs)
	obs = { ...
		'flag', 'olabels', 'observe', 'transforms', 'pdimension', 'description'; ...
		'n', 'n', @(a,~) abs(a).^2, [], [], 'atom number density'; ...
		'R', 'Re', @(a,~) real(a), [], [], 'real part of atomic field'; ...
		'I', 'Im', @(a,~) imag(a), [], [], 'imaginary part of atomic field'; ...
		'N', 'N', @(a,r) xint(abs(a).^2, r.dx, r), [], 1, 'total atom number'; ...
		'T', '<T>', @(a,r) xint(abs(a).^2, r.dx, r), [false true true], 1, 'total kinetic energy'; ...
		'x', '<x>', @(a,r) xint(r.x.*abs(a).^2, r.dx, r)./xint(abs(a).^2, r.dx, r), [], 1, 'x component of centre of mass'; ...
		'kx', '<kx>', @(a,r) xint(r.kx.*abs(a).^2, r.dx, r)./xint(abs(a).^2, r.dx, r), [false true true], 1, 'mean kx'; ...
		'Vx', 'V_x', @(a,r) xint(r.x.^2.*abs(a).^2, r.dx, r)./xint(abs(a).^2, r.dx, r), [], 1, 'variance along x axis'; ...
		'Vy', 'V_y', @(a,r) xint(r.y.^2.*abs(a).^2, r.dx, r)./xint(abs(a).^2, r.dx, r), [], 1, 'variance along y axis' ...
	};
	obs = cell2table(obs(2:end,2:end), 'VariableNames', obs(1,2:end), 'RowNames', obs(2:end,1));
end

if nargin == 0
	disp(obs(:,{'description'}))
	return
end

if nargin == 1 && istable(in)
	obs = [obs; in];
	return
end

if isfield(in, 'observe')
	n = numel(in.observe);
else
	n = 0;
end

for o = varargin
	o = o{1};
	if ischar(o)
		n = n + 1;
		ls = 0;
		if isfield(in, 'xlabels'), ls = numel(in.xlabels); end
		s = obs.olabels(o);
		if ls >= 2, s = strrep(s, 'x', in.xlabels{2}); end
		if ls >= 3, s = strrep(s, 'y', in.xlabels{3}); end
		if ls >= 4, s = strrep(s, 'z', in.xlabels{4}); end
		in.olabels(n) = s;	
		in.observe(n) = obs.observe(o);
		in.transforms(n) = obs.transforms(o);
		in.pdimension(n) = obs.pdimension(o);
	elseif isnumeric(o)
		in.images{n} = o;
	elseif isa(o, 'function_handle')
		in.compare{n} = o;
	end
end

end
function in = trap(in)
%TRAP	Initialise an XPDE structure for an atom trap
%
% Potential defaults to free atoms

if isfield(in.a, 'gamma'), in = make(in, 'dimension', 2); end
in = offer(in, 'ranges', repmat(nan, 1, in.dimension), 'dimension');
in = make(in, 'dimension', length(in.ranges), 'ranges');
if all(isfield(in.a, {'gamma', 'N'})) && isnan(in.ranges(2))
	in.ranges(2) = in.a.N;
end


assert(in.dimension > 1 && in.dimension < 5);
in = make(in, 'fields', 1);
in = offer(in, 'a', struct);

% remaining Lieb-Linniger deductions
in.a = make(in.a, 'g', 2*in.a.gamma, 'gamma');
offer(in.a, 'N', in.ranges(2), 'gamma');

% potential defaults to free atoms
in.a = offer(in.a, 'K', @(r) zeros(r.d.a));

LAPs = {@(D,r) D.x.^2, @(D,r) D.x.^2 + D.y.^2, @(D,r) D.x.^2 + D.y.^2 + D.z.^2};
in.a.LAP = LAPs{in.dimension - 1};
in.a.NL = @(a,~,r) (r.a.K(r) + r.a.g*abs(a).^2).*a;

in.linear = @(D,r) 0.5i*in.a.LAP(D,r);
in.da = @(a,w,r) -0.5i*in.a.NL(a,w,r);


end


function in = offer(in, field, value, varargin)
	% give in.field the default value, provided that has no value already and that in has the fields named in varargin 
	if ~isfield(in, field) && ...
		(isempty(varargin) || all(isfield(in, varargin)))
		in = setfield(in, field, value);
	end
end

function in = make(in, field, value, varargin)
	% FIXME comparing floats for equality
	if ~(isempty(varargin) || all(isfield(in, varargin))), return, end
	if isfield(in, field)
		assert(getfield(in, field) == value);
	else
		in = setfield(in, field, value);
	end
end
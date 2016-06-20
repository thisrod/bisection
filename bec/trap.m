function in = trap(in)
%TRAP	Initialise an XPDE structure for an atom trap
%
% in = TRAP(in)  takes an input structure with a minimal set of atom trap parameters, and fills in the rest.
%
% The fields of in.a are as follows.  Various atom trap functions compute these and add them as annotations, not just TRAP.
%
% op	equilibrium order parameter
% T	total kinetic energy of equilibrium order parameter
% K(r)	total trapping energy of equilibrium order parameter.  This
%		is a function, with the XPDE parameter stucture as an argument.
% R	total repulsion energy of equilibrium order parameter
%
% N	number of atoms in the field
% g		coefficient of |a|^2 in GPE
% gamma	Lieb-Linniger parameter
% bew	sound wave eigenfrequencies
% U	sound wave U modes
% V	sound wave V modes
% healing	healing length
%
%    See also: EQOP, BOGS

if isfield(in.a, 'gamma'), in = ensure(in, 'dimension', 2); end
in = offer(in, 'ranges', @() NaN(1, in.dimension), 'dimension');
in = ensure(in, 'dimension', @() length(in.ranges), 'ranges');
if all(isfield(in.a, {'gamma', 'N'})) && isnan(in.ranges(2))
	in.ranges(2) = in.a.N;
end


assert(in.dimension > 1 && in.dimension < 5);
in = ensure(in, 'fields', 1);
in = offer(in, 'a', struct);

% remaining Lieb-Linniger deductions
in.a = ensure(in.a, 'g', @() 2*in.a.gamma, 'gamma');
offer(in.a, 'N', @() in.ranges(2), 'gamma');

% potential defaults to free atoms
K = @(r) zeros(r.d.a);
in.a = offer(in.a, 'K', @() K);

LAPs = {@(D,r) D.x.^2, @(D,r) D.x.^2 + D.y.^2, @(D,r) D.x.^2 + D.y.^2 + D.z.^2};
in.a.LAP = LAPs{in.dimension - 1};
in.a.NL = @(a,~,r) (r.a.K(r) + r.a.g*abs(a).^2).*a;

in.linear = @(D,r) 0.5i*in.a.LAP(D,r);
in.da = @(a,w,r) -0.5i*in.a.NL(a,w,r);


end

% value can be a thunk.

function in = offer(varargin)
	in = munge(false, varargin{:});
end

function in = ensure(varargin)
	in = munge(true, varargin{:});
end

function in = munge(x, in, field, value, varargin)
	if ~(isempty(varargin) || all(isfield(in, varargin))), return, end
	if isa(value, 'function_handle'), value = value(); end
	if isfield(in, field)
		if x, assert(getfield(in, field) == value), end
	else
		in = setfield(in, field, value);
	end
end
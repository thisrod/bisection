function out = eqop(in, debug)
%EQOP     annotate input structure with equilibrium order parameter
%
%    in = EQOP(in) runs xsim to compute an equilibrium order parameter, and returns a structure annotated with op, T, K, R, and healing.
%
%    in = EQOP(in, 'debug') plots stuff that helps for debugging
%
%    See also: TRAP

debug = (nargin == 2 && strcmp(debug, 'debug'));

% data = {'T', 'K', 'R', 'healing'};
data = {'Re', 'Im'};	% stubs before observables are implemented

out = in;
in = order(in);
in.raw = true;

obs = data;
if debug, obs = [obs 'N' 'n' 5]; end
in = xinstrument(in, obs{:});

[~, in, rslt, raw] = xsim(in);
if debug, xgraph(rslt, in); end

rslt = rslt{1};

a = squeeze(raw{1,2}(1,1,end,:));  a = a(:);  out.a.op = a;
for o = data
	o = o{:};
	r = rslt{find(strcmp(in.olabels, o))};
	out.a = setfield(out.a, o, squeeze(r(1,1,end,:)));
end

% out.dV = in.dV;		% bogs used to need this

end
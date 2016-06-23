function out = eqop(in, debug)
%EQOP     annotate input structure with equilibrium order parameter
%
%    in = EQOP(in) runs xsim to compute an equilibrium order parameter, and returns a structure annotated with op, T, K, and R.
%
%    in = EQOP(in, 'debug') plots stuff that helps for debugging
%
%    See also: TRAP

debug = (nargin == 2 && strcmp(debug, 'debug'));

if in.ensembles(1) > 1
	warning(['Computing a ground state order parameter with an array of %d samples.  ' ...
		'This is likely to get messy.\n'], in.ensembles(1))
end

data = {'Teqm', 'Keqm', 'Reqm'};

out = in;
in = order(in);
in.raw = true;

obs = data;
if debug, obs = [obs 'N' 'n' 'nk']; end
in = xinstrument(in, obs{:});

[~, in, rslt, raw] = xsim(in);
if debug, xgraph(rslt, in); end

rslt = rslt{1};

a = squeeze(raw{1,2}(1,1,end,:));  a = a(:);  out.a.op = a;
for o = data
	o = o{:};
	r = rslt{find(strcmp(in.olabels, o))};
	out.a = setfield(out.a, o, squeeze(r(1,1,end,1)));
end

% out.dV = in.dV;		% bogs used to need this

end
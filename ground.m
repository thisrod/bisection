function in = ground(in)
%BDG     XPDE structure to draw Wigner samples of the ground state
%
%   out = GROUND(in) takes an XPDE input structure in of the form returned by
%   STATIC.  The ranges(1) and points(1) fields of this structure should be set by the caller, to values that are suitable for computing an equilibrium order parameter by imaginary time integration.  It computes the equilibrium order parameter and the sound wave modes.  The structure out is annotated with these, and its inital function generates samples of the Bogoliubov ground state according to
%   in.ensembles and in.a.T.  Finite temperature is NYI.

% in.a.T is the condensate temperature, expressed as a square wavenumber in dispersion units.

% Find equilibrium order parameter

gsop = order(in);
gsop.raw = true;
gsop = xinstrument(gsop, 'n', 'N', 'Kx');

[~, out, rslt, raw] = xsim(gsop);  rslt = rslt{1};
a = squeeze(raw{1,1,2}(:,:,end));  a = a(:);  in.a.op = a;
kix = find(strcmp(gsop.olabels, 'K^2'));  K = squeeze(rslt{kix}(1,end,:));

[in.a.bew, in.a.U, in.a.V] = buv(out,K,a);

in.noises = 2;  in.initial = @init;
% set ranges(1) and points(1) here for null propagation

end


function a = init(w,r)
	% strategy: project out the noise in each sound wave mode, add it back in with the right amplitude.  Consistency check: at high energy, the sound wave modes become the same as Hartree modes.
	
	% this method will work for a free condensate, where the sound wave modes are L^2 orthogonal.  Need to think about how an order parameter expands as components of trapped sound wave modes.
	
	U = r.a.U;  V = r.a.V;
	c = r.dV*sum(abs(U).^2);  c = 1./sqrt(c);
	Unorm = U.*repmat(c,r.points(2),1);
	b = Unorm'*([1 1i]*w/2).';
	a = U*b - V*conj(b);
end
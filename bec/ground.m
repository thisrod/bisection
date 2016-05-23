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
gsop = xinstrument(gsop, 'n', 'N', 'Kx', 'rshl');

[~, out, rslt, raw] = xsim(gsop);  rslt = rslt{1};
a = squeeze(raw{1,1,2}(:,:,end));  a = a(:);  in.a.op = a;
kix = find(strcmp(gsop.olabels, 'K^2'));  K = squeeze(rslt{kix}(1,end,:));
% Beware: xave fakes its answer to work around the grid length design flaw,
% so xave(uniform_value) is not equal to the uniform value!
in.a.healing = 1/sqrt(in.a.g*a(1)^2);

% modify buv so that the modes are scaled already
[in.a.bew, in.a.U, in.a.V] = buv(out,K,a);

if norm(in.a.U(:,1)) > 0.1*norm(in.a.op)
	warning(sprintf('The mean field approximation is dodgy: a normalised Bogoliubov mode has %.1e particles, but the order parameter has only %.1e particles.\n', out.dV*norm(in.a.U(:,1))^2, out.dV*norm(in.a.op)^2))
end

in.initial = @init;
% set ranges(1) and points(1) here for null propagation

end


function a = init(~,r)
	% new method.  Draw normal deviates and ignore w.
	
	modes = size(r.a.U,2);
	a = nan(r.nlattice,1);
	for i = 1:r.ensembles(1)
		w = randn(modes,2)*[1; 1i]/2;
		a((i-1)*r.nspace+1:i*r.nspace) = r.a.U*w - r.a.V*conj(w);
	end
end
function in = ground(in)
%BOGS     XPDE structure to draw Wigner samples of the ground state
%
%   out = BOGS(in) takes an XPDE input structure in of the form returned by
%   EQOP, for a trap with a constant potential as constructed by STATIC.  The ranges(1) and points(1) fields of this structure should be set by the caller, to values that are suitable for computing an equilibrium order parameter by imaginary time integration.  It computes the equilibrium order parameter and the sound wave modes.  The structure out is annotated with these, and its inital function generates samples of the Bogoliubov ground state according to
%   in.ensembles and in.a.T.  Finite temperature is NYI.

% in.a.T is the condensate temperature, expressed as a square wavenumber in dispersion units.

% [in.a.bew, in.a.U, in.a.V] = buv(in,in.a.Kludge,in.a.op);

% TODO: check that grid spacing is less than healing length, so that we're including all the real particles added by the Bogoliubov transformation.

in.grid = @buv;
in.initial = @init;
% set ranges(1) and points(1) here for null propagation

end


function a = init(~,r)
	% new method.  Draw normal deviates and ignore w.
	
	modes = size(r.a.U,2);
	a = nan(r.nlattice,1);
	for i = 1:r.ensembles(1)
		w = randn(modes,2)*[1; 1i]/2;
		a((i-1)*r.nspace+1:i*r.nspace) = r.a.op + r.a.U*w - r.a.V*conj(w);
	end
end
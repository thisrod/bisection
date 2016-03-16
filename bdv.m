function out = bdv(in,t)
%BDV     XPDE input structure for Bogoliubov ground states
% For now, set the field to the lagrangian of the input field.

out.name = 'Bogoliubov initial state';
for s = {'dimension', 'fields', 'ranges', 'c', 'points'}
	out.(s{1}) = in.(s{1});
end

out.points(1) = 2;  out.ranges(1) = 1;
out.transfer = @(w,r,a0,r0) tfr(t,w,r,a0,r0);

end

function a = tfr(t,w,r,a0,r0)
	% compute the Bogoliubov modes
	% this only works for one sample
	r.t = t;
	T = 1.44;  K = 16.85;  V = 19.7;
	% kludge, read from graphs.
	M = T + K + V;
	p = r.points;  g = r.ranges;
	Dxx = kron(eye(p(3)), ssd(p(2), g(2)));
	Dyy = kron(ssd(p(3), g(3)), eye(p(2)));
	LAP = Dxx + Dyy;
	NL = r.c.K(r) + 2*r.c.rpsn*abs(a0).^2 - M;
	NL = diag(NL(:));
	BdG = [-LAP+NL, -r.c.rpsn*diag(abs(a0).^2);
		r.c.rpsn*diag(abs(a0).^2), LAP-NL];
	tic; [ev,ew] = eig(BdG, 'vector'); toc

	save debug.mat
	a = LAP*a0.';
	a = a.';
end
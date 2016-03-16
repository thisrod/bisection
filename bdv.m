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

function b = tfr(t,w,r,a0,r0)
	% compute the Bogoliubov modes
	r.t = t;
	a = reshape(a0, [r.fields r.ensembles(1) r.points(2:end)]);
	a = squeeze(a(1,1,:,:));  a = a(:);

	p = r.points;  g = r.ranges;
	Dxx = kron(eye(p(3)), ssd(p(2), g(2)));
	Dyy = kron(ssd(p(3), g(3)), eye(p(2)));
	LAP = Dxx + Dyy;
	mu = 6.75 + 16.6 + 19.5;
	% kludge, read from graphs.
	size(r.c.K(r)), size(abs(a).^2)
	K = diag(r.c.K(r)) - mu + diag(2*r.c.rpsn*abs(a).^2);
	BdG = [-LAP+K, -r.c.rpsn*diag(a.^2);
		r.c.rpsn*diag(conj(a).^2), LAP-K];
	% project onto space orthogonal to a0
	U = qr([a eye(numel(a), numel(a)-1)]);  U = U(:, 2:end);
	U = kron(eye(2), U);
	tic; [ev,ew] = eig(U'*BdG*U, 'vector'); toc

	save debug.mat
	b = LAP*a;
end
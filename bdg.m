function out = bdg(in,t)
%BDG     XPDE input structure for Bogoliubov ground states
% For now, set the field to the lagrangian of the input field.

out.name = 'Bogoliubov initial state';
for s = {'dimension', 'fields', 'ranges', 'c', 'a', 'points'}
	if isfield(in, s{1}), out.(s{1}) = in.(s{1}); end
end

% TODO: get rid of t arguments, call static if K uses r.t

out.points(1) = 2;  out.ranges(1) = 1;
out.transfer = @(w,r,a0,r0) tfr(t,w,r,a0,r0);

end

function b = tfr(t,w,r,a0,r0)
	% compute Bogoliubov modes
	r.t = t;
	a = squeeze(a0);  a = a(:);

% Stub of 2D version
%	p = r.points;  g = r.ranges;
%	Dxx = kron(eye(p(3)), ssd(p(2), g(2)));
%	Dyy = kron(ssd(p(3), g(3)), eye(p(2)));
%	LAP = Dxx + Dyy;
%	mu = 6.75 + 16.6 + 19.5;
	mu = 1971.4285;
	% kludge, read from graph of V.
	LAP = ssd(r.points(2), r.ranges(2)/2);
	K = diag(r.a.K(r) - mu) + diag(2*r.a.g*abs(a).^2);
	BdG = [-LAP+K, -r.a.g*diag(a.^2);
		r.a.g*diag(conj(a).^2), LAP-K];
	% project onto space orthogonal to a0
	[U1,~] = qr([a eye(numel(a), numel(a)-1)]);  U1 = U1(:, 2:end);
	U = kron(eye(2), U1);
	tic; [ev,ew] = eig(U'*BdG*U, 'vector'); toc

	save debug.mat
	b = LAP*a;
end
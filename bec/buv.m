function r = buv(r)
%BUV     Add Bogoliubov sound wave modes to the XPDE parameter structure
%
%   r = BDG(r) takes an XPDE parameter structure, for an input structure of the form returned by
%   STATIC and EQOP.  It sets r.a.ew, r.a.U and r.a.V to the BdG frequencies and modes

a = r.a.op;
K = r.a.K(r);  K = K(1:r.nspace);
r = xgrid(r);	% get XPDE data

%   The missing mode is the ground state.
nm = r.nspace - 1;

mu = (r.a.Teqm + r.a.Keqm + r.a.Reqm)/r.a.N;

% form an orthonormal basis for the space orthogonal to a from
% columns 2:end of the Householder reflector that takes a to e1
opmode = a/norm(a);
R = eye(r.nspace);
v = R(:,1);  if sign(opmode(1)) < 0, v = -v; end
v = v + opmode;
R = R - 2*v*v'/(v'*v);  R = R(:, 2:end);

% Solve equation 4 of prl-78-1842

H0 = -ssd(r, 'lap') + diag(K) + r.a.g*diag(abs(a).^2) - mu*eye(r.nspace);
BdG = H0^2 + 2*r.a.g*diag(abs(a).^2)*H0;

[ev,ew] = eig(R'*BdG*R, 'vector');

if all(neareal(ew))
	ew = real(ew);
else
	error('Complex Bogoliubov mode frequencies')
end

% sort by increasing eigenvalue
[ew, i] = sort(ew);  ev = ev(:,i);

% get mode freqencies
ew = sqrt(ew);

% Fix numerical quirk where some modes are elliptically polarised
% N.B. this works so far, but might not be a general solution
BP = R*ev;  BM = H0*BP./repmat(ew',r.nspace,1);
U = (BP+BM)/2;  V = (BP-BM)/2;

keyboard

ixodd = find(sqrt(sum(imag(ev).^2)) > 0.1);
for i = 1:2:length(ixodd)
	ix = ixodd(i);
	B(:,ix:ix+1) = [real(B(:,ix)) imag(B(:,ix))];
	B(:,ix) = B(:,ix) / norm(B(:,ix));
	B(:,ix+1) = B(:,ix+1) / norm(B(:,ix+1));
end

% separate u and v modes, then normalise
modes = reshape(B,r.nspace,2,[]);
c = r.dV*sum(abs(modes(:,1,:)).^2 - abs(modes(:,2,:)).^2);
c = 1./sqrt(c);
modes = modes.*repmat(c,r.nspace,2,1);

U = squeeze(modes(:,1,:));  V = squeeze(modes(:,2,:));

assert(all(neareal(U(:))), 'Complex Bogoliubov mode functions')
assert(all(neareal(V(:))), 'Complex Bogoliubov mode functions')

if norm(U(:,1)) > 0.1*norm(a)
	warning(sprintf('The mean field approximation is dodgy: a normalised Bogoliubov mode has %.1e particles, but the order parameter has only %.1e particles.\n', r.dV*norm(U(:,1))^2, r.dV*norm(a)^2))
end

save buvdebug.mat ew U V
r.a.ew = ew;  r.a.U = U;  r.a.V = V;

end
	
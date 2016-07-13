function r = buv(r)
%BUV     Add Bogoliubov sound wave modes to the XPDE parameter structure
%
%   r = BDG(r) takes an XPDE parameter structure, for an input structure of the form returned by
%   STATIC and EQOP.  It sets r.a.ew, r.a.U and r.a.V to the BdG frequencies and modes

a = r.a.op;
K = r.a.K(r);  K = K(1:r.nspace);
r = xgrid(r);	% get XPDE data

%   There are nspace-1 modes, and the missing one is the ground state.

mu = (r.a.Teqm + r.a.Keqm + r.a.Reqm)/r.a.N;
M = mu*eye(r.nspace);
Bother = diag(r.a.g*abs(a).^2);
Bself = -ssd(r, 'lap') + diag(K) + 2*Bother;
BdG = [Bself-M, -Bother; Bother, -Bself+M];
% project onto space orthogonal to a0
% This could go wrong if, for some i, A(:,1:i) is almost rank-deficient.  Apparently a single Householder reflection is the right way to do it.  Cross that bridge if we come to it. 
[U1,~] = qr([a eye(numel(a), numel(a)-1)]);  U1 = U1(:, 2:end);
% FIXME the variable U is also used for the BdG modes
U = kron(eye(2), U1);
tic; [ev,ew] = eig(U'*BdG*U, 'vector'); toc

assert(all(neareal(ew)), 'Complex Bogoliubov mode frequencies')

% remove (u,v,e) (v,u,-e) degeneracy and sort by increasing eigenvalue
ev = ev(:, ew>0);  ew = ew(ew>0);
[ew, i] = sort(ew);  ev = ev(:,i);

% Fix numerical quirk where some modes are elliptically polarised
% N.B. this works so far, but might not be a general solution
B = U*ev;
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
	
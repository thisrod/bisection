function [ew, U, V] = buv(in, K, a)
%BUV     Find Bogoliubov sound wave modes
%
%   [ew, U, V] = BDG(in, a) takes an XPDE input structure in of the form returned by
%   STATIC; this must sample the trap density observable from xinstrument.  It also takes and the data value returned by xsim(order(in)).  It returns a vector of eigenvalues E, and matrices U and V with the corresponding mode functions as columns.

%   There are nspace-1 modes, and the missing one is the ground state.

nspace = prod(in.points(2:end));

mu = 1/in.a.healing^2;	% implementation restriction: free gas
M = mu*eye(nspace);
Bother = diag(in.a.g*abs(a).^2);
Bself = -ssd(in, 'lap') + diag(K) + 2*Bother;
BdG = [Bself-M, -Bother; Bother, -Bself+M];
% project onto space orthogonal to a0
% This could go wrong if, for some i, A(:,1:i) is almost rank-deficient.  Apparently a single Householder reflection is the right way to do it.  Cross that bridge if we come to it. 
[U1,~] = qr([a eye(numel(a), numel(a)-1)]);  U1 = U1(:, 2:end);
U = kron(eye(2), U1);
tic; [ev,ew] = eig(U'*BdG*U, 'vector'); toc

save debug1.mat

assert(norm(imag(ew)) / norm(ew) < 1e-5), ew = real(ew);

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
modes = reshape(B,nspace,2,[]);
c = in.dV*sum(abs(modes(:,1,:)).^2 - abs(modes(:,2,:)).^2);
c = 1./sqrt(c);
modes = modes.*repmat(c,nspace,2,1);

U = squeeze(modes(:,1,:));  V = squeeze(modes(:,2,:));
if norm(U(:,1)) > 0.1*norm(a)
	warning(sprintf('The mean field approximation is dodgy: a normalised Bogoliubov mode has %.1e particles, but the order parameter has only %.1e particles.\n', in.dV*norm(in.a.U(:,1))^2, in.dV*norm(in.a.op)^2))
end

end

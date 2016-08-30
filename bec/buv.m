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

BP = R*ev;  BM = H0*BP./repmat(ew',r.nspace,1);
U = (BP+BM)/2;  V = (BP-BM)/2;

assert(all(neareal(U(:))), 'Complex Bogoliubov mode functions')
assert(all(neareal(V(:))), 'Complex Bogoliubov mode functions')

% where modes are degeneate, impose conditions to diagonalise H

keys = round(ew, 10);		% identify ews that are within 1e-10

F = r.dV*(U'*U - V'*V);  G = r.dV*(U'*V-V'*U);
Un = r.dV*sum(abs(U).^2);  Vn = r.dV*sum(abs(V).^2);

UU = nan(size(U));  VV = nan(size(V));

for key = unique(keys)'
	ixs = find(keys == key);
	if numel(ixs) == 1
		UU(:,ixs) = U(:,ixs);  VV(:,ixs) = V(:,ixs);
		continue
	end
	assert(numel(ixs) == 2);
	j = ixs(1);  k = ixs(2);
	if r.dV*abs(G(j,k)) < 1e-6	% close enough to zero?
		m = 0;  l = -F(j,k)/F(k,k);
	else
		numerator = -F(j,k) + sqrt(F(j,k)^2-(Un(j)-Vn(j))*(Un(k)-Vn(k)));
		l = numerator/(Un(k)-Vn(k));
		m = numerator/(Un(j)-Vn(j));
	end

	% in the untrapped gas, mm will turn out zero
	ll(j) = l;  ll(k) = l;  mm(j) = m;  mm(k) = m;

	UU(:,j) = U(:,j) + l*U(:,k);  VV(:,j) = V(:,j) + l*V(:,k);  
	UU(:,k) = U(:,k) + m*U(:,j);  VV(:,k) = V(:,k) + m*V(:,j);  
end

% normalise modes to have boson commutation relations

UUn = r.dV*sum(abs(UU).^2);  VVn = r.dV*sum(abs(VV).^2);

c = UUn - VVn;
c = 1./sqrt(c);
UU = UU.*repmat(c,r.nspace,1);  VV = VV.*repmat(c,r.nspace,1);  

% keyboard	% run borth checks here

if norm(UU(:,1)) > 0.1*norm(a)
	warning(sprintf('The mean field approximation is dodgy: a normalised Bogoliubov mode has %.1e particles, but the order parameter has only %.1e particles.\n', r.dV*norm(UU(:,1))^2, r.dV*norm(a)^2))
end

r.a.ew = ew;  r.a.U = UU;  r.a.V = VV;

end
	
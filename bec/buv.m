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
W = [BP+BM; BP-BM]/2;

% where modes are degenerate, 
% and impose conditions to diagonalise H

WW = nan(size(W));
keys = round(ew, 10);		% identify ews that are within 1e-10

for key = unique(keys)'
	ixs = find(keys == key);
	
	% linearly combine re and im parts to get real modes
	w = W(:,ixs);
	[w,~,~] = svd([real(w), imag(w)], 'econ');
	w = w(:, 1:end/2);  u = w(1:end/2, :);  v = w(end/2+1:end, :);
	
	% set ll to the fraction of the pair that gets mixed into this mode
	if numel(ixs) == 1
		ll(ixs) = 0;
	elseif numel(ixs) == 2
		F = r.dV*(u'*u - v'*v);  G = r.dV*(u'*v-v'*u);
		if abs(G(1,2)) < 1e-8		% already orthogonal
			m = 0;  l = -F(1,2)/F(2,2);
		else
			ntr = -F(1,2) + sqrt(F(1,2)^2-F(1,1)*F(2,2));
			l = ntr/F(2,2);
			m = ntr/F(1,1);
		end
		ll(ixs) = [l m];  w = w*[1 m; l 1];
	else
		error 'Implementation restriction: can handle at most two degenerate modes'
	end
	
	WW(:,ixs) = w;
end

% normalise modes
UU = WW(1:end/2,:);  VV = WW(end/2+1:end, :);
c = r.dV*sum(abs(UU).^2-abs(VV).^2);
c = 1./sqrt(c);
WW = WW.*repmat(c, 2*r.nspace, 1);
UU = WW(1:end/2,:);  VV = WW(end/2+1:end, :);

% keyboard	% run borth checks here

if norm(UU(:,1)) > 0.1*norm(a)
	warning(sprintf('The mean field approximation is dodgy: a normalised Bogoliubov mode has %.1e particles, but the order parameter has only %.1e particles.\n', r.dV*norm(UU(:,1))^2, r.dV*norm(a)^2))
end

r.a.ew = ew;  r.a.U = UU;  r.a.V = VV;

end
	
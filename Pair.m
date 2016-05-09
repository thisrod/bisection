% shell to run order and buv for a one-dimensional free gas.

clear system gsop coherent 

system.name = 'free gas in one dimension, 70 modes';
system.a.gamma = 1e3;
system.a.N = 1e5;
system.points = [nan 70];
system.ensembles = [1 1 1];
system = trap(system);

gsop = order(system, 1e5);
gsop.ranges(1) = 0.5;
gsop.points(1) = 49;
gsop.steps = 30;
gsop.raw = true;
% gsop = xinstrument(gsop, 'n', 'N', @(~,in) in.a.N, 'T', 'K', 'V', @(~,in) in.a.g, 'g2');
gsop = xinstrument(gsop, 'n', 'N', 'Kx');

coherent = twop(static(system, 0));
coherent = xinstrument(coherent, 'ntw', 'g2tw');
% ground = bdg(system);


% Find equilibrium order parameter

[~, out, rslt, raw] = xsim(gsop);
rslt = rslt{1};
a = squeeze(raw{1,1,2}(:,:,end));  a = a(:);
kix = find(strcmp(gsop.olabels, 'K^2'));  K = squeeze(rslt{kix}(1,end,:));

% Find sound wave U and V modes

[ew, U, V] = buv(out,K,a);

% Plot sound wave stuff for validation

figure, plot(ew, '.k'), title('computed and expected BdG eigenvalues')
xlabel 'mode number', ylabel '\epsilon'
% w = ck, k is around n/2L, healing length is (2*gamma)^{-1/2} in LL normalisation
% FIXME include ordinary energy
n = 1:length(ew);  L = out.V;
k1 = 2*pi/L;  k = k1*n/2;  heal = 1/sqrt(2*out.a.gamma);
hold on, plot(n, k.*sqrt(heal.^-2+k.^2), '-k')
% plot(n, (n-2)*sqrt(out.a.gamma)/(4*L), '-k')

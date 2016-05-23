% shell to run order and buv for a one-dimensional free gas.

clear system gsop coherent 

system.name = 'free gas in one dimension, 70 modes';
system.a.gamma = 1e3;
system.a.N = 1e5;
system.points = [nan 70];
system.ensembles = [1 1 1];
system = trap(system);

gsop = static(system, 0);
gsop.ranges(1) = 0.5;
gsop.points(1) = 49;
gsop.steps = 30;
% gsop = xinstrument(gsop, 'n', 'N', @(~,in) in.a.N, 'T', 'K', 'V', @(~,in) in.a.g, 'g2');

coherent = twop(static(system, 0));
coherent = xinstrument(coherent, 'ntw', 'g2tw');
% ground = bdg(system);


% Find sound wave U and V modes

uvs = ground(gsop);

% Plot sound wave stuff for validation

figure, subplot 311
plot(uvs.a.bew, '.k'), title('computed and expected BdG eigenvalues')
ylabel '\epsilon'
% w = ck, k is around n/2L, healuvsg length is (2*gamma)^{-1/2} uvs LL normalisation
% FIXME include ordinary energy
n = 1:length(uvs.a.bew);
L = system.ranges(2);	% not quite, but good enough for now
k1 = 2*pi/L;  k = k1*n/2;  heal = 1/sqrt(2*system.a.gamma);
hold on, plot(n, k.*sqrt(heal.^-2+k.^2), '-k')
% plot(n, (n-2)*sqrt(out.a.gamma)/(4*L), '-k')

subplot 312, imagesc(uvs.a.U), title 'u modes', ylabel x
subplot 313, imagesc(uvs.a.V), title 'v modes', ylabel x
xlabel 'mode number'
colormap gray

% Plot sizes of modes

h = system.ranges(2) / system.points(2);
unms = h*sum(abs(uvs.a.U).^2);
vnms = h*sum(abs(uvs.a.V).^2);
figure, plot(n, unms, '.k', n, vnms, '.r'), legend u v
title 'Particle number of normalised sound wave modes'

% draw the initial state

uvs.ensembles = [70 2 1];
uvs = xinstrument(uvs, 'n');
xspde(uvs)
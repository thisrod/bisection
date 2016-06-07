% shell to run order and buv for a one-dimensional free gas.

clear system op coherent uvs wait

% Configure an untrapped 1D Bose gas with 1e5 atoms and Lieb-Linniger
% paramter 0.01.  This should be essentially coherent at zero
% temperature.
%
% The calculation uses the atomic separation as the unit of length,
% so that |\psi|=1 and the GPE right-hand side has the form 
% (D_11 + 2\gamma|\psi|^2)\psi.

system.name = 'free gas in one dimension, 70 modes';
system.a.gamma = 0.01;
system.a.N = 1e5;
system.points = [nan 70];
system.ensembles = [1 1 1];
system = trap(system);

% numerical parameters to find the inital order parameter by imaginary time

op = static(system, 0);
op.ranges(1) = 0.5;
op.points(1) = 49;
op.steps = 30;
op = eqop(op);

% simulate dynamics and look for oscillating correlations

wait = static(system, 0);
wait.ranges(1) = 12;
wait.points(1) = 49;
wait.steps = 30;
wait = xinstrument(wait, 'g2tw');

% check that g2 is 1 for a coherent order parameter

coherent = cogs(op);
coherent.ensembles = [90 1 2];
coherent = xinstrument(coherent, 'ntw', 'g2tw', @(t,in) 1);
wait.name = 'oscillating correlations from coherent initial state';
wait.ensembles = coherent.ensembles;
xspde({coherent, wait})

% configure Bogoliubov ground state sampling

uvs = bogs(op);

% Plot sound wave modes for validation

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
hold on, plot(n, 1./(2*sqrt(2)*k*heal), '-k');
title 'Particle number of normalised sound wave modes'

% sample a Bogoliubov ground state and sample g2

uvs.ensembles = [90 1 2];
uvs = xinstrument(uvs, 'N', 'Ntw', 'n', 'ntw', 'g2tw');

% run it

wait.name = 'stable correlations from bogoliubov initial state';
wait.ensembles = uvs.ensembles;
xspde({uvs, wait})


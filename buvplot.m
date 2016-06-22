load buvdebug.mat

figure, subplot 311
plot(ew, '.k'), title('computed and expected BdG eigenvalues')
ylabel '\epsilon'
% w = ck, k is around n/2L, healuvsg length is (2*gamma)^{-1/2} uvs LL normalisation
% FIXME include ordinary energy
n = 1:length(ew);
L = system.ranges(2);	% not quite, but good enough for now
k1 = 2*pi/L;  k = k1*n/2;  heal = 1/sqrt(2*system.a.gamma);
hold on, plot(n, k.*sqrt(heal.^-2+k.^2), '-k')
% plot(n, (n-2)*sqrt(out.a.gamma)/(4*L), '-k')

subplot 312, imagesc(U), title 'u modes', ylabel x
subplot 313, imagesc(V), title 'v modes', ylabel x
xlabel 'mode number'
colormap gray

% Plot sizes of modes

h = system.ranges(2) / system.points(2);
unms = h*sum(abs(U).^2);
vnms = h*sum(abs(V).^2);
figure, plot(n, unms, '.k', n, vnms, '.r'), legend u v
hold on, plot(n, 1./(2*sqrt(2)*k*heal), '-k');
title 'Particle number of normalised sound wave modes'


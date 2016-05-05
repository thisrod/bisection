% one dimensional BEC pair correlations

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
% gsop = xinstrument(gsop, 'n', 'N', @(~,in) in.a.N, 'T', 'K', 'V', @(~,in) in.a.g, 'g2');
gsop = xinstrument(gsop, 'N');

coherent = twop(static(system, 0));
coherent = xinstrument(coherent, 'ntw', 'g2tw');
ground = bdg(system);

% xspde(gsop)
[~, ~, rslt] = xspde({gsop, ground});
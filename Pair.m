% one dimensional BEC pair correlations

clear system gsop

system.name = 'free gas in one dimension, 70 modes';
system.a.gamma = 1e3;
system.a.N = 1e5;
system.points = [nan 70];
system.ensembles = 1e3;
system = trap(system);

gsop = order(system, 1e5);
gsop.ranges(1) = 0.5;
gsop.points(1) = 49;
gsop.steps = 30;
% gsop = xinstrument(gsop, 'n', 'N', @(~,in) in.a.N, 'T', 'K', 'V', @(~,in) in.a.g, 'g2');
gsop = xinstrument(gsop, 'N');

coherent = twop(static(system, 0));
coherent = xinstrument(coherent, 'N', @(~,in) in.a.N + in.points(2)/2, 'Re', 'Im');
ground = bdg(system, 0);

% xspde(gsop)
xspde({gsop, coherent})
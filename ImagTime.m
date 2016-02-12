function ImagTime()                        %%name of main function
vienna('obervables')
in = vienna('initial');
in = xinstrument(in, 'n', 5, 'K', @(t,in) 9, 'Vx');

% in.file = 'ImagTime.mat';

xspde(in)

end	% function BoseOne

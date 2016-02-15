function Kohn()                        %%name of main function
vienna('obervables')
icond = vienna('initial');
icond = xinstrument(icond, 'n', 5, 'K', @(t,in) 9, 'Vx');

kohn = vienna('static', 0);
in.name = 'dynamics of z-axis Kohn mode';
kohn.ranges(1) = 50;
kohn.transfer = @(w,r,a0,r0) exp(1i*0.01*r.x).*a0;
kohn = xinstrument(kohn, 'n', 10, 'x', 'kx');
kohn.steps = 150;

% in.file = 'ImagTime.mat';

xspde({icond, kohn})

end	% function BoseOne

function [e] = BoseOne()                        %%name of main function

global cfs
load cfs		% cols hold coefficients [x^4 1 x^2 y^2 z^2]

tmax = 5;

in.name =       'two dimensional Vienna trap, imaginary time initial state';
in.dimension =  3;                             % t,x,y
in.fields =     1; 
in.ranges =     [tmax,6,200];
% in.noises =     [2,2];                         %%xnoises, knoises per point
% in.ensembles =  [10,4,4];                      %%samples,ensembles,parallel
in.order = 0;
in.steps = 10*tmax;
in.step = @nrmstp;
in.graphs =     [1 0 0 0];
in.initial =    @(r,w,in) ones(size(r.x));
in.da  =        @Da;
in.linear =     @(D,in) 0.5*(D.x.^2 + D.y.^2); 
in.observe =    @(a,~,~) abs(a).^2;
gr.images =     [0];                     %%number of graphics images
% gr.pdimension = [4,1,1,1];                     %%maximum plot dimension
gr.olabels =    {'<|\psi|^2>'};
% gr.compares =   [1,1,1,1];                     %%comparisons?
% compare =    @Compare;                      %%Comparison handle

in2 = in;
in2.name =       'two dimensional Vienna trap splitting';
in2.da  =        @Db;
in2.linear =     @(D,in) 0.5*1i*(D.x.^2 + D.y.^2); 
in2.step = @xMP;
in2.ranges =     [17/1.368,6,200];
in2.steps = 50;
gr2 = gr;
gr.images =     [2];


e  =  xspde({in, in2},{gr, gr});                            %%Stochasic program
end                                            %%end of main function

function b = nrmstp(a,r,dt,dw,ft,in)
	b = xMP(a,r,dt,dw,ft,in);
	s = xint(abs(b).^2, in.dx, in);
	b = sqrt(7000/sum(s(:)))*b;
end

function da  =  Da(a,r,dt,dw)

global cfs

K = [r.x(:).^4 ones(size(r.x(:))) r.x(:).^2 r.y(:).^2]*cfs([1 2 3 5], 1);
K = reshape(K, size(r.x));
K = min(K, 100);		% trim unphysical part from quartic fit
da = -0.5*(K + 0.2255*abs(a).^2).*a*dt;

end

function da  =  Db(a,r,dt,dw)

global cfs
ts = (0:17)/1.368;
c = interp1(ts, cfs([1 2 3 5], :)', r.t)';
K = [r.x(:).^4 ones(size(r.x(:))) r.x(:).^2 r.y(:).^2]*c;
K = reshape(K, size(r.x));
K = min(K, 100);		% trim unphysical part from quartic fit
da = -0.5*1i*(K + 0.2255*abs(a).^2).*a*dt;

end

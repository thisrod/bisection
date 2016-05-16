function [r, rho] = tcent(x, y, z, K)
	% centre of mass of a thermal cloud
	[X,Y,Z] = ndgrid(x,y,z);
	rho = exp(-K*0.25^2);  rho = rho(:);
	r = [X(:) Y(:) Z(:)]'*rho/sum(rho);
end

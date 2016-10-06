function BdGExplicit()
% Sample a Bogoliubov ground state explicitly and compute the correlations.  This can be compared to the numerical diagonalisation in buv.m

% this converges very slowly
% note the healing length is around 4
% plan: keep sampling converged, increase whichever of N and R makes the bigger difference 

% N = 33, 70 points, 2000 samples, g2 = 9.44e-01, expected 8.90e-01
% 	g2raw = 9.64e-01, n^2 = 1.02e+00
% N = 33, 70 points, 20000 samples, g2 = 9.28e-01, expected 8.90e-01
% N = 33, 70 points, 200000 samples, g2 = 9.26e-01, expected 8.90e-01
% 	g2raw = 9.45e-01, n^2 = 1.02e+00
% N = 33, 210 points, 20000 samples, g2 = 9.13e-01, expected 8.90e-01
% N = 33, 210 points, 200000 samples, g2 = 9.14e-01, expected 8.90e-01
% 	g2raw = 9.33e-01, n^2 = 1.02e+00
% N = 100, 70 points, 20000 samples, g2 = 9.21e-01, expected 8.90e-01
% N = 100, 70 points, 200000 samples, g2 = 9.20e-01, expected 8.90e-01
% 	g2raw = 9.80e-01, n^2 = 1.06e+00
% N = 33, 700 points, 20000 samples, g2 = 9.60e-01, expected 8.90e-01
% 	g2raw = 9.76e-01, n^2 = 1.02e+00
% N = 33, 700 points, 200000 samples, g2 = 9.43e-01, expected 8.90e-01
% 	g2raw = 9.63e-01, n^2 = 1.02e+00
% N = 33, 700 points, 2000000 samples, g2 = 9.24e-01, expected 8.90e-01
% 	g2raw = 9.44e-01, n^2 = 1.02e+00

in.a.N = 33;
R = 700;
in.points = [1 R];
in.ranges = [0 in.a.N*(R-1)/R];	% correct for missing piece of circular grid
in.ensembles = [500 1 2];

% gms = 10.^(-4:1/3:-2/3);
gms = 3e-2;
gtwo = nan(2, length(gms));

for gm = gms
	system = in;
	system.a.gamma = gm;
	system = trap(system);
	system.a.healing = 1/sqrt(2*system.a.gamma);
	system.initial = @init;
	system = xinstrument(system, 'ntw', 'g2raw', 'Re', 'Im');
	system.name = sprintf('\\gamma = %.2e', gm);

	[~,~,out,~] = xspde(system);
%	[~,~,out,~] = xsim(system);
	out = out{1};
	gtwo(:, gms == gm) = [mean(out{2}(:,1,:,:))/mean(out{1}(:,1,:,:))^2 std(out{2}(:,1,:,:))];

	L = system.a.N;  n = system.points(2);
	k = 2*pi*(1:n/2)/L;  kk = system.a.healing*k;
	uu = ((kk+1./kk)./sqrt(kk.^2+2) + 1)/2;  vv = uu - 1;
	figure, plot(k, uu, '^k', k, vv, 'vk', 2*pi*[1 1]/system.a.healing, [0 100], '-r')
	title(sprintf('mode sizes, \\gamma = %.2e', gm)), ylim([0 3])
end

fprintf('N = %d, %d points, %d samples, g2 = %.2e, expected %.2e\n', ...
	 in.a.N, R, prod(in.ensembles), gtwo(1), 1-2*sqrt(gms)/pi)
fprintf('\tg2raw = %.2e, n^2 = %.2e\n', mean(out{2}(:,1,:,:)), mean(out{1}(:,1,:,:))^2)

% figure, ax = axes;
% errorbar(ax, gms, gtwo(1,:), gtwo(2,:), 'color', 'k', 'LineStyle', 'none');
% ax.XScale = 'log';
% hold on, plot(ax, gms, 1-2*sqrt(gms)/pi, '-k');

% @(t,in) 1-2*sqrt(system.a.gamma)/pi

end	% function BdGExplicit

function a = init(~,r)
	L = r.a.N;  n = r.nspace;  ens = r.ensembles(1);
	a = zeros(ens, n);
	for k = 2*pi*(1:n/2-1)/L	% two BdG modes per wave number
		z = [1 1i]*randn(2,2*ens)/2;  z = reshape(z,ens,2);
		kk = r.a.healing*k;
		uu = ((kk+1/kk)/sqrt(kk^2+2) + 1)/2;  vv = uu - 1;
		uu = sqrt(uu);  vv = sqrt(vv);
		a = a + (z(:,1)*uu-conj(z(:,1))*vv)*sin(k*r.xc{2});
		a = a + (z(:,2)*uu-conj(z(:,2))*vv)*cos(k*r.xc{2});
	end
	k = 2*pi*n/2/L;
	z = [1 1i]*randn(2,ens)/2;  z = reshape(z,ens,1);
	kk = r.a.healing*k;
	uu = ((kk+1/kk)/sqrt(kk^2+2) + 1)/2;  vv = uu - 1;
	uu = sqrt(uu);  vv = sqrt(vv);
	a = a + (z(:,1)*uu-conj(z(:,1))*vv)*(-1).^(1:n);
	a = sqrt(2/L)*a + ones(ens, n);	% order parameter
	a = reshape(a, r.d.a);
end % function init

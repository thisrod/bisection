% fit quadratic/quartic potentials to the Vienna trap data, for interpolation around the bottom of the trap.

% check units

prinax, [x, y, z, K] = loadk(13);  r0 = tcent(x,y,z,K);
x = x(26:end-25);  y = y(31:44);  K = K(26:end-25, 31:44, :);
% check which way coordinate transform goes
[X,Y,Z] = ndgrid(x,y,z);  r = [X(:)-r0(1) Y(:)-r0(2) Z(:)-r0(3)]*U';
wgt = exp(-K(:)/(2*2.88^2));
cfs = ([r(:,1).^4 r.^2].*repmat(wgt,1,4)) \ (K(:).*wgt);

Kfit = [r(:,1).^4 r.^2]*cfs;  Kfit = reshape(Kfit, size(K));
rr = sum((r').^2);

[~,i0] = min(abs(z-r0(3)));  i0 = 70;
clf
subplot(2,2,1), contour(y,x,K(:,:,i0), (0:0.2:2).^2), colormap([0 0 0])
hold on, axis equal
wgt = reshape(wgt, size(K));  contour(y,x,wgt(:,:,i0), 0:0.1:1, '-r')
r1 = r0+U*[0 0 z(i0)]';  plot(r1(2), r1(1), 'ok')
xlabel y, ylabel x
subplot(2,2,3),  contour(y,x,Kfit(:,:,i0), (0:0.2:2).^2), axis equal, colormap([0 0 0])
xlabel y, ylabel x

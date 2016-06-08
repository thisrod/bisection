function h = zplot(x, y, w, cmap)
%zplot Plot a function of a complex variable as a color image
%   zplot(x, y, w)
%   x and y are as for image.  w is a matrix of samples f(x+iy).
%   The modulus of f is plotted on the Argand plane as brightness, 
%   its argument as hue.  Colors might be excessively quantised.
%
%   TODO: describe cmap parameter
%
%   The data w are automatically reshaped to fit x and y.

if nargin < 4, cmap = phase(128); end
ff = linspace(0, 2*pi, size(cmap,1));
w = w(:);  f = mod(angle(w), 2*pi);
pixels = repmat(abs(w)/max(abs(w)), 1, 3).*interp1(ff, cmap, f);
pixels = reshape(pixels, numel(x), numel(y), 3);
image(x, y, pixels)
h = gca;
set(h, 'YDir', 'normal')
if nargout==0; clear h; end

end


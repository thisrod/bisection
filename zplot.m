function h = zplot(x, y, w)
%zplot Plot a function of a complex variable as a color image
%   zplot(x, y, w)
%   x and y are as for image.  w is a matrix of samples f(x+iy).
%   The modulus of f is plotted on the Argand plane as brightness, 
%   its argument as hue.  Colors might be excessively quantised.
%
%   The data w are automatically reshaped to fit x and y.

w = w(:);
pixels = hsv2rgb([0.5*(1+angle(w)/pi), ones(size(w)), abs(w)/max(abs(w))]);
pixels = reshape(pixels, numel(y), numel(x), 3);
h = image(x, y, pixels);
set(gca, 'YDir', 'normal')
if nargout==0; clear h; end

end


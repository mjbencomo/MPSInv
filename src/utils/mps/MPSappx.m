function [values,indices,weights] = MPSappx(grid,xc,q,s)
%MPSAPPX Narrow-support approximation of D^s delta(x-xc).
%
%   values = MPSappx(grid,xc,q,s)
%   [values,indices,weights] = MPSappx(grid,xc,q,s)
%
% grid may be a GridSpace1D object or a numeric uniform-grid vector. The
% approximation uses N=q+s consecutive points and satisfies
%
%   h*sum_j weights(j)*(x(indices(j))-xc)^k
%       = (-1)^s*s!*delta_{k,s},  k=0,...,q+s-1.
%
% q is the singular-source approximation order and s is the derivative
% order. The solve is performed in dimensionless coordinates to avoid the
% poor scaling caused by powers of h in the Vandermonde matrix.

arguments
    grid
    xc (1,1) double {mustBeFinite}
    q (1,1) double {mustBeInteger,mustBePositive}
    s (1,1) double {mustBeInteger,mustBeNonnegative}
end

if isa(grid,'GridSpace1D')
    x = grid.x.pts;
elseif isnumeric(grid) && isvector(grid)
    x = grid(:);
else
    error('MPSappx:InvalidGrid', ...
        'grid must be a GridSpace1D object or a numeric vector.');
end

x = x(:);
if numel(x) < 2
    error('MPSappx:GridTooSmall', ...
        'The grid must contain at least two points.');
end

dx = diff(x);
h = dx(1);
if h <= 0
    error('MPSappx:IncreasingGridRequired', ...
        'Grid points must be strictly increasing.');
end

scale = max([1;abs(x);abs(xc)]);
tolerance = 100*eps(scale);
if any(abs(dx-h) > tolerance)
    error('MPSappx:UniformGridRequired', ...
        'MPS approximations require a uniform grid.');
end

N = q+s;
if N > numel(x)
    error('MPSappx:InsufficientGridPoints', ...
        'The grid must contain at least q+s points.');
end

% Select the N grid points in [xc-N*h/2,xc+N*h/2), using a scaled
% tolerance to make the half-grid tie rule deterministic.
xLeft = xc-N*h/2;
j0 = find(x >= xLeft-tolerance,1,'first');
if isempty(j0) || j0+N-1 > numel(x)
    error('MPSappx:StencilOutsideGrid', ...
        'The complete q+s point MPS stencil must lie on the grid.');
end
indices = (j0:j0+N-1).';

% Dimensionless moment system. If weights=h^(-s-1)*c, then
% h*sum weights*(x-xc)^k = h^(k-s)*sum c*r^k.
r = (x(indices)-xc)/h;
A = zeros(N,N);
for k = 0:N-1
    A(k+1,:) = r.'.^k;
end
b = zeros(N,1);
b(s+1) = (-1)^s*factorial(s);
c = A\b;
weights = c/h^(s+1);

values = zeros(size(x));
values(indices) = weights;

% Preserve the orientation of a numeric input for backward compatibility.
if isnumeric(grid) && isrow(grid)
    values = values.';
end
end

% GridSpace represents a general spatial grid for a finite-difference problem.
% Provides a common interface for spatial grids of different dimensions.
%
% Properties:
%      dim - spatial dimension
%        N - array of number of grid points, for each dimension
%        h - array of grid sizes, for each dimension
%      dom - array of domains, for each dimension
%   numPts - total number of grid points over all dimensions

classdef (Abstract) GridSpace
    properties (SetAccess = protected)
        dim (1,1) double {mustBeInteger,mustBePositive} = 1
    end

    properties (Abstract, Dependent)
        N
        h
        dom
    end

    properties (Dependent)
        numPts
    end

    methods
        function value = get.numPts(obj)
            value = prod(obj.N);
        end
    end
end
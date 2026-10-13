% GridSpace represents a collection of points in physical space.
% Cartesian grids and scattered receiver sets share dim and numPts;
% spacing, coordinate grids, and mesh are specific to Cartesian grids.
%
% Properties:
%      dim - spatial dimension
%   numPts - total number of grid points over all dimensions

classdef (Abstract) GridSpace
    properties (SetAccess = protected)
        dim (1,1) double {mustBeInteger,mustBePositive} = 1
    end

    properties (Abstract, Dependent)
        numPts
    end
end

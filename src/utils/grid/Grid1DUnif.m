% Grid1DUnif represents a uniformly spaced one-dimensional grid.
% Extends Grid1D by storing the uniform grid spacing.
%
% Properties:
%   Inherited from Grid1D
%   h - uniform grid spacing

classdef (Abstract) Grid1DUnif < Grid1D
    properties (SetAccess = protected)
        h (1,1) double {mustBePositive} = 1
    end

    methods (Access = protected)
        function obj = Grid1DUnif(dom,pts,h)
            obj@Grid1D(dom,pts);
            obj.h = h;
        end
    end
end
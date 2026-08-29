% GridSpace1D represents a one-dimensional spatial grid.
%
% The coordinate grid must be either:
%   1. a uniformly spaced Grid1DUnif object, or
%   2. a Grid1D object containing exactly one point.
%
% A one-point grid is useful as an output or observation grid. Its grid
% spacing h is undefined, so it should not be used as a computational
% finite-difference grid.
%
% Properties:
%   x           - one-dimensional coordinate grid
%   N           - number of grid points
%   h           - uniform spacing; undefined for a one-point grid
%   dom         - spatial domain
%   isPointGrid - true when the grid contains exactly one point

classdef GridSpace1D < GridSpace
    properties (SetAccess = protected)
        x
    end

    properties (Dependent)
        N
        h
        dom
        isPointGrid
    end

    methods
        function obj = GridSpace1D(x)
            arguments
                x (1,1) Grid1D
            end

            % Permit uniform grids and singleton observation grids.
            if ~isa(x,'Grid1DUnif') && x.N ~= 1
                error('GridSpace1D:UnsupportedGrid', ...
                    ['The coordinate grid must be uniform or contain ' ...
                     'exactly one point.']);
            end

            obj.x = x;
            obj.dim = 1;
        end

        function value = get.N(obj)
            value = obj.x.N;
        end

        function value = get.h(obj)
            if obj.isPointGrid
                error('GridSpace1D:SpacingUndefined', ...
                    ['Grid spacing is undefined for a one-point grid. ' ...
                     'A one-point grid may be used as an output grid, ' ...
                     'but not as a computational finite-difference grid.']);
            end

            value = obj.x.h;
        end

        function value = get.dom(obj)
            value = obj.x.dom;
        end

        function value = get.isPointGrid(obj)
            value = obj.x.N == 1;
        end

        function X = mesh(obj)
            X = obj.x.pts;
        end
    end
end
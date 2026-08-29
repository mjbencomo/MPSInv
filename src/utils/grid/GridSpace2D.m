% GridSpace2D represents a Cartesian spatial grid in two dimensions.
%
% Each coordinate grid must be either:
%   1. a uniformly spaced Grid1DUnif object, or
%   2. a Grid1D object containing exactly one point.
%
% A singleton coordinate is useful for storing a solution along a line or
% at a single point. Grid spacing h is undefined when either coordinate is
% a singleton, so such grids should not be used as computational
% finite-difference grids.
%
% Examples:
%   Nx > 1, Ny > 1 : rectangular grid
%   Nx = 1, Ny > 1 : vertical line
%   Nx > 1, Ny = 1 : horizontal line
%   Nx = 1, Ny = 1 : single observation point

classdef GridSpace2D < GridSpace
    properties (SetAccess = protected)
        x
        y
    end

    properties (Dependent)
        N
        h
        dom
        singletonDimensions
        isLineGrid
        isPointGrid
    end

    methods
        function obj = GridSpace2D(x,y)
            arguments
                x (1,1) Grid1D
                y (1,1) Grid1D
            end

            obj.validateCoordinateGrid(x,'x');
            obj.validateCoordinateGrid(y,'y');

            obj.x = x;
            obj.y = y;
            obj.dim = 2;
        end

        function value = get.N(obj)
            value = [obj.x.N,obj.y.N];
        end

        function value = get.h(obj)
            if any(obj.singletonDimensions)
                error('GridSpace2D:SpacingUndefined', ...
                    ['Grid spacing is undefined when a coordinate grid ' ...
                     'contains only one point. Such grids may be used as ' ...
                     'output grids, but not as computational grids.']);
            end

            value = [obj.x.h,obj.y.h];
        end

        function value = get.dom(obj)
            value = [
                obj.x.dom
                obj.y.dom
                ];
        end

        function value = get.singletonDimensions(obj)
            value = [obj.x.N == 1,obj.y.N == 1];
        end

        function value = get.isLineGrid(obj)
            value = sum(obj.singletonDimensions) == 1;
        end

        function value = get.isPointGrid(obj)
            value = all(obj.singletonDimensions);
        end

        function [X,Y] = mesh(obj)
            [X,Y] = ndgrid(obj.x.pts,obj.y.pts);
        end
    end

    methods (Access = private)
        function validateCoordinateGrid(~,coordinate,name)
            if ~isa(coordinate,'Grid1DUnif') && coordinate.N ~= 1
                error('GridSpace2D:UnsupportedGrid', ...
                    ['The %s-coordinate grid must be uniform or contain ' ...
                     'exactly one point.'],name);
            end
        end
    end
end
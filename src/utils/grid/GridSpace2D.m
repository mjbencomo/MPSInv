% GridSpace2D represents a two-dimensional spatial grid.
% The spatial coordinates are specified by two Grid1DUnif objects.
%
% Properties:
%   x - Grid1DUnif
%   y - Grid1DUnif

classdef GridSpace2D < GridSpace
    properties (SetAccess = protected)
        x
        y
    end

    properties (Dependent)
        N
        h
        dom
    end

    methods
        function obj = GridSpace2D(x,y)
            arguments
                x (1,1) Grid1DUnif
                y (1,1) Grid1DUnif
            end

            obj.x = x;
            obj.y = y;
            obj.dim = 2;
        end

        function value = get.N(obj)
            value = [obj.x.N, obj.y.N];
        end

        function value = get.h(obj)
            value = [obj.x.h, obj.y.h];
        end

        function value = get.dom(obj)
            value = [
                obj.x.dom
                obj.y.dom
                ];
        end

        function [X,Y] = mesh(obj)
            [X,Y] = ndgrid(obj.x.pts,obj.y.pts);
        end
    end
end
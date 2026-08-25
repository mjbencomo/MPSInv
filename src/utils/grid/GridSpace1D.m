% GridSpace1D represents a one-dimensional spatial grid.
% The spatial coordinate is specified by a Grid1DUnif object.
%
% Properties:
%   x - Grid1DUnif

classdef GridSpace1D < GridSpace
    properties (SetAccess = protected)
        x
    end

    properties (Dependent)
        N
        h
        dom
    end

    methods
        function obj = GridSpace1D(x)
            arguments
                x (1,1) Grid1DUnif
            end

            obj.x = x;
            obj.dim = 1;
        end

        function value = get.N(obj)
            value = obj.x.N;
        end

        function value = get.h(obj)
            value = obj.x.h;
        end

        function value = get.dom(obj)
            value = obj.x.dom;
        end

        function X = mesh(obj)
            X = obj.x.pts;
        end
    end
end
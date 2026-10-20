% ReceiverSpace1D represents an ordered collection of receiver locations.
% points(j) is the coordinate of receiver j. Locations need not be sorted,
% uniformly spaced, or distinct: repeated points are separate channels.
% This is an observation space, not a finite-difference grid.

classdef ReceiverSpace1D < GridSpace
    properties (SetAccess = protected)
        points
    end

    properties (Dependent)
        numPts
    end

    methods
        function obj = ReceiverSpace1D(points)
            arguments
                points double {mustBeReal,mustBeFinite}
            end

            if ~isvector(points) || isempty(points)
                error('ReceiverSpace1D:InvalidPoints', ...
                    'Receiver coordinates must be a nonempty vector.');
            end

            obj.points = points(:);
            obj.dim = 1;
        end

        function value = get.numPts(obj)
            value = size(obj.points,1);
        end
    end
end

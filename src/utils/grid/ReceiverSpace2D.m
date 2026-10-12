% ReceiverSpace2D represents an ordered collection of receiver locations.
% points(j,:) = [x_j,y_j] contains the coordinates of receiver j.
% Rows need not form a Cartesian product; repeated rows are separate
% observation channels. This is not a finite-difference grid.

classdef ReceiverSpace2D < GridSpace
    properties (SetAccess = protected)
        points
    end

    properties (Dependent)
        numPts
    end

    methods
        function obj = ReceiverSpace2D(points)
            arguments
                points double {mustBeReal,mustBeFinite}
            end

            if isempty(points) || ~ismatrix(points) || size(points,2) ~= 2
                error('ReceiverSpace2D:InvalidPoints', ...
                    'Receiver coordinates must be a nonempty M-by-2 array.');
            end

            obj.points = points;
            obj.dim = 2;
        end

        function value = get.numPts(obj)
            value = size(obj.points,1);
        end
    end
end

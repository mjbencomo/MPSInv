% Grid1D represents a one-dimensional grid on a specified domain.
% 
% Properties:
%   dom - domain
%     N - Number of grid points
%   pts - grid points

classdef Grid1D
    properties (SetAccess = protected)
        dom (1,2) double = [0,1]
        N   (1,1) double {mustBeInteger, mustBePositive} = 1
        pts (:,1) double = 0
    end

    methods
        function obj = Grid1D(dom,pts)
            arguments
                dom (1,2) double
                pts double {mustBeVector}
            end
        
            pts = pts(:);

            % Validate domain
            if dom(1) >= dom(2)
                error('Grid1D:InvalidDomain', ...
                    'Domain must satisfy dom(1) < dom(2).');
            end

            % Require at least one grid point
            if isempty(pts)
                error('Grid1D:EmptyGrid', ...
                    'Grid must contain at least one point.');
            end

            % Check that all points lie in the domain
            if any(pts < dom(1) | pts > dom(2))
                error('Grid1D:PointsOutsideDomain', ...
                    'All grid points must lie inside the domain.');
            end

            % Check that points are strictly increasing
            if any(diff(pts) <= 0)
                error('Grid1D:InvalidPoints', ...
                    'Grid points must be strictly increasing.');
            end

            obj.dom = dom;
            obj.pts = pts;
            obj.N   = numel(pts);
        end
    end
end
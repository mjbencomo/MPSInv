classdef (Abstract) AcousticSolver < handle
    % AcousticSolver is the common base class for acoustic FD solvers.
    % Computational grids control time stepping; output grids control
    % where the returned solution is sampled and stored.

    properties (SetAccess = protected)
        computationalGrid
        computationalTimes
        outputGrid
        outputTimes
        dt
    end

    properties (Dependent)
        Nt
        NtOutput
    end

    methods
        function obj = AcousticSolver( ...
                computationalGrid,computationalTimes,outputGrid,outputTimes)
            arguments
                computationalGrid (1,1) GridSpace
                computationalTimes double {mustBeFinite}
                outputGrid (1,1) GridSpace
                outputTimes double {mustBeFinite}
            end

            if ~isvector(computationalTimes) || numel(computationalTimes) < 2
                error('AcousticSolver:InvalidComputationalTimes', ...
                    ['The computational time grid must contain at least ' ...
                     'two points.']);
            end

            computationalTimes = computationalTimes(:).';
            timeSteps = diff(computationalTimes);
            if any(timeSteps <= 0)
                error('AcousticSolver:InvalidComputationalTimes', ...
                    'Computational times must be strictly increasing.');
            end

            dt = timeSteps(1);
            tolerance = 100*eps(max(1,max(abs(computationalTimes))));
            if any(abs(timeSteps-dt) > tolerance)
                error('AcousticSolver:NonuniformComputationalTimes', ...
                    'This solver requires a uniform computational time grid.');
            end

            if isempty(outputTimes) || ~isvector(outputTimes)
                error('AcousticSolver:InvalidOutputTimes', ...
                    'Output times must be a nonempty vector.');
            end

            outputTimes = outputTimes(:).';
            if any(diff(outputTimes) <= 0)
                error('AcousticSolver:InvalidOutputTimes', ...
                    'Output times must be strictly increasing.');
            end

            timeTolerance = 100*eps(max(1,max(abs(computationalTimes))));
            if outputTimes(1) < computationalTimes(1)-timeTolerance || ...
                    outputTimes(end) > computationalTimes(end)+timeTolerance
                error('AcousticSolver:OutputTimesOutsideInterval', ...
                    ['Output times must lie within the computational ' ...
                     'time interval.']);
            end

            if outputGrid.dim ~= computationalGrid.dim
                error('AcousticSolver:OutputDimensionMismatch', ...
                    'Computational and output grids must have the same dimension.');
            end

            obj.computationalGrid = computationalGrid;
            obj.computationalTimes = computationalTimes;
            obj.outputGrid = outputGrid;
            obj.outputTimes = outputTimes;
            obj.dt = dt;
        end

        function value = get.Nt(obj)
            value = numel(obj.computationalTimes);
        end

        function value = get.NtOutput(obj)
            value = numel(obj.outputTimes);
        end
    end

    methods (Abstract)
        [pressure,velocity] = solve(obj,varargin)
    end
end

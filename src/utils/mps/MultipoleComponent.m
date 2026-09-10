classdef (Abstract) MultipoleComponent
    % MultipoleComponent describes the time-independent part of a
    % multipole source term.
    %
    % A component specifies source geometry and spatial discretization.
    % Concrete 1D and 2D subclasses add a time function to construct the
    % corresponding MultipoleTerm object used by the acoustic solvers.

    properties (SetAccess = private)
        approximationOrder (1,1) double = 4
    end

    methods
        function obj = MultipoleComponent(approximationOrder)
            arguments
                approximationOrder (1,1) double ...
                    {mustBeInteger,mustBePositive} = 4
            end

            obj.approximationOrder = approximationOrder;
        end
    end

    methods (Abstract)
        term = withTimeFunction(obj,timeFunction,options)
        % Combine this spatial component with a temporal function.
    end
end

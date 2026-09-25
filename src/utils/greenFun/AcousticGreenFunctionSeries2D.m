classdef AcousticGreenFunctionSeries2D < GreenFunctionSeries
    % AcousticGreenFunctionSeries2D is a collection of 2D acoustic kernels.

    methods
        function obj = AcousticGreenFunctionSeries2D( ...
                components,solver,options)
            arguments
                components
                solver
                options.OutputField (1,1) string ...
                    {mustBeMember(options.OutputField, ...
                    ["pressure","velocityX","velocityY"])} = "pressure"
            end

            components = ...
                AcousticGreenFunctionSeries2D.normalizeComponents( ...
                components);
            greenFunctions = cell(1,numel(components));
            for m = 1:numel(components)
                greenFunctions{m} = AcousticGreenFunction2D( ...
                    components{m},solver,OutputField=options.OutputField);
            end
            obj@GreenFunctionSeries(greenFunctions);
        end
    end

    methods (Static, Access = private)
        function components = normalizeComponents(components)
            if isempty(components)
                error('AcousticGreenFunctionSeries2D:EmptyComponents', ...
                    'components must contain at least one component.');
            elseif iscell(components)
                if any(~cellfun( ...
                        @(c) isa(c,'MultipoleComponent2D'),components))
                    error('AcousticGreenFunctionSeries2D:InvalidComponents', ...
                        'Every component must be a MultipoleComponent2D.');
                end
                components = reshape(components,1,[]);
            elseif isa(components,'MultipoleComponent2D')
                components = num2cell(reshape(components,1,[]));
            else
                error('AcousticGreenFunctionSeries2D:InvalidComponents', ...
                    'components must contain MultipoleComponent2D objects.');
            end
        end
    end
end

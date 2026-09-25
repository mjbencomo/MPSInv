classdef AcousticGreenFunctionSeries1D < GreenFunctionSeries
    % AcousticGreenFunctionSeries1D is a collection of 1D acoustic kernels.

    methods
        function obj = AcousticGreenFunctionSeries1D( ...
                components,solver,options)
            arguments
                components
                solver
                options.OutputField (1,1) string ...
                    {mustBeMember(options.OutputField, ...
                    ["pressure","velocity"])} = "pressure"
            end

            components = ...
                AcousticGreenFunctionSeries1D.normalizeComponents( ...
                components);
            greenFunctions = cell(1,numel(components));
            for m = 1:numel(components)
                greenFunctions{m} = AcousticGreenFunction1D( ...
                    components{m},solver,OutputField=options.OutputField);
            end
            obj@GreenFunctionSeries(greenFunctions);
        end
    end

    methods (Static, Access = private)
        function components = normalizeComponents(components)
            if isempty(components)
                error('AcousticGreenFunctionSeries1D:EmptyComponents', ...
                    'components must contain at least one component.');
            elseif iscell(components)
                if any(~cellfun( ...
                        @(c) isa(c,'MultipoleComponent1D'),components))
                    error('AcousticGreenFunctionSeries1D:InvalidComponents', ...
                        'Every component must be a MultipoleComponent1D.');
                end
                components = reshape(components,1,[]);
            elseif isa(components,'MultipoleComponent1D')
                components = num2cell(reshape(components,1,[]));
            else
                error('AcousticGreenFunctionSeries1D:InvalidComponents', ...
                    'components must contain MultipoleComponent1D objects.');
            end
        end
    end
end

classdef AcousticSource1D < handle
    % AcousticSource1D represents sums of separable 1D acoustic sources.
    %
    %   f_p(x,t) = sum_k g_{p,k}(x) w_{p,k}(t)
    %   f_v(x,t) = sum_k g_{v,k}(x) w_{v,k}(t)
    %
    % Smooth terms use spatial function handles. Localized terms store only
    % grid indices and weights, making them suitable for discrete multipoles.

    properties (Access = private)
        pressureTerms = struct('kind',{},'spatial',{},'timeFunction',{}, ...
            'grid',{},'indices',{},'weights',{})
        velocityTerms = struct('kind',{},'spatial',{},'timeFunction',{}, ...
            'grid',{},'indices',{},'weights',{})
    end

    properties
        pressureLabel (1,1) string = "Pressure source"
        pressureUnits (1,1) string = "Pa/s"
        velocityLabel (1,1) string = "Velocity source"
        velocityUnits (1,1) string = "m/s^2"
    end

    properties (Dependent)
        numberOfPressureTerms
        numberOfVelocityTerms
    end

    methods
        function obj = AcousticSource1D( ...
                pressureSpatialFunction,pressureTimeFunction, ...
                velocitySpatialFunction,velocityTimeFunction,options)
            % Construct a source, retaining the original single-term API.
            % AcousticSource1D() constructs the zero source.
            arguments
                pressureSpatialFunction (1,1) function_handle = ...
                    @(x) zeros(size(x))
                pressureTimeFunction (1,1) function_handle = ...
                    @(t) ones(size(t))
                velocitySpatialFunction (1,1) function_handle = ...
                    @(x) zeros(size(x))
                velocityTimeFunction (1,1) function_handle = ...
                    @(t) ones(size(t))
                options.PressureLabel (1,1) string = "Pressure source"
                options.PressureUnits (1,1) string = "Pa/s"
                options.VelocityLabel (1,1) string = "Velocity source"
                options.VelocityUnits (1,1) string = "m/s^2"
            end

            obj.pressureLabel = options.PressureLabel;
            obj.pressureUnits = options.PressureUnits;
            obj.velocityLabel = options.VelocityLabel;
            obj.velocityUnits = options.VelocityUnits;

            if nargin > 0
                obj.addPressureTerm(pressureSpatialFunction,pressureTimeFunction);
                obj.addVelocityTerm(velocitySpatialFunction,velocityTimeFunction);
            end
        end

        function value = get.numberOfPressureTerms(obj)
            value = numel(obj.pressureTerms);
        end

        function value = get.numberOfVelocityTerms(obj)
            value = numel(obj.velocityTerms);
        end

        function addPressureTerm(obj,spatialFunction,timeFunction)
            % Add a smooth separable pressure-source term.
            arguments
                obj
                spatialFunction (1,1) function_handle
                timeFunction (1,1) function_handle
            end
            obj.pressureTerms(end+1) = obj.smoothTerm( ...
                spatialFunction,timeFunction);
        end

        function addVelocityTerm(obj,spatialFunction,timeFunction)
            % Add a smooth separable velocity-source term.
            arguments
                obj
                spatialFunction (1,1) function_handle
                timeFunction (1,1) function_handle
            end
            obj.velocityTerms(end+1) = obj.smoothTerm( ...
                spatialFunction,timeFunction);
        end

        function addLocalizedPressureTerm(obj,grid,indices,weights,timeFunction)
            % Add a pressure term supported at selected grid indices.
            arguments
                obj
                grid (1,1) GridSpace1D
                indices double {mustBeInteger,mustBePositive}
                weights double {mustBeFinite}
                timeFunction (1,1) function_handle
            end
            obj.pressureTerms(end+1) = obj.localizedTerm( ...
                grid,indices,weights,timeFunction,"pressure");
        end

        function addLocalizedVelocityTerm(obj,grid,indices,weights,timeFunction)
            % Add a velocity term supported at selected grid indices.
            arguments
                obj
                grid (1,1) GridSpace1D
                indices double {mustBeInteger,mustBePositive}
                weights double {mustBeFinite}
                timeFunction (1,1) function_handle
            end
            obj.velocityTerms(end+1) = obj.localizedTerm( ...
                grid,indices,weights,timeFunction,"velocity");
        end

        function values = evaluatePressureSpatial(obj,x)
            % Sum smooth pressure spatial factors (localized terms omitted).
            arguments
                obj
                x double {mustBeFinite}
            end
            values = obj.evaluateSmoothSpatialSum( ...
                obj.pressureTerms,x,"pressure");
        end

        function value = evaluatePressureTime(obj,t)
            % Evaluate the time factor when pressure has exactly one term.
            arguments
                obj
                t (1,1) double {mustBeFinite}
            end
            value = obj.evaluateSingleTime(obj.pressureTerms,t,"pressure");
        end

        function values = evaluatePressure(obj,x,t)
            % Evaluate all smooth pressure terms at x and t.
            arguments
                obj
                x double {mustBeFinite}
                t (1,1) double {mustBeFinite}
            end
            values = obj.evaluateSmoothTerms(obj.pressureTerms,x,t,"pressure");
        end

        function values = evaluateVelocitySpatial(obj,x)
            % Sum smooth velocity spatial factors (localized terms omitted).
            arguments
                obj
                x double {mustBeFinite}
            end
            values = obj.evaluateSmoothSpatialSum( ...
                obj.velocityTerms,x,"velocity");
        end

        function value = evaluateVelocityTime(obj,t)
            % Evaluate the time factor when velocity has exactly one term.
            arguments
                obj
                t (1,1) double {mustBeFinite}
            end
            value = obj.evaluateSingleTime(obj.velocityTerms,t,"velocity");
        end

        function values = evaluateVelocity(obj,x,t)
            % Evaluate all smooth velocity terms at x and t.
            arguments
                obj
                x double {mustBeFinite}
                t (1,1) double {mustBeFinite}
            end
            values = obj.evaluateSmoothTerms(obj.velocityTerms,x,t,"velocity");
        end

        function G = samplePressureSpatial(obj,grid)
            % Sample summed spatial weights, treating all time factors as one.
            arguments
                obj
                grid (1,1) GridSpace1D
            end
            prepared = obj.preparePressure(grid);
            G = GridFunction(grid,Label=obj.pressureLabel,Units=obj.pressureUnits);
            G.setValues(obj.evaluatePreparedWithUnitTimes(prepared));
        end

        function G = sampleVelocitySpatial(obj,grid)
            % Sample summed spatial weights, treating all time factors as one.
            arguments
                obj
                grid (1,1) GridSpace1D
            end
            prepared = obj.prepareVelocity(grid);
            G = GridFunction(grid,Label=obj.velocityLabel,Units=obj.velocityUnits);
            G.setValues(obj.evaluatePreparedWithUnitTimes(prepared));
        end

        function prepared = preparePressure(obj,grid)
            % Cache all pressure spatial terms on grid before time stepping.
            arguments
                obj
                grid (1,1) GridSpace1D
            end
            prepared = obj.prepareTerms(obj.pressureTerms,grid,"pressure");
        end

        function prepared = prepareVelocity(obj,grid)
            % Cache all velocity spatial terms on grid before time stepping.
            arguments
                obj
                grid (1,1) GridSpace1D
            end
            prepared = obj.prepareTerms(obj.velocityTerms,grid,"velocity");
        end

        function values = evaluatePrepared(obj,prepared,t)
            % Evaluate a source returned by preparePressure/prepareVelocity.
            arguments
                obj
                prepared (1,1) struct
                t (1,1) double {mustBeFinite}
            end
            obj.validatePreparedSource(prepared);
            nTerms = numel(prepared.timeFunctions);
            if nTerms == 0
                values = zeros(prepared.nGridPoints,1);
                return
            end
            amplitudes = zeros(nTerms,1);
            for k = 1:nTerms
                amplitudes(k) = obj.evaluateTimeFactor( ...
                    prepared.timeFunctions{k},t,prepared.component);
            end
            % Solver state arrays and GridFunction values are full vectors.
            % Keep the cached spatial operator sparse, but return a full
            % forcing vector so callers see a consistent storage type.
            values = full(prepared.spatialWeights*amplitudes);
        end
    end

    methods (Static)
        function obj = zero()
            obj = AcousticSource1D();
        end

        function obj = combine(sources)
            % Combine an AcousticSource1D array or cell array.
            if isa(sources,'AcousticSource1D')
                sources = num2cell(sources);
            end
            if ~iscell(sources) || ...
                    any(~cellfun(@(s) isa(s,'AcousticSource1D'),sources))
                error('AcousticSource1D:InvalidSourceCollection', ...
                    'Input must contain only AcousticSource1D objects.');
            end
            obj = AcousticSource1D.zero();
            for k = 1:numel(sources)
                obj.pressureTerms = [obj.pressureTerms, ...
                    sources{k}.pressureTerms]; %#ok<AGROW>
                obj.velocityTerms = [obj.velocityTerms, ...
                    sources{k}.velocityTerms]; %#ok<AGROW>
            end
        end
    end

    methods (Access = private)
        function term = smoothTerm(~,spatialFunction,timeFunction)
            term = struct('kind',"smooth",'spatial',spatialFunction, ...
                'timeFunction',timeFunction,'grid',[], ...
                'indices',[],'weights',[]);
        end

        function term = localizedTerm(~,grid,indices,weights, ...
                timeFunction,componentName)
            indices = indices(:);
            weights = weights(:);
            if numel(indices) ~= numel(weights)
                error('AcousticSource1D:LocalizedSizeMismatch', ...
                    ['The %s localized indices and weights must have the ' ...
                     'same number of entries.'],componentName);
            end
            if any(indices > grid.N)
                error('AcousticSource1D:LocalizedIndexOutsideGrid', ...
                    'A %s source index lies outside its grid.',componentName);
            end
            if numel(unique(indices)) ~= numel(indices)
                error('AcousticSource1D:DuplicateLocalizedIndex', ...
                    'Localized source indices must be unique within a term.');
            end
            term = struct('kind',"localized",'spatial',[], ...
                'timeFunction',timeFunction,'grid',grid, ...
                'indices',indices,'weights',weights);
        end

        function prepared = prepareTerms(obj,terms,grid,componentName)
            nTerms = numel(terms);
            allLocalized = nTerms > 0 && ...
                all(arrayfun(@(term) term.kind == "localized",terms));
            if allLocalized
                spatialWeights = sparse(grid.N,nTerms);
            else
                spatialWeights = zeros(grid.N,nTerms);
            end
            timeFunctions = cell(1,nTerms);
            for k = 1:nTerms
                term = terms(k);
                if term.kind == "smooth"
                    spatialWeights(:,k) = obj.evaluateSpatialFactor( ...
                        term.spatial,grid.x.pts,componentName);
                else
                    if ~isequal(term.grid,grid)
                        error('AcousticSource1D:LocalizedGridMismatch', ...
                            ['A localized %s term must be used on the same ' ...
                             'grid on which it was constructed.'],componentName);
                    end
                    spatialWeights(term.indices,k) = term.weights;
                end
                timeFunctions{k} = term.timeFunction;
            end
            prepared = struct('spatialWeights',spatialWeights, ...
                'timeFunctions',{timeFunctions},'nGridPoints',grid.N, ...
                'component',componentName);
        end

        function values = evaluatePreparedWithUnitTimes(~,prepared)
            if isempty(prepared.timeFunctions)
                values = zeros(prepared.nGridPoints,1);
            else
                values = full(prepared.spatialWeights* ...
                    ones(numel(prepared.timeFunctions),1));
            end
        end

        function validatePreparedSource(~,prepared)
            required = {'spatialWeights','timeFunctions', ...
                'nGridPoints','component'};
            if ~all(isfield(prepared,required))
                error('AcousticSource1D:InvalidPreparedSource', ...
                    ['Input must be produced by preparePressure or ' ...
                     'prepareVelocity.']);
            end
        end

        function values = evaluateSmoothSpatialSum(obj,terms,x,componentName)
            values = zeros(size(x));
            for k = 1:numel(terms)
                if terms(k).kind == "smooth"
                    termValues = obj.evaluateSpatialFactor( ...
                        terms(k).spatial,x,componentName);
                    values = values+reshape(termValues,size(x));
                end
            end
        end

        function values = evaluateSmoothTerms(obj,terms,x,t,componentName)
            values = zeros(size(x));
            for k = 1:numel(terms)
                if terms(k).kind == "smooth"
                    spatial = obj.evaluateSpatialFactor( ...
                        terms(k).spatial,x,componentName);
                    temporal = obj.evaluateTimeFactor( ...
                        terms(k).timeFunction,t,componentName);
                    values = values+reshape(spatial,size(x)).*temporal;
                end
            end
        end

        function value = evaluateSingleTime(obj,terms,t,componentName)
            if isempty(terms)
                value = 1;
                return
            end
            if numel(terms) ~= 1
                error('AcousticSource1D:NonseparableComponent', ...
                    ['A component time factor is available only for one ' ...
                     'separable term. Use prepare and evaluatePrepared ' ...
                     'for a sum of terms.']);
            end
            value = obj.evaluateTimeFactor( ...
                terms(1).timeFunction,t,componentName);
        end

        function values = evaluateSpatialFactor(~,f,x,componentName)
            values = f(x);
            if ~isnumeric(values)
                error('AcousticSource1D:InvalidSpatialOutput', ...
                    'The %s spatial function must return numeric data.', ...
                    componentName);
            end
            if isscalar(values)
                values = repmat(values,size(x));
            elseif ~isequal(size(values),size(x))
                error('AcousticSource1D:InvalidSpatialSize', ...
                    ['The %s spatial function must return a scalar or an ' ...
                     'array having the same size as its coordinate input.'], ...
                    componentName);
            end
            if any(~isfinite(values),'all')
                error('AcousticSource1D:NonfiniteSpatialOutput', ...
                    'The %s spatial function returned a nonfinite value.', ...
                    componentName);
            end
            values = values(:);
        end

        function value = evaluateTimeFactor(~,f,t,componentName)
            value = f(t);
            if ~isnumeric(value) || ~isscalar(value) || ~isfinite(value)
                error('AcousticSource1D:InvalidTimeOutput', ...
                    ['The %s time function must return one finite numeric ' ...
                     'scalar for a scalar time input.'],componentName);
            end
        end
    end
end

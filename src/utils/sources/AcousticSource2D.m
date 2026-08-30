classdef AcousticSource2D < AcousticSource
    % AcousticSource2D represents sums of separable 2D acoustic sources.
    %
    %   f_p(x,y,t)  = sum_k g_{p,k}(x,y)  w_{p,k}(t)
    %   f_vx(x,y,t) = sum_k g_{vx,k}(x,y) w_{vx,k}(t)
    %   f_vy(x,y,t) = sum_k g_{vy,k}(x,y) w_{vy,k}(t)
    %
    % Smooth terms use spatial function handles of the form g(x,y).
    % Localized terms use MATLAB linear indices into a GridSpace2D array.

    properties (Access = private)
        pressureTerms = struct('kind',{},'spatial',{},'timeFunction',{}, ...
            'grid',{},'indices',{},'weights',{})
        velocityXTerms = struct('kind',{},'spatial',{},'timeFunction',{}, ...
            'grid',{},'indices',{},'weights',{})
        velocityYTerms = struct('kind',{},'spatial',{},'timeFunction',{}, ...
            'grid',{},'indices',{},'weights',{})
    end

    properties
        pressureLabel (1,1) string = "Pressure source"
        pressureUnits (1,1) string = "Pa/s"
        velocityXLabel (1,1) string = "x-velocity source"
        velocityXUnits (1,1) string = "m/s^2"
        velocityYLabel (1,1) string = "y-velocity source"
        velocityYUnits (1,1) string = "m/s^2"
    end

    properties (Dependent)
        numberOfPressureTerms
        numberOfVelocityXTerms
        numberOfVelocityYTerms
    end

    methods
        function obj = AcousticSource2D( ...
                pressureSpatialFunction,pressureTimeFunction, ...
                velocityXSpatialFunction,velocityXTimeFunction, ...
                velocityYSpatialFunction,velocityYTimeFunction,options)
            % Construct a source using an optional single term per component.
            % AcousticSource2D() constructs the zero source.
            arguments
                pressureSpatialFunction (1,1) function_handle = ...
                    @(x,y) zeros(size(x))
                pressureTimeFunction (1,1) function_handle = ...
                    @(t) ones(size(t))
                velocityXSpatialFunction (1,1) function_handle = ...
                    @(x,y) zeros(size(x))
                velocityXTimeFunction (1,1) function_handle = ...
                    @(t) ones(size(t))
                velocityYSpatialFunction (1,1) function_handle = ...
                    @(x,y) zeros(size(x))
                velocityYTimeFunction (1,1) function_handle = ...
                    @(t) ones(size(t))
                options.PressureLabel (1,1) string = "Pressure source"
                options.PressureUnits (1,1) string = "Pa/s"
                options.VelocityXLabel (1,1) string = "x-velocity source"
                options.VelocityXUnits (1,1) string = "m/s^2"
                options.VelocityYLabel (1,1) string = "y-velocity source"
                options.VelocityYUnits (1,1) string = "m/s^2"
            end

            obj.pressureLabel = options.PressureLabel;
            obj.pressureUnits = options.PressureUnits;
            obj.velocityXLabel = options.VelocityXLabel;
            obj.velocityXUnits = options.VelocityXUnits;
            obj.velocityYLabel = options.VelocityYLabel;
            obj.velocityYUnits = options.VelocityYUnits;

            if nargin > 0
                obj.addPressureTerm(pressureSpatialFunction,pressureTimeFunction);
                obj.addVelocityXTerm(velocityXSpatialFunction,velocityXTimeFunction);
                obj.addVelocityYTerm(velocityYSpatialFunction,velocityYTimeFunction);
            end
        end

        function value = get.numberOfPressureTerms(obj)
            value = numel(obj.pressureTerms);
        end

        function value = get.numberOfVelocityXTerms(obj)
            value = numel(obj.velocityXTerms);
        end

        function value = get.numberOfVelocityYTerms(obj)
            value = numel(obj.velocityYTerms);
        end

        function append(obj,other)
            % Append all terms from another 2D acoustic source.
            arguments
                obj
                other (1,1) AcousticSource2D
            end

            obj.pressureTerms = [obj.pressureTerms,other.pressureTerms];
            obj.velocityXTerms = [obj.velocityXTerms,other.velocityXTerms];
            obj.velocityYTerms = [obj.velocityYTerms,other.velocityYTerms];
        end

        function addPressureTerm(obj,spatialFunction,timeFunction)
            arguments
                obj
                spatialFunction (1,1) function_handle
                timeFunction (1,1) function_handle
            end
            obj.pressureTerms(end+1) = obj.smoothTerm(spatialFunction,timeFunction);
        end

        function addVelocityXTerm(obj,spatialFunction,timeFunction)
            arguments
                obj
                spatialFunction (1,1) function_handle
                timeFunction (1,1) function_handle
            end
            obj.velocityXTerms(end+1) = obj.smoothTerm(spatialFunction,timeFunction);
        end

        function addVelocityYTerm(obj,spatialFunction,timeFunction)
            arguments
                obj
                spatialFunction (1,1) function_handle
                timeFunction (1,1) function_handle
            end
            obj.velocityYTerms(end+1) = obj.smoothTerm(spatialFunction,timeFunction);
        end

        function addLocalizedPressureTerm(obj,grid,indices,weights,timeFunction)
            arguments
                obj
                grid (1,1) GridSpace2D
                indices double {mustBeInteger,mustBePositive}
                weights double {mustBeFinite}
                timeFunction (1,1) function_handle
            end
            obj.pressureTerms(end+1) = obj.localizedTerm( ...
                grid,indices,weights,timeFunction,"pressure");
        end

        function addLocalizedVelocityXTerm(obj,grid,indices,weights,timeFunction)
            arguments
                obj
                grid (1,1) GridSpace2D
                indices double {mustBeInteger,mustBePositive}
                weights double {mustBeFinite}
                timeFunction (1,1) function_handle
            end
            obj.velocityXTerms(end+1) = obj.localizedTerm( ...
                grid,indices,weights,timeFunction,"x-velocity");
        end

        function addLocalizedVelocityYTerm(obj,grid,indices,weights,timeFunction)
            arguments
                obj
                grid (1,1) GridSpace2D
                indices double {mustBeInteger,mustBePositive}
                weights double {mustBeFinite}
                timeFunction (1,1) function_handle
            end
            obj.velocityYTerms(end+1) = obj.localizedTerm( ...
                grid,indices,weights,timeFunction,"y-velocity");
        end

        function values = evaluatePressureSpatial(obj,x,y)
            obj.validateCoordinateArrays(x,y);
            values = obj.evaluateSmoothSpatialSum(obj.pressureTerms,x,y,"pressure");
        end

        function value = evaluatePressureTime(obj,t)
            arguments
                obj
                t (1,1) double {mustBeFinite}
            end
            value = obj.evaluateSingleTime(obj.pressureTerms,t,"pressure");
        end

        function values = evaluatePressure(obj,x,y,t)
            obj.validateCoordinateArrays(x,y);
            obj.validateScalarTime(t);
            values = obj.evaluateSmoothTerms(obj.pressureTerms,x,y,t,"pressure");
        end

        function values = evaluateVelocityXSpatial(obj,x,y)
            obj.validateCoordinateArrays(x,y);
            values = obj.evaluateSmoothSpatialSum( ...
                obj.velocityXTerms,x,y,"x-velocity");
        end

        function value = evaluateVelocityXTime(obj,t)
            arguments
                obj
                t (1,1) double {mustBeFinite}
            end
            value = obj.evaluateSingleTime(obj.velocityXTerms,t,"x-velocity");
        end

        function values = evaluateVelocityX(obj,x,y,t)
            obj.validateCoordinateArrays(x,y);
            obj.validateScalarTime(t);
            values = obj.evaluateSmoothTerms( ...
                obj.velocityXTerms,x,y,t,"x-velocity");
        end

        function values = evaluateVelocityYSpatial(obj,x,y)
            obj.validateCoordinateArrays(x,y);
            values = obj.evaluateSmoothSpatialSum( ...
                obj.velocityYTerms,x,y,"y-velocity");
        end

        function value = evaluateVelocityYTime(obj,t)
            arguments
                obj
                t (1,1) double {mustBeFinite}
            end
            value = obj.evaluateSingleTime(obj.velocityYTerms,t,"y-velocity");
        end

        function values = evaluateVelocityY(obj,x,y,t)
            obj.validateCoordinateArrays(x,y);
            obj.validateScalarTime(t);
            values = obj.evaluateSmoothTerms( ...
                obj.velocityYTerms,x,y,t,"y-velocity");
        end

        function G = samplePressureSpatial(obj,grid)
            arguments
                obj
                grid (1,1) GridSpace2D
            end
            prepared = obj.preparePressure(grid);
            G = GridFunction(grid,Label=obj.pressureLabel,Units=obj.pressureUnits);
            G.setValues(obj.evaluatePreparedWithUnitTimes(prepared));
        end

        function G = sampleVelocityXSpatial(obj,grid)
            arguments
                obj
                grid (1,1) GridSpace2D
            end
            prepared = obj.prepareVelocityX(grid);
            G = GridFunction(grid,Label=obj.velocityXLabel,Units=obj.velocityXUnits);
            G.setValues(obj.evaluatePreparedWithUnitTimes(prepared));
        end

        function G = sampleVelocityYSpatial(obj,grid)
            arguments
                obj
                grid (1,1) GridSpace2D
            end
            prepared = obj.prepareVelocityY(grid);
            G = GridFunction(grid,Label=obj.velocityYLabel,Units=obj.velocityYUnits);
            G.setValues(obj.evaluatePreparedWithUnitTimes(prepared));
        end

        function prepared = preparePressure(obj,grid)
            arguments
                obj
                grid (1,1) GridSpace2D
            end
            prepared = obj.prepareTerms(obj.pressureTerms,grid,"pressure");
        end

        function prepared = prepareVelocityX(obj,grid)
            arguments
                obj
                grid (1,1) GridSpace2D
            end
            prepared = obj.prepareTerms(obj.velocityXTerms,grid,"x-velocity");
        end

        function prepared = prepareVelocityY(obj,grid)
            arguments
                obj
                grid (1,1) GridSpace2D
            end
            prepared = obj.prepareTerms(obj.velocityYTerms,grid,"y-velocity");
        end

        function values = evaluatePrepared(obj,prepared,t)
            arguments
                obj
                prepared (1,1) struct
                t (1,1) double {mustBeFinite}
            end
            obj.validatePreparedSource(prepared);
            nTerms = numel(prepared.timeFunctions);
            if nTerms == 0
                values = zeros(prepared.gridSize);
                return
            end
            amplitudes = zeros(nTerms,1);
            for k = 1:nTerms
                amplitudes(k) = obj.evaluateTimeFactor( ...
                    prepared.timeFunctions{k},t,prepared.component);
            end
            values = reshape(full(prepared.spatialWeights*amplitudes), ...
                prepared.gridSize);
        end
    end

    methods (Static)
        function obj = zero()
            obj = AcousticSource2D();
        end

        function obj = combine(sources)
            if isa(sources,'AcousticSource2D')
                sources = num2cell(sources);
            end
            if ~iscell(sources) || ...
                    any(~cellfun(@(s) isa(s,'AcousticSource2D'),sources))
                error('AcousticSource2D:InvalidSourceCollection', ...
                    'Input must contain only AcousticSource2D objects.');
            end
            obj = AcousticSource2D.zero();
            for k = 1:numel(sources)
                obj.append(sources{k});
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
                error('AcousticSource2D:LocalizedSizeMismatch', ...
                    ['The %s localized indices and weights must have the ' ...
                     'same number of entries.'],componentName);
            end
            if any(indices > prod(grid.N))
                error('AcousticSource2D:LocalizedIndexOutsideGrid', ...
                    'A %s source index lies outside its grid.',componentName);
            end
            if numel(unique(indices)) ~= numel(indices)
                error('AcousticSource2D:DuplicateLocalizedIndex', ...
                    'Localized source indices must be unique within a term.');
            end
            term = struct('kind',"localized",'spatial',[], ...
                'timeFunction',timeFunction,'grid',grid, ...
                'indices',indices,'weights',weights);
        end

        function prepared = prepareTerms(obj,terms,grid,componentName)
            nTerms = numel(terms);
            nGridPoints = prod(grid.N);
            allLocalized = nTerms > 0 && ...
                all(arrayfun(@(term) term.kind == "localized",terms));
            if allLocalized
                spatialWeights = sparse(nGridPoints,nTerms);
            else
                spatialWeights = zeros(nGridPoints,nTerms);
            end
            timeFunctions = cell(1,nTerms);
            [X,Y] = grid.mesh();
            for k = 1:nTerms
                term = terms(k);
                if term.kind == "smooth"
                    spatialWeights(:,k) = obj.evaluateSpatialFactor( ...
                        term.spatial,X,Y,componentName);
                else
                    if ~isequal(term.grid,grid)
                        error('AcousticSource2D:LocalizedGridMismatch', ...
                            ['A localized %s term must be used on the same ' ...
                             'grid on which it was constructed.'],componentName);
                    end
                    spatialWeights(term.indices,k) = term.weights;
                end
                timeFunctions{k} = term.timeFunction;
            end
            prepared = struct('spatialWeights',spatialWeights, ...
                'timeFunctions',{timeFunctions},'gridSize',grid.N, ...
                'nGridPoints',nGridPoints,'component',componentName);
        end

        function values = evaluatePreparedWithUnitTimes(~,prepared)
            if isempty(prepared.timeFunctions)
                values = zeros(prepared.gridSize);
            else
                values = reshape(full(prepared.spatialWeights* ...
                    ones(numel(prepared.timeFunctions),1)),prepared.gridSize);
            end
        end

        function validatePreparedSource(~,prepared)
            required = {'spatialWeights','timeFunctions','gridSize', ...
                'nGridPoints','component'};
            if ~all(isfield(prepared,required))
                error('AcousticSource2D:InvalidPreparedSource', ...
                    ['Input must be produced by preparePressure, ' ...
                     'prepareVelocityX, or prepareVelocityY.']);
            end
        end

        function values = evaluateSmoothSpatialSum(obj,terms,x,y,componentName)
            values = zeros(size(x));
            for k = 1:numel(terms)
                if terms(k).kind == "smooth"
                    spatial = obj.evaluateSpatialFactor( ...
                        terms(k).spatial,x,y,componentName);
                    values = values+reshape(spatial,size(x));
                end
            end
        end

        function values = evaluateSmoothTerms(obj,terms,x,y,t,componentName)
            values = zeros(size(x));
            for k = 1:numel(terms)
                if terms(k).kind == "smooth"
                    spatial = obj.evaluateSpatialFactor( ...
                        terms(k).spatial,x,y,componentName);
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
                error('AcousticSource2D:NonseparableComponent', ...
                    ['A component time factor is available only for one ' ...
                     'separable term. Use prepare and evaluatePrepared ' ...
                     'for a sum of terms.']);
            end
            value = obj.evaluateTimeFactor( ...
                terms(1).timeFunction,t,componentName);
        end

        function values = evaluateSpatialFactor(~,f,x,y,componentName)
            values = f(x,y);
            if ~isnumeric(values)
                error('AcousticSource2D:InvalidSpatialOutput', ...
                    'The %s spatial function must return numeric data.', ...
                    componentName);
            end
            if isscalar(values)
                values = repmat(values,size(x));
            elseif ~isequal(size(values),size(x))
                error('AcousticSource2D:InvalidSpatialSize', ...
                    ['The %s spatial function must return a scalar or an ' ...
                     'array having the same size as its coordinate inputs.'], ...
                    componentName);
            end
            if any(~isfinite(values),'all')
                error('AcousticSource2D:NonfiniteSpatialOutput', ...
                    'The %s spatial function returned a nonfinite value.', ...
                    componentName);
            end
            values = values(:);
        end

        function value = evaluateTimeFactor(~,f,t,componentName)
            value = f(t);
            if ~isnumeric(value) || ~isscalar(value) || ~isfinite(value)
                error('AcousticSource2D:InvalidTimeOutput', ...
                    ['The %s time function must return one finite numeric ' ...
                     'scalar for a scalar time input.'],componentName);
            end
        end

        function validateCoordinateArrays(~,x,y)
            if ~isnumeric(x) || ~isnumeric(y) || ~isequal(size(x),size(y))
                error('AcousticSource2D:InvalidCoordinates', ...
                    'x and y must be numeric arrays having the same size.');
            end
            if any(~isfinite(x),'all') || any(~isfinite(y),'all')
                error('AcousticSource2D:InvalidCoordinates', ...
                    'x and y must contain only finite values.');
            end
        end

        function validateScalarTime(~,t)
            if ~isnumeric(t) || ~isscalar(t) || ~isfinite(t)
                error('AcousticSource2D:InvalidTimeInput', ...
                    't must be one finite numeric scalar.');
            end
        end
    end
end

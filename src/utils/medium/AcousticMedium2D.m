classdef AcousticMedium2D
    % AcousticMedium2D represents a continuous 2D acoustic medium.
    % Stores functions defining kappa(x,y) and beta(x,y), with factory
    % methods for common analytic, layered, inclusion, and tabulated models.

    properties (SetAccess = private)
        kappaFunction
        betaFunction
    end

    properties
        kappaLabel (1,1) string = "Bulk modulus"
        kappaUnits (1,1) string = "Pa"
        betaLabel (1,1) string = "Buoyancy"
        betaUnits (1,1) string = "m^3/kg"
        waveSpeedLabel (1,1) string = "Wave speed"
        waveSpeedUnits (1,1) string = "m/s"
    end

    %% Public methods
    methods
        %%%%%%%%
        function obj = AcousticMedium2D(kappaFunction,betaFunction,options)
            arguments
                kappaFunction (1,1) function_handle
                betaFunction  (1,1) function_handle
                options.KappaLabel (1,1) string = "Bulk modulus"
                options.KappaUnits (1,1) string = "Pa"
                options.BetaLabel (1,1) string = "Buoyancy"
                options.BetaUnits (1,1) string = "m^3/kg"
                options.WaveSpeedLabel (1,1) string = "Wave speed"
                options.WaveSpeedUnits (1,1) string = "m/s"
            end

            obj.kappaFunction = kappaFunction;
            obj.betaFunction = betaFunction;
            obj.kappaLabel = options.KappaLabel;
            obj.kappaUnits = options.KappaUnits;
            obj.betaLabel = options.BetaLabel;
            obj.betaUnits = options.BetaUnits;
            obj.waveSpeedLabel = options.WaveSpeedLabel;
            obj.waveSpeedUnits = options.WaveSpeedUnits;
        end

        %%%%%%%%
        function values = evaluateKappa(obj,x,y)
            % Evaluate kappa at numeric coordinate arrays x and y.
            values = obj.kappaFunction(x,y);
        end

        %%%%%%%%
        function values = evaluateBeta(obj,x,y)
            % Evaluate beta at numeric coordinate arrays x and y.
            values = obj.betaFunction(x,y);
        end

        %%%%%%%%
        function values = evaluateWaveSpeed(obj,x,y)
            % Evaluate c(x,y) = sqrt(kappa(x,y).*beta(x,y)).
            values = sqrt( ...
                obj.evaluateKappa(x,y).*obj.evaluateBeta(x,y));
        end
    end

    %% Static methods
    methods (Static)
        %%%%%%%%
        function obj = constant(options)
            arguments
                options.KappaValue (1,1) double ...
                    {mustBePositive,mustBeFinite} = 1
                options.BetaValue (1,1) double ...
                    {mustBePositive,mustBeFinite} = 1
            end

            kappa = AcousticMedium2D.makeConstantFunction( ...
                options.KappaValue);
            beta = AcousticMedium2D.makeConstantFunction( ...
                options.BetaValue);
            obj = AcousticMedium2D(kappa,beta);
        end

        %%%%%%%%
        function obj = linearX(options)
            arguments
                options.KappaRange (1,2) double ...
                    {mustBePositive,mustBeFinite} = [1,2]
                options.BetaRange (1,2) double ...
                    {mustBePositive,mustBeFinite} = [1,2]
                options.Domain (2,2) double ...
                    {mustBeFinite,mustBeValidDomain2D} = [0,1;0,1]
            end

            kappa = AcousticMedium2D.makeLinearXFunction( ...
                options.KappaRange,options.Domain(1,:));
            beta = AcousticMedium2D.makeLinearXFunction( ...
                options.BetaRange,options.Domain(1,:));
            obj = AcousticMedium2D(kappa,beta);
        end

        %%%%%%%%
        function obj = linearY(options)
            arguments
                options.KappaRange (1,2) double ...
                    {mustBePositive,mustBeFinite} = [1,2]
                options.BetaRange (1,2) double ...
                    {mustBePositive,mustBeFinite} = [1,2]
                options.Domain (2,2) double ...
                    {mustBeFinite,mustBeValidDomain2D} = [0,1;0,1]
            end

            kappa = AcousticMedium2D.makeLinearYFunction( ...
                options.KappaRange,options.Domain(2,:));
            beta = AcousticMedium2D.makeLinearYFunction( ...
                options.BetaRange,options.Domain(2,:));
            obj = AcousticMedium2D(kappa,beta);
        end

        %%%%%%%%
        function obj = smooth(options)
            arguments
                options.KappaRange (1,2) double ...
                    {mustBePositive,mustBeFinite,mustBeIncreasingRange} = [1,2]
                options.BetaRange (1,2) double ...
                    {mustBePositive,mustBeFinite,mustBeIncreasingRange} = [1,2]
                options.Domain (2,2) double ...
                    {mustBeFinite,mustBeValidDomain2D} = [0,1;0,1]
            end

            kappa = AcousticMedium2D.makeSmoothFunction( ...
                options.KappaRange,options.Domain);
            beta = AcousticMedium2D.makeSmoothFunction( ...
                options.BetaRange,options.Domain);
            obj = AcousticMedium2D(kappa,beta);
        end

        %%%%%%%%
        function obj = gaussianPeak(options)
            arguments
                options.KappaRange (1,2) double ...
                    {mustBePositive,mustBeFinite,mustBeIncreasingRange} = [1,2]
                options.BetaRange (1,2) double ...
                    {mustBePositive,mustBeFinite,mustBeIncreasingRange} = [1,2]
                options.Center (1,2) double {mustBeFinite} = [0.5,0.5]
                options.Width (1,2) double ...
                    {mustBePositive,mustBeFinite} = [1,1]
            end

            kappa = AcousticMedium2D.makeGaussianFunction( ...
                options.KappaRange,options.Center,options.Width,true);
            beta = AcousticMedium2D.makeGaussianFunction( ...
                options.BetaRange,options.Center,options.Width,true);
            obj = AcousticMedium2D(kappa,beta);
        end

        %%%%%%%%
        function obj = gaussianDip(options)
            arguments
                options.KappaRange (1,2) double ...
                    {mustBePositive,mustBeFinite,mustBeIncreasingRange} = [1,2]
                options.BetaRange (1,2) double ...
                    {mustBePositive,mustBeFinite,mustBeIncreasingRange} = [1,2]
                options.Center (1,2) double {mustBeFinite} = [0.5,0.5]
                options.Width (1,2) double ...
                    {mustBePositive,mustBeFinite} = [1,1]
            end

            kappa = AcousticMedium2D.makeGaussianFunction( ...
                options.KappaRange,options.Center,options.Width,false);
            beta = AcousticMedium2D.makeGaussianFunction( ...
                options.BetaRange,options.Center,options.Width,false);
            obj = AcousticMedium2D(kappa,beta);
        end

        %%%%%%%%
        function obj = layeredX(options)
            arguments
                options.KappaValues (1,:) double ...
                    {mustBePositive,mustBeFinite}
                options.BetaValues (1,:) double ...
                    {mustBePositive,mustBeFinite}
                options.RelativeWidths (1,:) double ...
                    {mustBePositive,mustBeFinite}
                options.Domain (2,2) double ...
                    {mustBeFinite,mustBeValidDomain2D} = [0,1;0,1]
            end

            AcousticMedium2D.validateLayerData( ...
                options.KappaValues,options.BetaValues, ...
                options.RelativeWidths);
            edges = AcousticMedium2D.makeLayerEdges( ...
                options.RelativeWidths,options.Domain(1,:));
            kappa = AcousticMedium2D.makeLayerXFunction( ...
                options.KappaValues,edges);
            beta = AcousticMedium2D.makeLayerXFunction( ...
                options.BetaValues,edges);
            obj = AcousticMedium2D(kappa,beta);
        end

        %%%%%%%%
        function obj = layeredY(options)
            arguments
                options.KappaValues (1,:) double ...
                    {mustBePositive,mustBeFinite}
                options.BetaValues (1,:) double ...
                    {mustBePositive,mustBeFinite}
                options.RelativeWidths (1,:) double ...
                    {mustBePositive,mustBeFinite}
                options.Domain (2,2) double ...
                    {mustBeFinite,mustBeValidDomain2D} = [0,1;0,1]
            end

            AcousticMedium2D.validateLayerData( ...
                options.KappaValues,options.BetaValues, ...
                options.RelativeWidths);
            edges = AcousticMedium2D.makeLayerEdges( ...
                options.RelativeWidths,options.Domain(2,:));
            kappa = AcousticMedium2D.makeLayerYFunction( ...
                options.KappaValues,edges);
            beta = AcousticMedium2D.makeLayerYFunction( ...
                options.BetaValues,edges);
            obj = AcousticMedium2D(kappa,beta);
        end

        %%%%%%%%
        function obj = circular(options)
            arguments
                options.KappaInside (1,1) double ...
                    {mustBePositive,mustBeFinite}
                options.KappaOutside (1,1) double ...
                    {mustBePositive,mustBeFinite}
                options.BetaInside (1,1) double ...
                    {mustBePositive,mustBeFinite}
                options.BetaOutside (1,1) double ...
                    {mustBePositive,mustBeFinite}
                options.Center (1,2) double {mustBeFinite} = [0.5,0.5]
                options.Radius (1,1) double ...
                    {mustBePositive,mustBeFinite} = 0.25
            end

            kappa = AcousticMedium2D.makeCircularFunction( ...
                options.KappaInside,options.KappaOutside, ...
                options.Center,options.Radius);
            beta = AcousticMedium2D.makeCircularFunction( ...
                options.BetaInside,options.BetaOutside, ...
                options.Center,options.Radius);
            obj = AcousticMedium2D(kappa,beta);
        end

        %%%%%%%%
        function obj = tabulated(options)
            arguments
                options.XCoordinates double ...
                    {mustBeVector,mustBeFinite}
                options.YCoordinates double ...
                    {mustBeVector,mustBeFinite}
                options.KappaValues double ...
                    {mustBePositive,mustBeFinite}
                options.BetaValues double ...
                    {mustBePositive,mustBeFinite}
                options.Method (1,1) string = "linear"
                options.ExtrapolationMethod (1,1) string = "none"
            end

            x = options.XCoordinates(:);
            y = options.YCoordinates(:);
            expectedSize = [numel(x),numel(y)];
            if ~isequal(size(options.KappaValues),expectedSize) || ...
                    ~isequal(size(options.BetaValues),expectedSize)
                error('AcousticMedium2D:InvalidTabulatedSize', ...
                    ['KappaValues and BetaValues must have size ' ...
                     '[numel(XCoordinates),numel(YCoordinates)].']);
            end
            if any(diff(x) <= 0) || any(diff(y) <= 0)
                error('AcousticMedium2D:InvalidCoordinates', ...
                    'Tabulated coordinates must be strictly increasing.');
            end

            kappaInterpolant = griddedInterpolant( ...
                {x,y},options.KappaValues, ...
                options.Method,options.ExtrapolationMethod);
            betaInterpolant = griddedInterpolant( ...
                {x,y},options.BetaValues, ...
                options.Method,options.ExtrapolationMethod);

            kappa = @(xq,yq) kappaInterpolant(xq,yq);
            beta = @(xq,yq) betaInterpolant(xq,yq);
            obj = AcousticMedium2D(kappa,beta);
        end
    end

    %% Private methods
    methods (Static, Access = private)
        %%%%%%%%
        function f = makeConstantFunction(value)
            f = @(x,y) value*ones(size(x));
        end
        
        %%%%%%%%
        function f = makeLinearXFunction(valueRange,xDomain)
            slope = diff(valueRange)/diff(xDomain);
            f = @(x,y) valueRange(1)+slope.*(x-xDomain(1));
        end

        %%%%%%%%
        function f = makeLinearYFunction(valueRange,yDomain)
            slope = diff(valueRange)/diff(yDomain);
            f = @(x,y) valueRange(1)+slope.*(y-yDomain(1));
        end

        %%%%%%%%
        function f = makeSmoothFunction(valueRange,domain)
            average = mean(valueRange);
            amplitude = diff(valueRange)/2;
            xLeft = domain(1,1);
            yLeft = domain(2,1);
            Lx = diff(domain(1,:));
            Ly = diff(domain(2,:));
            f = @(x,y) average+amplitude.* ...
                sin(2*pi*(x-xLeft)/Lx).*sin(2*pi*(y-yLeft)/Ly);
        end

        %%%%%%%%
        function f = makeGaussianFunction(valueRange,center,width,isPeak)
            gaussian = @(x,y) exp( ...
                -((x-center(1))./width(1)).^2 ...
                -((y-center(2))./width(2)).^2);
            if isPeak
                f = @(x,y) valueRange(1)+diff(valueRange).*gaussian(x,y);
            else
                f = @(x,y) valueRange(2)-diff(valueRange).*gaussian(x,y);
            end
        end

        %%%%%%%%
        function edges = makeLayerEdges(relativeWidths,domain)
            widths = relativeWidths/sum(relativeWidths)*diff(domain);
            edges = domain(1)+[0,cumsum(widths)];
            edges(end) = domain(2);
        end

        %%%%%%%%
        function f = makeLayerXFunction(values,edges)
            f = @(x,y) AcousticMedium2D.evaluateLayers(x,values,edges);
        end

        %%%%%%%%
        function f = makeLayerYFunction(values,edges)
            f = @(x,y) AcousticMedium2D.evaluateLayers(y,values,edges);
        end

        %%%%%%%%
        function result = evaluateLayers(coordinates,values,edges)
            layer = discretize(coordinates,edges);
            result = reshape(values(layer(:)),size(coordinates));
        end

        %%%%%%%%
        function f = makeCircularFunction(inside,outside,center,radius)
            f = @(x,y) outside+(inside-outside).* ...
                (((x-center(1)).^2+(y-center(2)).^2) <= radius^2);
        end

        %%%%%%%%
        function validateLayerData(kappaValues,betaValues,widths)
            n = numel(widths);
            if numel(kappaValues) ~= n || numel(betaValues) ~= n
                error('AcousticMedium2D:InvalidLayerData', ...
                    ['KappaValues, BetaValues, and RelativeWidths must ' ...
                     'contain the same number of entries.']);
            end
        end
    end
end


function mustBeIncreasingRange(value)
if value(1) > value(2)
    error('AcousticMedium2D:InvalidRange', ...
        'Range must have the form [minimum, maximum].');
end
end

function mustBeValidDomain2D(domain)
if any(domain(:,1) >= domain(:,2))
    error('AcousticMedium2D:InvalidDomain', ...
        'Each row of Domain must have the form [minimum, maximum].');
end
end

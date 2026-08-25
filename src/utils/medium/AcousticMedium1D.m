classdef AcousticMedium1D
    % AcousticMedium1D represents a continuous 1D acoustic medium.
    % Stores function handles defining kappa(x) and beta(x), with static
    % factory methods for constructing common material models.

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

    %%
    methods
        %%%%%%%%
        function obj = AcousticMedium1D(kappaFunction,betaFunction,options)
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
            obj.betaFunction  = betaFunction;
            obj.kappaLabel = options.KappaLabel;
            obj.kappaUnits = options.KappaUnits;
            obj.betaLabel = options.BetaLabel;
            obj.betaUnits = options.BetaUnits;
            obj.waveSpeedLabel = options.WaveSpeedLabel;
            obj.waveSpeedUnits = options.WaveSpeedUnits;
        end

        %%%%%%%%
        function values = evaluateKappa(obj,x)
            % Evaluate kappa at the specified spatial coordinates.
            % Input:
            %   x - Numeric array of spatial coordinates
            % Output:
            %   values - Array of kappa values with the same size as x
            values = obj.kappaFunction(x);
        end

        %%%%%%%%
        function values = evaluateBeta(obj,x)
            % Evaluate beta at the specified spatial coordinates.
            % Input:
            %   x - Numeric array of spatial coordinates
            % Output:
            %   values - Array of beta values with the same size as x
            values = obj.betaFunction(x);
        end

        %%%%%%%%
        function values = evaluateWaveSpeed(obj,x)
            % Evaluate c(x) = sqrt(kappa(x).*beta(x)).
            values = sqrt( ...
                obj.evaluateKappa(x).*obj.evaluateBeta(x));
        end
    end

    %%
    methods (Static)
        %%%%%%%%
        function obj = constant(options)
            arguments
                options.KappaValue (1,1) double ...
                    {mustBePositive,mustBeFinite} = 1
                options.BetaValue  (1,1) double ...
                    {mustBePositive,mustBeFinite} = 1
            end

            kappa = AcousticMedium1D.makeConstantFunction( ...
                options.KappaValue);

            beta = AcousticMedium1D.makeConstantFunction( ...
                options.BetaValue);

            obj = AcousticMedium1D(kappa,beta);
        end

        %%%%%%%%
        function obj = linear(options)
            arguments
                options.KappaRange (1,2) double ...
                    {mustBePositive,mustBeFinite} = [1,2]
                options.BetaRange (1,2) double ...
                    {mustBePositive,mustBeFinite} = [1,2]
                options.Domain (1,2) double ...
                    {mustBeFinite,mustBeIncreasingRange} = [0,1]
            end

            kappa = AcousticMedium1D.makeLinearFunction( ...
                options.KappaRange,options.Domain);

            beta = AcousticMedium1D.makeLinearFunction( ...
                options.BetaRange,options.Domain);

            obj = AcousticMedium1D(kappa,beta);
        end

        %%%%%%%%
        function obj = smooth(options)
            arguments
                options.KappaRange (1,2) double ...
                    {mustBePositive,mustBeFinite,mustBeIncreasingRange} = [1,2]
                options.BetaRange (1,2) double ...
                    {mustBePositive,mustBeFinite,mustBeIncreasingRange} = [1,2]
                options.Domain (1,2) double ...
                    {mustBeFinite,mustBeIncreasingRange} = [0,1]
            end

            kappa = AcousticMedium1D.makeSmoothFunction( ...
                options.KappaRange,options.Domain);

            beta = AcousticMedium1D.makeSmoothFunction( ...
                options.BetaRange,options.Domain);

            obj = AcousticMedium1D(kappa,beta);
        end

        %%%%%%%%
        function obj = gaussianPeak(options)
            arguments
                options.KappaRange (1,2) double ...
                    {mustBePositive,mustBeFinite,mustBeIncreasingRange} = [1,2]
                options.BetaRange (1,2) double ...
                    {mustBePositive,mustBeFinite,mustBeIncreasingRange} = [1,2]
                options.Center (1,1) double = 1.5
                options.Width  (1,1) double {mustBePositive} = 1
            end

            kappa = AcousticMedium1D.makeGaussianPeakFunction( ...
                options.KappaRange,options.Center,options.Width);

            beta = AcousticMedium1D.makeGaussianPeakFunction( ...
                options.BetaRange,options.Center,options.Width);

            obj = AcousticMedium1D(kappa,beta);
        end

        %%%%%%%%
        function obj = gaussianDip(options)
            arguments
                options.KappaRange (1,2) double ...
                    {mustBePositive,mustBeFinite,mustBeIncreasingRange} = [1,2]
                options.BetaRange (1,2) double ...
                    {mustBePositive,mustBeFinite,mustBeIncreasingRange} = [1,2]
                options.Center (1,1) double = 1.5
                options.Width  (1,1) double {mustBePositive} = 1
            end

            kappa = AcousticMedium1D.makeGaussianDipFunction( ...
                options.KappaRange,options.Center,options.Width);

            beta = AcousticMedium1D.makeGaussianDipFunction( ...
                options.BetaRange,options.Center,options.Width);

            obj = AcousticMedium1D(kappa,beta);
        end
    end

    %%
    methods (Static, Access = private)
        %%%%%%%%
        function f = makeConstantFunction(constValue)
            f = @(x) constValue*ones(size(x));
        end

        %%%%%%%%
        function f = makeLinearFunction(valueRange,domain)
            valueMin = valueRange(1);
            valueMax = valueRange(2);

            xLeft = domain(1);
            L     = domain(2) - domain(1);
            slope = (valueMax-valueMin)/L;

            f = @(x) valueMin + slope.*(x-xLeft);
        end

        %%%%%%%%
        function f = makeSmoothFunction(valueRange,domain)
            valueMin = valueRange(1);
            valueMax = valueRange(2);

            average   = 0.5*(valueMin + valueMax);
            amplitude = 0.5*(valueMax - valueMin);

            xLeft = domain(1);
            L     = domain(2) - domain(1);

            f = @(x) average + ...
                amplitude.*sin(2*pi*(x-xLeft)/L);
        end

        %%%%%%%%
        function f = makeGaussianPeakFunction(valueRange,center,width)
            valueMin = valueRange(1);
            valueMax = valueRange(2);

            f = @(x) valueMin + (valueMax-valueMin).* ...
                exp(-((x-center)./width).^2);
        end

        %%%%%%%%
        function f = makeGaussianDipFunction(valueRange,center,width)
            valueMin = valueRange(1);
            valueMax = valueRange(2);

            f = @(x) valueMax - (valueMax-valueMin).* ...
                exp(-((x-center)./width).^2);
        end
    end
end

function mustBeIncreasingRange(value)
if value(1) > value(2)
    error("AcousticMedium1D:InvalidRange", ...
        "Range must have the form [minimum, maximum].");
end
end

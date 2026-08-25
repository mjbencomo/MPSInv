classdef PML2D
    % PML2D stores parameters for four 2D absorbing layers.

    properties
        widthLeft (1,1) double {mustBeNonnegative,mustBeFinite} = 0
        widthRight (1,1) double {mustBeNonnegative,mustBeFinite} = 0
        widthBottom (1,1) double {mustBeNonnegative,mustBeFinite} = 0
        widthTop (1,1) double {mustBeNonnegative,mustBeFinite} = 0
        reflectionCoefficient (1,1) double ...
            {mustBePositive,mustBeFinite} = 1e-6
        polynomialDegree (1,1) double ...
            {mustBeInteger,mustBePositive} = 3
    end

    methods
        function obj = PML2D(options)
            arguments
                options.WidthLeft (1,1) double ...
                    {mustBeNonnegative,mustBeFinite} = 0
                options.WidthRight (1,1) double ...
                    {mustBeNonnegative,mustBeFinite} = 0
                options.WidthBottom (1,1) double ...
                    {mustBeNonnegative,mustBeFinite} = 0
                options.WidthTop (1,1) double ...
                    {mustBeNonnegative,mustBeFinite} = 0
                options.ReflectionCoefficient (1,1) double ...
                    {mustBePositive,mustBeFinite} = 1e-6
                options.PolynomialDegree (1,1) double ...
                    {mustBeInteger,mustBePositive} = 3
            end

            if options.ReflectionCoefficient >= 1
                error('PML2D:InvalidReflectionCoefficient', ...
                    'ReflectionCoefficient must be strictly less than one.');
            end

            obj.widthLeft = options.WidthLeft;
            obj.widthRight = options.WidthRight;
            obj.widthBottom = options.WidthBottom;
            obj.widthTop = options.WidthTop;
            obj.reflectionCoefficient = options.ReflectionCoefficient;
            obj.polynomialDegree = options.PolynomialDegree;
        end

        function [nLeft,nRight,nBottom,nTop,widths] = cellWidths(obj,h)
            % Round requested widths upward to whole grid cells.
            arguments
                obj
                h (1,2) double {mustBePositive,mustBeFinite}
            end
            nLeft = ceil(obj.widthLeft/h(1));
            nRight = ceil(obj.widthRight/h(1));
            nBottom = ceil(obj.widthBottom/h(2));
            nTop = ceil(obj.widthTop/h(2));
            widths = [h(1)*nLeft,h(1)*nRight, ...
                h(2)*nBottom,h(2)*nTop];
        end

        function sigma = dampingX(obj,x,xDomain,widths,cReference)
            arguments
                obj
                x double {mustBeFinite}
                xDomain (1,2) double {mustBeFinite}
                widths (1,2) double {mustBeNonnegative,mustBeFinite}
                cReference (1,1) double {mustBePositive,mustBeFinite}
            end
            sigma = obj.directionalDamping( ...
                x,xDomain,widths,cReference,"x");
        end

        function sigma = dampingY(obj,y,yDomain,widths,cReference)
            arguments
                obj
                y double {mustBeFinite}
                yDomain (1,2) double {mustBeFinite}
                widths (1,2) double {mustBeNonnegative,mustBeFinite}
                cReference (1,1) double {mustBePositive,mustBeFinite}
            end
            sigma = obj.directionalDamping( ...
                y,yDomain,widths,cReference,"y");
        end
    end

    methods (Access = private)
        function sigma = directionalDamping(obj,coordinates,domain, ...
                widths,cReference,direction)
            if domain(1) >= domain(2)
                error('PML2D:InvalidDomain', ...
                    '%s-domain endpoints must be strictly increasing.',direction);
            end

            sigma = zeros(size(coordinates));
            degree = obj.polynomialDegree;
            sigmaMax = zeros(1,2);
            active = widths > 0;
            sigmaMax(active) = -(degree+1)*cReference ...
                *log(obj.reflectionCoefficient)./(2*widths(active));

            coordinateScale = max([1;abs(coordinates(:)); ...
                abs(domain(:));widths(:)]);
            interfaceTolerance = 100*eps(coordinateScale);

            lower = coordinates < domain(1)-interfaceTolerance;
            if any(lower,'all')
                if widths(1) == 0
                    error('PML2D:ZeroLowerWidth', ...
                        ['Points lie below the lower %s-interface of a ' ...
                         'zero-width PML.'],direction);
                end
                eta = (domain(1)-coordinates(lower))/widths(1);
                sigma(lower) = sigmaMax(1)*eta.^degree;
            end

            upper = coordinates > domain(2)+interfaceTolerance;
            if any(upper,'all')
                if widths(2) == 0
                    error('PML2D:ZeroUpperWidth', ...
                        ['Points lie above the upper %s-interface of a ' ...
                         'zero-width PML.'],direction);
                end
                eta = (coordinates(upper)-domain(2))/widths(2);
                sigma(upper) = sigmaMax(2)*eta.^degree;
            end
        end
    end
end

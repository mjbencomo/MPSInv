classdef PML1D
    % PML1D stores parameters for left and right 1D absorbing layers.

    properties
        widthLeft (1,1) double {mustBeNonnegative,mustBeFinite} = 0
        widthRight (1,1) double {mustBeNonnegative,mustBeFinite} = 0
        reflectionCoefficient (1,1) double {mustBePositive,mustBeFinite} = 1e-6
        polynomialDegree (1,1) double {mustBeInteger,mustBePositive} = 3
    end

    methods
        function obj = PML1D(options)
            arguments
                options.WidthLeft (1,1) double ...
                    {mustBeNonnegative,mustBeFinite} = 0
                options.WidthRight (1,1) double ...
                    {mustBeNonnegative,mustBeFinite} = 0
                options.ReflectionCoefficient (1,1) double ...
                    {mustBePositive,mustBeFinite} = 1e-6
                options.PolynomialDegree (1,1) double ...
                    {mustBeInteger,mustBePositive} = 3
            end

            if options.ReflectionCoefficient >= 1
                error('PML1D:InvalidReflectionCoefficient', ...
                    'ReflectionCoefficient must be strictly less than one.');
            end

            obj.widthLeft = options.WidthLeft;
            obj.widthRight = options.WidthRight;
            obj.reflectionCoefficient = options.ReflectionCoefficient;
            obj.polynomialDegree = options.PolynomialDegree;
        end

        function [nLeft,nRight,widths] = cellWidths(obj,dx)
            % Round requested widths upward to whole grid cells.
            arguments
                obj
                dx (1,1) double {mustBePositive,mustBeFinite}
            end

            nLeft = ceil(obj.widthLeft/dx);
            nRight = ceil(obj.widthRight/dx);
            widths = dx*[nLeft,nRight];
        end

        function sigma = damping(obj,x,domain,widths,cReference)
            % Evaluate the polynomial damping profile at arbitrary points.
            arguments
                obj
                x double {mustBeFinite}
                domain (1,2) double {mustBeFinite}
                widths (1,2) double {mustBeNonnegative,mustBeFinite}
                cReference (1,1) double {mustBePositive,mustBeFinite}
            end

            if domain(1) >= domain(2)
                error('PML1D:InvalidDomain', ...
                    'Domain endpoints must be strictly increasing.');
            end

            sigma = zeros(size(x));
            degree = obj.polynomialDegree;

            % R is interpreted as the attenuation after propagation to the
            % outer boundary and back through the PML.
            sigmaMax = zeros(1,2);
            active = widths > 0;
            sigmaMax(active) = -(degree+1)*cReference ...
                *log(obj.reflectionCoefficient)./(2*widths(active));

            % Grid construction can place a nominal interface point a few
            % ulps outside the physical domain. Treat such points as lying
            % exactly on the interface so sigma is identically zero there.
            coordinateScale = max([1;abs(x(:));abs(domain(:));widths(:)]);
            interfaceTolerance = 100*eps(coordinateScale);

            left = x < domain(1)-interfaceTolerance;
            if any(left,'all')
                if widths(1) == 0
                    error('PML1D:ZeroLeftWidth', ...
                        'Points lie left of a zero-width PML.');
                end
                eta = (domain(1)-x(left))/widths(1);
                sigma(left) = sigmaMax(1)*eta.^degree;
            end

            right = x > domain(2)+interfaceTolerance;
            if any(right,'all')
                if widths(2) == 0
                    error('PML1D:ZeroRightWidth', ...
                        'Points lie right of a zero-width PML.');
                end
                eta = (x(right)-domain(2))/widths(2);
                sigma(right) = sigmaMax(2)*eta.^degree;
            end
        end
    end
end

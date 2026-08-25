% Grid1DUnifPrimal represents a uniform primal grid on a 1D domain.
% Grid points include both endpoints of the domain.

classdef Grid1DUnifPrimal < Grid1DUnif
    methods
        function obj = Grid1DUnifPrimal(dom,N)
            arguments
                dom (1,2) double
                N   (1,1) double ...
                    {mustBeInteger,mustBeGreaterThan(N,1)}
            end

            if dom(1) >= dom(2)
                error('Grid1DUnifPrimal:InvalidDomain', ...
                    'Domain must satisfy dom(1) < dom(2).');
            end

            h = (dom(2)-dom(1))/(N-1);
            pts = linspace(dom(1),dom(2),N).';
            obj@Grid1DUnif(dom,pts,h);
        end
    end
end
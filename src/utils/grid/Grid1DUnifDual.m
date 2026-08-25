% Grid1DUnifDual represents the dual grid associated with a uniform primal grid.
% Grid points are located midway between consecutive primal grid points.

classdef Grid1DUnifDual < Grid1DUnif
    methods
        function obj = Grid1DUnifDual(primal)
            arguments
                primal (1,1) Grid1DUnifPrimal
            end

            dom = primal.dom;
            h   = primal.h;
            pts = 0.5*( ...
                primal.pts(1:end-1) + ...
                primal.pts(2:end) );

            obj@Grid1DUnif(dom,pts,h);
        end
    end
end
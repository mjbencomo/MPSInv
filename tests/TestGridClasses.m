classdef TestGridClasses < matlab.unittest.TestCase
    % Tests for the one- and two-dimensional grid classes.

    properties
        PrimalGrid
    end

    methods (TestMethodSetup)
        function createPrimalGrid(testCase)
            testCase.PrimalGrid = Grid1DUnifPrimal([0,1],5);
        end
    end

    methods (Test)
        function primalGrid(testCase)
            g = testCase.PrimalGrid;

            testCase.verifyEqual(g.N,5);
            testCase.verifyEqual(g.h,0.25);
            testCase.verifyEqual(g.pts,[0;0.25;0.5;0.75;1], ...
                AbsTol=1e-14);
        end

        function dualGrid(testCase)
            gp = testCase.PrimalGrid;
            gd = Grid1DUnifDual(gp);

            testCase.verifyEqual(gd.N,4);
            testCase.verifyEqual(gd.h,gp.h);
            testCase.verifyEqual(gd.pts,[0.125;0.375;0.625;0.875], ...
                AbsTol=1e-14);
        end

        function gridSpace2D(testCase)
            x = testCase.PrimalGrid;
            y = Grid1DUnifPrimal([0,2],9);
            g = GridSpace2D(x,y);

            testCase.verifyEqual(g.dim,2);
            testCase.verifyEqual(g.N,[5,9]);
            testCase.verifyEqual(g.h,[0.25,0.25]);
            testCase.verifyEqual(g.dom,[0,1;0,2]);
        end

        function invalidDomain(testCase)
            testCase.verifyError( ...
                @() Grid1D([1,0],[0;0.5;1]), ...
                'Grid1D:InvalidDomain');
        end
    end
end

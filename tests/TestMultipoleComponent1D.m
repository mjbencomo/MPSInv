classdef TestMultipoleComponent1D < matlab.unittest.TestCase
    properties
        Grid
    end

    methods (TestMethodSetup)
        function createGrid(testCase)
            primal = Grid1DUnifPrimal([0,1],41);
            testCase.Grid = GridSpace1D(primal);
        end
    end

    methods (Test)
        function storesComponentData(testCase)
            component = MultipoleComponent1D( ...
                0.4,1,TargetField="velocity");

            testCase.verifyEqual(component.location,0.4);
            testCase.verifyEqual(component.derivativeOrder,1);
            testCase.verifyEqual(component.targetField,"velocity");
            testCase.verifyTrue(isa(component,'MultipoleComponent'));
        end

        function usesPressureAsDefaultTarget(testCase)
            component = MultipoleComponent1D(0.4,0);
            testCase.verifyEqual(component.targetField,"pressure");
        end

        function storesNoApproximationOrder(testCase)
            component = MultipoleComponent1D(0.4,0);
            testCase.verifyFalse(isprop(component,'approximationOrder'));
        end

        function returnsDescription(testCase)
            component = MultipoleComponent1D( ...
                0.4,1,TargetField="velocity");

            testCase.verifyEqual(component.description(), ...
                "velocity multipole, derivative order 1, at x = 0.4");
        end

        function stencilSatisfiesMomentConditions(testCase)
            xc = 0.437;
            s = 2;
            q = 4;
            component = MultipoleComponent1D(xc,s);

            [indices,weights,values] = ...
                component.createStencil(testCase.Grid,q);

            x = testCase.Grid.x.pts;
            h = testCase.Grid.h;
            testCase.verifyNumElements(indices,q+s);
            testCase.verifyNumElements(weights,q+s);
            testCase.verifyEqual(values(indices),weights);

            for k = 0:q+s-1
                actual = h*sum(weights.*(x(indices)-xc).^k);
                expected = 0;
                if k == s
                    expected = (-1)^s*factorial(s);
                end
                tolerance = 1e4*eps*max(1,abs(expected));
                testCase.verifyEqual(actual,expected, ...
                    AbsTol=tolerance);
            end
        end

        function acceptsOrderAtDiscretizationTime(testCase)
            component = MultipoleComponent1D(0.437,1);

            [indices2,~] = component.createStencil(testCase.Grid,2);
            [indices4,~] = component.createStencil(testCase.Grid,4);

            testCase.verifyNumElements(indices2,3);
            testCase.verifyNumElements(indices4,5);
        end
    end
end

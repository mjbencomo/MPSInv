classdef TestMultipoleComponent1D < matlab.unittest.TestCase
<<<<<<< HEAD
    methods (Test)
        function storesComponentDescription(testCase)
            component = MultipoleComponent1D( ...
                0.4,1, ...
                TargetField="velocity", ...
                ApproximationOrder=6);

            testCase.verifyTrue(isa(component,'MultipoleComponent'));
            testCase.verifyEqual(component.location,0.4);
            testCase.verifyEqual(component.derivativeOrder,1);
            testCase.verifyEqual(component.targetField,"velocity");
            testCase.verifyEqual(component.approximationOrder,6);
            testCase.verifyEqual(component.description(), ...
                "velocity multipole, derivative order 1, at x = 0.4");
        end

        function defaultsAreApplied(testCase)
            component = MultipoleComponent1D(0.4,0);

            testCase.verifyEqual(component.targetField,"pressure");
            testCase.verifyEqual(component.approximationOrder,4);
        end

        function constructsTermWithTimeFunction(testCase)
            component = MultipoleComponent1D( ...
                0.4,1,TargetField="velocity");
            term = component.withTimeFunction(@(t) 1+t,Amplitude=-2);

            testCase.verifyTrue(isa(term,'MultipoleTerm1D'));
            testCase.verifyEqual(term.location,component.location);
            testCase.verifyEqual( ...
                term.derivativeOrder,component.derivativeOrder);
            testCase.verifyEqual(term.targetField,component.targetField);
            testCase.verifyEqual(term.evaluateTime(0.5),-3);
        end

        function componentCanBeReused(testCase)
            component = MultipoleComponent1D(0.4,0);
            firstTerm = component.withTimeFunction(@(t) t);
            secondTerm = component.withTimeFunction(@(t) t.^2);

            testCase.verifyEqual(firstTerm.evaluateTime(0.5),0.5);
            testCase.verifyEqual(secondTerm.evaluateTime(0.5),0.25);
            testCase.verifyEqual(component.location,0.4);
        end

        function termUsesComponentApproximationOrder(testCase)
            primal = Grid1DUnifPrimal([0,1],41);
            pressureGrid = GridSpace1D(primal);
            velocityGrid = GridSpace1D(Grid1DUnifDual(primal));
            component = MultipoleComponent1D( ...
                0.45,1,ApproximationOrder=4);
            term = component.withTimeFunction(@(t) 1+t);

            source = term.discretize( ...
                pressureGrid,velocityGrid, ...
                ApproximationOrder=component.approximationOrder);

            testCase.verifyEqual(source.numberOfPressureTerms,1);
=======
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
>>>>>>> b0d5b1dcea97a019fdcf42ee218d06a571246d9f
        end
    end
end

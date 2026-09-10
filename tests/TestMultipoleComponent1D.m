classdef TestMultipoleComponent1D < matlab.unittest.TestCase
    methods (Test)
        function storesComponentDescription(testCase)
            component = MultipoleComponent1D( ...
                0.4,1, ...
                TargetField="velocity", ...
                ApproximationOrder=6, ...
                Label="Velocity dipole");

            testCase.verifyTrue(isa(component,'MultipoleComponent'));
            testCase.verifyEqual(component.location,0.4);
            testCase.verifyEqual(component.derivativeOrder,1);
            testCase.verifyEqual(component.targetField,"velocity");
            testCase.verifyEqual(component.approximationOrder,6);
            testCase.verifyEqual(component.label,"Velocity dipole");
        end

        function defaultsAreApplied(testCase)
            component = MultipoleComponent1D(0.4,0);

            testCase.verifyEqual(component.targetField,"pressure");
            testCase.verifyEqual(component.approximationOrder,4);
            testCase.verifyEqual(component.label,"");
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
        end
    end
end

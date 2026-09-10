classdef TestMultipoleComponent2D < matlab.unittest.TestCase
    methods (Test)
        function storesComponentDescription(testCase)
            component = MultipoleComponent2D( ...
                [0.4,0.6],[1,2], ...
                TargetField="velocityY", ...
                ApproximationOrder=6);

            testCase.verifyTrue(isa(component,'MultipoleComponent'));
            testCase.verifyEqual(component.location,[0.4,0.6]);
            testCase.verifyEqual(component.derivativeOrder,[1,2]);
            testCase.verifyEqual(component.targetField,"velocityY");
            testCase.verifyEqual(component.approximationOrder,6);
            testCase.verifyEqual(component.description(), ...
                ["velocityY multipole, derivative order [1,2], "+ ...
                 "at [x,y] = [0.4,0.6]"]);
        end

        function defaultsAreApplied(testCase)
            component = MultipoleComponent2D([0.4,0.6],[0,0]);

            testCase.verifyEqual(component.targetField,"pressure");
            testCase.verifyEqual(component.approximationOrder,4);
        end

        function constructsTermWithTimeFunction(testCase)
            component = MultipoleComponent2D( ...
                [0.4,0.6],[1,2],TargetField="velocityY");
            term = component.withTimeFunction(@(t) 1+t,Amplitude=-2);

            testCase.verifyTrue(isa(term,'MultipoleTerm2D'));
            testCase.verifyEqual(term.location,component.location);
            testCase.verifyEqual( ...
                term.derivativeOrder,component.derivativeOrder);
            testCase.verifyEqual(term.targetField,component.targetField);
            testCase.verifyEqual(term.evaluateTime(0.5),-3);
        end

        function componentCanBeReused(testCase)
            component = MultipoleComponent2D([0.4,0.6],[0,1]);
            firstTerm = component.withTimeFunction(@(t) t);
            secondTerm = component.withTimeFunction(@(t) t.^2);

            testCase.verifyEqual(firstTerm.evaluateTime(0.5),0.5);
            testCase.verifyEqual(secondTerm.evaluateTime(0.5),0.25);
            testCase.verifyEqual(component.location,[0.4,0.6]);
        end
    end
end

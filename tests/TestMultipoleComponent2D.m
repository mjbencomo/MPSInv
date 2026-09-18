classdef TestMultipoleComponent2D < matlab.unittest.TestCase
<<<<<<< HEAD
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
=======
    properties
        Grid
    end

    methods (TestMethodSetup)
        function createGrid(testCase)
            xPrimal = Grid1DUnifPrimal([0,1],41);
            yPrimal = Grid1DUnifPrimal([0,1],41);
            testCase.Grid = GridSpace2D(xPrimal,yPrimal);
        end
    end

    methods (Test)
        function storesComponentData(testCase)
            component = MultipoleComponent2D( ...
                [0.4,0.6],[1,2],TargetField="velocityY");

            testCase.verifyEqual(component.location,[0.4,0.6]);
            testCase.verifyEqual(component.derivativeOrder,[1,2]);
            testCase.verifyEqual(component.targetField,"velocityY");
            testCase.verifyTrue(isa(component,'MultipoleComponent'));
        end

        function usesPressureAsDefaultTarget(testCase)
            component = MultipoleComponent2D([0.4,0.6],[0,0]);
            testCase.verifyEqual(component.targetField,"pressure");
        end

        function storesNoApproximationOrder(testCase)
            component = MultipoleComponent2D([0.4,0.6],[0,0]);
            testCase.verifyFalse(isprop(component,'approximationOrder'));
        end

        function returnsDescription(testCase)
            component = MultipoleComponent2D( ...
                [0.4,0.6],[1,2],TargetField="velocityY");

            testCase.verifyEqual(component.description(), ...
                ["velocityY multipole, derivative order [1,2], " + ...
                 "at [x,y] = [0.4,0.6]"]);
        end

        function stencilSatisfiesTensorMomentConditions(testCase)
            q = 4;
            order = [1,2];
            location = [0.437,0.463];
            component = MultipoleComponent2D(location,order);

            [indices,weights,xIndices,yIndices] = ...
                component.createStencil(testCase.Grid,q);

            expectedSupport = prod(q+order);
            testCase.verifyNumElements(indices,expectedSupport);
            testCase.verifyNumElements(weights,expectedSupport);
            testCase.verifyNumElements(xIndices,q+order(1));
            testCase.verifyNumElements(yIndices,q+order(2));

            [ix,iy] = ind2sub(testCase.Grid.N,indices);
            x = testCase.Grid.x.pts(ix);
            y = testCase.Grid.y.pts(iy);
            hx = testCase.Grid.h(1);
            hy = testCase.Grid.h(2);
            rx = (x-location(1))/hx;
            ry = (y-location(2))/hy;
            scaledWeights = weights* ...
                hx^(order(1)+1)*hy^(order(2)+1);

            for kx = 0:q+order(1)-1
                for ky = 0:q+order(2)-1
                    actual = sum( ...
                        scaledWeights.*rx.^kx.*ry.^ky);
                    expectedX = 0;
                    expectedY = 0;
                    if kx == order(1)
                        expectedX = (-1)^order(1)* ...
                            factorial(order(1));
                    end
                    if ky == order(2)
                        expectedY = (-1)^order(2)* ...
                            factorial(order(2));
                    end
                    expected = expectedX*expectedY;
                    tolerance = 1e5*eps*max(1,abs(expected));
                    testCase.verifyEqual(actual,expected, ...
                        AbsTol=tolerance);
                end
            end
        end

        function acceptsOrderAtDiscretizationTime(testCase)
            component = MultipoleComponent2D( ...
                [0.437,0.463],[1,2]);

            [indices2,~] = component.createStencil(testCase.Grid,2);
            [indices4,~] = component.createStencil(testCase.Grid,4);

            testCase.verifyNumElements(indices2,(2+1)*(2+2));
            testCase.verifyNumElements(indices4,(4+1)*(4+2));
>>>>>>> b0d5b1dcea97a019fdcf42ee218d06a571246d9f
        end
    end
end

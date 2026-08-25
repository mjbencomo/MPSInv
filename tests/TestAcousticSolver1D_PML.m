classdef TestAcousticSolver1D_PML < matlab.unittest.TestCase
    properties
        grid
        medium
        pml
    end

    methods (TestMethodSetup)
        function createFixtures(testCase)
            testCase.grid = GridSpace1D(Grid1DUnifPrimal([0,1],41));
            testCase.medium = AcousticMedium1D.constant( ...
                KappaValue=1,BetaValue=1);
            testCase.pml = PML1D(WidthLeft=0.2,WidthRight=0.2, ...
                ReflectionCoefficient=1e-6,PolynomialDegree=3);
        end
    end

    methods (Test)
        function dampingVanishesInPhysicalDomain(testCase)
            solver = testCase.makeSolver(2,0:0.01:0.1,0);
            i = solver.physicalPressureIndices;
            testCase.verifyEqual(solver.sigmaP(i),zeros(numel(i),1));
            testCase.verifyGreaterThan(solver.sigmaP(1),0);
            testCase.verifyGreaterThan(solver.sigmaP(end),0);
        end

        function extendedGridUsesWholeCells(testCase)
            solver = testCase.makeSolver(2,[0,0.01],0);
            testCase.verifyEqual(solver.extendedGrid.h,testCase.grid.h, ...
                'AbsTol',10*eps);
            testCase.verifyEqual(solver.extendedGrid.dom,[-0.2,1.2], ...
                'AbsTol',10*eps);
        end

        function zeroDataRemainZeroForBothOrders(testCase)
            for order = [2,4]
                solver = testCase.makeSolver(order,0:0.01:0.1,[0,0.05,0.1]);
                [P,V] = solver.solve(@(x) 0*x,@(x) 0*x);
                testCase.verifyEqual(P.values,zeros(size(P.values)));
                testCase.verifyEqual(V.values,zeros(size(V.values)));
            end
        end

        function physicalInterfacesAreNotDirichletBoundaries(testCase)
            solver = testCase.makeSolver(2,[0,0.005],0);
            [P,~] = solver.solve(@(x) ones(size(x)),@(x) 0*x);
            testCase.verifyEqual(P.values,ones(testCase.grid.N,1), ...
                'AbsTol',1e-14);
        end

        function sourceIsRestrictedToPhysicalDomain(testCase)
            solver = testCase.makeSolver(2,[0,0.005],0.0025);
            source = AcousticSource1D(@(x) 2,@(t) 1);
            [P,~] = solver.solve(@(x) 0*x,@(x) 0*x,Source=source);
            testCase.verifyGreaterThan(max(abs(P.values)),0);
        end

        function sourceSuperpositionIsPreserved(testCase)
            times = 0:0.005:0.025;
            solver = testCase.makeSolver(2,times,times(end));
            source1 = AcousticSource1D(@(x) sin(pi*x),@(t) 1+t);
            source2 = AcousticSource1D(@(x) x.*(1-x),@(t) cos(t));
            combined = AcousticSource1D.combine({source1,source2});

            [P1,V1] = solver.solve(@(x) 0*x,@(x) 0*x,Source=source1);
            [P2,V2] = solver.solve(@(x) 0*x,@(x) 0*x,Source=source2);
            [P,V] = solver.solve(@(x) 0*x,@(x) 0*x,Source=combined);

            testCase.verifyEqual(P.values,P1.values+P2.values,AbsTol=1e-13);
            testCase.verifyEqual(V.values,V1.values+V2.values,AbsTol=1e-13);
        end

        function acceptsLocalizedTermsOnPhysicalGrids(testCase)
            times = 0:0.005:0.025;
            solver = testCase.makeSolver(2,times,times(end));
            physicalDualGrid = GridSpace1D( ...
                Grid1DUnifDual(testCase.grid.x));
            source = AcousticSource1D.zero();
            source.addLocalizedPressureTerm( ...
                testCase.grid,21,2,@(t) 1);
            source.addLocalizedVelocityTerm( ...
                physicalDualGrid,20,-1,@(t) 1);

            [P,V] = solver.solve(@(x) 0*x,@(x) 0*x,Source=source);

            testCase.verifyGreaterThan(max(abs(P.values),[],'all'),0);
            testCase.verifyGreaterThan(max(abs(V.values),[],'all'),0);
        end

        function rejectsLocalizedTermBuiltOnExtendedGrid(testCase)
            solver = testCase.makeSolver(2,[0,0.005],0.005);
            source = AcousticSource1D.zero();
            source.addLocalizedPressureTerm( ...
                solver.extendedGrid,3,1,@(t) 1);

            testCase.verifyError(@() solver.solve( ...
                @(x) 0*x,@(x) 0*x,Source=source), ...
                'AcousticSource1D:LocalizedGridMismatch');
        end

        function rejectsInvalidReflectionCoefficient(testCase)
            testCase.verifyError(@() PML1D(ReflectionCoefficient=1), ...
                'PML1D:InvalidReflectionCoefficient');
        end
    end

    methods (Access = private)
        function solver = makeSolver(testCase,order,times,outputTimes)
            solver = AcousticSolver1D_PML( ...
                testCase.grid,times,testCase.grid,outputTimes, ...
                testCase.medium,testCase.pml,SpatialOrder=order);
        end
    end
end

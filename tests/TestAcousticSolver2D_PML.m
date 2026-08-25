classdef TestAcousticSolver2D_PML < matlab.unittest.TestCase
    % Tests for the split-pressure 2D acoustic PML solver.

    properties
        grid
        medium
        pml
    end

    methods (TestMethodSetup)
        function createFixtures(testCase)
            testCase.grid = GridSpace2D( ...
                Grid1DUnifPrimal([0,1],41), ...
                Grid1DUnifPrimal([0,1],41));
            testCase.medium = AcousticMedium2D.constant( ...
                KappaValue=1,BetaValue=1);
            testCase.pml = PML2D( ...
                WidthLeft=0.2,WidthRight=0.2, ...
                WidthBottom=0.2,WidthTop=0.2, ...
                ReflectionCoefficient=1e-6,PolynomialDegree=3);
        end
    end

    methods (Test)
        function directionalDampingVanishesInPhysicalDomain(testCase)
            solver = testCase.makeSolver(2,[0,0.005],0);
            ix = solver.physicalPressureXIndices;
            iy = solver.physicalPressureYIndices;

            testCase.verifyEqual(solver.sigmaPX(ix,iy), ...
                zeros(numel(ix),numel(iy)));
            testCase.verifyEqual(solver.sigmaPY(ix,iy), ...
                zeros(numel(ix),numel(iy)));
            testCase.verifyGreaterThan(solver.sigmaPX(1,iy(20)),0);
            testCase.verifyGreaterThan(solver.sigmaPX(end,iy(20)),0);
            testCase.verifyGreaterThan(solver.sigmaPY(ix(20),1),0);
            testCase.verifyGreaterThan(solver.sigmaPY(ix(20),end),0);
        end

        function dampingProfilesAreDirectional(testCase)
            solver = testCase.makeSolver(2,[0,0.005],0);
            testCase.verifyEqual(solver.sigmaPX(:,1),solver.sigmaPX(:,end));
            testCase.verifyEqual(solver.sigmaPY(1,:),solver.sigmaPY(end,:));
        end

        function extendedGridUsesWholeCells(testCase)
            solver = testCase.makeSolver(2,[0,0.005],0);
            testCase.verifyEqual(solver.extendedGrid.h,testCase.grid.h, ...
                AbsTol=10*eps);
            testCase.verifyEqual(solver.extendedGrid.dom, ...
                [-0.2,1.2;-0.2,1.2],AbsTol=10*eps);
            testCase.verifyEqual(solver.extendedGrid.N,[57,57]);
        end

        function zeroDataRemainZeroForBothOrders(testCase)
            zero = @(x,y) zeros(size(x));
            for order = [2,4]
                solver = testCase.makeSolver( ...
                    order,0:0.005:0.02,[0,0.01,0.02]);
                [P,Vx,Vy] = solver.solve(zero,zero,zero);
                testCase.verifyEqual(P.values,zeros(size(P.values)));
                testCase.verifyEqual(Vx.values,zeros(size(Vx.values)));
                testCase.verifyEqual(Vy.values,zeros(size(Vy.values)));
            end
        end

        function physicalInterfacesAreNotDirichletBoundaries(testCase)
            solver = testCase.makeSolver(2,[0,0.005],0);
            zero = @(x,y) zeros(size(x));
            [P,~,~] = solver.solve(@(x,y) ones(size(x)),zero,zero);
            testCase.verifyEqual(P.values,ones(testCase.grid.N), ...
                AbsTol=1e-14);
        end

        function sourceSuperpositionIsPreserved(testCase)
            times = 0:0.005:0.02;
            solver = testCase.makeSolver(2,times,times(end));
            source1 = AcousticSource2D( ...
                @(x,y) sin(pi*x).*sin(pi*y),@(t) 1+t);
            source2 = AcousticSource2D( ...
                @(x,y) x.*(1-x).*y.*(1-y),@(t) cos(t));
            combined = AcousticSource2D.combine({source1,source2});
            zero = @(x,y) zeros(size(x));

            [P1,Vx1,Vy1] = solver.solve(zero,zero,zero,Source=source1);
            [P2,Vx2,Vy2] = solver.solve(zero,zero,zero,Source=source2);
            [P,Vx,Vy] = solver.solve(zero,zero,zero,Source=combined);

            testCase.verifyEqual(P.values,P1.values+P2.values,AbsTol=1e-13);
            testCase.verifyEqual(Vx.values,Vx1.values+Vx2.values,AbsTol=1e-13);
            testCase.verifyEqual(Vy.values,Vy1.values+Vy2.values,AbsTol=1e-13);
        end

        function acceptsLocalizedTermsOnPhysicalGrids(testCase)
            times = 0:0.005:0.02;
            solver = testCase.makeSolver(2,times,times(end));
            [~,Bx,By] = sampleMedium2D(testCase.grid,testCase.medium);
            source = AcousticSource2D.zero();
            pressureIndex = sub2ind(testCase.grid.N,21,21);
            source.addLocalizedPressureTerm( ...
                testCase.grid,pressureIndex,2,@(t) 1);
            source.addLocalizedVelocityXTerm(Bx.grid,20,-1,@(t) 1);
            source.addLocalizedVelocityYTerm(By.grid,20,1,@(t) 1);
            zero = @(x,y) zeros(size(x));

            [P,Vx,Vy] = solver.solve(zero,zero,zero,Source=source);
            testCase.verifyGreaterThan(max(abs(P.values),[],'all'),0);
            testCase.verifyGreaterThan(max(abs(Vx.values),[],'all'),0);
            testCase.verifyGreaterThan(max(abs(Vy.values),[],'all'),0);
        end

        function rejectsLocalizedTermBuiltOnExtendedGrid(testCase)
            solver = testCase.makeSolver(2,[0,0.005],0.005);
            source = AcousticSource2D.zero();
            source.addLocalizedPressureTerm( ...
                solver.extendedGrid,3,1,@(t) 1);
            zero = @(x,y) zeros(size(x));

            testCase.verifyError(@() solver.solve( ...
                zero,zero,zero,Source=source), ...
                'AcousticSource2D:LocalizedGridMismatch');
        end

        function usesOrderDependentCFLBound(testCase)
            % dx=dy=0.025; dt=0.016 gives CFL about 0.905.
            times = [0,0.016];
            testCase.verifyWarningFree(@() testCase.makeSolver(2,times,0));
            testCase.verifyError(@() testCase.makeSolver(4,times,0), ...
                'AcousticSolver2D_PML:CFLViolation');
        end

        function rejectsInvalidReflectionCoefficient(testCase)
            testCase.verifyError(@() PML2D(ReflectionCoefficient=1), ...
                'PML2D:InvalidReflectionCoefficient');
        end
    end

    methods (Access = private)
        function solver = makeSolver(testCase,order,times,outputTimes)
            solver = AcousticSolver2D_PML( ...
                testCase.grid,times,testCase.grid,outputTimes, ...
                testCase.medium,testCase.pml,SpatialOrder=order);
        end
    end
end

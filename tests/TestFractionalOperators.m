function tests = TestFractionalOperators
% Run with: results = runtests('TestFractionalOperators.m');
tests = functiontests(localfunctions);
end

function testAgainstDenseMatrices(testCase)
rng(17);
n = 16;
dt = 0.2;
u = randn(n, 3) + 1i*randn(n, 3);
for alpha = [0, 0.5, 1, 1.5, 2, 3.2]
    % Independent generalized binomial/product formulas.
    d = ones(n, 1);
    q = ones(n, 1);
    for k = 1:n-1
        d(k+1) = prod(((0:k-1) - alpha)./(1:k));
        q(k+1) = prod(((0:k-1) + alpha)./(1:k));
    end
    D = dt^(-alpha)*tril(toeplitz(d, [1, zeros(1,n-1)]));
    I = dt^alpha*tril(toeplitz(q, [1, zeros(1,n-1)]));
    verifyEqual(testCase, fractionalDerivative(u,dt,alpha), D*u, 'AbsTol',1e-10);
    verifyEqual(testCase, fractionalIntegral(u,dt,alpha), I*u, 'AbsTol',1e-10);
    verifyEqual(testCase, fractionalDerivative(u,dt,alpha,Transpose=true), D'*u, 'AbsTol',1e-10);
    verifyEqual(testCase, fractionalIntegral(u,dt,alpha,Transpose=true), I'*u, 'AbsTol',1e-10);
    recovered = fractionalIntegral(fractionalDerivative(u,dt,alpha),dt,alpha);
    verifyEqual(testCase, recovered, u, 'AbsTol',1e-9);
    recovered = fractionalDerivative(fractionalIntegral(u,dt,alpha),dt,alpha);
    verifyEqual(testCase, recovered, u, 'AbsTol',1e-9);
end
end

function testIntegerOrders(testCase)
u = [2; 3; -1; 4];
dt = 0.25;
verifyEqual(testCase, fractionalDerivative(u,dt,1), [u(1); diff(u)]/dt);
verifyEqual(testCase, fractionalIntegral(u,dt,1), dt*cumsum(u));
verifyEqual(testCase, fractionalIntegral(u,dt,2), dt^2*cumsum(cumsum(u)));
end

function testDimensionsAndEdges(testCase)
u = reshape(1:24, [4,3,2]);
for op = {@fractionalDerivative, @fractionalIntegral}
    f = op{1};
    v = f(u,0.1,0.5,Dimension=2);
    for j = 1:2
        for k = 1:4
            verifyEqual(testCase, v(k,:,j), f(u(k,:,j),0.1,0.5), 'AbsTol',1e-12);
        end
    end
    verifyEqual(testCase, f(u,0.1,0), double(u));
    verifyEqual(testCase, f(zeros(0,3),0.1,0.5), zeros(0,3));
    verifyEqual(testCase, f(2,0.25,0.5), 2*0.25^(0.5*(2*strcmp(func2str(f),'fractionalIntegral')-1)));
    row = 1:5;
    col = f(row.',0.1,0.5);
    verifyEqual(testCase, f(row,0.1,0.5), col.');
end
end

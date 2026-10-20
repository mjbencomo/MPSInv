function [x,info] = solver_CGLS(A,b,options)
%SOLVER_CGLS Solve a linear least-squares problem by CGLS.
%   X = SOLVER_CGLS(A,B) approximately solves
%
%       min_x ||A*x - B||_2,
%
%   where A is a LinearOp.
%
%   X = SOLVER_CGLS(A,B,Name=Value) specifies one or more options:
%
%       InitialGuess     Initial estimate. The default is the zero vector.
%       Tolerance        Relative tolerance for the normal residual.
%                        The default is 1e-6.
%       MaxIterations    Maximum number of CGLS iterations. The default is
%                        A.dimDomain.
%
%   [X,INFO] also returns convergence information in a structure with
%   fields:
%
%       flag                    0: converged
%                               1: maximum iterations reached
%                               2: numerical breakdown
%       converged               true when flag is zero
%       message                 termination description
%       iterations              number of completed iterations
%       dataResidualNorm        final ||B-A*X||_2
%       normalResidualNorm      final ||A'*(B-A*X)||_2
%       dataResidualHistory     data-residual norms at X0,...,X_k
%       normalResidualHistory   normal-residual norms at X0,...,X_k
%
%   Convergence occurs when
%
%       ||A'*(B-A*X_k)||_2 <= Tolerance*||A'*(B-A*X_0)||_2.
%
%   The implementation uses only A.fwd and A.adj. It does not form the
%   matrix representation of A or the normal-equations matrix A'*A.

    arguments
        A (1,1) LinearOp
        b (:,1) {mustBeNumeric,mustBeFinite}
        options.InitialGuess {mustBeNumeric,mustBeFinite} = []
        options.Tolerance (1,1) double ...
            {mustBeReal,mustBeFinite,mustBePositive} = 1e-6
        options.MaxIterations (1,1) double ...
            {mustBeInteger,mustBePositive} = A.dimDomain
    end

    if numel(b) ~= A.dimRange
        error('solver_CGLS:InvalidDataVector', ...
            'b must be a column vector with A.dimRange entries.');
    end

    if isempty(options.InitialGuess)
        x = zeros(A.dimDomain,1,'like',b);
    else
        x = options.InitialGuess;
        if ~iscolumn(x) || numel(x) ~= A.dimDomain
            error('solver_CGLS:InvalidInitialGuess', ...
                ['InitialGuess must be a column vector with ' ...
                 'A.dimDomain entries.']);
        end
    end

    maxit = options.MaxIterations;

    r = b - A.fwd(x);
    s = A.adj(r);
    p = s;

    sts = real(s'*s);
    normalResidual = sqrt(sts);
    stoppingThreshold = options.Tolerance*normalResidual;

    dataResHist = zeros(maxit + 1,1);
    normalResHist = zeros(maxit + 1,1);
    dataResHist(1) = norm(r);
    normalResHist(1) = normalResidual;

    if normalResidual == 0
        info = makeInfo(0,'Initial guess is a least-squares solution.', ...
            0,r,normalResidual,dataResHist(1),normalResHist(1), ...
            stoppingThreshold);
        return
    end

    flag = 1;
    message = 'Maximum number of iterations reached.';
    iter = 0;

    for k = 1:maxit
        q = A.fwd(p);
        qtq = real(q'*q);

        if ~isfinite(qtq) || qtq <= 0
            flag = 2;
            message = 'CGLS encountered a zero or nonfinite search image.';
            break
        end

        alpha = sts/qtq;
        x = x + alpha*p;
        r = r - alpha*q;
        s = A.adj(r);

        stsNew = real(s'*s);
        normalResidual = sqrt(stsNew);

        iter = k;
        dataResHist(k + 1) = norm(r);
        normalResHist(k + 1) = normalResidual;

        if ~isfinite(normalResidual)
            flag = 2;
            message = 'CGLS encountered a nonfinite normal residual.';
            break
        end

        if normalResidual <= stoppingThreshold
            flag = 0;
            message = 'Relative normal residual satisfies the tolerance.';
            break
        end

        beta = stsNew/sts;
        p = s + beta*p;
        sts = stsNew;
    end

    dataResHist = dataResHist(1:iter + 1);
    normalResHist = normalResHist(1:iter + 1);
    info = makeInfo(flag,message,iter,r,normalResidual, ...
        dataResHist,normalResHist,stoppingThreshold);
end

function info = makeInfo(flag,message,iterations,r,normalResidual, ...
        dataResHist,normalResHist,stoppingThreshold)
    info = struct( ...
        'flag',flag, ...
        'converged',flag == 0, ...
        'message',message, ...
        'iterations',iterations, ...
        'dataResidualNorm',norm(r), ...
        'normalResidualNorm',normalResidual, ...
        'dataResidualHistory',dataResHist, ...
        'normalResidualHistory',normalResHist, ...
        'stoppingThreshold',stoppingThreshold);
end

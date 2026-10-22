function [x,info] = solver_PCGLS(A,b,M,options)
%SOLVER_PCGLS Preconditioned CGLS (MPS notes, Algorithm 22).
%   [X,INFO] = SOLVER_PCGLS(A,B,M) minimizes ||A*X-B||_2 using A.fwd,
%   A.adj, and M.inv. No matrices or explicit inverses are assembled.
%
%   A is a LinearOp. M is a fixed Hermitian positive-definite InvLinearOp
%   of dimension A.dimDomain, approximating A'*A. M is the full normal-
%   equation preconditioner, not a factor L with M=L'*L. Positive
%   definiteness is a caller requirement; local checks cannot certify it.
%   The notes assume A has full column rank.
%
%   Name-value options:
%     InitialGuess    Column vector; [] means zero (default).
%     TrueSolution    Reference column vector; [] disables error tracking.
%     RelativeTolerance Relative normal-residual tolerance, default 1e-6.
%     AbsoluteTolerance Absolute normal-residual tolerance, default 0.
%                       Both tolerances are nonnegative and may be zero.
%                       Stop when norm(g_k) <= max(AbsoluteTolerance, ...
%                           RelativeTolerance*norm(g_0)), where g=A'*r.
%                       The initial guess is tested with the same rule.
%                       Data and preconditioned residuals are diagnostics.
%                       This replaces Algorithm 22's absolute tests and
%                       the former Tolerance option.
%     MaxIterations   Maximum completed updates, default A.dimDomain.
%
%   INFO fields:
%     flag                       0: tolerance met; 1: iteration limit;
%                                2: numerical/preconditioner breakdown
%     converged                  true when flag is zero
%     message                    termination description
%     iterations                 number of completed updates
%     stoppingThreshold          max(AbsTol,RelTol*norm(g_0))
%     dataResidualNorm           final norm(r)
%     normalResidualNorm         final norm(g)
%     preconditionedNormalResidualNorm
%                                final sqrt(real(g'*M.inv(g)))
%     dataResidualHistory        norm(r_k), k=0,...,iterations
%     normalResidualHistory      norm(g_k), k=0,...,iterations
%     preconditionedNormalResidualHistory
%                                sqrt(real(g_k'*M.inv(g_k)))
%     trueErrorNorm              final norm(X-TrueSolution), or []
%     trueErrorHistory           norm(X_k-TrueSolution), or []
%
%   All histories include the initial and final iterate. A preconditioned
%   residual norm is NaN if its inner product fails the local checks.
%   For complex inputs, adjoints and inner products are Hermitian.
%   TrueSolution never changes the updates or termination tests.

    arguments
        A (1,1) LinearOp
        b (:,1) {mustBeNumeric,mustBeFinite}
        M (1,1) InvLinearOp
        options.InitialGuess {mustBeNumeric,mustBeFinite} = []
        options.TrueSolution {mustBeNumeric,mustBeFinite} = []
        options.RelativeTolerance (1,1) double ...
            {mustBeReal,mustBeFinite,mustBeNonnegative} = 1e-6
        options.AbsoluteTolerance (1,1) double ...
            {mustBeReal,mustBeFinite,mustBeNonnegative} = 0
        options.MaxIterations (1,1) double ...
            {mustBeInteger,mustBePositive} = A.dimDomain
    end

    if numel(b) ~= A.dimRange
        error('solver_PCGLS:InvalidDataVector', ...
            'b must have A.dimRange entries.');
    end
    if M.dimDomain ~= A.dimDomain || M.dimRange ~= A.dimDomain
        error('solver_PCGLS:InvalidPreconditionerDimension', ...
            'M must map A.dimDomain-dimensional vectors to themselves.');
    end

    x = options.InitialGuess;
    if isempty(x)
        x = zeros(A.dimDomain,1,'like',b);
    elseif ~iscolumn(x) || numel(x) ~= A.dimDomain
        error('solver_PCGLS:InvalidInitialGuess', ...
            'InitialGuess must be a column vector of length A.dimDomain.');
    end
    xTrue = options.TrueSolution;
    trackTrueError = ~isempty(xTrue);
    if trackTrueError && (~iscolumn(xTrue) || numel(xTrue) ~= A.dimDomain)
        error('solver_PCGLS:InvalidTrueSolution', ...
            'TrueSolution must be a column vector of length A.dimDomain.');
    end

    maxit = options.MaxIterations;
    dataHist = zeros(maxit+1,1);
    normalHist = zeros(maxit+1,1);
    preHist = NaN(maxit+1,1);
    errorHist = [];
    if trackTrueError
        errorHist = zeros(maxit+1,1);
    end

    r = b - A.fwd(x);
    g = A.adj(r);
    stoppingThreshold = max(options.AbsoluteTolerance, ...
        options.RelativeTolerance*norm(g));
    flag = 1;
    message = 'Maximum number of iterations reached.';

    % At entry to pass k, x and r are the kth iterate and residual.
    % Testing before the next update also handles a stationary x0.
    for k = 0:maxit
        dataHist(k+1) = norm(r);
        normalHist(k+1) = norm(g);
        if trackTrueError
            errorHist(k+1) = norm(x-xTrue);
        end

        if ~isfinite(dataHist(k+1)) || ~isfinite(normalHist(k+1))
            flag = 2;
            message = 'PCGLS encountered a nonfinite residual.';
            break
        end

        w = M.inv(g);
        innerProduct = g'*w;
        rho = real(innerProduct);
        % Permit roundoff in a Hermitian inner product. This only checks
        % the direction encountered; it does not prove that M is HPD.
        innerScale = norm(g)*norm(w);
        imagTolerance = 100*eps(class(innerProduct))*innerScale;
        if ~isfinite(innerProduct) || ~isfinite(innerScale) || ...
                abs(imag(innerProduct)) > imagTolerance || rho < 0 || ...
                (rho == 0 && normalHist(k+1) > 0)
            flag = 2;
            message = ['Invalid preconditioned inner product: M must be ' ...
                       'fixed Hermitian positive definite.'];
            break
        end
        preHist(k+1) = sqrt(rho);

        if normalHist(k+1) <= stoppingThreshold
            flag = 0;
            message = 'Normal residual satisfies the relative/absolute tolerance.';
            break
        elseif k == maxit
            break
        end

        if k == 0
            z = w;
        else
            beta = rho/rhoPrevious;
            z = w + beta*z;
        end
        q = A.fwd(z);
        qtq = real(q'*q);
        if ~isfinite(qtq) || qtq <= 0
            flag = 2;
            message = 'PCGLS encountered a zero or nonfinite search image.';
            break
        end

        alpha = rho/qtq;
        x = x + alpha*z;
        r = r - alpha*q;
        g = A.adj(r);
        rhoPrevious = rho;
    end

    dataHist = dataHist(1:k+1);
    normalHist = normalHist(1:k+1);
    preHist = preHist(1:k+1);
    trueErrorNorm = [];
    if trackTrueError
        errorHist = errorHist(1:k+1);
        trueErrorNorm = errorHist(end);
    end
    info = struct( ...
        'flag',flag, ...
        'converged',flag == 0, ...
        'message',message, ...
        'iterations',k, ...
        'stoppingThreshold',stoppingThreshold, ...
        'dataResidualNorm',dataHist(end), ...
        'normalResidualNorm',normalHist(end), ...
        'preconditionedNormalResidualNorm',preHist(end), ...
        'dataResidualHistory',dataHist, ...
        'normalResidualHistory',normalHist, ...
        'preconditionedNormalResidualHistory',preHist, ...
        'trueErrorNorm',trueErrorNorm, ...
        'trueErrorHistory',errorHist);
end

function [x, iter, res_hist] = CGLS(A, AT, b, x0, tol, maxit)
% CGLS: Conjugate Gradient for Least Squares. This is taken from one of the
% books Bencomo gave me. 
% Solves min ||Ax - b||^2
%
% Inputs:
%   A     - function handle for forward operator
%   AT    - function handle for adjoint operator
%   b     - data vector
%   x0    - initial guess
%   tol   - stopping tolerance on ||s||
%   maxit - maximum iterations
%
% Outputs:
%   x         - reconstructed solution
%   iter      - iteration stopped on
%   r         - residual

 x = x0;

    r = b - A(x);
    s = AT(r);

    p = s;
    sts = s' * s;

    res_hist = zeros(maxit,1);

    for k = 1:maxit

        q = A(p);

        alpha = sts / (q' * q);
        x = x + (alpha * p);
        r = r - (alpha * q);

        s = AT(r);

        sts_new = s' * s;

        res_hist(k) = norm(s);   

        if res_hist(k) < tol
            res_hist = res_hist(1:k);
            iter = k;
            return;
        end

        beta = sts_new / sts;

        p = s + beta * p;

        sts = sts_new;
    end

    iter = maxit;
    res_hist = res_hist(1:maxit);
end
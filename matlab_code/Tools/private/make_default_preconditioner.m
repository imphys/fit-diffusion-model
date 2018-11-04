function [R ] = make_default_preconditioner( H, x)
% [R ] = make_default_preconditioner( H, x)
% Created default preconditioner for the fmin_fast optimization routine.
% Does a cholesky decomposition of H, with some adjustments when H is not
% positive definite.
% 
% Created by Dirk Poot, Erasmus MC
% 6-2-2013

[R, info] = chol(H);
if info>0
    diagH = full(diag(H));
    if any(~isfinite(diagH))
        maxlpcount = -1;
    else 
        maxlpcount = 4;
    end;
    lambda = .01* max(abs(diagH)); % set initial lambda to 1% of a very rough estimate of the maximum eigenvalue
    n = size(H,1);
    lpcnt = 0;
    while info>0 && lpcnt <= maxlpcount
        H = H + lambda * eye(n);
        [R, info] = chol(H);
        lambda = lambda * 4;
        lpcnt = lpcnt + 1;
    end;
    if info>0
        % we seem to be unable to find a preconditioner
        % Probably not all elements of H are finite
        % (or maybe non symmetric??)
        % total lambda =approx = 4^4*.01*max(abs(eig(H)) 
        %                      = 2.56 * max(abs(eig(H)) 
        % So even when H was negative definite, we should have solved
        % it before exitting the loop above.
        if any(isfinite(diagH))
            R = eye(n) * max(abs(diagH(isfinite(diagH)))) ;
        else
            R = eye(n);
        end;
    end;
end;


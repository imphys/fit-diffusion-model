function [R ] = make_advanced_preconditioner( H, x, minhess)
% [R ] = make_advanced_preconditioner( H, x)
% Creates an (slightly) advanced preconditioner, especially for the MCMC
% sampling routine, but is also applicable for the fmin_fast optimization routine.
% Does a cholesky decomposition of H, with some adjustments when H is not
% positive definite or badly conditioned.
% 
% Created by Dirk Poot, Erasmus MC
% 5-8-2014

if nargin<4 
    minhess = [];
end;
n= size(H,1);
if ~isempty( minhess )
%     diagH = full(diag(H));
%     adj = diagH<minhess;
%     if any(adj)
%         H = H + diag( adj .* (minhess-diagH) );
%     end;

    % make sure curvature is at least minhess =>
    % for all x with x(k)==1,  x'*H*x > minhess(k)
    % [V,E]=eig(H)
    % xtest = V*t   such that   Vk * t == 1
    % t'*V'*H*V*t == t'*V'*V*E*t    (V'*V == eye, since H is symmetric and real)
    %  == t'*E*t 
    %  minimize_t  sum_i (lambda_i * t_i^2) with sum_i Vk_i * t_i ==1
    % => t = tbase + null(Vk)*topt  
    %  minimize_topt  sum_i (lambda_i * (tbase_i+ sum_j null(Vk)_ij topt_j)^2)
    %  = minimize_topt sum_ijk (lambda_i * (tbase_i/n+ null(Vk)_ij topt_j) (tbase_i/n+ null(Vk)_ik topt_k))
    %  = minimize_topt sum_ijk (lambda_i * tbase_i/n tbase_i/n + 
    %                           lambda_i * tbase_i/n null(Vk)_ik topt_k +
    %                           lambda_i * tbase_i/n null(Vk)_ij topt_j +
    %                           lambda_i * null(Vk)_ik topt_k null(Vk)_ij topt_j)
    % =  minimize_topt  2*tbase' * E* null(Vk) topt
    %                  + topt' * null(Vk)'*E*null(Vk) * topt           
    % => topt = -(null(Vk)'*E*null(Vk)) \ (null(Vk)'*E* tbase)
    % Thus minimal x'*H*x with x(k)==1 is given by
    % ( tbase - null(Vk)*(null(Vk)'*E*null(Vk)) \ (null(Vk)'*E* tbase) )' *
    %  E * ( tbase - null(Vk)*(null(Vk)'*E*null(Vk)) \ (null(Vk)'*E* tbase) )
    % =  tbase'*E* tbase + 
    %  -2*(null(Vk)'*E* tbase)' * inv( null(Vk)'*E*null(Vk) ) *null(Vk)' * E *tbase
    %  + (null(Vk)'*E* tbase)' * inv( null(Vk)'*E*null(Vk) ) *null(Vk)' * E * null(Vk)*inv( null(Vk)'*E*null(Vk) ) * (null(Vk)'*E* tbase)
    % = tbase'*E* tbase - tbase' * E* null(Vk) * inv( null(Vk)'*E*null(Vk) ) * null(Vk)' * E *tbase
    %
    % That's quite complicated and expensive to compute (requires
    % eigenvalue decomposition and a loop over the parameters)
    % So let's try a slightly modified criterium:
    %
    % make sure curvature is at least minhess =>
    % for all x with x'*diag(minhess)*x==1,  x'*H*x > 1
    %  x = diag(1/sqrt(minhess))* y  =>
    % for all y with y'*y==1,  y'*diag(1./sqrt(minhess))*H*diag(1./sqrt(minhess))*y > 1
    % Hm = diag(1./sqrt(minhess))*H*diag(1./sqrt(minhess));
    % y'*Hm * y >1
    % [V,E]=eig(Hm)
    % y = V*t
    % t'*V'*Hm*V*t == t'*V'*V*E*t == t'*E*t
    % y'*y == t'*V'*V*t == t'*t ==1
    %
    %
    % Thats a lot easier: the lowest eigenvalue of the modified hessian
    % (Hm) should be larger than 1.
    % Could solve that in (mainly) 2 ways:
    % 1) compute explicit eigenvalue decomposition of Hm
    %    Add extra cost (to Hm) in directions of too small eigenvalues 
    % 2) approximate lowest eigenvalue
    %    add minhess * (1 - lowest eigenvalue) to H.
    % The first one leaves all large eigenvectors unaffected, which
    % potentially allows better steps, while the second one might be
    % computationally easier/faster, especially for large scale problems.
    
    % How to compute 1:
    % option a : using eig
    %   Hm =  diag(1./sqrt(minhess))*H*diag(1./sqrt(minhess))
    %   lambtest = min(eig(Hm));
    %   H = H + diag(minhess)*min(0, 1-lambtest )
    %
    % option b (?): using eig with B argument:
    %   [V,D]= eig( H, B) 
    %    => A*V = B* V *D
    %     ??
    % 
    
    if any(~isfinite(H))
        % for non finite H, minhess might actually be a better description.
        % Multiply minhess by (the arbitrary factor) 10 to compensate for the fact that minhess is a lower bound and make sure 
        % that the computed step will not be much to large. The hessian was not finite, so it's relatively likely that the gradient
        % is not that accurate as well.
        % For MCMC there is another reason to multiply minhess by a factor > 1:
        %    acceptprob =  P( xnew) / P( x ) * Q( x | xnew) / Q( xnew | x)
        %  When the hessian is not finite, probably the likelihood P( xnew ) is quite low
        %  This low likelihood should in that case not be compensated by a high Q( x | xnew).
        H = spdiags(minhess(:)*10,0,n,n);
    else
        recipminhess = spdiags(1./sqrt(minhess(:)),0,n,n);
        Hm =  recipminhess*H*recipminhess;
        lambtest = min(eig(Hm));
        if lambtest<1
        	% adjust diagonal of H so that min(eig(Hm)) == 1 :
            H(1:n+1:n*n) = H(1:n+1:n*n) + (1-lambtest)*minhess(:);
        end;
    end;
else
    if any(~isfinite(H))
        d = diag(H);
        dm = max(d(isfinite(d)))*ones(n,1);
        H = spdiags(dm(:),0,n,n);
    else
        e = eig(H);
        if min(e)<max(abs(e))*.01
            shift = max(abs(e))*.01-min(e);
            H(1:n+1:n*n) = H(1:n+1:n*n) + shift;
        end;
    end;
end;    

[R, info] = chol(H);

if info>0 
    error('should not reach, H should have been fixed before');
end;


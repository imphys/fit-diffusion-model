function [condest, maxeigv, mineigv, approxmaxeigv, approxmineigv] = inverse_free_condest( K , n , t, maxFevals, Minv)
% [condest, maxeigv, mineigv, approxmaxeigv, approxmineigv] = inverse_free_condest( K , n , t, maxFevals, Minv)
% Estimate condition number of symmetric positive definite matrix K.
% Doesn't use a decomposition of K and thus is still (relatively) efficient for very 
% large scale problems. Number of multiplications with K scales with condition number.
% (Thus: large condition number -> large computation time)
% INPUTS:
% K or Kmul: K: n x n matrix, symmetric positive definite matrix. 
%            Kmul, function that evaluates: Kmul(x) = K * x
% n : the size of k
% t : number of test vectors, default [2 2]
%     time scales linearly with this number, 
%     1 typically is sufficient but 2 might be better.
% maxFevals : default = inf; 
%             maximum number of hessian multiplications that this function
%             is allowed to perform (to constrain computation time)
%             This (nearly) directly limits the maximum condition number
%             that can be found. Especially the minimum eigenvalue might be
%             a substantial overbound.
% Minv: Optional; preconditioner for K. To be a good preconditioner Minv
%       should approximate inv(K). The conditioner number of 
%        sqrtm(Minv) * K * sqrtm(Minv) is computed. The algorithm only requires
%        multiplications with Minv. Minv can be a matrix or function, just
%        like K.
%
% OUTPUTS:
% condest : lower bound for the condition number
% maxeigv : lower bound for the maximum eigenvalue (typically within 1%)
% mineigv : upper bound for the lowest eigenvalue
% approxmaxeigv : approximate eigenvector(s) corresponding to the maximum eigenvalue
% approxmineigv : approximate eigenvector(s) corresponding to the lowest eigenvalue
%
% The returned condition number estimate is a strict lower bound on the 
% actual condition number. (i.e. the actual condition number is strictly higher (or equal))
% The approximation is typically accurate to within a factor 3, but is better for lower condition numbers.
% (As the appoximation accuracy is not formally proven, there might be exceptions with a 
%  worse approximation, although I wouldn't know how to construct them)
% This function assumes a continuous eigenvalue spectrum  in the sense that
% if the eigenvalue spectrum has substantial gaps, faster approximations
% are possible. If t is only 1 in rare cases a single high/low eigenvalue might
% be missed. 
% The highest eigenvalue is accurate to within a percent (or better)
% The relative error in the lowest eigenvalue might be substantial 
% (can be overestimated by a factor 2 to 3)
%
% Created by Dirk Poot, Erasmus MC
% 9-12-2011
if nargin<2
    n = size(K,1);
end;
if nargin<3
    t = [2 2];
end;
if nargin<4
    maxFevals = inf;
end;
if nargin<5
    Minv = @(x) x ; % default, identity preconditioner. Takes essiantially no computation time.
    hasM = false;
else
    hasM = true;
end;

% Algorithm that is used in this function :
% (fact:) test vector(s) are a sum of eigenvectors of K :
%  tstv = sum_i  v_i * c_i
% where c_i are the (scalar/row vector) weights of each eigenvector v_i of
% K. Each eigenvector v_i has associated eigenvalue lambda_i.
% Since we assume K is symmetric (hermitian) and positive definite, the
% eigenvectors are orthogonal and all eigenvalues positive.  
%
% We aim to find largest and smallest lambda_i's of K and approximate the
% accompagnying v_i. This is done with the following recursion on the test
% vector(s) tstv:
% tstv_j+1 = K* tstv_j - z_j * tstv_j 
%      == (sum_i lambda_i  * v_i * (c_j)_i) - (sum_i z_j * v_i* (c_j)_i )
%      == sum_i (lambda_i-z_j)* v_i * (c_j)_i 
%      == sum_i (Prod_j (lambda_i-z_j) ) * v_i * c0_i
% with tstv_0 = c0_i * v_i
% Since the v_i are orthonormal, the c0_i are standard normal distributed with the
% random initialization that is used for tstv_0.
% Note that the final formula has a polynomial in lambda. 
% In each iteration the next eigenvalue that is zero'd is the point within
% the current range of eigenvalues that thus far is is maximally amplified.
% This point is (approximately) found by searching the maximum of  the
% logarithm of the absolute of (Prod_j (lambda_i-z_j)) for all 
% points in between all already zero'd values (between the elements of
% sort(z)). By chosing this value numerical problems are avoided (no
% eigenvectors can dominate) and the extrema of 
% the interval of lambda's that exist can be most easily found. 
%
% The maximum eigenvalue found so far is computed by evaluating the
% 'effective eigenvalue' at each iteration: tstv_j'*K*tstv_j , where tstv
% is (implicately) made orthonormal.  This assumes that K is symmetric as
% otherwise the dot product might exceed the largest eigenvalue. (since the
% eigenvectors are then not  orthogonal). 
% However if 96% of the current 'effective eigenvalue' is larger than the
% maximum z, the next zero is at 99% of that value (instead of the maximum
% in between value). Iteration stops when the maximum contribution of
% lambda<max(lambda) is less than 1% of the norm of tstv.
%
% The minimum eigenvalue is found by zeroing the range between the
% 1.05*times the maximum eigenvalue (since it might be an underestimate)
% and the 1.2 x the current upperbound of the lowest eigenvalue.
%
% When a positive definite symmetric preconditioner is present the
% (implicit) recursion becomes: 
% tstv_j+1 = sqrtm(inv(M))*K*sqrtm(inv(M)) tstv_j - z_j * tstv_j  
% with tilde{tstv_j} = sqrtm(M) * tstv_j and multiplying both sides on the
% left with sqrtm(M) we get:
% tilde{tstv_j+1}  = K * inv(M) * tilde{tstv_j} - z_j * tilde{tstv_j}
% and the test relation tstv_j'*K*tstv_j  becomes:
%  tilde{tstv_j}'*inv(M) *K*inv(M) * tilde{tstv_j} 
% where tstv_j can be orthonormalized with the cholesky decomposition of 
%  tstv_j'*tstv_j == tilde{tstv_j}'*inv(M) * tilde{tstv_j} 
%

if isa(K,'function_handle')
    Kmul = K;
    clear K;
else
    Kmul = @(x) K*x;
end;
if isa(Minv,'function_handle')
    Minvmul = Minv;
    clear M
else
    Minvmul = @(x) Minv*x;
end;

numHMuls = 0; % number of multiplications that have been performed.
tstv = randn(n,t(1));
tstv = tstv * sqrtm(inv(tstv'*tstv)); % orthonormalize test vectors.
maxeigv = -inf;
zerodeigv = [];
it = 0;
nextzero = 0;
%%
while true
%%
    it = it+1;
    if numHMuls>=maxFevals
        break;
    end;
    Mtstv = Minvmul(tstv);
    KMtstv = Kmul(Mtstv);
    numHMuls = numHMuls+1;
%     tstev =  max(diag(tstvm'*tstv));
    if hasM
        [R]=chol(tstv'*Mtstv);
    else
        [dummy,R]=qr(tstv,0);
    end;
    tstev = max( eig( R'\ ((Mtstv'*KMtstv) / R )) ); % locate maximum eigenvalue in span of tstv. 
    tstv = KMtstv - tstv * nextzero;
    tstv = tstv * diag(1./sqrt(diag(tstv'*tstv))); % normalize. TODO: It might be possible to use the normalization constants.
    zerodeigv(end+1)= nextzero;

    if tstev > maxeigv
        maxeigv = tstev;
        maxeigv
    end;
    eigvsc = sum(log(abs(maxeigv-zerodeigv)));
    if max(zerodeigv)<.96*maxeigv
        nextzero = .99*maxeigv;
    else
        zerodeigv = sort(zerodeigv(:));
        cent = 0.5*(zerodeigv(1:end-1)+zerodeigv(2:end));
        centsc = sum(log(abs(bsxfun(@minus, zerodeigv, cent'))),1);
        [a, idx] = max(centsc);
        if a<eigvsc-4.6
            % maximum contribution of small eigenvalues is small enough (I think)
            disp('break');
            break;
        else
            % chop down the maximum.
            nextzero = cent(idx); 
        end;
    end;
%     nextzero 
end;
% plot(eigvalk, tstv)
approxmaxeigv = tstv;
it_eigmax = it;
%% we have found a good lower bound for the maximum eigenvalue.
% typically actual maximum is less than 1 % larger

% now bracket lowest eigenvalue in a similar way, but 
% convergence limit is much tighter, as the relative error in the lowest
% eigenvalue should not be larger than it's magnitude.
% note that convergence limit is specified in relative values.

tstv = randn(n,t(end));
tstv = tstv * sqrtm(inv(tstv'*tstv));
mineigv = inf;
approxmineigv = [];
zerodeigv = zeros(10,1);
links = zeros(2,size(zerodeigv,1)); % [next;prev]
cent = zeros(1,size(zerodeigv,1));
centsc = zeros(1,size(zerodeigv,1));
%nextidx =0; % signals zeroing before current lowest zero
it = 0;
nextzero = 1.05*maxeigv;
largeeigvsuprsc = 1.1; % if 0.5: mean influenced approximately constantly by each eigenvector in 'blocked part' of eigenvalue spectrum
                        % if 1: variance influenced approximately constantly by each eigenvector in 'blocked part' of eigenvalue spectrum
                        % therefore slighly larger, so variance contribution is even localized in the 'blocked part' of eigenvalue spectrum
%%
while true
%%
    it = it+1;
    if it>numel(zerodeigv)
        % expand matrices
        zerodeigv(2*it+10,:)=0;
        links(:,size(zerodeigv,1))=0;
        cent(:,size(zerodeigv,1)-1)=0;
        centsc(:,size(zerodeigv,1)-1)=0;
    end;
    
    % multiply vector
    if numHMuls>=maxFevals
        disp('maximum number of multiplications exhausted; approximation of minimum eigenvector might be bad');
        break;
    end;

    Mtstv = Minvmul(tstv);
    KMtstv = Kmul(Mtstv);
    
    numHMuls = numHMuls +1;
    % determine current 'eigenvalue'
    if 1
        % Old version, works for symmetric positive definite matrices
        if hasM
            [R]=chol(tstv'*Mtstv);
        else
            [dummy,R]=qr(tstv,0);
        end;
        tstev = (Mtstv'*KMtstv);
        mintstev = min( eig( R'\ (tstev / R )) ); % locate maximum eigenvalue in span of tstv. 
    else
        % new version. Try to get it also working for non symmetric positive matrices
        tstev = (KMtstv'*KMtstv);
        mintstev = sqrt( min( diag(tstev)  ) );
    end;
    tstv_prev = tstv;
    
    % zero requested eigenvalue and scale residual:
    tstv = KMtstv - tstv * nextzero;
    tstv_norm2 = diag(tstv'*tstv);
    tstv = tstv * diag(1./sqrt(tstv_norm2));
    % indicate that we zerod the eigenvalue:
    zerodeigv(it)= nextzero;
    % update eigvsc
    if mintstev < mineigv
        % new upper bound for lowest eigenvalue found. Recompute eigvsc
        mineigv = mintstev;
        approxmineigv = tstv_prev;
        mineigv
        if mineigv<0
            disp('Negative upper bound for lowest eigenvalue found.')
            disp('Inverse free condest can only estimate condition number of positive definite matrices.');
            disp('Searching is stopped and current test vectors and condest==inf are returned.');
            condest = inf;
            return;
        end;
        eigvsc = sum(log(abs(mineigv-zerodeigv(1:it))))+largeeigvsuprsc*log(mineigv);
    else
        eigvsc = eigvsc + log(abs(mineigv-nextzero));
    end;
    % update centsc
    centsc(1:it-2) = centsc(1:it-2) + log(abs(cent(1:it-2)-nextzero));
    if it==1
        links(:,it)=[1;1];
    else 
        % nextidx == k => new zero between zerodeigv(k) and zerodeigv(nxt)
        % (except when nextidx ==1, then between 0 and zerodeigv(nxt) )
        % cent( i - 1) == center between zerodeigv(i) and zerodeigv(links(1,i))
        nxt = links(1,nextidx);
        links(1,[nextidx it])=[it nxt];
        links(2,[it nxt])=[nextidx it];
        cent(it-1)= 0.5*(zerodeigv(it) + zerodeigv(nxt));
        centsc(it-1) = sum(log(abs(zerodeigv(1:it)-cent(it-1))))+largeeigvsuprsc*log(cent(it-1));
        if nextidx ~= 1
            cent(nextidx-1)= 0.5*(zerodeigv(it) + zerodeigv(nextidx));
            centsc(nextidx-1) = sum(log(abs(zerodeigv(1:it)-cent(nextidx-1))))+largeeigvsuprsc*log(cent(nextidx-1));
        end;            
    end;
    
        % compute standarddeviation of current eigenvalue spectrum w.r.t. mineigv
        % => sqrt(diag( tmp ) ) with 
        %   tmp = ( tstvm - mineigv * tstv_prev)'*( tstvm - mineigv * tstv_prev) 
        %         ( (tstvm - nextzero * tstv_prev) + (nextzero - mineigv )* tstv_prev)'*( (tstvm - nextzero * tstv_prev) + (nextzero - mineigv )* tstv_prev) 
        %       = tstv_norm2 + 2* (nextzero - mineigv ) * (tstvm - nextzero * tstv_prev)' * tstv_prev  + (nextzero - mineigv )^2 * tstv_prev'*tstv_prev
        %       = tstv_norm2 + 2* (nextzero - mineigv ) * tstvm' * tstv_prev  + ((nextzero - mineigv )^2 - 2* (nextzero - mineigv ) *nextzero)  * tstv_prev'*tstv_prev
        stdeig = sqrt( tstv_norm2 + 2*(nextzero - mineigv ) * diag(tstev) + ((nextzero - mineigv )^2 - 2* (nextzero - mineigv )) );
%         if 
%         end;

    % determine what to do next:
    if min(zerodeigv(1:it))>2*mineigv || it==1
        % if lowest zero was more than 2x current lowest eigenvalue, add zero just after current lowest eigenvalue.
        if it==1
            nextzero = 1.01*mineigv;
        else
            nextzero = 1.2*mineigv;
        end;
        nextidx = 1;
    else
        % chop down the worst contributing block:
        [a, idx] = max(centsc(1:it-1));
        nextidx = idx+1;
        if a<eigvsc-6%4.6
            % maximum contribution of large eigenvalues is small enough (I think)
            disp('break');
            break;
        else
            % chop down the maximum.
            nextzero = cent(idx); 
        end;
    end;
%     nextzero 
end;
% plot(eigvalk, tstv)
if isempty(approxmineigv)
    approxmineigv = tstv;
end;
it_eigmin = it;

condest = maxeigv / mineigv;
disp(['Condition number estimate:' num2str(condest) ' #Kmul e_max : ' num2str(it_eigmax) ' #Kmul e_min : ' num2str(it_eigmin)]);
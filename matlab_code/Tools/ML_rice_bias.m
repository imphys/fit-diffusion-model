function [bias, CRLB, ir, jr] = ML_rice_bias(fun, thetaML, sigma, n, parampriorfun,  gradPrior_linindex, hessPrior_linindex, d3Prior_linindex)
% [bias, CRLB] = ML_rice_bias(fun, thetaML, sigma [, n [, parampriorfun , gradPrior_linindex, hessPrior_linindex, d3Prior_linindex]])
% Computes the bias in a ML estimator with rice distributed data
% fun : intensity predicting function that is able to compute first and second order derivatives
%       (Uses the exact same interface as fit_MRI for intensity predicting functions)
%         [A, dAdpar, d2Adparidparj] = fun( thetaML(:, index_range) )
%          A      = nMRI x numel(index_range) predicted magnitudes
%          dAdpar = nMRI x numel(index_range) x size(par,1) 
%                           derivative of each A with respect to (w.r.t.) each parameter.
%                           ( = Jacobian of function for each column of theta)
%          d2AdparIdparJ = nMRI x numel(index_range) x (numel(I)) 
%                           second derivative (hessian) of each A w.r.t. each parameter combination.
%                           d2 A / d par(I(k)) d par(J(k)) 
%                           with [I,J] = find(triu(ones(size(par,1))));
% thetaML : parameters for which you want to compute the bias that an ML estimator
%           will introduce.
% sigma : default=1; noise level in the measurements
%         if sigma = 'lastInTheta' : the sigma is the last element of thetaML, 
%                   equal for all measurements of that voxel, and estimated with the ML estimator. 
%               (corresponds to fit_MRI : 'numPDFoptpar' = 1 )
% n     : default=1; number of repeats of the measurements. 
%           Repeat count of the entire set of measurements (as returned by fun)
%           The bias scales with 1/n.
%           This assumes all n repeats are used in a single ML estimation.
% parampriorfun : optional parameter prior function that was used in the 
%                 ML estimation of the parameters. The bias influence this
%                 has is also corrected. Typically you want parampriorfun 
%                 to cause a reduction of the bias so that ML_rice_bias is
%                 able to find a better correction for the remaining bias.
%                 See fit_MRI : opts.parameterPrior for more information.
%  gradPrior_linindex, hessPrior_linindex, d3Prior_linindex : 
%                 linear index of the non-zero elements of the gradient
%                 hessian and third derivative of the parameter prior
%                 function. The hessian and third derivative are assumed to
%                 be fully symmetric and the indices should contain only 1
%                 copy of each value: 
%                 e.g.: hess_linindex = find( triu( ones(size(thetaML,1)) ) )
%                       (or a subset of that, or use tril)
%
% OUTPUTS:
% bias : size(thetaML), bias in thetaML
%        To correct for the bias in thetaML subtract it as in:
%        thetaML_corrected = thetaML - bias;
% CRLB : Cramer Rao lower bound of thetaML.
%        Requesting this output is almost for free since the CRLB of
%        thetaML is computed within the bias computation.
%
% Created by Dirk Poot, Erasmus MC, 21-3-2012

if nargin<3 || isempty(sigma)
    sigma = 1;
end;
if nargin<4 || isempty(n)
    n = 1;
end;
hasparamprior = nargin>=5 && ~isempty(parampriorfun);

sigmaInTheta = isequal( sigma ,'lastInTheta');

sztht = size(thetaML);
szsigma = size(sigma);
spatialsigma = isequal(sztht(2:end),szsigma(2:end)) && ~sigmaInTheta;
sigmasel = sigma; % default fill of sigma
numtr = prod(sztht(2:end));
npar = sztht(1) - sigmaInTheta;
mat = true(npar,npar);
[ i , j ] = find(mat);
selmask = logical(triu(mat));
[ ir , jr ] = find(selmask);
rep2 =zeros(npar,npar);
rep2(ir+(jr-1)*npar)=1:numel(ir);
rep2(jr+(ir-1)*npar)=1:numel(ir);
rep = ones(1,numel(i));
% see 'some_formulas.pdf' ('some_formulas.tex') for a derivation of the equations
maxnuminblock = 100;
bias = zeros(size(thetaML));
if nargout>1
    selmaske = logical(triu(ones( sztht(1) )));
    CRLB = zeros( [ nnz(selmaske), sztht(2:end) ] );
end;
if hasparamprior
% parse: gradPrior_linindex, hessPrior_linindex, d3Prior_linindex
    [tmp_i, tmp_j] = ind2sub([sztht(1) sztht(1)], hessPrior_linindex);
    tmp_msk = find(tmp_i~=tmp_j);
    hessPrior_assignindex = sub2ind( [sztht(1) sztht(1)], [tmp_i;tmp_i(tmp_msk)], [tmp_j;tmp_j(tmp_msk)] );
    hessPrior_readindex = [ (1:numel(hessPrior_linindex))';tmp_msk];
    
    prior_d3_k = zeros(sztht(1),sztht(1),sztht(1));
    tmp = cell(1,3);
    [tmp{:}] = ind2sub(sztht(1)*[1 1 1], d3Prior_linindex);
    prms = perms(1:3);
    for k=1:size(prms ,1)
        tmp_lin = sub2ind( sztht(1)*[1 1 1], tmp{prms(k,:)}); 
        prior_d3_k(tmp_lin) = 1:numel(d3Prior_linindex);
    end;
    d3Prior_assignindex = find(prior_d3_k);
    d3Prior_readindex = prior_d3_k( d3Prior_assignindex );
end;
stindx = 1;
% ws = warning('query');
progressbar('start',numtr,[],'MinTimeInterval',1);
while stindx<=numtr
    edindx = min(numtr, stindx + maxnuminblock);
    numtrinblock = edindx-stindx+1;
    thetasel = thetaML(:,stindx:edindx);
    if hasparamprior
        [prior_f, prior_g, prior_h, prior_d3] = parampriorfun( thetasel ) ;
    end;    
    if sigmaInTheta
        sigmasel = thetasel(end,:);
        thetasel = thetasel(1:end-1,:);
    elseif spatialsigma
        sigmasel = sigma(:,stindx:edindx);
    end;
    [A, dAdpar, d2Adpar] = fun( thetasel ); % first and second derivative required.
    dAdpar = permute(dAdpar,[1 3 2]);
    d2Adpar = permute(d2Adpar,[1 3 2]);
    
    if hasparamprior
        if sigmaInTheta
            [expectd2IA, expectd2IAs, expectd2Iss, expectd3IA , expectd3IsAA, expectd3IAAs, expectd3IAss, expectd3IsAs, expectd3Isss , expectd3bIAAA, expectd3bIAAs, expectd3bIAss, expectd3bIsss] = riceExpectFuncs_s( A , sigmasel) ;
        else
            [expectd2IA, expectd3IA , expectd3bIAAA] = riceExpectFuncs_s( A , sigmasel, [1 4 10]) ;
        end;
    else
        if sigmaInTheta
            [expectd2IA, expectd2IAs, expectd2Iss, expectd3IA , expectd3IsAA, expectd3IAAs, expectd3IAss, expectd3IsAs, expectd3Isss ] = riceExpectFuncs_s( A , sigmasel) ;
        else
            [expectd2IA, expectd3IA ] = riceExpectFuncs( A , sigmasel) ;
        end;
    end;
    ws2 = warning('off','MATLAB:nearlySingularMatrix');
    ws3 = warning('off','MATLAB:illConditionedMatrix');
    ws4 = warning('off','MATLAB:singularMatrix');
    if hasparamprior
        for k=1:numtrinblock
            I = reshape( expectd2IA(:,k)'* (dAdpar(:,i,k).*dAdpar(:,j,k)), npar,npar);
            if sigmaInTheta
                Its = expectd2IAs(:,k)'*dAdpar(:,:,k);
                Iss = sum(expectd2Iss(:,k));
                I = [I Its';Its Iss]; %#ok<AGROW>
            end;
            Ip_grad = I;
            Ip_grad(gradPrior_linindex,gradPrior_linindex) =  Ip_grad(gradPrior_linindex,gradPrior_linindex) + 0* prior_g(:,k) * prior_g(:,k)' ;
            I(hessPrior_assignindex) = I(hessPrior_assignindex) - prior_h(hessPrior_readindex,k);
%             Ip_grad = I; % DEBUG: should they be different?
            invI = inv(I);
            if nargout>1
                CRLB(:,k+stindx-1) = invI(selmaske);
            end;        

            d2Adpar_f = reshape(d2Adpar(:,rep2,k),[size(A,1) npar npar]);

            sum3_i = dAdpar(:,:,k)*invI(1:npar,:); % sum first derivative of S with invI, so can be sum over i,j,k

            invI2 = invI * Ip_grad * invI; % sum over m and n
            sum_kmn = dAdpar(:,:,k) * invI2(1:npar,:);
            
            % Compute J, part : dS dS dS * (- df^3 /f^2 + df d2f / f)
            sum_ko = sum( sum3_i(:,1:npar) .* dAdpar(:,:,k) ,2); % sum over k and o (where sum3_i is already summed over one of them)
            bias_J1 = sum( bsxfun(@times, sum_ko .* (- expectd3bIAAA(:,k) - 2*expectd3IA(:,k)) , sum3_i) , 1);
            % Compute J, part : dS/dj dS/didk * (df^2 /f )
            sum_ko2 = sum( bsxfun(@times, reshape(sum3_i(:,1:npar), [size(A,1) 1 npar]),  d2Adpar_f) , 3);
            bias_J2 = (expectd2IA(:,k)'* sum_ko2) * invI(1:npar,:); % sum over l and i

            % Compute K, part : dS dS dS * (2 df^3 /f^2 - 3 df d2f / f)
            sum_mnko = sum( sum_kmn(:,1:npar) .* dAdpar(:,:,k) ,2);
            bias_K1 = sum( bsxfun(@times, sum_mnko .* (2 * expectd3bIAAA(:,k) + 6*expectd3IA(:,k)) , sum3_i));
            % Compute K, part : A * (df^2 /f )
            %    with A = d2S/dido dS/dk + d2S/dodk dS/di + d2S/didk dS/do 
            Aterm = bsxfun(@times, reshape( d2Adpar_f , [size(A,1) npar*npar])*reshape(invI2(1:npar,1:npar),[],1) , sum3_i) + ... % sum d2S/dodk dS/di over mnoki (so only j and l are left)
                2 * sum( bsxfun(@times, reshape(dAdpar(:,:,k) * invI2(1:npar,1:npar), [size(A,1) 1 npar]),  d2Adpar_f ) ,3 ) * invI(1:npar,:); % sum d2S/didk dS/do+ d2S/dido dS/dk  over mnoki
            bias_K2 = - (expectd2IA(:,k)'* Aterm) ;
            
            % explicity create f''' and reduce that; we could do a bit more
            % efficient, but thats not easy and this computation is not
            % the most time consuming (I assume)
            prior_d3_k( d3Prior_assignindex ) = prior_d3( d3Prior_readindex );
            biasf3 = ((reshape(prior_d3_k,[sztht(1), sztht(1)*sztht(1)])*invI2(:)))'* invI; % sum over mnok and i
            
            if sigmaInTheta
                % we already computed with j==sigma (as this does not influence the integrals)
                % so compute parts with i == sigma; o==sigma; k==sigma
                %                       i&o == sigma; i&k==sigma;o&k==sigma
                %                       i&o&k==sigma
                % Compute J part : dS dS dS * (- df^3 /f^2 + df d2f / f)
%                 bias_J1 = bias_J1 + (sum_ko'*(- expectd3bIAAs(:,k) - 2*expectd3IAAs(:,k))) * invI(end,:)... % with i => sigma
%                          + 2* sum( bsxfun(@times, sum3_i(:,end) .*(- expectd3bIAAs(:,k) - expectd3IAAs(:,k)-expectd3IsAA(:,k)) , sum3_i) , 1)...  % with k (or o) => sigma
%                          + 2* ( sum3_i(:,end)'*(- expectd3bIAss(:,k) - expectd3IAss(:,k)-expectd3IsAs(:,k)) ) * invI(end,:) ...; %  with i and k (or o) => sigma
%                          +  (- expectd3bIAss(:,k) - 2*expectd3IsAs(:,k))'* sum3_i *  invI(end,end)...%  with k and o => sigma
%                          +  sum(- expectd3bIsss(:,k) - 2*expectd3Isss(:,k) ,1 )* invI(end,end) * invI(end,:); % with i and o and k == sigma
                bias_J1 = bias_J1 ...
                        + sum( bsxfun(@times, sum3_i(:,end) .*(- expectd3bIAAs(:,k) - 2* expectd3IAAs(:,k)) , sum3_i) , 1)...  % (1) with k => sigma
                        + (sum_ko'*(- expectd3bIAAs(:,k) - 2*expectd3IAAs(:,k))) * invI(end,:)... % (2) with i => sigma
                        + ( sum3_i(:,end)'*(- expectd3bIAss(:,k) - 2*expectd3IAss(:,k)) ) * invI(end,:) ... % (3) with i and k => sigma
                        + sum( bsxfun(@times, sum3_i(:,end) .*(- expectd3bIAAs(:,k) - 2* expectd3IsAA(:,k) ) , sum3_i) , 1)...  % (4) with o => sigma
                        +  (- expectd3bIAss(:,k) - 2*expectd3IsAs(:,k))'* sum3_i *  invI(end,end) ...% (5)  with k and o => sigma
                        + ( sum3_i(:,end)'*(- expectd3bIAss(:,k) - 2*expectd3IsAs(:,k)) ) * invI(end,:) ... % (6) with i and o => sigma
                        +  sum(- expectd3bIsss(:,k) - 2*expectd3Isss(:,k) ,1 )* invI(end,end) * invI(end,:); % (7) with i and o and k == sigma
                % Compute J, part : dS/dj dS/didk * (df^2 /f )
                bias_J2 = bias_J2  ... %0:  with i==sigma or k ==sigma : second derivative wrt S is zero.
                         + invI(end,1:npar) * reshape( expectd2IAs(:,k)'* reshape(d2Adpar_f, [], npar*npar), npar, npar) * invI(1:npar,:); % (4) with o == sigma
%                 bias_J2 = bias_J2  ... %0:  with i==sigma or k ==sigma : second derivative wrt S is zero.
%                          + invI(end,1:npar) * reshape( expectd2IAs(:,k)'* reshape(d2Adpar_f, [], npar*npar), npar, npar) * invI(1:npar,:); % with o == sigma
                          % with higher derivatives, second derivative wrt S is zero
                          
                % Compute K, part : dS dS dS * (2 df^3 /f^2 - 3 df d2f / f)
%                 bias_K1 = bias_K1 + (sum_mnko'*(2* expectd3bIAAs(:,k) + 6 *expectd3IAAs(:,k))) * invI(end,:)... % with i => sigma
%                          + 2* sum( bsxfun(@times, sum_kmn(:,end) .*(2* expectd3bIAAs(:,k) + 3* expectd3IAAs(:,k) +3 * expectd3IsAA(:,k)) , sum3_i) , 1)...  % with k (or o) => sigma
%                          + 2* ( sum_kmn(:,end)'*(2* expectd3bIAss(:,k) +3* expectd3IAss(:,k) +3*expectd3IsAs(:,k)) ) * invI(end,:) ...; %  with i and k (or o) => sigma
%                          +  (2* expectd3bIAss(:,k) + 6*expectd3IsAs(:,k))'* sum3_i *  invI2(end,end)...%  with k and o => sigma
%                          +  sum(2* expectd3bIsss(:,k) + 6*expectd3Isss(:,k) ,1 )* invI2(end,end) * invI(end,:); % with i and o and k == sigma
                bias_K1 = bias_K1 ...
                         + sum( bsxfun(@times, sum_kmn(:,end) .*(2* expectd3bIAAs(:,k) + 4* expectd3IAAs(:,k) +2 * expectd3IsAA(:,k)) , sum3_i) , 1)...  % (1) with k  => sigma
                         + (sum_mnko'*(2* expectd3bIAAs(:,k) + 4 *expectd3IAAs(:,k) + 2 * expectd3IsAA(:,k))) * invI(end,:)... % (2) with i => sigma
...%                          + sum( sum_kmn(:,end) .*(2* expectd3bIAss(:,k) + 2* expectd3IAss(:,k) +2 * expectd3IsAs(:,k))  , 1) * invI(end,:) ...  % (3) with k and i => sigma
                         + sum( sum_kmn(:,end) .*(2* expectd3bIAss(:,k) + 2* expectd3IAss(:,k) +4 * expectd3IsAs(:,k))  , 1) * invI(end,:) ...  % (3) with k and i => sigma
                         + sum( bsxfun(@times, sum_kmn(:,end) .*(2* expectd3bIAAs(:,k) + 4* expectd3IAAs(:,k) +2 * expectd3IsAA(:,k)) , sum3_i) , 1)...  % (4) with o => sigma
                         + (2* expectd3bIAss(:,k) + 4*expectd3IsAs(:,k) + 2*expectd3IAss(:,k))'* sum3_i *  invI2(end,end)...%  (5) with k and o => sigma
                         + sum_kmn(:,end)'*(2* expectd3bIAss(:,k) + 4* expectd3IsAs(:,k) +2*expectd3IAss(:,k))  * invI(end,:) ...; % (6) with i and o => sigma
                         + sum(2* expectd3bIsss(:,k) + 6*expectd3Isss(:,k) ,1 )* invI2(end,end) * invI(end,:); % (7) with i and o and k == sigma
                % Compute K, part : dS/dj dS/didk * (df^2 /f )
                bias_K2 = bias_K2 ...
                         - invI2(end,1:npar) * reshape( expectd2IAs(:,k)'* reshape(d2Adpar_f, [], npar*npar), npar, npar) * invI(1:npar,:)...  % (1) with k == sigma
                         - (( expectd2IAs(:,k)'* reshape(d2Adpar_f, [], npar*npar)) * reshape(invI2(1:npar,1:npar),[],1)) * invI(end,:) ... % (2) with i==sigma 
                         - invI2(end,1:npar) * reshape( expectd2IAs(:,k)'* reshape(d2Adpar_f, [], npar*npar), npar, npar) * invI(1:npar,:); % (4) with k == sigma
%                 bias_K2 = bias_K2 + (( expectd2IAs(:,k)'* reshape(d2Adpar_f, [], npar*npar)) * reshape(invI2(1:npar,1:npar),[],1)) * invI(end,:) ... % with i==sigma 
%                          + 2* invI2(end,1:npar) * reshape( expectd2IAs(:,k)'* reshape(d2Adpar_f, [], npar*npar), npar, npar) * invI(1:npar,:); % with k (or o) == sigma
                          % with higher derivatives, second derivative wrt S is zero
                          

                % validate with fully expanded code:
                npnt = size(A,1);
                dAdpark = dAdpar(:,:,k);
                %            l      o     i      k
                refJ = zeros(npnt, npar+1,npar+1,npar+1);
                refK = zeros(npnt, npar+1,npar+1,npar+1);
                
                refJ(:,1:npar,1:npar,1:npar) = bsxfun(@times, bsxfun(@times, bsxfun(@times, -expectd3bIAAA(:,k), dAdpark), reshape(dAdpark,[npnt 1 npar 1])), reshape(dAdpark,[npnt 1 1 npar])) ...
                                              +bsxfun(@times, bsxfun(@times, bsxfun(@times, -2*expectd3IA(:,k) , dAdpark), reshape(dAdpark,[npnt 1 npar 1])), reshape(dAdpark,[npnt 1 1 npar])) ...
                                              +bsxfun(@times,                bsxfun(@times, expectd2IA(:,k),dAdpark),                                 reshape(d2Adpar_f,[npnt 1 npar npar]));             
                % k->s
                refJ(:,1:npar,1:npar,npar+1) = bsxfun(@times,               bsxfun(@times, -expectd3bIAAs(:,k),  dAdpark), reshape(dAdpark,[npnt 1 npar 1])) ...
                                              +bsxfun(@times, bsxfun(@times, -2*expectd3IAAs(:,k),dAdpark), reshape(dAdpark,[npnt 1 npar 1]));
                % i->s
                refJ(:,1:npar,npar+1,1:npar) = bsxfun(@times,               bsxfun(@times, -expectd3bIAAs(:,k),  dAdpark), reshape(dAdpark,[npnt 1 1 npar])) ...
                                              +bsxfun(@times, bsxfun(@times,-2*expectd3IAAs(:,k),dAdpark), reshape(dAdpark,[npnt 1 1 npar]));
                % k&i->s
                refJ(:,1:npar,npar+1,npar+1) = bsxfun(@times, -expectd3bIAss(:,k),  dAdpark ) ...
                                              +bsxfun(@times, -2*expectd3IAss(:,k), dAdpark);
                %o->s                         
                refJ(:,npar+1,1:npar,1:npar) = bsxfun(@times,               bsxfun(@times, -expectd3bIAAs(:,k),  reshape(dAdpark,[npnt 1 npar 1])), reshape(dAdpark,[npnt 1 1 npar])) ...
                                              +bsxfun(@times, bsxfun(@times,-2*expectd3IsAA(:,k),reshape(dAdpark,[npnt 1 npar 1])), reshape(dAdpark,[npnt 1 1 npar])) ...
                                              +bsxfun(@times, expectd2IAs(:,k), reshape(d2Adpar_f,[npnt 1 npar npar])) ;
                %o&k->s                         
                refJ(:,npar+1,1:npar,npar+1) = bsxfun(@times, -expectd3bIAss(:,k),  reshape(dAdpark ,[npnt 1 npar 1])) ...
                                              +bsxfun(@times, -2*expectd3IsAs(:,k),reshape(dAdpark,[npnt 1 npar 1]));
                %o&i->s
                refJ(:,npar+1,npar+1,1:npar) = bsxfun(@times, -expectd3bIAss(:,k),  reshape(dAdpark ,[npnt 1 1 npar])) ...
                                              +bsxfun(@times, -2*expectd3IsAs(:,k),  reshape(dAdpark ,[npnt 1 1 npar]) );
                %o&i&k->s
                refJ(:,npar+1,npar+1,npar+1) = -expectd3bIsss(:,k) ...
                                               -2*expectd3Isss(:,k);
                
                refK(:,1:npar,1:npar,1:npar) = bsxfun(@times, bsxfun(@times, bsxfun(@times, 2*expectd3bIAAA(:,k), dAdpark), reshape(dAdpark,[npnt 1 npar 1])), reshape(dAdpark,[npnt 1 1 npar])) ...
                                              +bsxfun(@times, bsxfun(@times, bsxfun(@times, 6*expectd3IA(:,k)   , dAdpark), reshape(dAdpark,[npnt 1 npar 1])), reshape(dAdpark,[npnt 1 1 npar])) ...
                                              +bsxfun(@times,                bsxfun(@times, -expectd2IA(:,k),dAdpark),                                   reshape(d2Adpar_f,[npnt 1 npar npar])) ...
                                              +bsxfun(@times,                bsxfun(@times, -expectd2IA(:,k),reshape(d2Adpar_f,[npnt npar 1 npar])),   reshape(dAdpark,[npnt 1 npar 1])                   ) ...
                                              +bsxfun(@times,                bsxfun(@times, -expectd2IA(:,k),reshape(d2Adpar_f,[npnt npar npar 1])),   reshape(dAdpark,[npnt 1 1 npar])                   ) ;
                % k->s
                refK(:,1:npar,1:npar,npar+1) = bsxfun(@times,               bsxfun(@times, 2*expectd3bIAAs(:,k),  dAdpark), reshape(dAdpark,[npnt 1 npar 1])) ...
                                              +bsxfun(@times, bsxfun(@times, 4*expectd3IAAs(:,k),dAdpark), reshape(dAdpark,[npnt 1 npar 1]))...
                                              +bsxfun(@times, bsxfun(@times, 2*expectd3IsAA(:,k),dAdpark), reshape(dAdpark,[npnt 1 npar 1]))...
                                              +bsxfun(@times, -expectd2IAs(:,k), reshape(d2Adpar_f,[npnt npar npar 1])) ;
                % i->s
                refK(:,1:npar,npar+1,1:npar) = permute(refK(:,1:npar,1:npar,npar+1),[1 2 4 3]);
                %o->s 
                refK(:,npar+1,1:npar,1:npar) = permute(refK(:,1:npar,1:npar,npar+1),[1 4 2 3]);
                
                % i&k->s
                refK(:,1:npar,npar+1,npar+1) = bsxfun(@times, 2*expectd3bIAss(:,k),  dAdpark)  ...
                                              +bsxfun(@times, 2*expectd3IAss(:,k),dAdpark)...
                                              +bsxfun(@times, 4*expectd3IsAs(:,k),dAdpark);
                % o&k->s                         
                refK(:,npar+1,1:npar,npar+1) = permute(refK(:,1:npar,npar+1,npar+1),[1 3 2 4]);
                % o&i->s
                refK(:,npar+1,npar+1,1:npar) = permute(refK(:,1:npar,npar+1,npar+1),[1 3 4 2]);
                %o&i&k->s
                refK(:,npar+1,npar+1,npar+1) = 2*expectd3bIsss(:,k) ...
                                              +6*expectd3Isss(:,k);
                
                
                bias_refJ = sum( reshape( reshape( permute(  refJ   ,[1 3 2 4]),[npnt*(npar+1) , (npar+1)*(npar+1)]) * invI(:)  , [npnt (npar+1)])* invI ,1);
                refK(end+1,:,:,:) = reshape(prior_d3_k,[1 npar+1 npar+1, npar+1]);
                bias_refK = sum( reshape( reshape( permute(  refK   ,[1 3 2 4]),[(npnt+1)*(npar+1) , (npar+1)*(npar+1)]) * invI2(:) , [(npnt+1) (npar+1)])* invI ,1);
                bias_ref = prior_g(:,k)'* invI(gradPrior_linindex,:) + bias_refJ + .5*bias_refK;
%                 bias_ref_nosigma = sum( reshape( reshape( permute(ref(:,1:npar,1:npar,1:npar),[1 3 2 4]),[npnt*(npar) , (npar)*(npar)])*reshape(invI(1:npar,1:npar),[],1) , [npnt (npar)])* invI(1:npar,:) ,1);
                
%                 disp( ( bias_d3part + bias_d2part) - bias_ref )
            end;
            
            % compute bias vector including sigma part; but ignoring any
            % sigma part in i o k 
            bias_k_unscaled = prior_g(:,k)'* invI(gradPrior_linindex,:) + bias_J1 + bias_J2 + .5*(bias_K1+ bias_K2 + biasf3);

%             if sigmaInTheta
%                 disp( bias_k_unscaled - bias_ref )
%                 bias_k_unscaled = bias_ref;
%             end;
            bias(:,k+stindx-1) = bias_k_unscaled(:)/n;
        end;
    else 
        % not hasparamprior
        for k=1:numtrinblock
            I = reshape( expectd2IA(:,k)'* (dAdpar(:,i,k).*dAdpar(:,j,k)), npar,npar);
            if sigmaInTheta
                Its = expectd2IAs(:,k)'*dAdpar(:,:,k);
                Iss = sum(expectd2Iss(:,k));
                I = [I Its';Its Iss]; %#ok<AGROW>
            end;
            invI = inv(I);
            if nargout>1
                CRLB(:,k+stindx-1) = invI(selmaske);
            end;        
            
            d2Adpar_f = reshape(d2Adpar(:,rep2,k),[size(A,1) npar npar]);
            
            sum3_i = dAdpar(:,:,k)*invI(1:npar,:); % sum first derivative of S with invI, so can be sum over i,j,k
            sum3_ko = sum(dAdpar(:,:,k).*sum3_i(:,1:npar), 2); % summed over entire invI, with first derivative: summed over j,k for 3rd derivative integral
            
            bias_d3part = (sum3_ko.*expectd3IA(:,k))'* sum3_i; 
            
            Aterm = -.5* bsxfun(@times, reshape( d2Adpar_f , [size(A,1) npar*npar])*reshape(invI(1:npar,1:npar),[],1) , sum3_i) ; % sum d2S/dodk dS/di over ok and i (so only j and l are left)
                   % + sum( bsxfun(@times,reshape(sum3_i(:,1:npar), [size(A,1) 1 npar]),d2Adpar_f ) ,3 ) * invI(1:npar,:); % sum d2S/didk dS/do over ko and i + sum d2S/dido dS/dk  over ok and i
            bias_d2part = expectd2IA(:,k)'* Aterm;
%             sum2_jk = reshape( sum( bsxfun(@times, sum3_i , reshape(d2Adpar(:,rep2,k),[size(A,1) npar npar])) ,2), [size(A,1) npar]);

%             sumA_jk = 2*sum2_jk + bsxfun(@times, dAdpar(:,:,k), d2Adpar(:,rep2,k)*invI(:));

%             if sigmaInTheta
%                 invIts = invI(1:npar,end)';
%                 invIss = invI(end,end);
%                 invIe = invI(1:npar,:);
%                 invItse = invI(end,:);
%                 invI = invI(1:npar,1:npar); % theta x theta part
%             end;
%             bias_k_unscaled = sum( bsxfun(@times, sum3_i , sum3_jk.*expectd3IA(:,k) ), 1) ...
%                   + (expectd2IA(:,k)' * sum2_jk(:,:) - .5*( expectd2IA(:,k)'*sumA_jk(:,:)) ) *invI ;
            if sigmaInTheta
                % we already computed with j==sigma (as this does not influence the integrals)
                % so compute parts with i == sigma; o==sigma; k==sigma
                %                       i&o == sigma; i&k==sigma;o&k==sigma
                %                       i&o&k==sigma
                % Compute J part : dS dS dS * (- df^3 /f^2 + df d2f / f)
%                 bias_d3part = bias_d3part + (sum3_ko' * expectd3IsAA(:,k) ) * invI(end,:) ... % with i => sigma
%                          +  sum( bsxfun(@times, sum3_i(:,end) .*((1-1)*expectd3IsAA(:,k) + 2* expectd3IAAs(:,k) ) , sum3_i) , 1)...  % with k or o => sigma
%                          +  (sum3_i(:,end)'*(expectd3IAss(:,k)+expectd3IsAs(:,k)) ) * invI(end,:) ...; %  with i and k (or o) => sigma
%                          +  (expectd3IAss(:,k)'* sum3_i) *  invI(end,end)...%  with k and o => sigma
%                          +  sum( expectd3Isss(:,k) ,1 )* invI(end,end) * invI(end,:); % with i and o and k == sigma
                     
                bias_d3part = bias_d3part ...
                         +  sum( bsxfun(@times, sum3_i(:,end) .* expectd3IsAA(:,k)  , sum3_i) , 1)...  % (1) with k  => sigma
                         +  (expectd3IsAA(:,k)'*sum3_ko )*invI(end,:) ...  % (2) with  i => sigma
                         +  (sum3_i(:,end)'*( 2*expectd3IsAs(:,k) - expectd3IAss(:,k) ) ) * invI(end,:) ...; %  (3) with i and k => sigma
                         +  sum( bsxfun(@times, sum3_i(:,end) .* (2*expectd3IAAs(:,k)-expectd3IsAA(:,k)) , sum3_i) ,1)  ...; % (4) with o => sigma
                         +  (expectd3IAss(:,k)'* sum3_i) *  invI(end,end)...% (5)  with k and o => sigma
                         +  (expectd3IAss(:,k)'* sum3_i(:,end)) *  invI(end,:)...% (6) with i and o => sigma
                         +  sum( expectd3Isss(:,k) ,1 )* invI(end,end) * invI(end,:); % (7) with i and o and k == sigma
                     
                bias_d2part = bias_d2part ...
                         -.5* (invI(end,1:npar)* reshape( expectd2IAs(:,k)'* reshape(d2Adpar_f, [], npar*npar), npar, npar) ) * invI(1:npar,:) ... % (1) with k => sigma 
                         -.5* (( expectd2IAs(:,k)'* reshape(d2Adpar_f, [], npar*npar)) * reshape(invI(1:npar,1:npar),[],1)) * invI(end,:) ... % (2) with i => sigma 
                         +.5* (invI(end,1:npar)* reshape( expectd2IAs(:,k)'* reshape(d2Adpar_f, [], npar*npar), npar, npar) ) * invI(1:npar,:) ; % (3) with o==sigma 
                          % with higher derivatives, second derivative wrt S is zero     
                     
%                 tmpts = -.5*(expectd2IAs(:,k)'* d2Adpar(:,:,k));
%                 tmpts = reshape(tmpts(rep2),[npar ,npar]) ;
%                 F112 = reshape( sum( bsxfun(@times, dAdpar(:,i,k).*dAdpar(:,j,k), expectd3IAAs(:,k)) ,1), [npar npar]) + tmpts ;
%                 F211 = reshape( sum( bsxfun(@times, dAdpar(:,i,k).*dAdpar(:,j,k), expectd3IsAA(:,k)) ,1), [npar npar]) - tmpts ;
%                 F122 = expectd3IAss(:,k)'*dAdpar(:,:,k);
%                 F212 = expectd3IsAs(:,k)'*dAdpar(:,:,k);
%                 F222 = sum( expectd3Isss(:,k) );
%                 % * performes sum: j&k   i&keep s            j&k   i&keep s
%                 biasthtadj = invIts * F112 * invIe + (F112(:)'*invI(:)) * invItse ...
%                             +invIts * F211 * invIe ...  % maybe * 2
%                             +invIts * F122'* invItse ...  % maybe * 2
%                             +invIss * F212 * invIe + ((invIts * F212') * invItse)...
%                             +invIss * F222 * invItse;
%                 bias_k_unscaled = [bias_k_unscaled 0]+biasthtadj;

                % validate with fully expanded code:
                npnt = size(A,1);
                dAdpark = dAdpar(:,:,k);
                %            l      o     i      k
                ref = zeros(npnt, npar+1,npar+1,npar+1);
                
                ref(:,1:npar,1:npar,1:npar) = bsxfun(@times, bsxfun(@times, bsxfun(@times, expectd3IA(:,k), dAdpark), reshape(dAdpark,[npnt 1 npar 1])), reshape(dAdpark,[npnt 1 1 npar])) ...
                                             +bsxfun(@times,                bsxfun(@times, .5*expectd2IA(:,k),dAdpark),                                 reshape(d2Adpar_f,[npnt 1 npar npar])) ...
                                             +bsxfun(@times,                bsxfun(@times,-.5*expectd2IA(:,k),reshape(d2Adpar_f,[npnt npar 1 npar])),   reshape(dAdpark,[npnt 1 npar 1])                   ) ...
                                             +bsxfun(@times,                bsxfun(@times,-.5*expectd2IA(:,k),reshape(d2Adpar_f,[npnt npar npar 1])),   reshape(dAdpark,[npnt 1 1 npar])                   );
                % k->s
                ref(:,1:npar,1:npar,npar+1) = bsxfun(@times, bsxfun(@times, expectd3IsAA(:,k),dAdpark), reshape(dAdpark,[npnt 1 npar 1]))...
                                             +bsxfun(@times, -.5*expectd2IAs(:,k), reshape(d2Adpar_f,[npnt npar npar 1])) ;
                % i->s
                ref(:,1:npar,npar+1,1:npar) = bsxfun(@times, bsxfun(@times, expectd3IsAA(:,k),dAdpark), reshape(dAdpark,[npnt 1 1 npar]))...
                                             +bsxfun(@times, -.5*expectd2IAs(:,k), reshape(d2Adpar_f,[npnt npar 1 npar])) ;
                % k&i->s
                ref(:,1:npar,npar+1,npar+1) = bsxfun(@times, -expectd3IAss(:,k), dAdpark) ...
                                             +bsxfun(@times, 2*expectd3IsAs(:,k), dAdpark) ;
                %o->s                         
                ref(:,npar+1,1:npar,1:npar) = bsxfun(@times, bsxfun(@times, 2*expectd3IAAs(:,k),reshape(dAdpark,[npnt 1 npar 1])), reshape(dAdpark,[npnt 1 1 npar])) ...
                                             +bsxfun(@times, bsxfun(@times, - expectd3IsAA(:,k),reshape(dAdpark,[npnt 1 npar 1])), reshape(dAdpark,[npnt 1 1 npar])) ...
                                             +bsxfun(@times, .5*expectd2IAs(:,k), reshape(d2Adpar_f,[npnt 1 npar npar])) ;
                %o&k->s                         
                ref(:,npar+1,1:npar,npar+1) = bsxfun(@times,   expectd3IAss(:,k),reshape(dAdpark,[npnt 1 npar 1]));
                %o&i->s
                ref(:,npar+1,npar+1,1:npar) = bsxfun(@times,   expectd3IAss(:,k),reshape(dAdpark,[npnt 1 1 npar]));
                %o&i&k->s
                ref(:,npar+1,npar+1,npar+1) =  expectd3Isss(:,k);
                
                bias_ref = sum( reshape( reshape( permute(ref,[1 3 2 4]),[npnt*(npar+1) , (npar+1)*(npar+1)])*invI(:) , [npnt (npar+1)])* invI ,1);
                bias_ref_nosigma = sum( reshape( reshape( permute(ref(:,1:npar,1:npar,1:npar),[1 3 2 4]),[npnt*(npar) , (npar)*(npar)])*reshape(invI(1:npar,1:npar),[],1) , [npnt (npar)])* invI(1:npar,:) ,1);
                
%                 disp( ( bias_d3part + bias_d2part) - bias_ref )
                bias_d3part = bias_ref;bias_d2part =0;
            end;
            bias_k_unscaled = bias_d3part + bias_d2part;
            bias(:,k+stindx-1) = bias_k_unscaled(:)/n;
        end;
    end;
    warning(ws4);
    warning(ws3);
    warning(ws2)
    stindx = edindx+1;
    progressbar(edindx);
end;
progressbar('ready');
